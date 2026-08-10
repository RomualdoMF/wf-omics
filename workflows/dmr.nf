// Orchestrates the two opt-in DMR (Differentially Methylated Region) comparisons
// built on top of ont-methylDMR-kit's DSS-based calling (see modules/local/dmr.nf
// for the actual processes and the rationale for hooking in pre-concat_bedmethyl):
//
//   - dmr_haplotype_compare: haplotype 1 vs haplotype 2 of the same sample.
//     Only meaningful with --phased. Uses workflow mod's modkit_H1_per_chr /
//     modkit_H2_per_chr emits directly.
//   - dmr_sample_compare: the sample vs an external reference, given as either a
//     BAM/CRAM (piled up here, reusing workflows/methyl.nf's own sample_probs/
//     modkit processes) or an already-built bedmethyl (split by chromosome here
//     instead). Uses workflow mod's sample_bedmethyl_per_chr emit (the
//     haplotype-agnostic/combined side, regardless of --phased).
//
// Both run DSS for 5mC and 5hmC (the two modified bases workflows/methyl.nf's
// own modkit invocation produces). Annotation (gencode/imprinted) and the HTML
// report always run once DMRs are called -- same as upstream ont-methylDMR-kit,
// where those aren't optional either. Only plotting (modbamtools/methylartist)
// is toggled, via dmr_haplotype_compare_args / dmr_sample_compare_args -- a
// literal subset of upstream's own CLI flag vocabulary (--plot,
// --methylartist_only, --gene_list <path>, --imprinted).

include {
    prep_dmr_bedmethyl;
    split_bedmethyl_by_chr;
    list_reference_chromosomes;
    dmr_calling;
    aggregate_dmrs;
    download_dmr_annotations;
    annotate_dmrs;
    prepare_gene_list;
    report_dmrs;
    plot_dmr_modbamtools;
    plot_dmr_methylartist;
    plot_phased_dmr_modbamtools;
    plot_phased_dmr_methylartist;
} from '../modules/local/dmr.nf'

include { sample_probs; modkit } from './methyl.nf'

// The modified bases workflows/methyl.nf's modkit pileups for (see workflow mod
// in that file, --modified-bases 5mC 5hmC) -- "m"/"h" are modkit's own bedMethyl
// mod-code characters.
def DMR_MOD_TYPES = ["m", "h"]

// Parses the literal subset of upstream ont-methylDMR-kit CLI flags this
// integration understands out of a raw args list
// (params.dmr_haplotype_compare_args / params.dmr_sample_compare_args).
// Annotation + report are unconditional (see header comment); only plotting is
// gated by these.
def dmr_stage_flags(List args) {
    def gene_list_idx = args.indexOf("--gene_list")
    return [
        plot: args.contains("--plot"),
        methylartist_only: args.contains("--methylartist_only"),
        imprinted: args.contains("--imprinted"),
        gene_list: (gene_list_idx >= 0 && gene_list_idx + 1 < args.size()) ? args[gene_list_idx + 1] : "",
    ]
}


workflow dmr_haplotype_compare {
    take:
        modkit_H1_per_chr    // tuple(meta[alias,sq,...], group='1', bedmethyl.gz) -- workflow mod emit
        modkit_H2_per_chr    // tuple(meta, group='2', bedmethyl.gz)
        phased_bam           // tuple(name, bam, bai) -- the haplotagged BAM, for plotting
        reference             // ref_channel: tuple(fasta, fai, cache, REF_PATH)
        args                  // List<String>, params.dmr_haplotype_compare_args
    main:
        def flags = dmr_stage_flags(args)
        def comparison_label = "haplotype_compare"

        h1_keyed = modkit_H1_per_chr | map { meta, group, bed -> tuple(meta.sq, meta.alias, bed) }
        h2_keyed = modkit_H2_per_chr | map { meta, group, bed -> tuple(meta.sq, bed) }
        // [sq, alias, bed_h1, bed_h2]
        paired = h1_keyed | join(h2_keyed, by: 0)
        alias_ch = paired | map { sq, alias, bed1, bed2 -> alias } | first()

        annotations_ready = download_dmr_annotations().ready
        gene_list = prepare_gene_list(annotations_ready, flags.imprinted, flags.gene_list).gene_list
        // gencode GTF used for plotting comes from the same annotations tarball as
        // annotate_dmrs' BED -- resolved only once the download is confirmed ready.
        gencode_annotation = annotations_ready | map { ready ->
            tuple(
                file("${params.dmr_annotations_dir}/annotations/gencode.v46.GRCh38.annotation.sorted.gtf.gz"),
                file("${params.dmr_annotations_dir}/annotations/gencode.v46.GRCh38.annotation.sorted.gtf.gz.tbi")
            )
        }

        dmr_tables = Channel.empty()
        annotated_tables = Channel.empty()
        reports = Channel.empty()
        plots = Channel.empty()

        for (mod_char in DMR_MOD_TYPES) {
            prepped1 = prep_dmr_bedmethyl(paired.map { sq, alias, bed1, bed2 -> tuple(sq, bed1) }, mod_char)
            prepped2 = prep_dmr_bedmethyl(paired.map { sq, alias, bed1, bed2 -> tuple(sq, bed2) }, mod_char)
            chr_pairs = prepped1.dss_bed | join(prepped2.dss_bed, by: 0)

            dmr_out = dmr_calling(chr_pairs)
            agg_key = alias_ch | map { alias -> tuple(alias, comparison_label, mod_char) }
            agg_out = aggregate_dmrs(
                agg_key,
                dmr_out.chr_dmrs | map { chr, f -> f } | collect,
                dmr_out.status_log | collect
            )
            dmr_tables = dmr_tables | mix(agg_out.dmr_table)

            ann_out = annotate_dmrs(agg_out.dmr_table, annotations_ready, flags.imprinted)
            annotated_tables = annotated_tables | mix(ann_out.annotated)

            rep_out = report_dmrs(ann_out.annotated, ann_out.annotation_log)
            reports = reports | mix(rep_out.report)

            if (flags.plot) {
                if (!flags.methylartist_only) {
                    mb = plot_phased_dmr_modbamtools(ann_out.annotated, gencode_annotation.collect(), gene_list.collect(), phased_bam.collect())
                    plots = plots | mix(mb.dmr_plots)
                }
                ma = plot_phased_dmr_methylartist(ann_out.annotated, gencode_annotation.collect(), gene_list.collect(), reference.collect(), phased_bam.collect())
                plots = plots | mix(ma.dmr_plots)
            }
        }
    emit:
        dmr_table = dmr_tables
        annotated = annotated_tables
        report = reports
        plots = plots
}


workflow dmr_sample_compare {
    take:
        sample_bedmethyl_per_chr  // tuple(meta, group='*', bedmethyl.gz) -- workflow mod emit (combined side)
        reference_path             // String: params.dmr_sample_compare_reference (BAM/CRAM or bedmethyl)
        chromosome_codes           // ArrayList, e.g. ["chr1","1",...,"chrX","X"]
        reference                   // ref_channel: tuple(fasta, fai, cache, REF_PATH)
        sample_bam                  // tuple(name, bam, bai) -- the sample's own whole-genome BAM, for plotting
        args                         // List<String>, params.dmr_sample_compare_args
    main:
        def flags = dmr_stage_flags(args)
        def comparison_label = "sample_compare"
        def lower_ref = reference_path.toLowerCase()
        def is_bam = lower_ref.endsWith('.bam') || lower_ref.endsWith('.cram')

        sample_keyed = sample_bedmethyl_per_chr | map { meta, group, bed -> tuple(meta.sq, meta.alias, bed) }
        alias_ch = sample_keyed | map { sq, alias, bed -> alias } | first()

        // Reference side, raw (unfiltered by mod code) per-chromosome bedmethyl --
        // converges both branches to the same shape (tuple(chr, raw_bedmethyl))
        // consumed by prep_dmr_bedmethyl below, same as prep_dmr_bedmethyl does
        // for the sample side.
        def reference_available_for_plotting = is_bam
        if (is_bam) {
            def idx_ext = lower_ref.endsWith('.cram') ? '.crai' : '.bai'
            ref_bam_ch = Channel.fromPath(reference_path, checkIfExists: true)
                | map { it -> tuple(['alias': it.baseName], it, file("${it}${idx_ext}", checkIfExists: true)) }

            ref_probs = sample_probs(
                ref_bam_ch | map { meta, xam, xai -> tuple(xam, xai, meta) },
                reference.collect()
            )
            ref_chroms = list_reference_chromosomes(
                ref_bam_ch | map { meta, xam, xai -> tuple(xam, xai) },
                chromosome_codes
            ).chroms | splitText() | map { it.trim() } | filter { it }

            ref_modkit_input = ref_chroms
                | combine(ref_bam_ch)
                | combine(ref_probs.probs)
                | map { sq, meta, xam, xai, probs -> tuple(meta + [sq: sq, probs: probs], xam, xai) }

            reference_keyed = modkit(ref_modkit_input, reference.collect(), '--combine-strands --cpg').modkit
                | map { meta, group, bed -> tuple(meta.sq, bed) }
        } else {
            reference_keyed = split_bedmethyl_by_chr(
                Channel.fromPath(reference_path, checkIfExists: true),
                chromosome_codes
            ).chr_beds
                | flatten
                | map { f -> tuple(f.name.replaceFirst(/^chr_split__/, '').replaceFirst(/\.bed$/, ''), f) }
        }

        annotations_ready = download_dmr_annotations().ready
        gene_list = prepare_gene_list(annotations_ready, flags.imprinted, flags.gene_list).gene_list
        gencode_annotation = annotations_ready | map { ready ->
            tuple(
                file("${params.dmr_annotations_dir}/annotations/gencode.v46.GRCh38.annotation.sorted.gtf.gz"),
                file("${params.dmr_annotations_dir}/annotations/gencode.v46.GRCh38.annotation.sorted.gtf.gz.tbi")
            )
        }

        dmr_tables = Channel.empty()
        annotated_tables = Channel.empty()
        reports = Channel.empty()
        plots = Channel.empty()

        if (flags.plot && !reference_available_for_plotting) {
            log.warn "dmr_sample_compare: --plot requested but ${reference_path} is a bedmethyl, not a BAM/CRAM -- skipping plots (no reads to plot for the reference side)."
        }

        for (mod_char in DMR_MOD_TYPES) {
            sample_prepped = prep_dmr_bedmethyl(sample_keyed.map { sq, alias, bed -> tuple(sq, bed) }, mod_char)
            ref_prepped = prep_dmr_bedmethyl(reference_keyed, mod_char)
            chr_pairs = sample_prepped.dss_bed | join(ref_prepped.dss_bed, by: 0)

            dmr_out = dmr_calling(chr_pairs)
            agg_key = alias_ch | map { alias -> tuple(alias, comparison_label, mod_char) }
            agg_out = aggregate_dmrs(
                agg_key,
                dmr_out.chr_dmrs | map { chr, f -> f } | collect,
                dmr_out.status_log | collect
            )
            dmr_tables = dmr_tables | mix(agg_out.dmr_table)

            ann_out = annotate_dmrs(agg_out.dmr_table, annotations_ready, flags.imprinted)
            annotated_tables = annotated_tables | mix(ann_out.annotated)

            rep_out = report_dmrs(ann_out.annotated, ann_out.annotation_log)
            reports = reports | mix(rep_out.report)

            if (flags.plot && reference_available_for_plotting) {
                // ref_bam_ch was built above, inside the `if (is_bam)` branch --
                // in scope here since reference_available_for_plotting == is_bam.
                ref_bam_for_plot = ref_bam_ch | map { meta, xam, xai -> tuple(meta.alias, xam, xai) }
                if (!flags.methylartist_only) {
                    mb = plot_dmr_modbamtools(ann_out.annotated, gencode_annotation.collect(), gene_list.collect(), sample_bam.collect(), ref_bam_for_plot.collect())
                    plots = plots | mix(mb.dmr_plots)
                }
                ma = plot_dmr_methylartist(ann_out.annotated, gencode_annotation.collect(), gene_list.collect(), reference.collect(), sample_bam.collect(), ref_bam_for_plot.collect())
                plots = plots | mix(ma.dmr_plots)
            }
        }
    emit:
        dmr_table = dmr_tables
        annotated = annotated_tables
        report = reports
        plots = plots
}
