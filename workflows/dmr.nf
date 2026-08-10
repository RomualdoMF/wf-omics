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
//
// NOTE on the 5mC/5hmC ("m"/"h") duplication below: DSL2 forbids invoking the
// same process more than once from within a single workflow scope (calling it
// from a `for (mod_char in ["m","h"])` loop hits "Process has already been
// used"), so instead of looping, every process used per mod type is included
// twice under a `_m`/`_h` alias (see the include block) and each stage's logic
// is written out twice, once per mod type, rather than factored into a loop.

include {
    prep_dmr_bedmethyl as prep_dmr_bedmethyl_side1_m;
    prep_dmr_bedmethyl as prep_dmr_bedmethyl_side1_h;
    prep_dmr_bedmethyl as prep_dmr_bedmethyl_side2_m;
    prep_dmr_bedmethyl as prep_dmr_bedmethyl_side2_h;
    split_bedmethyl_by_chr;
    list_reference_chromosomes;
    dmr_calling as dmr_calling_m;
    dmr_calling as dmr_calling_h;
    aggregate_dmrs as aggregate_dmrs_m;
    aggregate_dmrs as aggregate_dmrs_h;
    download_dmr_annotations;
    annotate_dmrs as annotate_dmrs_m;
    annotate_dmrs as annotate_dmrs_h;
    prepare_gene_list;
    report_dmrs as report_dmrs_m;
    report_dmrs as report_dmrs_h;
    plot_dmr_modbamtools as plot_dmr_modbamtools_m;
    plot_dmr_modbamtools as plot_dmr_modbamtools_h;
    plot_dmr_methylartist as plot_dmr_methylartist_m;
    plot_dmr_methylartist as plot_dmr_methylartist_h;
    plot_phased_dmr_modbamtools as plot_phased_dmr_modbamtools_m;
    plot_phased_dmr_modbamtools as plot_phased_dmr_modbamtools_h;
    plot_phased_dmr_methylartist as plot_phased_dmr_methylartist_m;
    plot_phased_dmr_methylartist as plot_phased_dmr_methylartist_h;
} from '../modules/local/dmr.nf'

include { sample_probs; modkit } from './methyl.nf'

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

        // ---- 5mC ("m") ----
        prepped1_m = prep_dmr_bedmethyl_side1_m(paired.map { sq, alias, bed1, bed2 -> tuple(sq, bed1) }, "m")
        prepped2_m = prep_dmr_bedmethyl_side2_m(paired.map { sq, alias, bed1, bed2 -> tuple(sq, bed2) }, "m")
        chr_pairs_m = prepped1_m.dss_bed | join(prepped2_m.dss_bed, by: 0)

        dmr_out_m = dmr_calling_m(chr_pairs_m)
        agg_key_m = alias_ch | map { alias -> tuple(alias, comparison_label, "m") }
        agg_out_m = aggregate_dmrs_m(
            agg_key_m,
            dmr_out_m.chr_dmrs | map { chr, f -> f } | collect,
            dmr_out_m.status_log | collect
        )
        dmr_tables = dmr_tables | mix(agg_out_m.dmr_table)

        ann_out_m = annotate_dmrs_m(agg_out_m.dmr_table, annotations_ready, flags.imprinted)
        annotated_tables = annotated_tables | mix(ann_out_m.annotated)

        rep_out_m = report_dmrs_m(ann_out_m.annotated, ann_out_m.annotation_log)
        reports = reports | mix(rep_out_m.report)

        if (flags.plot) {
            if (!flags.methylartist_only) {
                mb_m = plot_phased_dmr_modbamtools_m(ann_out_m.annotated, gencode_annotation.collect(), gene_list.collect(), phased_bam.collect())
                plots = plots | mix(mb_m.dmr_plots)
            }
            ma_m = plot_phased_dmr_methylartist_m(ann_out_m.annotated, gencode_annotation.collect(), gene_list.collect(), reference.collect(), phased_bam.collect())
            plots = plots | mix(ma_m.dmr_plots)
        }

        // ---- 5hmC ("h") ----
        prepped1_h = prep_dmr_bedmethyl_side1_h(paired.map { sq, alias, bed1, bed2 -> tuple(sq, bed1) }, "h")
        prepped2_h = prep_dmr_bedmethyl_side2_h(paired.map { sq, alias, bed1, bed2 -> tuple(sq, bed2) }, "h")
        chr_pairs_h = prepped1_h.dss_bed | join(prepped2_h.dss_bed, by: 0)

        dmr_out_h = dmr_calling_h(chr_pairs_h)
        agg_key_h = alias_ch | map { alias -> tuple(alias, comparison_label, "h") }
        agg_out_h = aggregate_dmrs_h(
            agg_key_h,
            dmr_out_h.chr_dmrs | map { chr, f -> f } | collect,
            dmr_out_h.status_log | collect
        )
        dmr_tables = dmr_tables | mix(agg_out_h.dmr_table)

        ann_out_h = annotate_dmrs_h(agg_out_h.dmr_table, annotations_ready, flags.imprinted)
        annotated_tables = annotated_tables | mix(ann_out_h.annotated)

        rep_out_h = report_dmrs_h(ann_out_h.annotated, ann_out_h.annotation_log)
        reports = reports | mix(rep_out_h.report)

        if (flags.plot) {
            if (!flags.methylartist_only) {
                mb_h = plot_phased_dmr_modbamtools_h(ann_out_h.annotated, gencode_annotation.collect(), gene_list.collect(), phased_bam.collect())
                plots = plots | mix(mb_h.dmr_plots)
            }
            ma_h = plot_phased_dmr_methylartist_h(ann_out_h.annotated, gencode_annotation.collect(), gene_list.collect(), reference.collect(), phased_bam.collect())
            plots = plots | mix(ma_h.dmr_plots)
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
        // consumed by prep_dmr_bedmethyl_side2_* below, same as
        // prep_dmr_bedmethyl_side1_* does for the sample side.
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

        // ---- 5mC ("m") ----
        sample_prepped_m = prep_dmr_bedmethyl_side1_m(sample_keyed.map { sq, alias, bed -> tuple(sq, bed) }, "m")
        ref_prepped_m = prep_dmr_bedmethyl_side2_m(reference_keyed, "m")
        chr_pairs_m = sample_prepped_m.dss_bed | join(ref_prepped_m.dss_bed, by: 0)

        dmr_out_m = dmr_calling_m(chr_pairs_m)
        agg_key_m = alias_ch | map { alias -> tuple(alias, comparison_label, "m") }
        agg_out_m = aggregate_dmrs_m(
            agg_key_m,
            dmr_out_m.chr_dmrs | map { chr, f -> f } | collect,
            dmr_out_m.status_log | collect
        )
        dmr_tables = dmr_tables | mix(agg_out_m.dmr_table)

        ann_out_m = annotate_dmrs_m(agg_out_m.dmr_table, annotations_ready, flags.imprinted)
        annotated_tables = annotated_tables | mix(ann_out_m.annotated)

        rep_out_m = report_dmrs_m(ann_out_m.annotated, ann_out_m.annotation_log)
        reports = reports | mix(rep_out_m.report)

        if (flags.plot && reference_available_for_plotting) {
            ref_bam_for_plot = ref_bam_ch | map { meta, xam, xai -> tuple(meta.alias, xam, xai) }
            if (!flags.methylartist_only) {
                mb_m = plot_dmr_modbamtools_m(ann_out_m.annotated, gencode_annotation.collect(), gene_list.collect(), sample_bam.collect(), ref_bam_for_plot.collect())
                plots = plots | mix(mb_m.dmr_plots)
            }
            ma_m = plot_dmr_methylartist_m(ann_out_m.annotated, gencode_annotation.collect(), gene_list.collect(), reference.collect(), sample_bam.collect(), ref_bam_for_plot.collect())
            plots = plots | mix(ma_m.dmr_plots)
        }

        // ---- 5hmC ("h") ----
        sample_prepped_h = prep_dmr_bedmethyl_side1_h(sample_keyed.map { sq, alias, bed -> tuple(sq, bed) }, "h")
        ref_prepped_h = prep_dmr_bedmethyl_side2_h(reference_keyed, "h")
        chr_pairs_h = sample_prepped_h.dss_bed | join(ref_prepped_h.dss_bed, by: 0)

        dmr_out_h = dmr_calling_h(chr_pairs_h)
        agg_key_h = alias_ch | map { alias -> tuple(alias, comparison_label, "h") }
        agg_out_h = aggregate_dmrs_h(
            agg_key_h,
            dmr_out_h.chr_dmrs | map { chr, f -> f } | collect,
            dmr_out_h.status_log | collect
        )
        dmr_tables = dmr_tables | mix(agg_out_h.dmr_table)

        ann_out_h = annotate_dmrs_h(agg_out_h.dmr_table, annotations_ready, flags.imprinted)
        annotated_tables = annotated_tables | mix(ann_out_h.annotated)

        rep_out_h = report_dmrs_h(ann_out_h.annotated, ann_out_h.annotation_log)
        reports = reports | mix(rep_out_h.report)

        if (flags.plot && reference_available_for_plotting) {
            ref_bam_for_plot_h = ref_bam_ch | map { meta, xam, xai -> tuple(meta.alias, xam, xai) }
            if (!flags.methylartist_only) {
                mb_h = plot_dmr_modbamtools_h(ann_out_h.annotated, gencode_annotation.collect(), gene_list.collect(), sample_bam.collect(), ref_bam_for_plot_h.collect())
                plots = plots | mix(mb_h.dmr_plots)
            }
            ma_h = plot_dmr_methylartist_h(ann_out_h.annotated, gencode_annotation.collect(), gene_list.collect(), reference.collect(), sample_bam.collect(), ref_bam_for_plot_h.collect())
            plots = plots | mix(ma_h.dmr_plots)
        }
    emit:
        dmr_table = dmr_tables
        annotated = annotated_tables
        report = reports
        plots = plots
}
