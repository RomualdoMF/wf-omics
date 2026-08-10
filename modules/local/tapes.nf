// TAPES (https://github.com/RomualdoMF/tapes, a fork of a-xavier/tapes patched for
// dbNSFP5.x column compatibility -- see docker/tapes/README or romualdomf/tapes:fixed on
// Docker Hub) -- ACMG/AMP criteria (PVS1..BP7) + pathogenicity probability classification
// on top of the VEP CSQ annotation already produced by modules/local/vep.nf, including the
// dbNSFP5.3.1a in-silico predictor / ClinVar / gnomAD columns wired in via
// params.vep_plugin_args (see nextflow.config, the "dbNSFP" entry).
//
// Runs per-contig against the same per-contig annotated SNP VCFs already produced by
// annotate_snp_vcf in main.nf (the `annotations` channel, *before* concat_snp_vcfs merges
// them into one whole-sample VCF) -- ACMG classification is purely per-variant, so there's
// no need to wait for the merge, same reasoning as modules/local/tldr.nf running on
// per-contig BAMs instead of the final whole-genome one.

process run_tapes {
    label "tapes"
    cpus 1
    memory 4.GB
    input:
        tuple val(xam_meta), path(vcf), path(tbi)
        val(genome)
    output:
        tuple val(xam_meta), path("*.tapes.txt"), optional: true, emit: table
    script:
        // TAPES' bundled acmg_db (romualdomf/tapes:fixed, /opt/tapes/acmg_db) only has
        // hg19/hg38 reference data -- same guard as run_vep in modules/local/vep.nf.
        """
        if [[ "${genome}" != "hg38" ]] && [[ "${genome}" != "hg19" ]]; then
            echo "skipping TAPES ACMG classification for build '${genome}' (hg19/hg38 only)"
        else
            out_base=\$(basename ${vcf} .vcf.gz)
            # TAPES' own VCF reader opens the file directly (no gzip support), see
            # src/vep_process.py::vep_process_vcf.
            gunzip -c ${vcf} > input.vcf
            python3 /opt/tapes/tapes.py sort \
                -i input.vcf -o \${out_base}.tapes.txt --tab --acmg \
                --assembly ${genome} --acmg_db /opt/tapes/acmg_db
        fi
        """
}


process merge_tapes {
    // merge the per-contig TAPES tables into one, keeping only the first header --
    // same awk idiom as modules/local/tldr.nf::merge_tldr.
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        path(tables)
        val(xam_meta)
    output:
        path("*.wf_snp.tapes.txt")
    script:
        """
        awk 'FNR==1 && NR!=1 {next} {print}' ${tables} > ${xam_meta.alias}.wf_snp.tapes.txt
        """
}


workflow tapes_classify {
    take:
        annotated_snp_contigs  // tuple(xam_meta, vcf.gz, vcf.gz.tbi) per contig -- the
                                // `annotations` channel from main.nf's SNP annotation block
        genome                 // "hg38" / "hg19" / other
    main:
        per_contig_tables = run_tapes(annotated_snp_contigs, genome).table

        // same pattern as modules/local/tldr.nf::workflow tldr -- xam_meta is per-contig
        // here too (carries sq:/id:), so dedupe down to just meta.alias for the merged
        // output name.
        alias = annotated_snp_contigs
            | map { meta, vcf, tbi -> ['alias': meta.alias] }
            | unique

        merged_table = merge_tapes(
            per_contig_tables.map { meta, table -> table }.collect(),
            alias
        )
    emit:
        table = merged_table
}
