// CNV classification ensemble: runs ISV + ClassifyCNV on top of AnnotSV's
// .wf_cnv.annotsv.tsv, merges everything into .wf_cnv.annotated.tsv, and
// writes the same columns back into the original .wf_cnv.vcf.gz as INFO
// fields, producing .wf_cnv.annotated.vcf.gz.
//
// Split into two processes/containers on purpose:
//   - run_cnv_ensemble (label "cnv_ensemble", bin/cnv_ensemble_classify.py):
//     needs ISV + ClassifyCNV, hence the custom romualdomf/cnv-ensemble
//     image (docker/cnv_ensemble/Dockerfile) -- no download needed at
//     runtime, unlike VEP/AnnotSV.
//   - annotate_vcf_with_tsv (label "wf_common",
//     bin/annotate_vcf_from_annotsv.py): only needs pysam, which the
//     pipeline's existing ontresearch/wf-common image already has (used
//     elsewhere for reporting) -- so this reuses that image instead of
//     duplicating pysam into cnv-ensemble. It's the same process the SV path
//     uses (see modules/local/annotsv.nf), just fed the ensemble-merged TSV
//     here instead of AnnotSV's raw one.
//
// Optional feature, off by default (params.cnv_ensemble). Only meaningful
// when params.annotsv is also on, since it consumes AnnotSV's own output
// directly -- see workflows/wf-human-cnv.nf for how the two are chained.

include {
    annotate_vcf_with_tsv
} from './annotsv.nf'

process run_cnv_ensemble {
    label "cnv_ensemble"
    cpus 2
    memory 4.GB
    input:
        tuple val(xam_meta), path("input.annotsv.tsv")
        val(genome)
    output:
        tuple val(xam_meta), path("${xam_meta.alias}.wf_cnv.annotated.tsv"), emit: annotated_tsv, optional: true
    script:
        def tsv_name = "${xam_meta.alias}.wf_cnv.annotated.tsv"
        """
        if [[ "${genome}" == "hg38" ]] || [[ "${genome}" == "hg19" ]]; then
            cnv_ensemble_classify.py \
                --annotsv-tsv input.annotsv.tsv \
                --classifycnv-dir \${CLASSIFYCNV_DIR} \
                --genome-build ${genome} \
                --output-tsv ${tsv_name}
        fi
        """
}


workflow cnv_ensemble {
    take:
        annotsv_tuple  // tuple(xam_meta, <sample>.wf_cnv.annotsv.tsv)
        cnv_vcf_tuple  // tuple(xam_meta, <sample>.wf_cnv.vcf.gz, .tbi) -- the VCF AnnotSV was run on
        genome         // "hg38" / "hg19" / other
    main:
        tsv_result = run_cnv_ensemble(annotsv_tuple, genome).annotated_tsv
        vcf_result = annotate_vcf_with_tsv(tsv_result, cnv_vcf_tuple, "cnv").annotated_vcf
    emit:
        annotated_tsv = tsv_result
        annotated_vcf = vcf_result
}
