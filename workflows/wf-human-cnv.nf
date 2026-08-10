include {
    callCNV;
    getVersions;
    add_snp_tools_to_versions;
    bgzip_and_index_vcf;
    makeReport
} from "../modules/local/wf-human-cnv.nf"

include {
    mosdepth
} from "../modules/local/common.nf"

include {
    annotate_vcf
} from "../modules/local/vep.nf"

include {
    annotsv
} from "../modules/local/annotsv.nf"

include {
    cnv_ensemble
} from "../modules/local/cnv_ensemble.nf"

workflow cnv {
    take:
        bam
        ref
        clair3_vcf
        bed
        genome_build
        workflow_params
    main:
        // get mosdepth results for window size 1000
        mosdepth(bam, bed, ref, "1000", false, "bed")
        mosdepth_stats = mosdepth.out.mosdepth_tuple.map{ meta, bed, dist, threshold -> [bed, dist, threshold]}
        mosdepth_summary = mosdepth.out.summary
        if (params.depth_intervals) {
            mosdepth_perbase = mosdepth.out.perbase
        } else {
            mosdepth_perbase = Channel.from("$projectDir/data/OPTIONAL_FILE")
        }

        mosdepth_all = mosdepth_stats.concat(mosdepth_summary).concat(mosdepth_perbase).collect()
        
        cnvs = callCNV(clair3_vcf, mosdepth_all, ref, genome_build)
        spectre_vcf = cnvs.spectre_vcf
        spectre_vcf_bgzipped = bgzip_and_index_vcf(spectre_vcf)
        spectre_bed = cnvs.spectre_bed
        spectre_karyotype = cnvs.spectre_karyotype

        // check if SnpEff annotations have been requested
        if (!params.annotation) {
            spectre_final_vcf = spectre_vcf_bgzipped
            annotsv_tsv = Channel.empty()
            cnv_ensemble_tsv = Channel.empty()
            cnv_ensemble_vcf = Channel.empty()
        }
        else {
            // append '*' to indicate that annotation should be performed on all chr at once
            vcf_for_annotation = spectre_vcf_bgzipped.map{ it << '*' }
            spectre_final_vcf = annotate_vcf(vcf_for_annotation, genome_build, "cnv").annot_vcf

            // optionally rank/annotate the CNVs further with AnnotSV, and (on
            // top of that) run the ISV/ClassifyCNV ensemble on AnnotSV's output
            if (params.annotsv) {
                annotsv_result = annotsv(spectre_final_vcf, genome_build, "cnv").annotsv_tsv
                annotsv_tsv = annotsv_result.map{ meta, tsv -> tsv }

                if (params.cnv_ensemble) {
                    ensemble_result = cnv_ensemble(annotsv_result, spectre_final_vcf, genome_build)
                    cnv_ensemble_tsv = ensemble_result.annotated_tsv.map{ meta, tsv -> tsv }
                    cnv_ensemble_vcf = ensemble_result.annotated_vcf.map{ meta, vcf, tbi -> [vcf, tbi] }
                } else {
                    cnv_ensemble_tsv = Channel.empty()
                    cnv_ensemble_vcf = Channel.empty()
                }
            } else {
                annotsv_tsv = Channel.empty()
                cnv_ensemble_tsv = Channel.empty()
                cnv_ensemble_vcf = Channel.empty()
            }
        }

        software_versions_tmp = getVersions()
        software_versions = add_snp_tools_to_versions(software_versions_tmp)
        if (params.output_report){
            report = makeReport(software_versions.collect(), workflow_params, spectre_bed, spectre_karyotype, genome_build)
        } else {
            report = Channel.empty()
        }

    emit:
        output = spectre_final_vcf.map{ meta, vcf, tbi -> [vcf, tbi]}.concat(report, annotsv_tsv, cnv_ensemble_tsv, cnv_ensemble_vcf)
        cnv_vcf = spectre_final_vcf
}