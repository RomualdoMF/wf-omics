// AnnotSV (https://github.com/lgmgeo/AnnotSV) annotation/ranking for structural
// variants and CNVs. Runs on the already VEP-annotated .wf_sv.vcf.gz (from
// wf-human-sv) and .wf_cnv.vcf.gz (from wf-human-cnv). Optional feature, off by
// default (params.annotsv).
//
// Like the VEP cache, the human annotation data (~5GB, Annotations_Human_*.tar.gz)
// is not baked into the container: it's downloaded once into
// params.annotsv_annotations_dir (idempotent, mkdir-lock, same pattern as
// modules/local/vep.nf) and bind-mounted at runtime (see base.config,
// withLabel:annotsv containerOptions). The container itself
// (docker/annotsv/Dockerfile) only has AnnotSV's code/tools.
//
// Extra AnnotSV flags can be set per SV-type via params.annotsv_sv_custom_args
// / params.annotsv_cnv_custom_args (raw token lists, e.g. ["-hpo", "HP:0001156"]).
// NOTE: AnnotSV's own `-vcf 1` (native VCF output, via the bundled
// variantconvert tool) doesn't work -- it produced a header-only VCF with no
// #CHROM line or records, in every variantconvertMode tried, even after
// fixing an unrelated pandas-version crash in the same code path. Rather than
// keep chasing that, variantconvert (and its Python stack) has been removed
// from the AnnotSV image entirely (see docker/annotsv/Dockerfile) -- passing
// "-vcf"/"1" in either custom_args list now just does nothing, since
// variantconvert isn't there to invoke. Both SV and CNV get an annotated VCF
// a different way instead: annotate_vcf_with_tsv below writes AnnotSV's own
// TSV columns back into the original VCF as INFO fields directly (via
// bin/annotate_vcf_from_annotsv.py), no variantconvert involved. The CNV path
// runs that on top of the fuller ISV/ClassifyCNV-merged TSV instead of
// AnnotSV's raw one -- see modules/local/cnv_ensemble.nf.

file(params.annotsv_annotations_dir).mkdirs()


process download_annotsv_annotations {
    // ensures params.annotsv_annotations_dir/Annotations_Human exists, downloading
    // it from the AnnotSV authors' server if missing. Same idempotent/lockable
    // approach as download_vep_cache in modules/local/vep.nf.
    label "annotsv"
    cpus 1
    memory 2.GB
    output:
        val(true), emit: ready
    script:
        def target = "${params.annotsv_annotations_dir}/Annotations_Human"
        def lock = "${params.annotsv_annotations_dir}/.download.lock"
        def tarball = "${params.annotsv_annotations_dir}/Annotations_Human.tar.gz"
        def url = "https://www.lbgi.fr/~geoffroy/Annotations/" +
            "Annotations_Human_${params.annotsv_human_annotations_version}.tar.gz"
        """
        if [ ! -d "${target}" ]; then
            until mkdir "${lock}" 2>/dev/null; do
                [ -d "${target}" ] && break
                sleep 5
            done
            if [ ! -d "${target}" ]; then
                curl -L -C - -o "${tarball}" "${url}"
                tar -xzf "${tarball}" -C "${params.annotsv_annotations_dir}"
                rm -f "${tarball}"
            fi
            rmdir "${lock}" 2>/dev/null || true
        fi
        """
}


process run_annotsv {
    // ranks/annotates structural variants; produces AnnotSV's native TSV output.
    // Skipped for genomes other than hg19/hg38 (AnnotSV also supports CHM13/mm9/
    // mm10, not wired up here since the rest of the pipeline only handles
    // hg19/hg38 -- see run_vep in modules/local/vep.nf for the same convention).
    // no publishDir here: the TSV flows into the `report`/artifacts channel in
    // workflows/wf-human-sv.nf and wf-human-cnv.nf, published once via
    // output_sv/output_cnv (same convention as the rest of those outputs, e.g.
    // the VEP-annotated VCF).
    label "annotsv"
    cpus 2
    memory 8.GB
    input:
        tuple val(xam_meta), path("input.vcf.gz"), path("input.vcf.gz.tbi")
        val(genome)
        val(output_label)  // e.g. "sv", "cnv" -- becomes "<alias>.wf_<label>.annotsv.tsv"
        val(annotations_ready)
    output:
        tuple val(xam_meta), path("${xam_meta.alias}.wf_${output_label}.annotsv.tsv"), emit: annotsv_tsv, optional: true
    script:
        def build = genome == 'hg19' ? 'GRCh37' : 'GRCh38'
        def out_name = "${xam_meta.alias}.wf_${output_label}.annotsv.tsv"
        // per-label extra AnnotSV flags (raw tokens, appended as-is) -- e.g.
        // params.annotsv_sv_custom_args = ["-vcf", "1"]
        def custom_args = (output_label == 'sv' ? params.annotsv_sv_custom_args : params.annotsv_cnv_custom_args).join(' ')
        """
        if [[ "${genome}" == "hg38" ]] || [[ "${genome}" == "hg19" ]]; then
            AnnotSV \
                -SVinputFile input.vcf.gz \
                -outputDir . \
                -outputFile ${out_name} \
                -genomeBuild ${build} \
                -annotationsDir ${params.annotsv_annotations_dir} \
                -overwrite 1 \
                ${custom_args}
        fi
        """
}


process annotate_vcf_with_tsv {
    // writes an AnnotSV-derived TSV's columns back into the VCF it was
    // derived from as INFO fields -- see bin/annotate_vcf_from_annotsv.py.
    // Called two ways:
    //   - workflows/wf-human-sv.nf: directly on AnnotSV's own raw
    //     <sample>.wf_sv.annotsv.tsv (no ensemble step for SV).
    //   - modules/local/cnv_ensemble.nf: on the ISV/ClassifyCNV-merged
    //     <sample>.wf_cnv.annotated.tsv instead (a superset of AnnotSV's own
    //     columns) -- annotate_vcf_from_annotsv.py auto-detects which shape
    //     of TSV it was given.
    // label "wf_common" (not "annotsv"): only needs pysam, which the
    // pipeline's existing ontresearch/wf-common image already has (used
    // elsewhere for reporting) -- reuses that instead of adding pysam to the
    // annotsv image.
    // no publishDir here: the output flows into the `report`/`output`
    // channel of whichever workflow calls this, published once via
    // output_sv/output_cnv like the rest of those outputs.
    label "wf_common"
    cpus 1
    memory 2.GB
    input:
        tuple val(xam_meta), path("input.annotated.tsv")
        tuple val(xam_meta2), path("input.vcf.gz"), path("input.vcf.gz.tbi")
        val(output_label)  // "sv" or "cnv" -- becomes "<alias>.wf_<label>.annotated.vcf.gz"
    output:
        tuple val(xam_meta), path("${xam_meta.alias}.wf_${output_label}.annotated.vcf.gz"), path("${xam_meta.alias}.wf_${output_label}.annotated.vcf.gz.tbi"), emit: annotated_vcf
    script:
        def vcf_name = "${xam_meta.alias}.wf_${output_label}.annotated.vcf.gz"
        """
        annotate_vcf_from_annotsv.py \
            --annotated-tsv input.annotated.tsv \
            --input-vcf input.vcf.gz \
            --output-vcf ${vcf_name}
        """
}


workflow annotsv {
    take:
        vcf_tuple     // tuple(xam_meta, vcf.gz, vcf.gz.tbi) -- .wf_sv.vcf.gz or .wf_cnv.vcf.gz
        genome        // "hg38" / "hg19" / other
        output_label  // "sv" or "cnv"
    main:
        ready = download_annotsv_annotations().ready
        result = run_annotsv(vcf_tuple, genome, output_label, ready).annotsv_tsv
    emit:
        annotsv_tsv = result
}
