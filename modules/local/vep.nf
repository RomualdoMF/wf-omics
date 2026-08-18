// VEP (Ensembl Variant Effect Predictor) annotation, replacing the original
// SnpEff (functional consequences) + SnpSift (ClinVar annotate/filter) steps.
//
// Design notes:
// - The `ensemblorg/ensembl-vep` image does not ship bcftools, so contig
//   splitting (see prepare_annotation_input) is done in the default pipeline
//   container (which already carries bcftools/tabix, see concat_vcfs in
//   modules/local/common.nf), and only the actual `vep` call runs inside the
//   vep_annotation-labelled container.
// - Clinical significance is VEP's own standard CLIN_SIG annotation
//   (--check_existing, see run_vep), sourced entirely from the already
//   downloaded cache -- no separate ClinVar VCF, no extra download. There is
//   no separate ClinVar-only output VCF anymore -- the full annotated VCF is
//   what gets published and fed to the report, which reads CLIN_SIG/SIFT
//   straight out of CSQ (see bin/workflow_glue/report_snp.py).
// - `annotate_vcf` is a workflow with the exact same name, inputs and output
//   as the process it replaces in modules/local/common.nf, so main.nf /
//   wf-human-cnv.nf / wf-human-sv.nf don't need to change how they call it.
// - VEP cache/plugins are large, shared, external directories -- not something
//   Nextflow's normal per-task file staging is meant for. They are bind-mounted
//   directly into the container at the same absolute path on host and guest
//   (see base.config, withLabel:vep_annotation containerOptions), and
//   downloaded in place (idempotently, with a simple mkdir-lock to stay safe
//   if snp/cnv/sv annotation all request them at once) rather than being
//   tracked as regular Nextflow outputs.

// make sure the target directories exist before anything tries to bind-mount them
file(params.vep_cache_dir).mkdirs()
file(params.vep_plugins_dir).mkdirs()


process download_vep_cache {
    // ensures params.vep_cache_dir contains the requested species/assembly/version,
    // downloading it via VEP's own installer if it doesn't. Safe to call more
    // than once concurrently (snp/cnv/sv annotation each call this): only the
    // first caller to grab the mkdir-lock actually downloads, the rest wait.
    label "vep_annotation"
    cpus 2
    memory 4.GB
    input:
        val(assembly)
    output:
        val(true), emit: ready
    script:
        def target = "${params.vep_cache_dir}/${params.vep_species}/${params.vep_cache_version}_${assembly}"
        def lock = "${params.vep_cache_dir}/.download_${params.vep_species}_${params.vep_cache_version}_${assembly}.lock"
        """
        if [ ! -d "${target}" ]; then
            until mkdir "${lock}" 2>/dev/null; do
                [ -d "${target}" ] && break
                sleep 5
            done
            if [ ! -d "${target}" ]; then
                perl /opt/vep/src/ensembl-vep/INSTALL.pl --AUTO c \
                    --SPECIES ${params.vep_species} --ASSEMBLY ${assembly} \
                    --CACHE_VERSION ${params.vep_cache_version} \
                    --CACHEDIR ${params.vep_cache_dir} \
                    --NO_HTSLIB --NO_UPDATE --QUIET
            fi
            rmdir "${lock}" 2>/dev/null || true
        fi
        """
}


process download_vep_plugins {
    // ensures params.vep_plugins_dir is populated, same idempotent/lockable
    // approach as download_vep_cache above.
    label "vep_annotation"
    cpus 1
    memory 2.GB
    output:
        val(true), emit: ready
    script:
        // the lock has to live *inside* vep_plugins_dir: that's the only path
        // guaranteed to be bind-mounted into the container (see base.config),
        // a sibling/parent path wouldn't be visible in there. Hidden entries
        // (dotfiles, including the lock itself) are excluded from the
        // "already populated" check below so the lock doesn't fool it.
        def lock = "${params.vep_plugins_dir}/.download.lock"
        """
        has_content() { [ -n "\$(ls -A ${params.vep_plugins_dir} 2>/dev/null | grep -v '^\\.')" ]; }
        if ! has_content; then
            until mkdir "${lock}" 2>/dev/null; do
                has_content && break
                sleep 5
            done
            if ! has_content; then
                perl /opt/vep/src/ensembl-vep/INSTALL.pl --AUTO p \
                    --PLUGINSDIR ${params.vep_plugins_dir} \
                    --PLUGINS ${params.vep_install_plugins} \
                    --NO_HTSLIB --NO_UPDATE --QUIET
            fi
            rmdir "${lock}" 2>/dev/null || true
        fi
        """
}


process prepare_annotation_input {
    // split the input VCF down to a single contig (SNP path, annotated
    // per-contig for parallelism) or pass the whole file through unchanged
    // (SV/CNV path, marked with contig == '*').
    cpus 1
    memory 2.GB
    input:
        tuple val(xam_meta), path("input.vcf.gz"), path("input.vcf.gz.tbi"), val(contig)
        val(output_label)
    output:
        tuple val(xam_meta), path("prepared.vcf.gz"), env(FULL_OUTPUT_LABEL), emit: prepared
    script:
        """
        if [ "${contig}" == '*' ]; then
            cp input.vcf.gz prepared.vcf.gz
            FULL_OUTPUT_LABEL="${output_label}"
        else
            bcftools view -r ${contig} input.vcf.gz | bgzip > prepared.vcf.gz
            FULL_OUTPUT_LABEL="${output_label}.${contig}"
        fi
        """
}


process run_vep {
    // runs functional consequence + ClinVar annotation with VEP, replacing
    // `snpEff ann` + `SnpSift annotate`. Clinical significance and SIFT come
    // from VEP's own standard annotation (--check_existing, --sift), sourced
    // entirely from the already-downloaded cache -- no separate ClinVar VCF or
    // extra download involved. See bin/workflow_glue/report_snp.py for how
    // the report reads the resulting CLIN_SIG/SIFT fields back out of CSQ.
    // Falls back to a plain passthrough for genomes other than hg19/hg38,
    // same behaviour as the SnpEff step it replaces.
    label "vep_annotation"
    cpus 4
    memory 8.GB
    input:
        tuple val(xam_meta), path("prepared.vcf.gz"), val(output_label)
        val(genome)
        val(cache_ready)
        val(plugins_ready)
    output:
        tuple val(xam_meta), path("${xam_meta.alias}.wf_${output_label}.vcf.gz"), path("${xam_meta.alias}.wf_${output_label}.vcf.gz.tbi"), emit: annot_vcf
    script:
        def assembly = genome == 'hg19' ? 'GRCh37' : 'GRCh38'
        def custom_flags = params.vep_custom_args.collect { "--custom ${it}" }.join(' ')
        def plugin_flags = params.vep_plugin_args.collect { "--plugin ${it}" }.join(' ')
        def out_name = "${xam_meta.alias}.wf_${output_label}.vcf.gz"
        """
        # --output_file ${out_name} (not STDOUT piped via shell `>`): without --fasta, some
        # plugins (seen with NMD/UTRAnnotator/CADD-on-indels in particular) make VEP print a
        # "no FASTA file specified" warning straight to its own STDOUT before --warning_file
        # even takes effect -- when that STDOUT is a shell pipe into the compressed VCF file
        # (the old `--output_file STDOUT ... > ${out_name}` form), the warning text corrupts
        # the bgzip stream and tabix then fails with "tbx_index_build failed". Writing to a
        # real named output file sidesteps this entirely: VEP's own bgzip writer never shares
        # a stream with terminal chatter, whatever prints it and regardless of --warning_file.
        if [[ "${genome}" != "hg38" ]] && [[ "${genome}" != "hg19" ]]; then
            cp prepared.vcf.gz ${out_name}
        else
            vep --input_file prepared.vcf.gz --output_file ${out_name} --vcf --compress_output bgzip \
                --cache --offline \
                --dir_cache ${params.vep_cache_dir} --dir_plugins ${params.vep_plugins_dir} \
                --species ${params.vep_species} --assembly ${assembly} \
                --cache_version ${params.vep_cache_version} \
                --fork ${task.cpus} --no_stats --force_overwrite \
                --sift b --check_existing \
                --warning_file STDERR \
                ${custom_flags} ${plugin_flags} ${params.vep_extra_args} \
                > vep.log 2>&1
        fi
        tabix -p vcf ${out_name}
        """
}


workflow annotate_vcf {
    take:
        vcf_contig_tuple  // tuple(xam_meta, vcf.gz, vcf.gz.tbi, contig)
        genome            // "hg38" / "hg19" / other
        output_label      // e.g. "snp", "sv", "cnv"
    main:
        assembly = genome.map { it == 'hg19' ? 'GRCh37' : 'GRCh38' }
        cache_ready = download_vep_cache(assembly).ready
        plugins_ready = download_vep_plugins().ready
        prepared = prepare_annotation_input(vcf_contig_tuple, output_label).prepared
        final_out = run_vep(prepared, genome, cache_ready, plugins_ready).annot_vcf
    emit:
        annot_vcf = final_out
}
