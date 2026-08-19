// tldr (https://github.com/adamewing/tldr) -- identifies and annotates
// transposable-element-mediated insertions in long-read alignments. Runs
// per-contig, against the same intermediate haplotagged BAMs straglr's
// call_str already consumes (workflows/wf-human-snp.nf, the `str_bams`
// emit) -- i.e. the phase where alignment has already happened but the
// per-contig BAMs haven't been merged back into one whole-genome BAM yet.
// See docker/tldr/Dockerfile for what's bundled in the image and why
// (minimap2/samtools/mafft/exonerate/tabix are all real runtime
// dependencies of tldr itself, not just of the upstream alignment).

process call_tldr {
    // 8GB (fixed, no retry scaling) was fine for smaller contigs, but larger ones
    // (e.g. chr4, ~220k read clusters) got OOM-killed (exit 137) with no retry --
    // no errorStrategy was set here at all, so Nextflow's default ('terminate')
    // aborted the whole run on the very first failure. Scale with task.attempt.
    label "tldr"
    cpus { params.tldr_procs }
    memory { MemoryScaling.forAttempt(MemoryScaling.SERIES_16, task.attempt, params.max_memory) }
    errorStrategy {task.exitStatus in [137, 140] ? 'retry' : 'finish'}
    maxRetries { MemoryScaling.retriesNeeded(MemoryScaling.SERIES_16, params.max_memory) }
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        tuple val(xam_meta.sq), path("*.table.txt"), optional: true
    script:
        def chr = xam_meta.sq
        def outbase = "${xam_meta.alias}.${chr}"
        // -b/--bams, -r/--ref, -o/--outbase and -c/--chroms are not exposed
        // as params (see nextflow.config): the BAM/ref are wired from the
        // pipeline's own channels, outbase is a deterministic per-contig
        // name so merge_tldr can concatenate the results afterwards, and
        // --chroms would be redundant -- this per-contig BAM already only
        // contains reads for ${chr}.
        def nonref_opt = params.tldr_nonref ? "-n ${params.tldr_nonref}" : ""
        def max_cluster_opt = params.tldr_max_cluster_size ? "--max_cluster_size ${params.tldr_max_cluster_size}" : ""
        def use_pickles_opt = params.tldr_use_pickles ? "--use_pickles ${params.tldr_use_pickles}" : ""
        """
        n_reads=\$(samtools view -c ${xam})
        if [[ "\$n_reads" -gt 0 ]]; then
            tldr \
                -b ${xam} \
                -e ${params.tldr_elts} \
                -r ${ref} \
                -p ${params.tldr_procs} \
                -m ${params.tldr_minreads} \
                --embed_minreads ${params.tldr_embed_minreads} \
                -o ${outbase} \
                --max_te_len ${params.tldr_max_te_len} \
                --min_te_len ${params.tldr_min_te_len} \
                --min_alt_frac ${params.tldr_min_alt_frac} \
                --min_alt_depth ${params.tldr_min_alt_depth} \
                --min_total_depth_frac ${params.tldr_min_total_depth_frac} \
                ${max_cluster_opt} \
                --wiggle ${params.tldr_wiggle} \
                --flanksize ${params.tldr_flanksize} \
                ${nonref_opt} \
                ${params.tldr_color_consensus ? "--color_consensus" : ""} \
                ${params.tldr_detail_output ? "--detail_output" : ""} \
                --extend_consensus ${params.tldr_extend_consensus} \
                ${params.tldr_trdcol ? "--trdcol" : ""} \
                ${params.tldr_keep_pickles ? "--keep_pickles" : ""} \
                ${use_pickles_opt}
        else
            echo "no reads on ${chr}, skipping"
        fi
        """
}


process merge_tldr {
    // merge the per-contig tables into one, keeping only the first header
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        path(tables)
        val(xam_meta)
    output:
        path("*.wf_tldr.table.txt")
    script:
        """
        awk 'FNR==1 && NR!=1 {next} {print}' ${tables} > ${xam_meta.alias}.wf_tldr.table.txt
        """
}


process getVersions {
    label "tldr"
    cpus 1
    // This process finishes in well under a second, which races Nextflow's own
    // background resource-trace collector (nxf_mem_watch, polling /proc/$pid inside
    // the container) against the task's FD teardown -- observed failing with
    // inconsistent exit codes (1, 141/SIGPIPE) despite versions.txt always coming out
    // correct, reproducible even standalone outside the pipeline. Not something
    // fixable from here (it's in Nextflow's own generated wrapper, not this script) --
    // the `sleep 1` gives the watcher time to poll at least once before the process
    // exits, and retry is a safety net for whatever timing slips through anyway.
    // Purely diagnostic (feeds the software-versions report, nothing downstream
    // depends on it) -- if it's still hitting the race after a couple retries,
    // give up gracefully rather than failing an entire multi-hour run over a
    // missing version string.
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2
    output:
        path "versions.txt"
    script:
        """
        tldr -v | sed 's/ /,/' >> versions.txt
        samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
        minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
        mafft --version 2>&1 | head -n 1 | sed 's/^/mafft,/' >> versions.txt
        exonerate --version | head -n 1 | sed 's/exonerate from exonerate version/exonerate,/' >> versions.txt
        sleep 1
        """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_tldr {
    // publish inputs to output directory
    label "tldr"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    script:
    """
    echo "Writing output files"
    """
}


workflow tldr {
    take:
        bam_channel   // per-contig (xam, xam_idx, xam_meta) -- e.g. clair_vcf.str_bams
        ref_channel
    main:
        // turn ref channel into value channel so it can be used more than once
        ref_as_value = ref_channel.collect()

        per_contig_tables = call_tldr(bam_channel, ref_as_value)

        software_versions = getVersions()

        // bam_channel.xam_meta is per contig (containing sq: and id:), so
        // just use meta.alias for the merged output name -- same pattern as
        // workflows/wf-human-str.nf's concat_str_vcfs/merge_tsv.
        alias = bam_channel
            | map { xam, xai, meta -> ['alias': meta.alias] }
            | unique

        merged_table = merge_tldr(
            per_contig_tables.map { sq, table -> table }.collect(),
            alias
        )

    emit:
        output = merged_table.concat(software_versions).flatten()
        table = merged_table
}
