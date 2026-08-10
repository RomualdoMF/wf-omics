// DMR (Differentially Methylated Region) calling, adapted from
// https://github.com/NyagaM/ont-methylDMR-kit (commit dcbc3da31db19f477194cdd33ebc3001554f888e).
//
// Two opt-in comparisons (see workflows/dmr.nf):
//   - dmr_haplotype_compare: haplotype 1 vs haplotype 2 of the same sample (--phased only).
//   - dmr_sample_compare: the sample vs an external reference (BAM/CRAM or bedmethyl).
//
// Unlike the upstream tool -- which always re-splits its input bedmethyl by
// chromosome before running DSS, even when starting from a pre-built bedmethyl --
// this integration hooks in *before* workflows/methyl.nf's concat_bedmethyl merges
// modkit's already-per-chromosome output into one file. That per-chromosome
// parallelism is reused directly for prep_dmr_bedmethyl/dmr_calling; only an
// external, already-merged reference bedmethyl (dmr_sample_compare with a
// .bed/.bedmethyl reference) needs the upstream-style split_bedmethyl_by_chr.
//
// DSS parameters (delta/p.threshold/minlen/minCG/dis.merge) are hardcoded, same as
// upstream -- not exposed as params (see dmr_calling below).

// Annotation data (gencode v46 promoter/exon/intron BED + imprinted gene list) is not
// vendored into this repo (~73MB tarball, committed directly in the upstream repo with
// no separate stable download URL) -- it's fetched once from a pinned upstream commit
// into params.dmr_annotations_dir, same idempotent/mkdir-lock pattern as
// download_vep_cache/download_annotsv_annotations.
file(params.dmr_annotations_dir).mkdirs()


process list_reference_chromosomes {
    // For dmr_sample_compare when the reference is a BAM/CRAM: intersects the
    // chromosomes actually present in that BAM's header with chromosome_codes,
    // so modkit only gets asked to pileup regions that exist (a BAM given by the
    // user won't necessarily have gone through the same bam_flagstats-based
    // contig filtering the rest of the pipeline uses for the primary sample).
    // No label -- uses the default container (has samtools), same convention as
    // modules/local/common.nf's getGenome.
    cpus 1
    memory 2.GB
    input:
        tuple path(xam), path(xam_idx)
        val chromosome_codes
    output:
        path("chroms.txt"), emit: chroms
    script:
        def allow = chromosome_codes.join('\n')
        """
        cat <<'ALLOWED' > allowed.txt
${allow}
ALLOWED
        samtools view -H ${xam} --no-PG | grep '^@SQ' | sed -nE 's,.*SN:([^[:space:]]*).*,\\1,p' | sort -u > present.txt
        comm -12 <(sort allowed.txt) present.txt > chroms.txt
        """
}


process download_dmr_annotations {
    // ensures params.dmr_annotations_dir/annotations exists, downloading it from a
    // pinned ont-methylDMR-kit commit if missing. Only invoked when annotate/report
    // is actually requested (see dmr_stage_flags in workflows/dmr.nf).
    label "dmr_analysis"
    cpus 1
    memory 2.GB
    output:
        val(true), emit: ready
    script:
        def target = "${params.dmr_annotations_dir}/annotations"
        def lock = "${params.dmr_annotations_dir}/.download.lock"
        def tarball = "${params.dmr_annotations_dir}/annotations.tar.gz"
        // pinned to the commit analysed/integrated in this session
        def url = "https://raw.githubusercontent.com/NyagaM/ont-methylDMR-kit/" +
            "dcbc3da31db19f477194cdd33ebc3001554f888e/annotations.tar.gz"
        """
        if [ ! -d "${target}" ]; then
            until mkdir "${lock}" 2>/dev/null; do
                [ -d "${target}" ] && break
                sleep 5
            done
            if [ ! -d "${target}" ]; then
                curl -L -C - -o "${tarball}" "${url}"
                mkdir -p "${target}"
                tar -xzf "${tarball}" -C "${target}"
                rm -f "${tarball}"
            fi
            rmdir "${lock}" 2>/dev/null || true
        fi
        """
}


process prep_dmr_bedmethyl {
    // Filters an already per-chromosome, bgzip-compressed bedmethyl (as produced by
    // workflows/methyl.nf's modkit/modkit_phase, pre-concat_bedmethyl) for one
    // modified-base code, and reformats it into the 4-column table DSS expects
    // (chr, pos, N=valid coverage, X=modified-read count). No chromosome splitting
    // needed here -- the input is already scoped to one chromosome.
    label "dmr_analysis"
    cpus 1
    memory 2.GB
    input:
        tuple val(chr), path(bedmethyl_gz)
        val(mod_char)  // "m" (5mC) or "h" (5hmC)
    output:
        tuple val(chr), path("${chr}__${mod_char}.dss.tsv"), emit: dss_bed
    script:
        """
        echo -e "chr\\tpos\\tN\\tX" > "${chr}__${mod_char}.dss.tsv"
        zcat -f ${bedmethyl_gz} | \\
            awk -F'\\t' -v mc="${mod_char}" 'BEGIN{OFS="\\t"} \$4 == mc && \$5 >= 5 {print \$1, \$3, \$5, \$12}' \\
            >> "${chr}__${mod_char}.dss.tsv"
        """
}


process split_bedmethyl_by_chr {
    // Only used for dmr_sample_compare when dmr_sample_compare_reference is an
    // already-merged bedmethyl (not a BAM/CRAM we pileup ourselves). Splits the
    // raw bedmethyl rows by chromosome only (no mod-code/coverage filtering --
    // that's prep_dmr_bedmethyl's job, applied uniformly afterwards for both the
    // BAM-derived and bedmethyl-derived reference paths), restricted to the same
    // chromosome_codes list the rest of the pipeline uses so the per-chromosome
    // join with the sample side lines up. Runs once (not once per mod type),
    // since 5mC/5hmC rows share the same split.
    label "dmr_analysis"
    cpus 2
    memory 4.GB
    input:
        path bedmethyl
        val chromosome_codes  // ArrayList, e.g. ["chr1","1",...,"chrX","X"]
    output:
        path "chr_split__*.bed", emit: chr_beds
    script:
        def allow = chromosome_codes.join(',')
        """
        zcat -f ${bedmethyl} | \\
            awk -F'\\t' -v allow="${allow}" '
                BEGIN {
                    OFS = "\\t"
                    n = split(allow, arr, ",")
                    for (i = 1; i <= n; i++) allowed[arr[i]] = 1
                }
                (\$1 in allowed) {
                    print > ("chr_split__" \$1 ".bed")
                }'
        """
}


process dmr_calling {
    // DSS/bsseq DMR calling for one chromosome, one pair of prepped bedmethyl
    // tables. Parameters match upstream ont-methylDMR-kit exactly (not exposed as
    // params -- same decision upstream made, see its README).
    label "dmr_analysis"
    cpus 4
    memory 4.GB
    time 1.h
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    maxRetries 1
    maxForks 10
    input:
        tuple val(chr), path(bed1), path(bed2)
    output:
        tuple val(chr), path("dmrs_${chr}.bed"), emit: chr_dmrs
        path("dmr_status_${chr}.log"), emit: status_log
        path("debug_${chr}"), emit: debug_output, optional: true
    script:
        """
        mkdir -p debug_${chr}

        n_sites_1=\$(tail -n +2 ${bed1} | wc -l)
        n_sites_2=\$(tail -n +2 ${bed2} | wc -l)
        echo "Chromosome ${chr}: Sample1 has \${n_sites_1} sites, Sample2 has \${n_sites_2} sites" > dmr_status_${chr}.log

        if [ \${n_sites_1} -lt 10 ] || [ \${n_sites_2} -lt 10 ]; then
            echo "Insufficient sites for DMR analysis on chromosome ${chr}" >> dmr_status_${chr}.log
            echo -e "chr\\tstart\\tend\\tlength\\tnCG\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_${chr}.bed
            exit 0
        fi

        cat <<EOF > dmr_analysis_${chr}.R
        library(DSS)
        require(bsseq)

        dat1 = read.table("${bed1}", header=TRUE)
        dat2 = read.table("${bed2}", header=TRUE)

        BSobj = makeBSseqData(list(dat1, dat2), c("C1","C2"))
        dmlTest = DMLtest(BSobj, group1="C1", group2="C2", smoothing=TRUE, ncores=${task.cpus})

        dmrs = callDMR(dmlTest,
                       delta=0.10,
                       p.threshold=0.01,
                       minlen=100,
                       minCG=10,
                       dis.merge=100,
                       pct.sig=0.5)
        write.table(dmrs, file="dmrs_${chr}.tsv", row.names=FALSE, quote=FALSE, sep="\\t")
EOF

        Rscript dmr_analysis_${chr}.R > r_output_${chr}.log 2>&1

        if [[ -f dmrs_${chr}.tsv && -s dmrs_${chr}.tsv ]]; then
            echo -e "chr\\tstart\\tend\\tlength\\tnCG\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_${chr}.bed
            tail -n +2 dmrs_${chr}.tsv >> dmrs_${chr}.bed
            n_dmrs=\$(tail -n +2 dmrs_${chr}.bed | wc -l)
            echo "Found \${n_dmrs} DMRs for chromosome ${chr}" >> dmr_status_${chr}.log
        else
            echo "No significant DMRs found for chromosome ${chr}" >> dmr_status_${chr}.log
            echo -e "chr\\tstart\\tend\\tlength\\tnCG\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > dmrs_${chr}.bed
        fi

        mv dmr_analysis_${chr}.R r_output_${chr}.log debug_${chr}/
        rm -f dmrs_${chr}.tsv
        """
}


process aggregate_dmrs {
    // Concatenates the per-chromosome DMR beds/logs from dmr_calling into one
    // final table + status log, sorted by chr/pos.
    label "dmr_analysis"
    cpus 1
    memory 2.GB
    input:
        tuple val(alias), val(comparison_label), val(mod_char)  // e.g. ("SAMPLE01", "haplotype_compare", "m")
        path(chr_dmr_files)
        path(status_logs)
    output:
        tuple val(alias), val(comparison_label), val(mod_char), path("${alias}.wf_mods.${comparison_label}.${mod_char == 'h' ? '5hmC' : '5mC'}.dmr_table.tsv"), emit: dmr_table
        path("${alias}.wf_mods.${comparison_label}.${mod_char == 'h' ? '5hmC' : '5mC'}.dmr_status.log"), emit: status_log
    script:
        def out_label = (mod_char == 'h' ? '5hmC' : '5mC')
        def table_name = "${alias}.wf_mods.${comparison_label}.${out_label}.dmr_table.tsv"
        def log_name = "${alias}.wf_mods.${comparison_label}.${out_label}.dmr_status.log"
        """
        echo -e "chr\\tstart\\tend\\tlength\\tnCG\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat" > "${table_name}"
        for dmr_file in dmrs_*.bed; do
            [ -s "\${dmr_file}" ] || continue
            tail -n +2 "\${dmr_file}" >> "${table_name}"
        done

        if [ \$(tail -n +2 "${table_name}" | wc -l) -gt 0 ]; then
            head -1 "${table_name}" > sorted.tmp
            tail -n +2 "${table_name}" | sort -k1,1V -k2,2n >> sorted.tmp
            mv sorted.tmp "${table_name}"
        fi

        echo "DMR Analysis Summary (${comparison_label}, ${out_label})" > "${log_name}"
        echo "===================" >> "${log_name}"
        cat dmr_status_*.log >> "${log_name}"
        total_dmrs=\$(tail -n +2 "${table_name}" | wc -l)
        echo "" >> "${log_name}"
        echo "Total DMRs across all chromosomes: \${total_dmrs}" >> "${log_name}"
        """
}


process annotate_dmrs {
    // bedtools intersect against the gencode promoter/exon/intron BED, same as
    // upstream annotate_dmrs.nf, with an optional imprinted-gene filter.
    label "dmr_analysis"
    cpus 1
    memory 2.GB
    input:
        tuple val(alias), val(comparison_label), val(mod_char), path(dmr_table)
        val(annotations_ready)
        val(imprinted_only)
    output:
        tuple val(alias), val(comparison_label), val(mod_char), path("${alias}.wf_mods.${comparison_label}.${mod_char == 'h' ? '5hmC' : '5mC'}.dmr_annotated.tsv"), emit: annotated
        path("${alias}.wf_mods.${comparison_label}.${mod_char == 'h' ? '5hmC' : '5mC'}.dmr_annotation.log"), emit: annotation_log
    script:
        def out_label = (mod_char == 'h' ? '5hmC' : '5mC')
        def annotation_bed = "${params.dmr_annotations_dir}/annotations/gencode.v46.annotation.exon-promoters-introns.sorted.bed"
        def imprinted_tsv = "${params.dmr_annotations_dir}/annotations/imprinted_genes.tsv"
        def out_name = "${alias}.wf_mods.${comparison_label}.${out_label}.dmr_annotated.tsv"
        def log_name = "${alias}.wf_mods.${comparison_label}.${out_label}.dmr_annotation.log"
        """
        echo -e 'chr\\tstart\\tend\\tlength\\tnCG\\tmeanMethy1\\tmeanMethy2\\tdiff.Methy\\tareaStat\\tannotation_chr\\tannotation_start\\tannotation_end\\tstrand\\tannotation\\tbiotype\\tgene' > annotation_header.txt

        tail -n +2 ${dmr_table} > dmrs_clean.bed
        cp -L ${annotation_bed} annotation_clean.bed

        bedtools intersect -a dmrs_clean.bed -b annotation_clean.bed -wa -wb > dmrs_annotated.tmp
        (cat annotation_header.txt; sort dmrs_annotated.tmp | uniq) > "${out_name}"

        touch "${log_name}"
        total_dmrs=\$(wc -l < dmrs_clean.bed)
        annotated_dmrs=\$(sort dmrs_annotated.tmp | uniq | wc -l)
        echo "Total DMRs: \${total_dmrs}" >> "${log_name}"
        echo "DMR Annotations: \${annotated_dmrs}" >> "${log_name}"

        if [ "${imprinted_only}" = "true" ] && [ -f "${imprinted_tsv}" ]; then
            tail -n +2 "${imprinted_tsv}" | cut -f1 | sort | uniq > imprinted_genes_list.txt
            awk -F'\\t' 'NR==FNR{genes[\$1]; next} FNR==1 || \$16 in genes' \\
                imprinted_genes_list.txt "${out_name}" > "${out_name}.imprinted.tmp"
            mv "${out_name}.imprinted.tmp" "${out_name}"
            imprinted_dmrs=\$(tail -n +2 "${out_name}" | cut -f1-3 | sort | uniq | wc -l)
            echo "DMRs overlapping imprinted genes: \${imprinted_dmrs}" >> "${log_name}"
        fi

        rm -f dmrs_annotated.tmp annotation_header.txt dmrs_clean.bed annotation_clean.bed imprinted_genes_list.txt
        """
}


process prepare_gene_list {
    // Resolves the gene list used to filter which DMRs get plotted (see
    // dmr_stage_flags in workflows/dmr.nf): an explicit --gene_list path, the
    // imprinted-gene list (--imprinted, derived from the same downloaded
    // annotations used by annotate_dmrs), or an empty file (plot everything --
    // the plotting processes treat a zero-size gene list as "no filter").
    label "dmr_analysis"
    cpus 1
    memory 1.GB
    input:
        val(annotations_ready)
        val(imprinted_only)
        val(user_gene_list)  // String path, or "" if not provided
    output:
        path("gene_list.txt"), emit: gene_list
    script:
        def imprinted_tsv = "${params.dmr_annotations_dir}/annotations/imprinted_genes.tsv"
        """
        if [ -n "${user_gene_list}" ]; then
            cp "${user_gene_list}" gene_list.txt
        elif [ "${imprinted_only}" = "true" ] && [ -f "${imprinted_tsv}" ]; then
            tail -n +2 "${imprinted_tsv}" | cut -f1 | sort -u > gene_list.txt
        else
            touch gene_list.txt
        fi
        """
}


process report_dmrs {
    // pandas/matplotlib HTML summary report, adapted from upstream
    // generate_dmr_report.nf. Needs matplotlib (not available in the pipeline's
    // existing wf_common image), so it uses ont-methylDMR-kit's own report image
    // rather than adding matplotlib to wf_common.
    label "dmr_report"
    cpus 1
    memory 2.GB
    input:
        tuple val(alias), val(comparison_label), val(mod_char), path(annotated_tsv)
        path(annotation_log)
    output:
        path("${alias}.wf_mods.${comparison_label}.${mod_char == 'h' ? '5hmC' : '5mC'}.dmr_report.html"), emit: report
    script:
        def out_label = (mod_char == 'h' ? '5hmC' : '5mC')
        def out_name = "${alias}.wf_mods.${comparison_label}.${out_label}.dmr_report.html"
        """
#!/usr/bin/env python
import os
os.environ['MPLCONFIGDIR'] = os.getcwd()
import pandas as pd
import matplotlib.pyplot as plt
import base64
import io

summary_text = ""
if os.path.exists("${annotation_log}"):
    with open("${annotation_log}", 'r') as f:
        summary_text = f.read().replace('\\n', '<br>')

try:
    df = pd.read_csv("${annotated_tsv}", sep='\\t')
    plots = []

    if 'annotation' in df.columns:
        df['annotation'] = df['annotation'].replace(['promoter_plus', 'promoter_minus'], 'promoters')
        fig, ax = plt.subplots(figsize=(6, 4))
        counts = df['annotation'].value_counts().head(10).sort_values(ascending=True)
        color_palette = ['#73a3d4', '#ff9d5c', '#66c37f', '#ff7e79', '#b194e6', '#c8a17e', '#f78ce0', '#b1b1b1', '#f2f17b', '#92e4e1']
        plot_colors = color_palette[:len(counts)]
        bars = ax.barh(counts.index, counts.values, color=plot_colors)
        ax.set_xlabel('Count')
        ax.set_title('DMRs by Genomic Annotation')
        for bar, val in zip(bars, counts.values):
            ax.text(bar.get_width()+0.5, bar.get_y()+bar.get_height()/2, str(val), va='center')
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        buf.seek(0)
        plots.append(base64.b64encode(buf.read()).decode())
        plt.close()

    if 'biotype' in df.columns and not df['biotype'].dropna().empty:
        fig, ax = plt.subplots(figsize=(8, 5), subplot_kw=dict(aspect="equal"))
        biotype_counts = df['biotype'].value_counts()
        threshold = biotype_counts.sum() * 0.02
        other_count = biotype_counts[biotype_counts < threshold].sum()
        biotype_counts_main = biotype_counts[biotype_counts >= threshold]
        if other_count > 0:
            biotype_counts_main.loc['Other'] = other_count
        color_palette = ['#4682B4', '#FF8C00', '#2E8B57', '#5F9EA0', '#FFA500', '#3CB371', '#6495ED', '#FFD700']
        colors = color_palette[:len(biotype_counts_main)]
        wedges, texts, autotexts = ax.pie(biotype_counts_main, autopct='%1.1f%%', startangle=90,
                                          pctdistance=0.85, colors=colors,
                                          wedgeprops=dict(width=0.4, edgecolor='w'))
        ax.legend(wedges, biotype_counts_main.index, title="Biotypes", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
        plt.setp(autotexts, size=10, weight="bold")
        ax.set_title("Distribution of DMRs by Gene Biotype")
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        buf.seek(0)
        plots.append(base64.b64encode(buf.read()).decode())
        plt.close()

    display_cols = ['chr', 'start', 'end', 'length', 'nCG', 'meanMethy1',
                    'meanMethy2', 'diff.Methy', 'annotation', 'biotype', 'gene']
    display_cols = [c for c in display_cols if c in df.columns]
    df_display = df[display_cols].drop_duplicates()
    for col in ('meanMethy1', 'meanMethy2', 'diff.Methy'):
        if col in df_display.columns:
            df_display[col] = df_display[col].round(3)

    table_html = df_display.to_html(index=False, table_id="dmr_table", classes="display", escape=False)
    plots_html = '<div class="plot-container">'
    for plot_b64 in plots:
        plots_html += f'<div class="plot-box"><img src="data:image/png;base64,{plot_b64}" alt="summary plot"></div>'
    plots_html += '</div>'
except Exception as e:
    table_html = f"<p>Error processing data: {str(e)}</p>"
    plots_html = ""

html = f'''<!DOCTYPE html>
<html>
<head>
<title>DMR Report - ${alias} (${comparison_label}, ${out_label})</title>
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
<style>
    body {{font-family:"Segoe UI", Arial, sans-serif; margin:20px; background:#f5f5f5;}}
    .page-container {{max-width: 1600px; margin: auto;}}
    h1 {{color:#333; text-align:center; border-bottom: 1px solid #ccc; padding-bottom: 10px; margin-bottom: 30px;}}
    h2 {{color:#333; border-bottom: 1px solid #eee; padding-bottom: 8px; margin-top: 0; margin-bottom: 20px; font-size: 1.5em;}}
    .content-box {{background: white; padding: 25px; margin-bottom: 30px; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.05); border: 1px solid #e9ecef;}}
    .summary-text {{background:#e8f4f8; padding:15px; border-radius:5px; font-family:monospace;}}
    .plot-container {{display:flex; flex-wrap:wrap; justify-content:center; gap:30px;}}
    .plot-box {{padding:15px; border:1px solid #ddd; border-radius:8px; background-color:#fdfdfd; flex:1; min-width:400px; max-width:48%; text-align:center;}}
    .plot-box img {{max-width:100%; height:auto;}}
    table.dataTable {{width:100%!important; border-collapse: collapse;}}
    table.dataTable th, table.dataTable td {{ padding: 12px 15px; text-align: left; border-bottom: 1px solid #dddddd; border-left: none; border-right: none;}}
    table.dataTable thead th {{ background-color: #4CAF50; color: white; border-bottom: 2px solid #45a049; }}
    table.dataTable tbody tr:hover {{ background-color: #f1f1f1; }}
</style>
</head>
<body>
<div class="page-container">
    <h1>Differentially Methylated Regions (DMRs) Summary<br><small>${alias} &mdash; ${comparison_label} &mdash; ${out_label}</small></h1>
    <div class="content-box"><h2>Analysis Summary</h2><div class="summary-text">{summary_text if summary_text else "No summary available"}</div></div>
    <div class="content-box"><h2>Summary Plots</h2>{plots_html}</div>
    <div class="content-box"><h2>DMR Table</h2>{table_html}</div>
</div>
<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
<script>
    \$(document).ready(function() {{
        \$('#dmr_table').DataTable({{"pageLength": 25, "order": [[3, "desc"]]}});
    }});
</script>
</body>
</html>'''

with open("${out_name}", "w") as f:
    f.write(html)
"""
}


process plot_dmr_modbamtools {
    // Plots each significant DMR region from two BAMs (H1/H2, or sample/reference).
    // gene_list_path may be an empty file (plot everything) or a real gene list.
    label "dmr_plot_modbamtools"
    cpus 2
    memory 4.GB
    input:
        tuple val(alias), val(comparison_label), val(mod_char), path(annotated_tsv)
        tuple path(gtf_file), path(gtf_tbi_file)
        path gene_list_path
        tuple val(name1), path(bam1), path(bai1)
        tuple val(name2), path(bam2), path(bai2)
    output:
        path("*.html"), emit: dmr_plots, optional: true
        path("plot_status.log"), emit: status_log
    script:
        """
        genes_of_interest="genes_of_interest.txt"
        if [[ -s "${gene_list_path}" ]]; then
            cp "${gene_list_path}" "\${genes_of_interest}"
        else
            touch "\${genes_of_interest}"
        fi
        touch plot_status.log
        echo "0" > plot_count.tmp

        awk 'NR > 1 { key = \$1 FS \$2 FS \$3 FS \$NF; if (!seen[key]++) print }' ${annotated_tsv} | \\
        while IFS=\$'\\t' read -r chr start end length nCG meanMethy1 meanMethy2 diffMethy areaStat annotation_chr annotation_start annotation_end strand annotation biotype gene; do
            if [[ ! -s "\${genes_of_interest}" ]] || grep -q -w "\${gene}" "\${genes_of_interest}"; then
                region="\${chr}:\${start}-\${end}"
                output_prefix="${alias}.wf_mods.${comparison_label}.\${chr}_\${start}-\${end}_\${gene}"
                if modbamtools plot -r \${region} -g ${gtf_file} \\
                    -s ${bam1.baseName},${bam2.baseName} -p \${output_prefix} ${bam1} ${bam2} -o ./; then
                    echo \$(( \$(cat plot_count.tmp) + 1 )) > plot_count.tmp
                    echo "Plotted \${region} (\${gene})" >> plot_status.log
                else
                    echo "Error plotting \${region} (\${gene})" >> plot_status.log
                fi
            fi
        done
        echo "Total plots generated: \$(cat plot_count.tmp)" >> plot_status.log
        rm -f plot_count.tmp
        """
}


process plot_dmr_methylartist {
    label "dmr_plot_methylartist"
    cpus 2
    memory 4.GB
    input:
        tuple val(alias), val(comparison_label), val(mod_char), path(annotated_tsv)
        tuple path(gtf_file), path(gtf_tbi_file)
        path gene_list_path
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        tuple val(name1), path(bam1), path(bai1)
        tuple val(name2), path(bam2), path(bai2)
    output:
        path("*.png"), emit: dmr_plots, optional: true
        path("plot_status.log"), emit: status_log
    script:
        def mod_flag = mod_char == 'h' ? 'h' : 'm'
        """
        genes_of_interest="genes_of_interest.txt"
        if [[ -s "${gene_list_path}" ]]; then
            cp "${gene_list_path}" "\${genes_of_interest}"
        else
            touch "\${genes_of_interest}"
        fi
        touch plot_status.log
        echo "0" > plot_count.tmp

        awk 'NR > 1 { key = \$1 FS \$2 FS \$3 FS \$NF; if (!seen[key]++) print }' ${annotated_tsv} | \\
        while IFS=\$'\\t' read -r chr start end length nCG meanMethy1 meanMethy2 diffMethy areaStat annotation_chr annotation_start annotation_end strand annotation biotype gene; do
            if [[ ! -s "\${genes_of_interest}" ]] || grep -q -w "\${gene}" "\${genes_of_interest}"; then
                region="\${chr}:\${start}-\${end}"
                output_prefix="${alias}.wf_mods.${comparison_label}.\${chr}_\${start}-\${end}_\${gene}"
                if methylartist locus --interval \${region} --gtf ${gtf_file} \\
                    --bams ${bam1},${bam2} --ref ${ref} --motif CG --mods ${mod_flag} \\
                    --outfile \${output_prefix} --labelgenes --nomask; then
                    echo \$(( \$(cat plot_count.tmp) + 1 )) > plot_count.tmp
                    echo "Plotted \${region} (\${gene})" >> plot_status.log
                else
                    echo "Error plotting \${region} (\${gene})" >> plot_status.log
                fi
            fi
        done
        echo "Total methylartist plots generated: \$(cat plot_count.tmp)" >> plot_status.log
        rm -f plot_count.tmp
        """
}


process plot_phased_dmr_modbamtools {
    // Haplotype-compare variant: a single haplotagged BAM, plotted with
    // modbamtools' own -hp (per-haplotype) mode instead of two separate BAM
    // tracks -- there's only one physical BAM here (H1/H2 are HP tags within
    // it), unlike dmr_sample_compare where two real BAMs exist.
    label "dmr_plot_modbamtools"
    cpus 2
    memory 4.GB
    input:
        tuple val(alias), val(comparison_label), val(mod_char), path(annotated_tsv)
        tuple path(gtf_file), path(gtf_tbi_file)
        path gene_list_path
        tuple val(name), path(bam), path(bai)
    output:
        path("*.html"), emit: dmr_plots, optional: true
        path("plot_status.log"), emit: status_log
    script:
        """
        genes_of_interest="genes_of_interest.txt"
        if [[ -s "${gene_list_path}" ]]; then
            cp "${gene_list_path}" "\${genes_of_interest}"
        else
            touch "\${genes_of_interest}"
        fi
        touch plot_status.log
        echo "0" > plot_count.tmp

        awk 'NR > 1 { key = \$1 FS \$2 FS \$3 FS \$NF; if (!seen[key]++) print }' ${annotated_tsv} | \\
        while IFS=\$'\\t' read -r chr start end length nCG meanMethy1 meanMethy2 diffMethy areaStat annotation_chr annotation_start annotation_end strand annotation biotype gene; do
            if [[ ! -s "\${genes_of_interest}" ]] || grep -q -w "\${gene}" "\${genes_of_interest}"; then
                region="\${chr}:\${start}-\${end}"
                output_prefix="${alias}.wf_mods.${comparison_label}.\${chr}_\${start}-\${end}_\${gene}"
                if modbamtools plot -r \${region} -g ${gtf_file} \\
                    -s ${bam.baseName} -p \${output_prefix} -hp ${bam} -o ./; then
                    echo \$(( \$(cat plot_count.tmp) + 1 )) > plot_count.tmp
                    echo "Plotted phased \${region} (\${gene})" >> plot_status.log
                else
                    echo "Error plotting phased \${region} (\${gene})" >> plot_status.log
                fi
            fi
        done
        echo "Total phased plots generated: \$(cat plot_count.tmp)" >> plot_status.log
        rm -f plot_count.tmp
        """
}


process plot_phased_dmr_methylartist {
    label "dmr_plot_methylartist"
    cpus 2
    memory 4.GB
    input:
        tuple val(alias), val(comparison_label), val(mod_char), path(annotated_tsv)
        tuple path(gtf_file), path(gtf_tbi_file)
        path gene_list_path
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        tuple val(name), path(bam), path(bai)
    output:
        path("*.png"), emit: dmr_plots, optional: true
        path("plot_status.log"), emit: status_log
    script:
        def mod_flag = mod_char == 'h' ? 'h' : 'm'
        """
        genes_of_interest="genes_of_interest.txt"
        if [[ -s "${gene_list_path}" ]]; then
            cp "${gene_list_path}" "\${genes_of_interest}"
        else
            touch "\${genes_of_interest}"
        fi
        touch plot_status.log
        echo "0" > plot_count.tmp

        awk 'NR > 1 { key = \$1 FS \$2 FS \$3 FS \$NF; if (!seen[key]++) print }' ${annotated_tsv} | \\
        while IFS=\$'\\t' read -r chr start end length nCG meanMethy1 meanMethy2 diffMethy areaStat annotation_chr annotation_start annotation_end strand annotation biotype gene; do
            if [[ ! -s "\${genes_of_interest}" ]] || grep -q -w "\${gene}" "\${genes_of_interest}"; then
                region="\${chr}:\${start}-\${end}"
                output_prefix="${alias}.wf_mods.${comparison_label}.\${chr}_\${start}-\${end}_\${gene}"
                if methylartist locus --interval \${region} --gtf ${gtf_file} \\
                    --bams ${bam} --ref ${ref} --motif CG --mods ${mod_flag} \\
                    --outfile \${output_prefix} --labelgenes --nomask --phased --ignore_ps; then
                    echo \$(( \$(cat plot_count.tmp) + 1 )) > plot_count.tmp
                    echo "Plotted phased \${region} (\${gene})" >> plot_status.log
                else
                    echo "Error plotting phased \${region} (\${gene})" >> plot_status.log
                fi
            fi
        done
        echo "Total methylartist phased plots generated: \$(cat plot_count.tmp)" >> plot_status.log
        rm -f plot_count.tmp
        """
}
