# OMICS (UNDER CONSTRUCTION)

PRE-ALFA

All-in-one workflow from base calling to ACMG/AMP classification and reports for Single-Nucleotide Variant (SNV), Structural Variant (SV), Copy Number Variation (CNV), Short Tandem Repeats (STR) and Differential Metilation (DMR) and .

## Introduction

The core of this workflow is based on https://github.com/epi2me-labs/wf-human-variation 

This repository contains a [nextflow](https://www.nextflow.io/) workflow for analysing variation in human genomic data. Specifically this workflow can perform the following:

* diploid variant calling
* single nucleotide and small indels variants (SNVs) calling
* structural variants (SVs) calling
* modified bases and differentially methylated regions (DMRs) calling
* copy number variants (CNVs) calling
* short tandem repeats (STRs) expansion genotyping
* transposable elements insertions (TEIs) identification
* pathogenicity prediction (ACMG/AMP classification) of SNVs with TAPES
* pathogenicity prediction (ACMG/AMP classification) of SVs with AnnotSV
* pathogenicity prediction (ACMG/AMP classification) of CNVs with AnnotSV/ClassifyCNV/ISV
* pathogenicity prediction of STRs with STRchive
* graphical reports for each type of result

<figure>
<img src="docs/images/wf-omics.drawio.svg" alt="wf-omics overview schematic."/>
<figcaption>Schematic depicting wf-omics workflow.</figcaption>
</figure>

The tools embedded in individual sub-workflows within wf-human-variation are specifically designed for use with whole-genome Oxford Nanopore Technologies sequencing data. While 20x average coverage is the absolute minimum requirement for the workflow to run, we recommend an average coverage above 30x to ensure optimal performance. Usage below the minimum coverage may cause the workflow to terminate with an error, or yield unexpected outcomes.




## Compute requirements

Recommended requirements:

+ CPUs = 32
+ Memory = 128GB

Minimum requirements:

+ CPUs = 16
+ Memory = 32GB

Approximate run time: Variable depending on whether it is targeted sequencing or whole genome sequencing, as well as coverage and the individual analyses requested. For instance, a 90X human sample run (options: `--snp --sv --mod --str --cnv --phased --sex XY`) takes less than 8h with recommended resources.

