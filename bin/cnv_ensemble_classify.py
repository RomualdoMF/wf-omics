#!/usr/bin/env python3
"""Run an ensemble of CNV annotation/classification tools on top of AnnotSV output.

Takes a `<sample>.wf_cnv.annotsv.tsv` (produced by AnnotSV, see
modules/local/annotsv.nf) and runs it through:
  - ISV        https://github.com/tsladecek/isv_package  (annotate + predict)
  - ClassifyCNV https://github.com/Genotek/ClassifyCNV

then merges the AnnotSV columns with each tool's output columns (each set
prefixed with the tool's name) into a single `<sample>.wf_cnv.annotated.tsv`.
The first four (unprefixed) columns of the output are the canonical
coordinates -- Chr, Start, Stop, Type -- taken from AnnotSV; every other
column is prefixed AnnotSV_ / ISV_ / ClassifyCNV_ and rows are joined back on
those four canonical columns (not on row order), so a tool dropping a CNV
(e.g. ClassifyCNV skips alt-contig / zero-length CNVs) just leaves NaNs for
that tool's columns rather than desynchronising the table.

Writing those same columns back into the original .wf_cnv.vcf.gz as INFO
fields (producing .wf_cnv.annotated.vcf.gz) is a separate script,
bin/annotate_cnv_vcf.py -- deliberately not done here, see that script's
docstring for why.

X-CNV (https://github.com/kbvstmd/XCNV) was intentionally left out: it only
supports GRCh37/hg19 (this pipeline's CNVs are GRCh38) and needs ~3GB of data
from a single non-CDN server -- out of scope for now.

Wired into the pipeline as modules/local/cnv_ensemble.nf (process
run_cnv_ensemble, label "cnv_ensemble"), which runs this script inside the
romualdomf/cnv-ensemble container (docker/cnv_ensemble/Dockerfile) -- that
image already has ISV and ClassifyCNV installed, with ClassifyCNV's checkout
at $CLASSIFYCNV_DIR (baked in, no external data download needed). This file
lives in bin/ (not scripts/) for that reason: Nextflow puts bin/ on PATH for
every process regardless of container, so it's called bare as
`cnv_ensemble_classify.py`, not via an absolute path.

Can also be run standalone, outside the pipeline/container, e.g. for
debugging -- see Setup/Usage below.

--------------------------------------------------------------------------
Setup (only needed to run this standalone, outside the pipeline container)
--------------------------------------------------------------------------
1. ISV (pip package, works with modern numpy/pandas/scikit-learn despite its
   own pinned requirements being years out of date -- see _patch_sklearn_gb_losses
   below for the one compatibility shim that's actually needed):
       pip install isv --no-deps
       pip install numpy pandas xgboost shap numba plotly sklearn-json scikit-learn

2. ClassifyCNV (not on PyPI; clone it and pass --classifycnv-dir):
       git clone https://github.com/Genotek/ClassifyCNV.git
   Needs bedtools >=2.27.1 on PATH.

--------------------------------------------------------------------------
Usage
--------------------------------------------------------------------------
    cnv_ensemble_classify.py \
        --annotsv-tsv OMICS_09.wf_cnv.annotsv.tsv \
        --classifycnv-dir /path/to/ClassifyCNV \
        --genome-build hg38 \
        --output-tsv OMICS_09.wf_cnv.annotated.tsv
"""
import argparse
import subprocess
import sys
import tempfile
import types
from pathlib import Path

import pandas as pd

CANONICAL_COLUMNS = ['Chr', 'Start', 'Stop', 'Type']


def _patch_sklearn_gb_losses():
    """isv depends on the unmaintained `sklearn_json` package (last released
    2021), which does `from sklearn.ensemble import ..., _gb_losses` at
    import time -- a private submodule removed in modern scikit-learn. That
    import is only ever used by sklearn_json to deserialize a
    GradientBoostingClassifier; isv's actual bundled model is XGBoost-based,
    so that code path never runs -- a dummy stand-in module is enough to let
    the (otherwise-unused) top-level import succeed.
    """
    import sklearn.ensemble as _ens
    if hasattr(_ens, '_gb_losses'):
        return
    shim = types.ModuleType('sklearn.ensemble._gb_losses')
    for name in ('BinomialDeviance', 'ExponentialLoss', 'MultinomialDeviance'):
        setattr(shim, name, type(name, (), {}))
    sys.modules['sklearn.ensemble._gb_losses'] = shim
    _ens._gb_losses = shim


_patch_sklearn_gb_losses()
import isv  # noqa: E402 -- must be imported after the shim above


def load_annotsv_cnvs(annotsv_tsv):
    """Read an AnnotSV output TSV, keep one row per CNV.

    AnnotSV emits one "full" row per structural variant plus additional
    "split" rows (one per overlapping gene) -- only the "full" rows are one
    row per CNV, which is what every downstream tool here expects.
    """
    df = pd.read_csv(annotsv_tsv, sep='\t', dtype=str)
    full = df[df['Annotation_mode'] == 'full'].copy()

    # AnnotSV strips the "chr" prefix internally (SV_chrom is e.g. "4", not
    # "chr4"); put it back for consistency with the pipeline's own VCF
    # contig naming (chr1, chr2, ... -- see the wf_cnv.vcf.gz header).
    full['Chr'] = 'chr' + full['SV_chrom'].str.replace(r'^chr', '', regex=True)
    full['Start'] = full['SV_start']
    full['Stop'] = full['SV_end']
    full['Type'] = full['SV_type']

    other_cols = [c for c in df.columns if c not in ('SV_chrom', 'SV_start', 'SV_end', 'SV_type')]
    full = full.rename(columns={c: f'AnnotSV_{c}' for c in other_cols})
    ordered = CANONICAL_COLUMNS + [f'AnnotSV_{c}' for c in other_cols]
    return full[ordered].reset_index(drop=True)


def run_isv(cnvs):
    """Run isv.annotate() + isv.predict() on the CNV list.

    Returns a DataFrame with Chr/Start/Stop/Type plus ISV_-prefixed columns
    (the annotate() feature columns and the predict() pathogenicity score).
    """
    # isv's numba-jitted internals expect real ints for start/end, not the
    # strings load_annotsv_cnvs() carries everything else as.
    cnv_list = [
        [chrom, int(start), int(stop), cnv_type]
        for chrom, start, stop, cnv_type in cnvs[CANONICAL_COLUMNS].values.tolist()
    ]
    annotated = isv.annotate(cnv_list)
    scores = isv.predict(annotated, proba=True)

    result = annotated.rename(columns={
        'chrom': 'Chr', 'start': 'Start', 'end': 'Stop', 'cnv_type': 'Type',
    })
    result['Start'] = result['Start'].astype(str)
    result['Stop'] = result['Stop'].astype(str)
    result['pathogenicity_score'] = scores

    other_cols = [c for c in result.columns if c not in CANONICAL_COLUMNS]
    result = result.rename(columns={c: f'ISV_{c}' for c in other_cols})
    return result[CANONICAL_COLUMNS + [f'ISV_{c}' for c in other_cols]]


def run_classifycnv(cnvs, classifycnv_dir, genome_build):
    """Run ClassifyCNV (subprocess) on the CNV list.

    Returns a DataFrame with Chr/Start/Stop/Type plus ClassifyCNV_-prefixed
    columns (Scoresheet.txt, minus its own VariantID/coordinate columns).
    """
    classifycnv_dir = Path(classifycnv_dir)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        bed_path = tmp / 'input.bed'
        cnvs[CANONICAL_COLUMNS].to_csv(bed_path, sep='\t', header=False, index=False)
        outdir = tmp / 'results'
        subprocess.run(
            [sys.executable, str(classifycnv_dir / 'ClassifyCNV.py'),
             '--infile', str(bed_path),
             '--GenomeBuild', genome_build,
             '--outdir', str(outdir)],
            check=True, cwd=classifycnv_dir,
        )
        result = pd.read_csv(outdir / 'Scoresheet.txt', sep='\t', dtype=str)

    result = result.rename(columns={'Chromosome': 'Chr', 'End': 'Stop'})
    other_cols = [c for c in result.columns if c not in CANONICAL_COLUMNS + ['VariantID']]
    result = result.rename(columns={c: f'ClassifyCNV_{c}' for c in other_cols})
    return result[CANONICAL_COLUMNS + [f'ClassifyCNV_{c}' for c in other_cols]]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '--annotsv-tsv', required=True, type=Path,
        help='AnnotSV output TSV (<sample>.wf_cnv.annotsv.tsv)')
    parser.add_argument(
        '--classifycnv-dir', required=True, type=Path,
        help='Path to a cloned Genotek/ClassifyCNV checkout')
    parser.add_argument(
        '--genome-build', default='hg38', choices=['hg19', 'hg38'],
        help='Genome build of the input coordinates (default: hg38)')
    parser.add_argument(
        '--output-tsv', required=True, type=Path,
        help='Output path (<sample>.wf_cnv.annotated.tsv)')
    args = parser.parse_args()

    cnvs = load_annotsv_cnvs(args.annotsv_tsv)
    print(f'{len(cnvs)} CNV(s) loaded from {args.annotsv_tsv}', file=sys.stderr)

    isv_result = run_isv(cnvs)
    print(f'ISV: annotated {len(isv_result)} CNV(s)', file=sys.stderr)

    classifycnv_result = run_classifycnv(cnvs, args.classifycnv_dir, args.genome_build)
    print(f'ClassifyCNV: classified {len(classifycnv_result)} CNV(s)', file=sys.stderr)

    merged = cnvs.merge(isv_result, on=CANONICAL_COLUMNS, how='left')
    merged = merged.merge(classifycnv_result, on=CANONICAL_COLUMNS, how='left')

    merged.to_csv(args.output_tsv, sep='\t', index=False)
    print(f'Wrote {len(merged)} row(s), {len(merged.columns)} column(s) to {args.output_tsv}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
