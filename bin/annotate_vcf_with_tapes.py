#!/usr/bin/env python3
"""Write TAPES' Probability_Path / Prediction_ACMG_tapes columns back into the SNP
VCF as INFO fields.

TAPES (see modules/local/tapes.nf, bin/annotate_vcf_from_annotsv.py for the
equivalent CNV/SV-side script) reads the VEP-annotated SNP VCF directly. It has
a decompose-multiallelics code path (tapes.py's src/t_func.py::decompose_vcf,
gated by test_if_decomposed sampling the input's first 2000 lines), but on a
real run against this pipeline's per-contig VCFs it never actually decomposes --
confirmed empirically: every one of a real output's multiallelic ALT lists
(58,847 of them, matching the source VCF's multiallelic record count exactly)
comes back from TAPES exactly as comma-joined as the VCF's own ALT field, not
split into one row per allele. So each TAPES row corresponds to one whole VCF
record (all of its ALT alleles together), matched by (Chr, Start, Ref, Alt) with
Alt taken as the full comma-joined list -- not decomposed per-allele -- and Chr
*without* the "chr" prefix the VCF's own CHROM values carry. Both injected INFO
fields are therefore Number=1 (one value per record), not Number=A.

The TAPES table is ~150 columns and can run into the millions of rows for a
whole-genome SNP set -- unlike bin/annotate_vcf_from_annotsv.py (which loads
every column via csv.DictReader for AnnotSV's much smaller CNV/SV row counts),
this only ever keeps the two target columns per row to stay memory-reasonable.

Usage:
    annotate_vcf_with_tapes.py \
        --tapes-tsv OMICS_09.wf_snp.tapes.txt \
        --input-vcf OMICS_09.wf_snp.vcf.gz \
        --output-vcf OMICS_09.wf_snp.annotated.vcf.gz
"""
import argparse
import csv
import sys
from pathlib import Path

import pysam

# same rationale as bin/annotate_vcf_from_annotsv.py -- TAPES rows can exceed the
# csv module's default 131072-byte field limit.
csv.field_size_limit(sys.maxsize)

REQUIRED_COLUMNS = ('Chr', 'Start', 'Ref', 'Alt', 'Probability_Path', 'Prediction_ACMG_tapes')


def load_tapes_lookup(tsv_path):
    """Returns {(chrom, pos, ref, alt): (Probability_Path, Prediction_ACMG_tapes)},
    chrom carrying the same "chr" prefix as the VCF's own CHROM column, and alt
    exactly as TAPES wrote it -- the full comma-joined ALT list for multiallelic
    records, not decomposed (see module docstring).
    """
    lookup = {}
    with open(tsv_path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        try:
            idx = {col: header.index(col) for col in REQUIRED_COLUMNS}
        except ValueError as exc:
            sys.exit(f'{tsv_path}: missing expected column ({exc})')

        row_count = 0
        for row in reader:
            row_count += 1
            chrom = row[idx['Chr']]
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom
            key = (chrom, row[idx['Start']], row[idx['Ref']], row[idx['Alt']])
            lookup[key] = (row[idx['Probability_Path']], row[idx['Prediction_ACMG_tapes']])
    print(f'{row_count} row(s) loaded from {tsv_path}', file=sys.stderr)
    return lookup


def annotate_vcf(lookup, input_vcf, output_vcf):
    vcf_in = pysam.VariantFile(input_vcf)
    vcf_in.header.info.add(
        'TAPES_Probability_Path', 1, 'Float',
        'Probability_Path column from TAPES ACMG classification')
    vcf_in.header.info.add(
        'TAPES_Prediction_ACMG', 1, 'String',
        'Prediction_ACMG_tapes column from TAPES ACMG classification')

    vcf_out = pysam.VariantFile(output_vcf, 'wz', header=vcf_in.header)

    record_count = 0
    matched_count = 0
    for record in vcf_in:
        record_count += 1
        alt = ','.join(record.alts or ())
        row = lookup.get((record.chrom, str(record.pos), record.ref, alt))
        if row is not None:
            matched_count += 1
            record.info['TAPES_Probability_Path'] = float(row[0])
            # VCF INFO values can't contain raw whitespace (not a delimiter, but
            # several parsers still choke on it) -- TAPES' own values like
            # "Likely Benign" need the space swapped out.
            record.info['TAPES_Prediction_ACMG'] = row[1].replace(' ', '_')
        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()
    pysam.tabix_index(output_vcf, preset='vcf', force=True)
    return record_count, matched_count


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--tapes-tsv', required=True, type=Path)
    parser.add_argument('--input-vcf', required=True, type=Path)
    parser.add_argument('--output-vcf', required=True, type=Path)
    args = parser.parse_args()

    lookup = load_tapes_lookup(args.tapes_tsv)
    record_count, matched_count = annotate_vcf(
        lookup, str(args.input_vcf), str(args.output_vcf))
    print(
        f'{record_count} VCF record(s), {matched_count} matched to a TAPES row '
        f'-> {args.output_vcf}',
        file=sys.stderr)


if __name__ == '__main__':
    main()
