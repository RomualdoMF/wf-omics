#!/usr/bin/env python3
"""Write an AnnotSV-derived TSV's columns back into a VCF as INFO fields.

Accepts two shapes of input TSV, auto-detected by its header:

  - AnnotSV's own raw output (<sample>.wf_sv.annotsv.tsv /
    <sample>.wf_cnv.annotsv.tsv, see modules/local/annotsv.nf run_annotsv) --
    columns like SV_chrom/SV_start/SV_end/SV_type plus everything else
    AnnotSV computed. Only the "full" rows (one per variant; AnnotSV also
    emits "split" rows, one per overlapping gene) are used. Every column
    except the coordinate/ID ones is prefixed AnnotSV_ in the output.

    Matched to VCF records by (CHROM, ID) -- AnnotSV's own TSV carries the
    source VCF's ID verbatim (e.g. "Sniffles2.DEL.364CS0"), which is unique
    per record and sidesteps a real mismatch: AnnotSV reports translocations
    as SV_type=TRA, but Sniffles2's VCF calls the same records SVTYPE=BND --
    matching on (CHROM, POS, SVTYPE) as before meant "TRA" != "BND" always
    failed to match, so every BND record silently got zero AnnotSV
    annotation (confirmed on a real run: 0 of 124 BND records annotated,
    vs ~62% for DEL). ID-based matching doesn't care what either side calls
    the SV type.

  - An already-merged ensemble TSV (<sample>.wf_cnv.annotated.tsv, produced
    by bin/cnv_ensemble_classify.py) -- already has Chr/Start/Stop/Type as
    its first four columns and every other column already prefixed
    (AnnotSV_ / ISV_ / ClassifyCNV_), so it's used as-is. This shape has no
    ID column (ClassifyCNV/ISV don't carry it through), so it keeps matching
    on (CHROM, POS, SVTYPE) as before -- CNV callers don't have an
    equivalent BND-vs-TRA naming clash, so this was never the bug.

Either way, every remaining column is added as its own INFO field
(sanitizing the column name into a valid VCF INFO ID, and the value into one
that doesn't break the INFO field's own `;`/`=`/whitespace-delimited
syntax). Records with no matching TSV row (e.g. an alt-contig variant
ClassifyCNV skipped, or one AnnotSV didn't annotate) are passed through
unannotated.

This is what both the SV and CNV paths use to build their own
.wf_sv.annotated.vcf.gz / .wf_cnv.annotated.vcf.gz -- instead of AnnotSV's
native `-vcf 1` output (via its bundled variantconvert tool), which produced
a header-only VCF with no records in every mode tried, and has since been
removed from the pipeline's AnnotSV image entirely (see
docker/annotsv/Dockerfile). Only dependency is pysam -- no reason to pull in
pandas for what's fundamentally a row-matching join.

Wired into the pipeline as the annotate_vcf_with_tsv process in
modules/local/annotsv.nf (label "wf_common", so it reuses the pipeline's
existing ontresearch/wf-common image -- that already has pysam installed for
other reporting steps, no need for a dedicated container). Lives in bin/ (not
scripts/) so Nextflow's automatic bin/-on-PATH convention picks it up.

Usage:
    annotate_vcf_from_annotsv.py \
        --annotated-tsv OMICS_09.wf_sv.annotsv.tsv \
        --input-vcf OMICS_09.wf_sv.vcf.gz \
        --output-vcf OMICS_09.wf_sv.annotated.vcf.gz
"""
import argparse
import csv
import re
import sys
from pathlib import Path

import pysam

# AnnotSV rows for large/complex SVs can pack a huge Gene_name / overlapping-gene
# list into a single field (seen on real SV calls spanning many genes), well past
# the csv module's default 131072-byte field limit -- raise it before any TSV read.
csv.field_size_limit(sys.maxsize)

# ensemble (CNV) TSV shape: no ID column, keeps matching on these four.
CANONICAL_COLUMNS = ['Chr', 'Start', 'Stop', 'Type']
# raw AnnotSV TSV shape: matched on these two instead (see module docstring).
RAW_MATCH_COLUMNS = ['Chr', 'ID']
# raw AnnotSV columns never re-emitted as their own AnnotSV_* INFO field:
#   - SV_chrom/ID: consumed to build the match key (Chr/ID above), redundant
#     with the VCF's own CHROM/ID.
#   - SV_start/SV_end: redundant with the VCF's own POS/END/CIPOS-derived
#     coordinates.
#   - INFO/FORMAT: AnnotSV copies the source VCF's raw INFO/FORMAT strings
#     into these columns verbatim -- re-emitting them as AnnotSV_INFO/
#     AnnotSV_FORMAT would nest a whole INFO-field's worth of ';'/'='-delimited
#     text inside a single sanitized INFO value, unreadable and redundant with
#     the output VCF's own INFO/FORMAT columns.
# SV_type is kept (not in this list): AnnotSV's own classification can
# legitimately differ from the VCF's SVTYPE (e.g. TRA vs BND), so it's worth
# surfacing as AnnotSV_SV_type rather than silently dropped.
RAW_ANNOTSV_CONSUMED_COLUMNS = ('SV_chrom', 'SV_start', 'SV_end', 'ID', 'INFO', 'FORMAT')

_INFO_ID_RE = re.compile(r'[^0-9A-Za-z_.]')
_WHITESPACE_RE = re.compile(r'\s+')


def _sanitize_info_id(name):
    """VCF INFO IDs may only contain [0-9A-Za-z_.]."""
    return _INFO_ID_RE.sub('_', name)


def _sanitize_info_value(value):
    """VCF INFO values may not contain ';', '=', or whitespace (unescaped)."""
    value = value.replace(';', ',').replace('=', ':')
    return _WHITESPACE_RE.sub('_', value.strip())


def _load_raw_annotsv_rows(reader):
    """Normalise AnnotSV's own raw TSV output into the canonical shape,
    keeping only "full" rows (one per variant -- AnnotSV also emits "split"
    rows, one per overlapping gene, which would each try to claim the same
    VCF record).
    """
    rows = [row for row in reader if row.get('Annotation_mode') == 'full']
    other_cols = [c for c in reader.fieldnames if c not in RAW_ANNOTSV_CONSUMED_COLUMNS]

    normalized = []
    for row in rows:
        new_row = {
            'Chr': 'chr' + re.sub(r'^chr', '', row['SV_chrom']),
            'ID': row['ID'],
        }
        for col in other_cols:
            new_row[f'AnnotSV_{col}'] = row[col]
        normalized.append(new_row)

    annotation_cols = [f'AnnotSV_{c}' for c in other_cols]
    return normalized, annotation_cols, RAW_MATCH_COLUMNS


def load_annotated_tsv(tsv_path):
    """Returns (rows, annotation_cols, match_keys): rows keyed by match_keys
    (either ['Chr', 'ID'] for AnnotSV's own raw output, or the older
    ['Chr', 'Start', 'Stop', 'Type'] for the merged ensemble TSV -- see
    module docstring for why), and annotation_cols listing every other
    column to inject into the VCF.
    """
    with open(tsv_path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if 'Chr' in reader.fieldnames:
            # already-canonical ensemble TSV (bin/cnv_ensemble_classify.py output)
            rows = list(reader)
            annotation_cols = [c for c in reader.fieldnames if c not in CANONICAL_COLUMNS]
            return rows, annotation_cols, CANONICAL_COLUMNS
        return _load_raw_annotsv_rows(reader)


def _vcf_record_key(record, match_keys):
    if match_keys == RAW_MATCH_COLUMNS:
        return (record.chrom, record.id)
    return (record.chrom, str(record.pos), record.info.get('SVTYPE'))


def annotate_vcf(rows, annotation_cols, match_keys, input_vcf, output_vcf):
    lookup = {
        tuple(row[k] for k in match_keys): row
        for row in rows
    }

    vcf_in = pysam.VariantFile(input_vcf)
    for col in annotation_cols:
        vcf_in.header.info.add(
            _sanitize_info_id(col), 1, 'String',
            f'{col} column from AnnotSV annotation')

    vcf_out = pysam.VariantFile(output_vcf, 'wz', header=vcf_in.header)

    annotated_count = 0
    for record in vcf_in:
        key = _vcf_record_key(record, match_keys)
        row = lookup.get(key)
        if row is not None:
            annotated_count += 1
            for col in annotation_cols:
                value = row.get(col, '')
                if value in ('', 'nan', None):
                    continue
                record.info[_sanitize_info_id(col)] = _sanitize_info_value(value)
        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()
    pysam.tabix_index(output_vcf, preset='vcf', force=True)
    return annotated_count


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '--annotated-tsv', required=True, type=Path,
        help='AnnotSV raw TSV or a merged ensemble TSV (auto-detected)')
    parser.add_argument(
        '--input-vcf', required=True, type=Path,
        help='Original VCF the TSV was derived from (.wf_sv.vcf.gz / .wf_cnv.vcf.gz)')
    parser.add_argument(
        '--output-vcf', required=True, type=Path,
        help='Output path (<sample>.wf_sv.annotated.vcf.gz / .wf_cnv.annotated.vcf.gz)')
    args = parser.parse_args()

    rows, annotation_cols, match_keys = load_annotated_tsv(args.annotated_tsv)
    print(
        f'{len(rows)} row(s), {len(annotation_cols)} column(s) loaded from '
        f'{args.annotated_tsv}, matching on {match_keys}',
        file=sys.stderr)

    annotated_count = annotate_vcf(
        rows, annotation_cols, match_keys, str(args.input_vcf), str(args.output_vcf))
    print(f'Annotated {annotated_count} record(s) in {args.output_vcf}', file=sys.stderr)


if __name__ == '__main__':
    main()
