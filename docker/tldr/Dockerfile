# tldr (https://github.com/adamewing/tldr) -- identifies and annotates
# transposable-element-mediated insertions in long-read alignments.
#
# Meant to run against the per-contig, pre-merge haplotagged BAMs already
# produced during SNP calling/phasing (see workflows/wf-human-snp.nf,
# `haplotagged_ctg_bams` / the `str_bams` emit -- same channel straglr's
# call_str already consumes), not the final whole-genome BAM.
#
# NOTE on dependencies -- checked tldr's own source (tldr/tldr) directly,
# don't assume from its tldr.yml conda recipe alone:
#   - minimap2 IS required, even though the input BAM is already aligned:
#     tldr's own algorithm extracts soft-clipped candidate reads, builds a
#     local consensus with mafft, then re-maps that consensus against both
#     the TE reference library and the genome via `subprocess`-invoked
#     minimap2 to confirm/annotate the insertion. This is internal to tldr's
#     own logic, unrelated to how the input BAM was produced upstream.
#   - samtools (CLI binary) IS required: tldr shells out to it directly
#     (`samtools faidx`, `index`, `view`, `sort`, `mpileup`).
#   - tabix (the standalone CLI binary) is not actually shelled out to by
#     tldr's own code -- the only tabix-file access in its source is
#     `pysam.Tabixfile(...)` (for the optional `--nonref` annotation file),
#     which uses pysam's own bundled htslib bindings. Included anyway (a few
#     hundred KB) so a .tbi can be built inside this container too, e.g. for
#     preparing a `--nonref` file, rather than only ever reading one.
#   - exonerate and mafft are also required (subprocess calls for the
#     assembly/alignment steps of tldr's own insertion-calling algorithm).
#
# All five (samtools, minimap2, mafft, exonerate, tabix) are in Debian
# bookworm's main archive at versions matching or newer than tldr.yml's own
# pins (samtools 1.16.1 vs >=1.2 required, minimap2 2.24 exact match, mafft
# 7.505 vs 7.490, exonerate 2.4.0 exact match, tabix 1.16) -- apt is used
# instead of conda here, consistent with the other custom images in this repo
# (docker/annotsv, docker/cnv_ensemble), and avoids pulling in a conda/mamba
# base image just for five packages that Debian already ships.
#
# tldr's own setup.py pins are loose lower bounds only (cython>=0.29.0,
# pysam>=0.8.1, bx-python>=0.5.0, scikit-learn>=0.20.0, numpy>=1.9.0) -- no
# upper bounds, so no repeat of the isv/cnv_ensemble pinned-dependency mess
# (see docker/cnv_ensemble/Dockerfile) is expected; modern versions are
# installed via pip with no extra pinning.
#
# Published at romualdomf/tldr (hardcoded in base.config, withLabel:tldr --
# same reasoning as annotsv/cnv_ensemble there: Docker Hub pull rate limits
# are keyed to whoever pulls, not the publishing account).
#   docker build -t romualdomf/tldr:1.2.2 docker/tldr
#   docker push romualdomf/tldr:1.2.2

FROM debian:bookworm-slim

ARG TLDR_VERSION=1.2.2

RUN apt-get update && apt-get install -y --no-install-recommends \
        git ca-certificates bash coreutils \
        samtools minimap2 mafft exonerate tabix \
        python3 python3-pip \
        procps \
    && rm -rf /var/lib/apt/lists/*
# procps (ps) isn't a dependency of tldr itself -- Nextflow execs `ps` inside
# the container to collect task resource metrics, and errors out without it
# (bit us before building the AnnotSV and cnv-ensemble images).

# Debian bookworm's pip refuses plain "pip install" (PEP 668,
# externally-managed-environment) without this.
ENV PIP_BREAK_SYSTEM_PACKAGES=1

# tldr's setup.py checks minimap2/samtools/mafft/exonerate are already on
# PATH at install time (all installed via apt above), then installs its own
# loosely-pinned Python deps (cython, pysam, bx-python, scikit-learn, numpy)
# via pip automatically -- no need to list them separately here.
RUN git clone --branch "${TLDR_VERSION}" --depth 1 \
        https://github.com/adamewing/tldr.git /opt/tldr \
    && pip install --no-cache-dir /opt/tldr

# tldr's own bundled human TE reference library (no external ~GB-scale
# download needed, unlike VEP's cache or AnnotSV's human annotations).
ENV TLDR_REF_DIR=/opt/tldr/ref

WORKDIR /data
