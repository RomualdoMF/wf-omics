# AnnotSV (https://github.com/lgmgeo/AnnotSV) container for wf-omics.
#
# Bakes in AnnotSV itself (Tcl/Tk code only, no variantconvert -- see below)
# -- NOT the human annotation data (Annotations_Human_*.tar.gz is ~5GB). That
# data is downloaded
# separately, on the host, into params.annotsv_annotations_dir and bind-mounted
# into the container at runtime (see base.config, withLabel:annotsv
# containerOptions, and modules/local/annotsv.nf download_annotsv_annotations).
#
# Published at romualdomf/annotsv (hardcoded in base.config, withLabel:annotsv --
# see the note there for why: Docker Hub pull rate limits are keyed to the
# puller, not the publisher, so there's no benefit to making this configurable).
# If you rebuild with a new ANNOTSV_VERSION, bump params.wf.annotsv_tag in
# base.config to match.
#   docker build -t romualdomf/annotsv:3.5.10 docker/annotsv
#   docker push romualdomf/annotsv:3.5.10

FROM debian:bookworm-slim

ARG ANNOTSV_VERSION=3.5.10

# dependency list per AnnotSV's own bioconda recipe (run requirements):
# https://github.com/bioconda/bioconda-recipes/blob/master/recipes/annotsv/meta.yaml
# -- minus python3/pip and its variantconvert-only packages (pyfaidx, natsort,
# pandas, tqdm, polars, appdirs), see the note below. AnnotSV's own core
# (Tcl/Tk) annotation engine doesn't touch Python at all.
# (procps is not an AnnotSV dependency -- Nextflow itself execs `ps` inside
# the container to collect task resource metrics, and errors out without it)
RUN apt-get update && apt-get install -y --no-install-recommends \
        make git curl unzip ca-certificates coreutils bash \
        tcl8.6 tk8.6 \
        bedtools bcftools \
        openjdk-17-jre-headless \
        procps \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf "$(command -v tclsh8.6)" /usr/bin/tclsh
# (that last symlink: Debian's tcl8.6 package only ships tclsh8.6; AnnotSV's
# Makefile shells out to a plain "tclsh" to detect the Tcl version, and
# silently gets an empty result -- breaking its install paths -- if missing)

# in-place install (PREFIX=.): $ANNOTSV/bin/AnnotSV is used directly from the
# cloned tree, no separate copy step -- this is the exact flow documented in
# share/doc/AnnotSV/quickstart.md. Deliberately NOT running
# `make install-human-annotation` / `install-mouse-annotation` here.
#
# No python3/pip in this image on purpose: AnnotSV's `make install` also
# unconditionally runs an `install-variantconvert` sub-target (clones
# SamuelNicaise/variantconvert, needed for -vcf 1 native VCF output), which
# has no flag to skip it and no `pip3`/`pip`/`python -m pip` to try without
# Python installed -- its own Makefile already treats that failure as
# non-fatal (falls back through all three, then just removes its own
# "installed" flag file) so `make install` still succeeds, just without
# variantconvert. That's fine here: variantconvert's -vcf 1 output turned out
# to be broken anyway (header-only VCF, no records, in every mode tried), so
# the pipeline builds its own annotated VCF a different way instead --
# writing AnnotSV's TSV output back into the VCF as INFO fields directly, see
# modules/local/annotsv.nf (annotate_vcf_with_tsv) and
# bin/annotate_vcf_from_annotsv.py. The variantconvert source still gets
# cloned (harmless, a few MB) since that clone step itself doesn't depend on
# Python -- just cleaned up below to not carry it around uselessly.
RUN git clone --branch "v${ANNOTSV_VERSION}" --depth 1 \
        https://github.com/lgmgeo/AnnotSV.git /opt/AnnotSV \
    && cd /opt/AnnotSV \
    && make PREFIX=. install \
    && rm -rf share/python3/variantconvert

ENV ANNOTSV=/opt/AnnotSV
ENV PATH="${ANNOTSV}/bin:${PATH}"

WORKDIR /data
