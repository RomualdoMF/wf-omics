# ISV (https://github.com/tsladecek/isv_package) + ClassifyCNV
# (https://github.com/Genotek/ClassifyCNV) combined into one container for
# wf-omics's CNV classification ensemble (see bin/cnv_ensemble_classify.py).
#
# Combined into a single image rather than one-per-tool because: ClassifyCNV
# has no Python package dependencies at all (stdlib + bedtools only -- checked
# its source directly), nothing it needs conflicts with ISV's stack, and the
# two are always run together as a single logical "CNV ensemble" step in the
# script above -- see wf-omics.mdc, "Ensemble de classificação de CNV", for
# the fuller reasoning. (AnnotSV stays in its own separate image: it's used
# independently and bind-mounts its own ~20GB annotation directory, which has
# nothing to do with these two.)
#
# Published at romualdomf/cnv-ensemble (hardcoded in base.config,
# withLabel:cnv_ensemble -- same reasoning as withLabel:annotsv there: Docker
# Hub pull rate limits are keyed to whoever pulls, not to the publishing
# account, so there's no benefit to making the namespace a separate param).
#   docker build -t romualdomf/cnv-ensemble:1.0.0 docker/cnv_ensemble
#   docker push romualdomf/cnv-ensemble:1.0.0

FROM debian:bookworm-slim

ARG CLASSIFYCNV_VERSION=v1.1.1

RUN apt-get update && apt-get install -y --no-install-recommends \
        git curl ca-certificates bash coreutils \
        bedtools \
        python3 python3-pip \
        procps \
    && rm -rf /var/lib/apt/lists/*
# procps (ps) isn't a dependency of either tool -- Nextflow itself execs `ps`
# inside the container to collect task resource metrics, and errors out
# without it (bit us once already building the AnnotSV image).

# Debian bookworm's pip refuses plain "pip install" (PEP 668,
# externally-managed-environment) without this.
ENV PIP_BREAK_SYSTEM_PACKAGES=1

# isv's own pinned requirements (pandas<2, numpy<2, xgboost<2, an abandoned
# sklearn_json release, ...) are years out of date and pandas<2 doesn't even
# have a prebuilt wheel for this Python version -- but the actual code works
# fine against current versions. Install isv itself without its pins, then
# supply modern versions of everything it needs ourselves.
#
# NOTE: this doesn't fix isv's dependency `sklearn_json` importing a private
# sklearn symbol (`sklearn.ensemble._gb_losses`) that no longer exists in
# current scikit-learn -- that's patched at *runtime* by
# bin/cnv_ensemble_classify.py (_patch_sklearn_gb_losses), not here, so
# anything importing `isv` inside this container needs that shim applied
# first.
#
# (No pysam here: writing the ensemble's columns back into the CNV VCF as
# INFO fields is a separate script, bin/annotate_cnv_vcf.py, run in the
# pipeline's existing ontresearch/wf-common image instead -- that already
# has pysam installed, so there's no reason to duplicate it in this image.)
RUN pip install --no-cache-dir --no-deps isv==0.3.16 \
    && pip install --no-cache-dir \
        numpy pandas xgboost shap numba plotly sklearn-json scikit-learn

# ClassifyCNV ships its reference databases inside the repo itself -- no
# separate data download needed (unlike AnnotSV's human annotations).
RUN git clone --branch "${CLASSIFYCNV_VERSION}" --depth 1 \
        https://github.com/Genotek/ClassifyCNV.git /opt/ClassifyCNV

ENV CLASSIFYCNV_DIR=/opt/ClassifyCNV

WORKDIR /data
