#!/usr/bin/env bash
# Generate nextflow_schema with updated basecaller enumerations
#
# This script uses `nextflow config` to obtain the basecaller container,
# creates JSON arrays of the models using the container's list-models script
# and injects them with jq to create nextflow_schema.json.new.
set -euo pipefail
LC_ALL=C

TARGET=$1

get_artifactory "${WHALEFISH_ARTIFACTORY_METADATA_ROOT}/dorado/latest/models.json"

# Get basecaller_cfg with a valid Clair3 model
# CW-4166 Note we skip including khz models which are synonymous with the model name without the distinction ... for now
awk 'NR>1 && $1 == "dorado" && $3 != "-" {print $2}' data/clair3_models.tsv | grep -v khz_ | sort | uniq > valid_clair3_dorado_models.ls
awk 'NR>1 && $1 == "dorado" && $3 == "-" {print $2}' data/clair3_models.tsv | grep -v khz_ | sort | uniq > invalid_clair3_dorado_models.ls
awk 'NR>1 && $1 != "dorado" && $3 != "-" {print $2}' data/clair3_models.tsv | sort | uniq > valid_clair3_nondorado_models.ls

# Get basecaller_cfg with a valid Dorado model in models.json (should be sorted but lets be sure)
jq -r '.simplex[] | select(startswith("rna") | not)' models.json | sort > simplex_models.ls
comm -12 valid_clair3_dorado_models.ls simplex_models.ls > basecalling_models.ls # basecalling models WITH a clair3 model

cat valid_clair3_dorado_models.ls invalid_clair3_dorado_models.ls | sort > mentioned_dorado_models.ls
comm -13 mentioned_dorado_models.ls simplex_models.ls > unmentioned_dorado_models.ls
if [[ $(wc -l < unmentioned_dorado_models.ls) -ne 0 ]]; then
    echo "There are Dorado models in the container that are not mapped to a Clair3 model:"
    sed 's,^,* ,' unmentioned_dorado_models.ls
    exit 1
fi
comm -23 valid_clair3_dorado_models.ls simplex_models.ls | cat - valid_clair3_nondorado_models.ls | sort > clair3_only_models.ls

# Convert model lists to JSON arrays
SIMPLEX_MODELS=$(cat simplex_models.ls | sed '$a\custom' | cat - clair3_only_models.ls | jq -Rn '[inputs]')

jq \
    -j \
    --indent 4 \
    --argjson simplex_models "${SIMPLEX_MODELS}" \
    '(.definitions.advanced_options.properties.override_basecaller_cfg.enum) = $simplex_models' \
    ${TARGET}/nextflow_schema.json > ${TARGET}/nextflow_schema.json.new

echo "----------------------------------------------------------------------------------------"
if ! diff -b nextflow_schema.json nextflow_schema.json.new; then
    echo "----------------------------------------------------------------------------------------"
    echo "Model schema requires updating."
    echo "Inspect the diff before amending."
    echo "To amend nextflow_schema.json, copy the diff from the logs and run patch -p0 <DIFF>"
    exit 1
else
    echo "nextflow_schema.json matches available models for latest dorado container."
fi
