#!/usr/bin/env bash

set -Eeuo pipefail


# =============================================================================
# Arguments
# =============================================================================

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <parameters.yaml> <conda-environment>"
    echo
    echo "Example:"
    echo "  $0 params/params01.yaml rachis-qiime2-2026.7"
    exit 1
fi

MYPARAMS="$1"
ENV="$2"

if [[ ! -f "$MYPARAMS" ]]; then
    echo "Error: parameters file not found: $MYPARAMS" >&2
    exit 1
fi


# =============================================================================
# Project configuration
# =============================================================================

PROJECT_DIR="$(pwd)"

BASEURL="https://raw.githubusercontent.com/lauromoraes/microbiom/main/nb-templates"

STEPSDIR="${PROJECT_DIR}/nb-templates"


# =============================================================================
# Read experiment name
# =============================================================================

EXPERIMENT="$(
    awk '
        /^[[:space:]]*experiment_name[[:space:]]*:/ {
            sub(/^[^:]*:[[:space:]]*/, "")
            gsub(/["'\''"]/, "")
            print
            exit
        }
    ' "$MYPARAMS"
)"

if [[ -z "$EXPERIMENT" ]]; then
    echo "Error: experiment_name was not found in $MYPARAMS" >&2
    exit 1
fi


# =============================================================================
# Experiment directories
# =============================================================================

EXPERIMENTFOLDER="${PROJECT_DIR}/experiments/${EXPERIMENT}"
EXECUTEDDIR="${EXPERIMENTFOLDER}/nb-executed-steps"

mkdir -p \
    "$STEPSDIR" \
    "$EXPERIMENTFOLDER" \
    "$EXECUTEDDIR"


# =============================================================================
# Temporary environment
# =============================================================================
#
# PIPELINE_TEMP_DIR can be defined externally, for example:
#
#   export PIPELINE_TEMP_DIR=/mnt/data/tmp/microbiom
#
# If it is not defined, /tmp/microbiom is used.
#
# Using an external environment variable keeps machine-specific paths
# separate from the project source code.
# =============================================================================

TEMP_BASE="${PIPELINE_TEMP_DIR:-/tmp/microbiom}"
PIPELINE_TMP="${TEMP_BASE}/${EXPERIMENT}"

export TMPDIR="${PIPELINE_TMP}/qiime2"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

export JOBLIB_TEMP_FOLDER="${PIPELINE_TMP}/joblib"

export QIIME_CACHE="${PIPELINE_TMP}/qiime2-cache"

mkdir -p \
    "$TMPDIR" \
    "$JOBLIB_TEMP_FOLDER" \
    "$QIIME_CACHE"


# =============================================================================
# Temporary environment cleanup
# =============================================================================
#
# Temporary files are removed only when the pipeline completes successfully.
#
# If the pipeline fails, temporary files are preserved to facilitate
# debugging.
# =============================================================================

cleanup_temp_environment() {
    local exit_code=$?

    echo

    if [[ $exit_code -eq 0 ]]; then

        echo "Cleaning temporary environment..."

        if [[ -n "${PIPELINE_TMP:-}" && -d "$PIPELINE_TMP" ]]; then
            rm -rf -- "$PIPELINE_TMP"
            echo "Removed temporary directory:"
            echo "  $PIPELINE_TMP"
        fi

    else

        echo "Pipeline failed with exit code: $exit_code" >&2
        echo >&2
        echo "Temporary files were preserved for debugging:" >&2
        echo "  $PIPELINE_TMP" >&2

    fi

    return "$exit_code"
}

trap cleanup_temp_environment EXIT


# =============================================================================
# Pipeline steps
# =============================================================================

STEPS=(
    # "step-prepare-data"
    # "step-quality-control"
    # "step-rarefaction-analysis"
    # "step-metataxonomy"
    # "step-diversity-analysis"
    # "step-abundance-analysis"
    # "step-lefse-analysis"
    # "step-picrust2-analysis"
    # "step-report"
    "step-ancom"
)


# =============================================================================
# Pipeline information
# =============================================================================

echo
echo "============================================================"
echo "Microbiome pipeline"
echo "============================================================"
echo
echo "Experiment name : $EXPERIMENT"
echo "Project directory: $PROJECT_DIR"
echo "Parameters file : $MYPARAMS"
echo
echo "Temporary environment:"
echo "  TMPDIR             : $TMPDIR"
echo "  JOBLIB_TEMP_FOLDER : $JOBLIB_TEMP_FOLDER"
echo "  QIIME_CACHE        : $QIIME_CACHE"
echo


# =============================================================================
# Save experiment parameters
# =============================================================================

cp -f \
    "$MYPARAMS" \
    "${EXPERIMENTFOLDER}/$(basename "$MYPARAMS")"


# =============================================================================
# Download utils.py
# =============================================================================

UTILS_FILE="${STEPSDIR}/utils.py"

if [[ ! -f "$UTILS_FILE" ]]; then

    echo "Downloading utils.py..."

    wget \
        --quiet \
        --show-progress \
        "${BASEURL}/utils.py" \
        -O "$UTILS_FILE"

fi


# Keep copies of utils.py alongside the executed notebooks and project root.

cp -f \
    "$UTILS_FILE" \
    "${EXECUTEDDIR}/utils.py"

cp -f \
    "$UTILS_FILE" \
    "${PROJECT_DIR}/utils.py"


# =============================================================================
# Conda configuration
# =============================================================================

CONDA_BASE="${CONDA_BASE:-/opt/miniconda}"
CONDA_SH="${CONDA_BASE}/etc/profile.d/conda.sh"

if [[ ! -f "$CONDA_SH" ]]; then
    echo "Error: Conda initialization file not found:" >&2
    echo "  $CONDA_SH" >&2
    exit 1
fi

# shellcheck disable=SC1091
source "$CONDA_SH"


# =============================================================================
# Validate Conda environment
# =============================================================================

if ! conda env list \
    | awk '{print $1}' \
    | grep -Fxq "$ENV"
then
    echo "Error: Conda environment not found: $ENV" >&2
    echo >&2

    conda env list

    exit 1
fi


# =============================================================================
# Activate Conda environment
# =============================================================================

conda activate "$ENV"

echo "Conda environment:"
echo "  ${CONDA_DEFAULT_ENV:-unknown}"
echo
echo "Python executable:"
echo "  $(command -v python)"
echo
echo "QIIME executable:"
echo "  $(command -v qiime || true)"
echo


# =============================================================================
# Validate QIIME 2 environment
# =============================================================================

if ! command -v qiime >/dev/null 2>&1; then
    echo "Error: qiime was not found after activating environment '$ENV'." >&2
    exit 1
fi


python - <<'PY'
import sys
import qiime2

print(f"Python version : {sys.version.split()[0]}")
print(f"QIIME 2 version: {qiime2.__version__}")
PY


# =============================================================================
# Validate temporary directories
# =============================================================================

echo
echo "Temporary filesystem status:"
echo

df -h \
    "$TMPDIR" \
    "$JOBLIB_TEMP_FOLDER" \
    "$QIIME_CACHE"

echo


# =============================================================================
# Refresh QIIME 2 plugin cache
# =============================================================================

echo "Refreshing QIIME 2 plugin cache..."

qiime dev refresh-cache


# =============================================================================
# Validate Papermill
# =============================================================================

if python -c 'import papermill' >/dev/null 2>&1; then

    PAPERMILL_COMMAND=(
        python
        -m
        papermill
    )

else

    echo "Error: papermill is not installed in Conda environment '$ENV'." >&2
    echo >&2
    echo "Check with:" >&2
    echo "  conda run -n '$ENV' python -c 'import papermill'" >&2
    echo >&2
    echo "Install, if compatible with the environment:" >&2
    echo "  conda install -n '$ENV' -c conda-forge papermill ipykernel" >&2

    exit 1

fi


# =============================================================================
# Execute pipeline
# =============================================================================

TOTAL_STEPS="${#STEPS[@]}"

for i in "${!STEPS[@]}"; do

    STEP="${STEPS[$i]}"

    STEP_NUMBER=$((i + 1))

    echo
    echo "============================================================"
    echo "Pipeline step ${STEP_NUMBER}/${TOTAL_STEPS}: ${STEP}"
    echo "============================================================"
    echo

    STEPFILE="${STEPSDIR}/${STEP}.ipynb"

    EXECUTEDFILE="${EXECUTEDDIR}/${STEP}-${EXPERIMENT}.ipynb"


    # -------------------------------------------------------------------------
    # Download notebook template
    # -------------------------------------------------------------------------

    if [[ ! -f "$STEPFILE" ]]; then

        echo "Downloading notebook:"
        echo "  ${BASEURL}/${STEP}.ipynb"

        wget \
            --quiet \
            --show-progress \
            "${BASEURL}/${STEP}.ipynb" \
            -O "$STEPFILE"

    fi


    # -------------------------------------------------------------------------
    # Validate notebook
    # -------------------------------------------------------------------------

    if [[ ! -s "$STEPFILE" ]]; then
        echo "Error: notebook is empty or invalid:" >&2
        echo "  $STEPFILE" >&2
        exit 1
    fi


    # -------------------------------------------------------------------------
    # Execute notebook
    # -------------------------------------------------------------------------

    echo "Input notebook:"
    echo "  $STEPFILE"
    echo

    echo "Output notebook:"
    echo "  $EXECUTEDFILE"
    echo

    "${PAPERMILL_COMMAND[@]}" \
        "$STEPFILE" \
        "$EXECUTEDFILE" \
        --parameters_file "$MYPARAMS"

done


# =============================================================================
# Pipeline completed
# =============================================================================

echo
echo "============================================================"
echo "Pipeline completed successfully"
echo "============================================================"
echo
echo "Experiment:"
echo "  $EXPERIMENT"
echo
echo "Experiment directory:"
echo "  $EXPERIMENTFOLDER"
echo
echo "Executed notebooks:"
echo "  $EXECUTEDDIR"
echo
