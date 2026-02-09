#! /bin/bash

MYPARAMS=$1; # Parameters file path
ENV=$2; # Virtual Environment name or source path

EXPERIMENT=$(grep 'experiment_name' ${MYPARAMS} | awk '{print $2}' | tr -d "\"'")
echo "Experiment name: ${EXPERIMENT}";
EXPERIMENTFOLDER="${PWD}/experiments/${EXPERIMENT}";

BASEURL="https://raw.githubusercontent.com/lauromoraes/microbiom/main/nb-templates";

# Define pipeline steps
STEPS=(
    "step-prepare-data"
    "step-quality-control"
    "step-rarefaction-analysis"
    "step-metataxonomy"
    "step-diversity-analysis"
    "step-abundance-analysis"
    "step-lefse-analysis"
    "step-picrust2-analysis"
    );

STEPSDIR="${PWD}/nb-templates"
if ! [ -d "$STEPSDIR" ]; then
  echo "Creating directory for notebook templates: ${STEPSDIR}";
  mkdir -p ${STEPSDIR};
fi

if ! [ -d "$EXPERIMENTFOLDER" ]; then
  echo "Creating directory for experiment artifacts: ${EXPERIMENTFOLDER}";
  mkdir -p ${EXPERIMENTFOLDER};
  cp ${1} ${EXPERIMENTFOLDER}; # Copy pamameters YAML file to experiment
fi

EXECUTEDDIR="${EXPERIMENTFOLDER}/nb-executed-steps"
if ! [ -d "$EXECUTEDDIR" ]; then
  echo "Creating directory for executed notebooks: ${EXECUTEDDIR}";
  mkdir -p ${EXECUTEDDIR};
fi

# Download the utils.py file
if ! [ -f "${STEPSDIR}/utils.py" ]; then
  echo "Downloading: ${STEPSDIR}/utils.py";
  wget "${BASEURL}/utils.py" -O "${STEPSDIR}/utils.py";
  cp "${BASEURL}/utils.py" "${EXECUTEDDIR}/utils.py"
  cp "${BASEURL}/utils.py" "."
fi



echo "Processing parameters from: ${MYPARAMS}";

# Activate virtual environment with all dependencies
conda init bash
source ~/anaconda3/etc/profile.d/conda.sh;
conda activate ${ENV};
qiime dev refresh-cache;

# Execute each step
for i in "${!STEPS[@]}"; do
	echo "====== Executing Pipeline STEP $((i+1)): ${STEPS[i]} ======";

	# Define paths
	STEPFILE="${STEPSDIR}/${STEPS[i]}.ipynb";
	EXECUTEDFILE="${EXECUTEDDIR}/${STEPS[i]}-${EXPERIMENT}.ipynb";

	# Download notebook if it does not exist
	if ! [ -f "$STEPFILE" ]; then
		echo "... Downloading file: ${STEPFILE} ...";
		wget "${BASEURL}/${STEPS[i]}.ipynb" -O "${STEPSDIR}/${STEPS[i]}.ipynb";
	fi

	# Execute notebook
	echo ">>> Executing STEP file: ${STEPFILE} <<<";
	papermill "${STEPFILE}" "${EXECUTEDFILE}" -f "${MYPARAMS}";
done

# Deactivate the virtual environment
conda deactivate;
