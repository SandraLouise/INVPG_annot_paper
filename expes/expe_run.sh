#! /bin/bash

###################################################################
## Authors: Sandra Romain, Claire Lemaitre
##
## manages the automation, launches the pipeline (with SBATCH) for a given dataset and with the graph building tools given by the user.
##
###################################################################

##################################################
# Environment and parameters
##################################################

## Input arguments
DIR_INPUT_FILES=$1
EXPE_ID=$2
PATH_OUTPUT=$3
PATH_TRUTH=$4
shift 4
declare -a TOOLS=($@)

## Script
SCRIPT_PIPELINE="run_pipeline.sh"

ENV_BEDTOOLS=". envbedtools-2.27.1.sh"

##################################################
# Run pipeline for each tool given
##################################################

for tool in "${TOOLS[@]}"; do

    input_file="${DIR_INPUT_FILES}/${tool}_input_${EXPE_ID}.txt"

    sbatch --job-name=${tool}_${EXPE_ID} $SCRIPT_PIPELINE $input_file $PATH_OUTPUT $PATH_TRUTH $EXPE_ID $tool

done