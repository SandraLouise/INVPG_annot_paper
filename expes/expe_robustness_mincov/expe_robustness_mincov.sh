#! /bin/bash

###################################################################
## Authors: Sandra Romain, Claire Lemaitre
##
## calls the previous script for a list of tools and a list of min-cov values
##
###################################################################

##################################################
# Environment and parameters
##################################################

## Input arguments
PATH_INPUT=$1
EXPE_ID=$2
PATH_TRUTH=$3 
PATH_OUTPUT=$4
shift 4
declare -a TOOLS=($@)

## Script
SCRIPT_MINCOV="run_annot_eval_mincov.sh"


##################################################
# Run pipeline for each tool given
##################################################

for tool in "${TOOLS[@]}"; do

	for mincov in 0.1 0.25 0.5 0.75 0.9 0.95; do
	    sbatch --job-name=${tool}_${EXPE_ID}_${mincov} $SCRIPT_MINCOV $PATH_INPUT $EXPE_ID $tool $PATH_TRUTH $PATH_OUTPUT $mincov
	done

done