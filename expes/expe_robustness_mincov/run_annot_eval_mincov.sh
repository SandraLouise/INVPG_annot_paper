#! /bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G

###################################################################
## Authors: Sandra Romain, Claire Lemaitre
##
## runs the annotation and evaluation for a given tool, expe and a given min-cov value.
##
###################################################################

## Input arguments
PATH_INPUT=$1
EXPE_ID=$2
TOOL=$3
PATH_TRUTH=$4
PATH_OUTPUT=$5  # should be different from PATH_INPUT (otherwise merged stat files will be overwritten)
MINCOV=$6


## Environments
. /local/env/envconda.sh
ENV_INVPGANNOT="invpgannot"
ENV_BEDTOOLS=". envbedtools-2.27.1.sh"

## Scripts
SCRIPT_EVALUTATION="../eval_and_stats/graph_annot_statistics.sh"


## Input files
INPUT_ID="${EXPE_ID}_${TOOL}"

GRAPH="${PATH_INPUT}/graphs/${INPUT_ID}.gfa"
BUBBLES="${PATH_INPUT}/bubbles/${INPUT_ID}.vcf"


## Output files
OUTPUT_ID="${EXPE_ID}_${TOOL}_${MINCOV}"

OUTPUT_ANNOT="${PATH_OUTPUT}/inversions/${OUTPUT_ID}.bed"
OUTPUT_STATS="${PATH_OUTPUT}/inversions/${OUTPUT_ID}.stats"
OUTPUT_EVAL="${PATH_OUTPUT}/evaluation/${OUTPUT_ID}.eval"
OUTPUT_MERGED_EVAL="${PATH_OUTPUT}/evaluation/${EXPE_ID}.eval"
OUTPUT_MERGED_INTERSECT="${PATH_OUTPUT}/evaluation/${EXPE_ID}.intersect"

mkdir -p ${PATH_OUTPUT}
mkdir -p ${PATH_OUTPUT}/inversions
mkdir -p ${PATH_OUTPUT}/evaluation
touch $OUTPUT_MERGED_EVAL
touch $OUTPUT_MERGED_INTERSECT


#--------------------------------------------#
# 3. Annotate bubbles
#--------------------------------------------#
conda activate $ENV_INVPGANNOT

invpg -v $BUBBLES -g $GRAPH -o $OUTPUT_ANNOT -m $MINCOV -d 10

conda deactivate

#--------------------------------------------#
# 4. Evaluate annotation
#--------------------------------------------#

$SCRIPT_EVALUTATION $GRAPH $PATH_TRUTH $OUTPUT_ANNOT $OUTPUT_STATS $OUTPUT_ID > $OUTPUT_EVAL
grep -v '^PG' $OUTPUT_EVAL | awk -v OFS='\t' -v mincov="$MINCOV" '{ print $1,$2,mincov,$3,$4 }' >> $OUTPUT_MERGED_EVAL

