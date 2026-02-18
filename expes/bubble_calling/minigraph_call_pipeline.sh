#!/bin/bash
#SBATCH --job-name=gfatools_pipeline
#SBATCH --cpus-per-task=8

###################################################################
## Author: Sandra Romain
##
## This script produces a VCF describing the variant bubbles from a minigraph 
## pangenome.
## This pipeline uses gfatools bubble, minigraph call and scripts 
## from minigraph git repository.
##
## Arguments (in their respective order in the command):
##     1) path to output directory
##     2) full path to minigraph graph (rGFA)
##     3) full path to reference sequence (FASTA)
##     4...n) full path to haplotype(s) sequence (FASTA)
##
###################################################################

# Input arguments
OUTPUT_BUBBLES=$1
GRAPH=$2
INPUT_FILE=$3

declare -a ARR_HAPS=($(cat $INPUT_FILE))
REF="${ARR_HAPS[0]}"
HAP_COUNT=$(( ${#ARR_HAPS[@]} - 1 ))

# Environment and parameters
ENV="minigraph/"
MGUTILS="minigraph/misc/mgutils.js"
MGUTILS_ES6="minigraph/misc/mgutils-es6.js"
THREADS=8

# Output parameters
DIR_INTERMEDIARY=$(echo $OUTPUT_BUBBLES | sed "s/.vcf/_$(date +%s)/g")
mkdir $DIR_INTERMEDIARY
OUT_GFATOOLS="${DIR_INTERMEDIARY}/gfatools.bed"
OUT_PATH_REF="${DIR_INTERMEDIARY}/mgcall.0ref.bed"
OUT_MERGE_BED="${DIR_INTERMEDIARY}/merge.bed"

# Prepare environment
. /local/env/envconda.sh
conda activate $ENV

# List graph bubbles
gfatools bubble $GRAPH > $OUT_GFATOOLS

# Extract reference path in bubbles
minigraph -cxasm --call -t8 $GRAPH $REF > $OUT_PATH_REF

# Extract haplotype(s) path in bubbles
for i in $(seq 1 $HAP_COUNT); do

    OUT_PATH_HAP="${DIR_INTERMEDIARY}/mgcall.hap${i}.bed"
    HAP="${ARR_HAPS[$i]}"

    minigraph -cxasm --call -t8 $GRAPH $HAP > $OUT_PATH_HAP

done

# Merge per-sample call
ls ${DIR_INTERMEDIARY}/mgcall.*.bed > ${DIR_INTERMEDIARY}/samples.txt
paste ${DIR_INTERMEDIARY}/mgcall.*.bed | $MGUTILS merge -s ${DIR_INTERMEDIARY}/samples.txt - > $OUT_MERGE_BED

# Convert to VCF
$MGUTILS_ES6 merge2vcf -r0 $OUT_MERGE_BED > $OUTPUT_BUBBLES

# End
conda deactivate
rm -r $DIR_INTERMEDIARY