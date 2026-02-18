#! /bin/bash

###################################################################
## Author: Sandra Romain
##
## This script computes all graph and inversion annotation statistics
## and recall needed to reproduce the tables and figures in the
## INVPG-annot paper.
##
###################################################################

# Input arguments
GRAPH_GFA=$1
TRUTH_BED=$2
ANNOT_BED=$3
ANNOT_STAT=$4
FOUND_MISSING=$5

# Loading env
ENV_BEDTOOLS=". envbedtools-2.27.1.sh"
ENV_PYTHON=". envpython-3.9.5.sh"
SCRIPT_REDUNDANCY="redundancy_stats.py"
SCRIPT_FOUND_MISSING="found_missing_inversions.py"

###################################################################

# Expe factors
expe="$(basename -- $GRAPH_GFA .gfa)"
IFS='_' read -r -a expe_array <<< "$expe"
TOOL="${expe_array[2]}"
DIV="${expe_array[1]}"
HAP_NB="${expe_array[0]}"

echo -e "PG\tDIV\tStat_name\tStat_value"

# Graph statistics
graph_size="$( grep '^S' $GRAPH_GFA | awk '{sum+=length($3)} END{print sum}' )"
node_number="$( grep -c '^S' $GRAPH_GFA )"
edge_number="$( grep -c '^L' $GRAPH_GFA )"

echo -e "${TOOL}\t${DIV}\tGraph_size\t${graph_size}"
echo -e "${TOOL}\t${DIV}\tNode_number\t${node_number}"
echo -e "${TOOL}\t${DIV}\tEdge_number\t${edge_number}"

# Annotation statistics
annot_stats=$( cat $ANNOT_STAT )
IFS=$'\n' read -d "\034" -r -a annot_stats_array <<<"${annot_stats}\034"
for stat in "${annot_stats_array[@]}"; do
    echo -e "${TOOL}\t${DIV}\t${stat}"
done

# Annotation evaluation
$ENV_BEDTOOLS

recall="$( bedtools intersect -a $TRUTH_BED -b $ANNOT_BED -f 0.5 -r -c | awk '{ if ($NF !=0) print }' | wc -l )"
falsepos="$( bedtools intersect -a $ANNOT_BED -b $TRUTH_BED -v | wc -l )"
echo -e "${TOOL}\t${DIV}\tRecall\t${recall}"
echo -e "${TOOL}\t${DIV}\tFalse_positives\t${falsepos}"

ALL_INTERSECT="${ANNOT_BED}.all_intersect"
bedtools intersect -a $TRUTH_BED -b $ANNOT_BED -wao > $ALL_INTERSECT
$ENV_PYTHON

# Intersect output for figures
RECALL_INTERSECT="${ANNOT_BED}.recall_intersect"
bedtools intersect -a $TRUTH_BED -b $ANNOT_BED -wao -f 0.5 -r | awk -v OFS='\t' -v pg="$TOOL" -v div="$DIV" '{print pg,div,$2,$3,$5,$6,$7}' > $RECALL_INTERSECT

redundancy_stats=$( python3 $SCRIPT_REDUNDANCY $ALL_INTERSECT $RECALL_INTERSECT )
IFS=$'\n' read -d "\034" -r -a redundancy_stats_array <<<"${redundancy_stats}\034"
for stat in "${redundancy_stats_array[@]}"; do
    echo -e "${TOOL}\t${DIV}\t${stat}"
done

# Found VS missing inversions
python3 $SCRIPT_FOUND_MISSING $TRUTH_BED $ALL_INTERSECT $RECALL_INTERSECT > $FOUND_MISSING

# Remove intermediary files
rm $ALL_INTERSECT