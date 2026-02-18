#! /bin/bash

#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --constraint avx2

###################################################################
## Authors: Sandra Romain, Claire Lemaitre
##
## for each graph building tool, contains the pipeline
##
###################################################################

##################################################
# Environment and parameters
##################################################

## Input arguments
PATH_INPUT_FILE=$1
PATH_OUTPUT=$2
PATH_TRUTH=$3
EXPE_ID=$4
TOOL=$5

## Environments
ENV_VG="vg1.65.0"
ENV_INVPGANNOT="invpgannot"
ENV_BEDTOOLS=". envbedtools-2.27.1.sh"

## Scripts
SCRIPT_CACTUS_BUILD="../graph_construction_scripts/build_cactus_graph.sh"
SCRIPT_MINIGRAPH_BUILD="../graph_construction_scripts/build_minigraph_graph.sh"
SCRIPT_MGC_BUILD="../graph_construction_scripts/build_mgc_graph_v2.9.9.sh"
SCRIPT_PGGB_BUILD="../graph_construction_scripts/build_pggb_v0.7.4_graph.sh"
SCRIPT_MINIGRAPH_CALL="../bubble_calling/minigraph_call_pipeline.sh"
SCRIPT_EVALUTATION="../eval_and_stats/graph_annot_statistics.sh"

## Output parameters
OUTPUT_ID="${EXPE_ID}_${TOOL}"

OUTPUT_GRAPH="${PATH_OUTPUT}/graphs/${OUTPUT_ID}.gfa"
OUTPUT_BUBBLES="${PATH_OUTPUT}/bubbles/${OUTPUT_ID}.vcf"
OUTPUT_ANNOT="${PATH_OUTPUT}/inversions/${OUTPUT_ID}.bed"
OUTPUT_STATS="${PATH_OUTPUT}/inversions/${OUTPUT_ID}.stats"
OUTPUT_EVAL="${PATH_OUTPUT}/evaluation/${OUTPUT_ID}.eval"
FOUND_MISSING="${PATH_OUTPUT}/evaluation/${OUTPUT_ID}.found_missing.txt"

mkdir -p ${PATH_OUTPUT}/graphs
mkdir -p ${PATH_OUTPUT}/bubbles
mkdir -p ${PATH_OUTPUT}/inversions
mkdir -p ${PATH_OUTPUT}/evaluation

## Other
TRUTH_REF_ID=$(grep -v '^#' $PATH_TRUTH | head -n 1 | awk '{ print $1 }')

##################################################
# DIFFERING STEPS BETWEEN TOOLS
##################################################

if [ $TOOL = "minigraph" ]; then

    #--------------------------------------------#
    # 1a. Build MINIGRAPH graph
    #--------------------------------------------#
    $SCRIPT_MINIGRAPH_BUILD $PATH_INPUT_FILE $OUTPUT_GRAPH

    #--------------------------------------------#
    # 2a. Call bubbles in MINIGRAPH graph
    #--------------------------------------------#
    $SCRIPT_MINIGRAPH_CALL $OUTPUT_BUBBLES $OUTPUT_GRAPH $PATH_INPUT_FILE

    # rename reference in bubble VCF
    ref_id="$(grep -v '^#' $OUTPUT_BUBBLES | head -n 1 | awk '{ print $1 }')"
    sed -i "s/${ref_id}/${TRUTH_REF_ID}/g" $OUTPUT_BUBBLES

else
    #--------------------------------------------#
    # 1b. Build CACTUS graph
    #--------------------------------------------#
    if [ $TOOL = "cactus" ]; then
        $SCRIPT_CACTUS_BUILD $PATH_INPUT_FILE $OUTPUT_GRAPH

        # for bubble calling
        ref_path_tag="REFERENCE"
    
    #--------------------------------------------#
    # 1c. Build MINIGRAPH-CACTUS graph
    #--------------------------------------------#
    elif [ $TOOL = "mgc" ]; then
        $SCRIPT_MGC_BUILD $PATH_INPUT_FILE $OUTPUT_GRAPH

        # for bubble calling
        ref_path_tag="REFERENCE"
    
    #--------------------------------------------#
    # 1d. Build PGGB graph (v0.7.4)
    #--------------------------------------------#
    elif [ $TOOL = "pggb" ]; then
        $SCRIPT_PGGB_BUILD $PATH_INPUT_FILE $OUTPUT_GRAPH

        # for bubble calling
        ref_path_tag=$(grep '^>' $(head -n 1 $PATH_INPUT_FILE) | sed 's/>//g' | awk '{ split($1,header,"#"); print header[1] }')

    fi 

    #--------------------------------------------#
    # 2bcde. Call bubbles
    #--------------------------------------------#
    conda activate $ENV_VG

    ref_path=$(vg paths -Mx $OUTPUT_GRAPH | grep "${ref_path_tag}" | awk '{ print $1 }')
    vg deconstruct -p $ref_path -a $OUTPUT_GRAPH > $OUTPUT_BUBBLES

    # rename reference in bubble VCF in case different from truth file (for evaluation)
    vcf_ref_id=$(grep -v '^#' $OUTPUT_BUBBLES | head -n 1 | awk '{print $1}')
    sed -i "s/${vcf_ref_id}/${TRUTH_REF_ID}/g" $OUTPUT_BUBBLES

    conda deactivate

fi

##################################################
# COMMON STEPS BETWEEN TOOLS
##################################################

#--------------------------------------------#
# 3. Annotate bubbles
#--------------------------------------------#

conda activate $ENV_INVPGANNOT

invpg -v $OUTPUT_BUBBLES -g $OUTPUT_GRAPH -o $OUTPUT_ANNOT -m 0.5 -d 10

conda deactivate

#--------------------------------------------#
# 4. Evaluate annotation
#--------------------------------------------#

$SCRIPT_EVALUTATION $OUTPUT_GRAPH $PATH_TRUTH $OUTPUT_ANNOT $OUTPUT_STATS $FOUND_MISSING > $OUTPUT_EVAL

