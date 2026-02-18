#! /bin/bash
#SBATCH --cpus-per-task=16

t=16
INPUT_FILE=$1
OUT_GFA=$2

###################################################################
## Authors: Sandra Romain, Claire Lemaitre
##
## This script builds a pangenome graph from assemblies in separated files
## + INPUT_FILE must be the path to a file containing all the haplotype file paths (each on each line)
## + Reference genome must be on the first line
## + Works whatever the number of genomes/haplotypes
##
###################################################################

GENOME_LIST=$(tr '\n' ' ' < $INPUT_FILE)

#ENV
conda activate minigraph

minigraph -cxggs -t$t $GENOME_LIST > $OUT_GFA

conda deactivate
