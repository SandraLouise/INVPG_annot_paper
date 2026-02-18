#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --constraint avx2

# Parameters
# ---------
INPUT_FILE=$1
FINAL_GFA=$2

###################################################################
## Authors: Sandra Romain, Claire Lemaitre
##
## This script builds a pangenome graph from assemblies in separated files
## + assemblies must be in .fasta/.fa format
## + one file per specimen, can contain multiple headers
##
###################################################################

#============================================================
# Define work environnment
#============================================================
. /local/env/envconda.sh
ENV_CACTUS="cactus/cactus_env/bin/activate"
ENV_HALVG="cactus-bin-v2.6.7/venv_cactus-v2.6.7/bin/activate"
#Note: using this old env, because contains hal2vg command  (warning old version of vg : 1.50, but here only for vg convert)

# prepare working directory + copying the input file in it  (this dir will be deleted at the end)
RESDIR="TODELETE_Cactus_$(date +%s)"
mkdir -p $RESDIR
inputFile="pipeline.txt"
cp $INPUT_FILE $RESDIR/$inputFile


# Place into suitable directory
cd $RESDIR

#Getting the reference name : first field in the second line of input file
refGenome=$(sed -n '2p' $inputFile | cut -f1)

# Run Progressive Cactus
source $ENV_CACTUS
cactus --binariesMode singularity ./js $inputFile cactus.hal
deactivate

# Convert HAL to GFA
source $ENV_HALVG
hal2vg --noAncestors --refGenomes $refGenome cactus.hal > cactus.pg
vg convert --gfa-out -W cactus.pg > cactus.gfa
deactivate

cd ..
mv $RESDIR/cactus.gfa $FINAL_GFA

rm -rf $RESDIR
