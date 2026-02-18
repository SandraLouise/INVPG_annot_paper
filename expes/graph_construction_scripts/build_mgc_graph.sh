#!/bin/bash
#SBATCH --job-name=MGCsim
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --constraint avx2

###################################################################
## Authors: Siegfried Dubois, Sandra Romain, Claire Lemaitre
##
## This script builds a pangenome graph from assemblies in separated files
## + assemblies must be in .fasta/.fa format
## + one file per specimen, can contain multiple headers
##
###################################################################

INPUT_FILE=$1
# typically : mgc_input.txt
FINAL_GFA=$2

# Working dirs
OUTPATH="TODELETE_MGC_$(date +%s)"  # dir with all intermediary files of Minigraph-Cactus (will be deleted at the end of the script)
TOIL=$OUTPATH"/jobstore/"
TEMP=$OUTPATH"/tempfiles/"

# Loading env
. /local/env/envconda.sh
ENV_CACTUS="cactus_2_8_2/cactus/cactus_env/bin/activate"
ENV_VG="vg1.61.0"

source $ENV_CACTUS

# Creating dirs
mkdir -p $OUTPATH $TOIL $TEMP

# Copy input file, because mgc re-writes it
INPIPELINE=$OUTPATH/pipeline.txt
cp $INPUT_FILE $INPIPELINE

# VARS
CLIP=0
FILTER=0

# Getting reference name (first line in file before \t)
echo "$(head -n 1 $INPIPELINE)" | cut -d$'\t' -f1 > $TEMP"tempfile.txt"
NAME_REF=`cat $TEMP"tempfile.txt"`
JB=$TOIL".js_0"
[ -d $JB ] && rm -r $JB


mkdir $OUTPATH"/CACTUS"
OUT=$OUTPATH"/CACTUS/graph"
mkdir $OUTPATH"/CACTUS/workdir"
mkdir $OUTPATH"/CACTUS/outdir"
cactus-minigraph $JB $INPIPELINE $OUT.gfa --reference $NAME_REF --binariesMode singularity #--configFile cactus_config.xml
cactus-graphmap $JB $INPIPELINE $OUT.gfa $OUT.paf  --reference $NAME_REF --outputFasta $OUT.sv.gfa.fa.gz --binariesMode singularity #--configFile cactus_config.xml
cactus-align $JB $INPIPELINE $OUT.paf $OUT.hal --pangenome --outGFA --outVG --reference $NAME_REF --workDir $OUTPATH"/CACTUS/workdir" --binariesMode singularity #--configFile cactus_config.xml
cactus-graphmap-join $JB --vg $OUT.vg --outDir $OUTPATH"/CACTUS/outdir" --outName "final" --reference $NAME_REF --clip $CLIP --filter $FILTER --binariesMode singularity #--configFile cactus_config.xml
if [ -d $JB ]
then
    rm -r $JB
fi
GRAPH=$OUTPATH"/CACTUS/outdir/final.full.gfa"
gzip -d $GRAPH".gz"

conda deactivate

# Convert the graph in GFA1.1 to GFA1.0 format with all P-lines
conda activate $ENV_VG
vg convert -g -f -W $GRAPH > $FINAL_GFA

conda deactivate

rm -rf $OUTPATH

