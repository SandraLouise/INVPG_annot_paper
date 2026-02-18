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
SEQS=$OUTPATH"/sequences/"

# Loading env
. /local/env/envconda.sh
ENV_CACTUS="apptainer run ../cactus_v2.9.9.sif"
ENV_VG="vg1.61.0"

# Creating dirs
mkdir -p $OUTPATH $TOIL $TEMP $SEQS

# Copy input file, because mgc re-writes it
INPIPELINE=$OUTPATH/pipeline.txt

# Create a pipeline file with relative paths and copy files to output
WD="$(pwd)"
REPATH=$(cat <<END
import subprocess,os,pathlib

with open("$INPUT_FILE",'r',encoding='utf-8') as creader:
    with open("$INPIPELINE",'w',encoding='utf-8') as cwriter:
        for line in creader:
            os.system("cp" + " " + line.strip().split('\t')[1] + " " + "$SEQS")
            cwriter.write(line.split('\t')[0]+"\t"+str(os.path.relpath("$SEQS"+str(pathlib.Path(line.strip().split('\t')[1]).name), "$WD")+"\n"))
END
)

FILE="$(python3 -c "$REPATH")"

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
$ENV_CACTUS cactus-minigraph $JB $INPIPELINE $OUT.gfa --reference $NAME_REF #--configFile cactus_config.xml
$ENV_CACTUS cactus-graphmap $JB $INPIPELINE $OUT.gfa $OUT.paf  --reference $NAME_REF --outputFasta $OUT.sv.gfa.fa.gz #--configFile cactus_config.xml
$ENV_CACTUS cactus-align $JB $INPIPELINE $OUT.paf $OUT.hal --pangenome --outGFA --outVG --reference $NAME_REF --workDir $OUTPATH"/CACTUS/workdir" #--configFile cactus_config.xml
$ENV_CACTUS cactus-graphmap-join $JB --vg $OUT.vg --outDir $OUTPATH"/CACTUS/outdir" --outName "final" --reference $NAME_REF --clip $CLIP --filter $FILTER #--configFile cactus_config.xml
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