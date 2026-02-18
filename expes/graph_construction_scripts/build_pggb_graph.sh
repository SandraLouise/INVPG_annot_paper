#!/bin/bash
#SBATCH --job-name=PGGBsim
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

###################################################################
## Authors: Siegfried Dubois, Sandra Romain, Claire Lemaitre
##
## This script builds a pangenome graph from assemblies in a single file
## + assemblies must be in .fasta/.fa format
## + one file for all genomes, can contain multiple headers
##
###################################################################

# Parameters 
INPUT_FILE=$1 #input file containing all genome file paths (one per line)
FINAL_GFA=$2

# Loading env
ENV_PGGB="pggb"
ENV_SAMTOOLS=". envsamtools-1.15.sh"

# Creating output dir
OUTPUT="TODELETE_PGGB_$(date +%s)"
mkdir -p $OUTPUT

# PGGB requires all haplotypes to be combined into a single file + indexed with samtools
source $ENV_SAMTOOLS
MERGED_FA=$OUTPUT/merged_haplotypes.fa      #file with all merged haplotypes
COUNT=0   # number of genomes
> $MERGED_FA  # Ensure file starts empty 

while IFS= read -r line; do
  if [[ -f "$line" ]]; then
    cat "$line" >> $MERGED_FA
    ((COUNT++))
  else
    echo "Warning: File not found - $line" >&2
  fi
done < $INPUT_FILE

echo "Concatenated $file_count files."

samtools faidx $MERGED_FA
conda deactivate


# Creating the variation graph
conda activate $ENV_PGGB
pggb -i $MERGED_FA -o $OUTPUT -n $COUNT -t 8 -p 90 -s 5k 

conda deactivate

# rename final gfa 
mv $MERGED_FA.*.smooth.final.gfa $FINAL_GFA

rm -rf $OUTPUT
