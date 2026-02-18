#! /bin/bash

###################################################################
## Authors: Sandra Romain
##
## Run RagTag[^1] for each query haplotype
## [^1]: Alonge, M., Lebeigle, L., Kirsche, M., Jenike, K., Ou, S., Aganezov, S., ... & Soyk, S. (2022). Automated assembly scaffolding using RagTag elevates a new tomato system for high-throughput genome editing. Genome biology, 23(1), 258.
###################################################################

conda activate utils

ref=$1 #reference genome (fasta)
query=$2 #query genome (fasta)
exclude=$3 #txt file containing ref chromosomes to exclude, one per line
t=$4 #number of threads
outputdir=$5

ragtag.py scaffold $ref $query -e $exclude -t $t -o $outputdir

conda deactivate
