#! /bin/bash

###################################################################
## Authors: Sandra Romain
##
## Run RagTag[^1] for each query haplotype
## [^1]: Alonge, M., Lebeigle, L., Kirsche, M., Jenike, K., Ou, S., Aganezov, S., ... & Soyk, S. (2022). Automated assembly scaffolding using RagTag elevates a new tomato system for high-throughput genome editing. Genome biology, 23(1), 258.
###################################################################

REFERENCE="../NGSdatasets/HumanGenomeSVConsortium/GRCh38_chrom_as_vcf.fna"
HAPDIR=$1
OUTDIR=$2

SCRIPT_RAGTAG="../data/HUMAN/02_run_ragtag_scaffold.sh"

# Filter chromosomes to use from reference genomes
cd $OUTDIR
chrom_to_exclude="ragtag_exclude_extra.txt"
#grep '^>' $REFERENCE | awk '{ print $1 }' | sed -r 's/>//g' | grep -F '_' > $chrom_to_exclude
#echo 'chrEBV' >> $chrom_to_exclude
echo 'chrmitochondrion' > $chrom_to_exclude

# Run Ragtag on each given haplotype
for path_to_hap in $(ls ${HAPDIR}/*.*.fa); do

    arr_path=(${path_to_hap//"/"/ })
    hap_file=${arr_path[6]}
    arr_file=(${hap_file//".f"/ })
    hap_id=${arr_file[0]}

    sbatch --cpus-per-task=16 $SCRIPT_RAGTAG $REFERENCE $path_to_hap $chrom_to_exclude 16 $hap_id

done
