#! /bin/bash

###################################################################
## Authors: Sandra Romain
##
## Extract chr7 contigs for each haplotype
##
###################################################################

#============================================================
# Define work environnment
#============================================================
ENV_SAMTOOLS="envsamtools-1.15.sh"

#------------------------------------------------------------
# Human data
#------------------------------------------------------------
SELECT_CHR=$1
FULLHAP_DIR=$2
RAGTAG_DIR=$3
CHR_DIR=$4

#============================================================

. $ENV_SAMTOOLS

mkdir $CHR_DIR
cd $RAGTAG_DIR

for hap_dir in $(ls -d */); do

    hap_id=${hap_dir::-1}
    hap_contigs=${hap_id}_${SELECT_CHR}_contigs.txt

    grep "^${SELECT_CHR}_RagTag" ${hap_id}/ragtag.scaffold.agp | grep 'W' | awk '{ print $6 }' > $hap_contigs

    samtools faidx ${FULLHAP_DIR}/${hap_id}.fa -r $hap_contigs > ${CHR_DIR}/${hap_id}_${SELECT_CHR}.fa

done
