#! /bin/bash

###################################################################
## Authors: Sandra Romain
###################################################################

# Input
refFA=$1
invVCF=$2

# Fixed input & env
nHAP=8
invPCT=50
simSCRIPT='../data/sim_data/scripts/00_simulate_inv_in_fasta.py'

# Simulate haplotypes
for i in $(seq 1 $nHAP); do

    outSeqID="simulated_hap_50pctINV_${i}"
    outFile="simulated_hap_50pctINV_${i}.fa"

    touch $outFile
    python3 $simSCRIPT $refFA $invVCF $outSeqID $invPCT

done
