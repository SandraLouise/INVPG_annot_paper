#! /bin/python3



def import_fasta(file):
    seq = ""

    with open(file, "r") as f:
        for l in f:
            if l[0] != ">":
                seq = seq + l.rstrip()
    
    return seq

def import_vcf(file):
    d_snp = {}
    d_vcf = {}

    with open(file, "r") as f:

        for l in f:
            if l[0] != "#":

                parsed_l = l.rstrip().split("\t")
                pos = int(parsed_l[1])
                alt = parsed_l[4]
               
                #-------------------------------------------------------
                # Add snp to dict  
                #-------------------------------------------------------              
                d_snp[pos] = alt
                d_vcf[pos] = l

    return d_snp, d_vcf

"""=========================================================
      /!\ Works only when working with ONE chromosome.
(the fasta file given must contain only the chr of interest) 
                    + VCF must be sorted
========================================================="""

import sys
import random

inFASTA, inVCF, outFA_name, pctSNP, seed = sys.argv[1:]

random.seed(int(seed))

inSEQ = import_fasta(inFASTA)
dSNP, dVCF = import_vcf(inVCF)

expectedLEN = len(inSEQ)

#===========================================================
# 1. Select SNP to add
#===========================================================

full_snp_list = list(dSNP.keys())
pctSNP = int(pctSNP)

n_snp = min(len(full_snp_list), int(len(full_snp_list) * (pctSNP/100)))

selected_snp_list = sorted(random.sample(full_snp_list, n_snp))

#===========================================================
# 2. Add snp in sequence
#===========================================================

list_char = list(inSEQ)

for pos in selected_snp_list:
    list_char[pos] = dSNP[pos]
    
mutSEQ = "".join(list_char)

if len(inSEQ) != len(mutSEQ):
    print("wrong length for mutSEQ")
    exit

#===========================================================
# 3. Print SEQ with SNPs
#===========================================================


# Output fasta
with open(outFA_name, "w") as fasta:
    fasta.write(f">simulated\n")
    fasta.write(mutSEQ)

# Output associated VCF
outVCF = f"{outFA_name}.vcf"
with open(outVCF, "w") as vcf:
    for snp in selected_snp_list:
        vcf.write(dVCF[snp])