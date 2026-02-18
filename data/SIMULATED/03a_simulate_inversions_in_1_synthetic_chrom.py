#! /bin/python3

def get_complementary(seq):
    
    d = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}

    c_seq = ""
    for char in reversed(seq):
        c_seq = c_seq + d[char]
    
    return c_seq

def import_fasta(file):
    seq = ""

    with open(file, "r") as f:
        for l in f:
            if l[0] != ">":
                seq = seq + l.rstrip()
    
    return seq

def import_vcf(file):
    d_inv = {}
    d_vcf = {}

    with open(file, "r") as f:

        prev_inv_s, prev_inv_e = 0, 0

        for l in f:
            if l[0] != "#":

                parsed_l = l.rstrip().split("\t")
                pos = int(parsed_l[1])
                end = int(parsed_l[7].split("END=")[1].split(";")[0])

                #-------------------------------------------------------
                # Check if overlaps previous SV
                #-------------------------------------------------------

                prev_s = prev_inv_s
                prev_e = prev_inv_e

                if prev_s < pos < prev_e:
                    continue
                
                #-------------------------------------------------------
                # Add INV to dict  
                #-------------------------------------------------------              
                if "INV" in parsed_l[2]:
                    d_inv[pos] = end
                    d_vcf[pos] = l
                    prev_inv_s, prev_inv_e = pos, end

    return d_inv, d_vcf

"""=========================================================
      /!\ Works only when working with ONE chromosome.
(the fasta file given must contain only the chr of interest) 
                    + VCF must be sorted
========================================================="""

import sys
import random

inFASTA, inVCF, outFA_name, pctINV, outFA_header, seed = sys.argv[1:]

random.seed(int(seed))

inSEQ = import_fasta(inFASTA)
dINV, dVCF = import_vcf(inVCF)

expectedLEN = len(inSEQ)

#===========================================================
# 1. Select INV to add
#===========================================================

full_inv_list = list(dINV.keys())
pctINV = int(pctINV)

n_inv = min(len(full_inv_list), int(len(full_inv_list) * (pctINV/100)))

selected_inv_list = sorted(random.sample(full_inv_list, n_inv))

#===========================================================
# 2. Add INV in sequence
#===========================================================

inv_starts = sorted(selected_inv_list)
seq_wINV = inSEQ[0:inv_starts[0]]

for i in range(len(inv_starts)):
    start = inv_starts[i]
    end = dINV[start]

    if i < len(inv_starts) - 1:
        nextStart = inv_starts[i+1]
    else:
        nextStart = len(inSEQ)

    #-------------------------------------------------------
    # Add inverted seq between start and end
    #-------------------------------------------------------
    seq_wINV = seq_wINV + get_complementary(inSEQ[start:end])

    #-------------------------------------------------------
    # Add following seq until next inv or seq end
    #-------------------------------------------------------
    seq_wINV = seq_wINV + inSEQ[end:nextStart]

if len(inSEQ) != len(seq_wINV):
    print("wrong length for seq_wINV")
    exit

#===========================================================
# 3. Print SEQ with SVs
#===========================================================

# If only INVs
outSEQ = seq_wINV

# Output fasta
outFA = f"{outFA_name}.fa"
with open(outFA, "w") as fasta:
    fasta.write(f">{outFA_header}\n")
    fasta.write(f"{outSEQ}\n")

# Output associated VCF
outVCF = f"{outFA_name}.vcf"
with open(outVCF, "w") as vcf:
    for inv in inv_starts:
        vcf.write(dVCF[inv])
