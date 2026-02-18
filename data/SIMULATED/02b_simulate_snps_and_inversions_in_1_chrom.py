#! /bin/python3


"""=========================================================
      Idea here is to select a % of SNPs and INV from a VCF file 
      to simulate in an output haplotype sequence
      BUT : SNPs inside inversions are selected if and only if the corresponding inversion is selected
      NOTE : here, the level of div inside inversions is twice the level outside
========================================================="""

#note: another possibility is to select 50% of SNPs inside inversions as fixed and the remaining unfixed. To obtain that, modify lines 133 and 149 : instead of putting the SNP in the dicoInv2snp randomly choose between dicoInv2snp and outside with proba 0.5.

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

def import_inv_vcf(file):
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

def import_snp_vcf(file):
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

inFASTA, inSNP, inINV, outFA_name, pctINV, outFA_header, allFixed, seed = sys.argv[1:]

#allFixed should be set to 0 to get half the SNPs inside inversions fixed ; or 1 for all SNPs fixed
allFixed=int(allFixed)

random.seed(int(seed))

inSEQ = import_fasta(inFASTA)
dSNP, dsnpVCF = import_snp_vcf(inSNP)
dINV, dinvVCF = import_inv_vcf(inINV)

expectedLEN = len(inSEQ)


#===========================================================
# 0. Prepare two SNP lists : inside and outside INVs
#===========================================================

list_inv=sorted(dINV.keys())
iInv = 0
currStart = list_inv[iInv]
currEnd = dINV[list_inv[iInv]]
lastEnd = dINV[list_inv[-1]]
dicoInv2snp = {}

#init dico
for posinv in list_inv:
	dicoInv2snp[posinv] = []
outside = []
for pos in sorted(dSNP.keys()):
	if pos < currStart or pos > lastEnd:
		outside.append(pos)
	else:
		if pos == currStart:
			continue
			# on ne fait rien
		elif pos < currEnd:
			isFixed = 1
			if not allFixed:
				isFixed = random.randint(0,1)
			if isFixed:
				dicoInv2snp[list_inv[iInv]].append(pos)
			else:
				outside.append(pos)
		else:
			keep = True
			if pos == currEnd:  # on update l'inversion mais on ne répertorie nulle part ce snp
				keep = False
			# change currStart and currEnd ; may need to skip inversions
			iInv += 1
			currStart = list_inv[iInv]
			currEnd = dINV[list_inv[iInv]]
			while pos >= currEnd and iInv < len(list_inv)-1:
				iInv += 1
				currStart = list_inv[iInv]
				currEnd = dINV[list_inv[iInv]]
			if pos < currStart and keep :
				outside.append(pos)
			elif pos > currStart and pos < currEnd:
				isFixed = 1
				if not allFixed:
					isFixed = random.randint(0,1)
				if isFixed:
					dicoInv2snp[list_inv[iInv]].append(pos)
				else:
					outside.append(pos)
			elif pos > lastEnd:
				outside.append(pos)

#===========================================================
# 1. Select INV to add
#===========================================================

full_inv_list = list(dINV.keys())
pctINV = int(pctINV)

n_inv = min(len(full_inv_list), int(len(full_inv_list) * (pctINV/100)))

selected_inv_list = sorted(random.sample(full_inv_list, n_inv))


#===========================================================
# 2. Select SNPs to add : all the SNPs inside inv + 50% of the SNPs outside INV
#===========================================================
# other possibilities (not implemented) : 
#  - to get the same diversity level inside/outside inv : count number of fixed ones + sample the remaining randomly...
#  - to get a mix of fixed and unfixed SNPs inside inv : in step 0, do not select all snps inside inv but half of them... : modify lines 133 and 149 : instead of putting the SNP in the dicoInv2snp randomly choose between dicoInv2snp and outside with proba 0.5.

# Snps in selected inv
snps_in_invs = []
for i in selected_inv_list:
	snps_in_invs+=dicoInv2snp[i]

#SNPs outside invs
pctSNP = 50
n_snp = min(len(outside), int(len(outside) * (pctSNP/100)))
snps_outside = random.sample(outside, n_snp)

#print(len(snps_in_invs))
#print(len(snps_outside))

selected_snp_list = sorted(snps_in_invs + snps_outside)

#===========================================================
# 3. Add SNPs in sequence
#===========================================================

list_char = list(inSEQ)

for pos in selected_snp_list:
    list_char[pos] = dSNP[pos]
    
mutSEQ = "".join(list_char)

if len(inSEQ) != len(mutSEQ):
    print("wrong length for mutSEQ")
    exit


#===========================================================
# 4. Add INV in sequence
#===========================================================

inv_starts = sorted(selected_inv_list)
seq_wINV = mutSEQ[0:inv_starts[0]]

for i in range(len(inv_starts)):
    start = inv_starts[i]
    end = dINV[start]

    if i < len(inv_starts) - 1:
        nextStart = inv_starts[i+1]
    else:
        nextStart = len(mutSEQ)

    #-------------------------------------------------------
    # Add inverted seq between start and end
    #-------------------------------------------------------
    seq_wINV = seq_wINV + get_complementary(mutSEQ[start:end])

    #-------------------------------------------------------
    # Add following seq until next inv or seq end
    #-------------------------------------------------------
    seq_wINV = seq_wINV + mutSEQ[end:nextStart]

if len(mutSEQ) != len(seq_wINV):
    print("wrong length for seq_wINV")
    exit

#===========================================================
# 5. Print SEQ with SVs + VCFs
#===========================================================

# If only INVs
outSEQ = seq_wINV

# Output fasta
outFA = f"{outFA_name}.fa"
with open(outFA, "w") as fasta:
    fasta.write(f">{outFA_header}\n")
    fasta.write(outSEQ)

# Output associated VCF for inversions
outVCF = f"{outFA_name}.vcf"
with open(outVCF, "w") as vcf:
    for inv in inv_starts:
        vcf.write(dinvVCF[inv])

# Output associated VCF for SNPs
outVCF2 = f"{outFA_name}_snps.vcf"
with open(outVCF2, "w") as vcf:
	for snp in selected_snp_list:
		vcf.write(dsnpVCF[snp])