#! /bin/python

import sys
import random
from datetime import date

refFASTA, nbINV, outVCF, seed = sys.argv[1:]

# Setting seed for reproducibility
random.seed(int(seed))

#===========================================================
# Get ref chrom length from fasta
#===========================================================
refLEN = 0

with open(refFASTA, "r") as file:
    for line in file:

        if line[0] == ">":
            refIDbed = line.rstrip()[1:]
            refID = refIDbed.split("#")[0]  # restricting the name in the vcf before the first #, to ease the manipulation in R
        
        else:
            refLEN += len(line.rstrip())

#===========================================================
# Simulate random SVs
#===========================================================
numINV = int(nbINV)
minINV, maxINV = 50, 1000000     # INV size range : 50 bp - 1 Mb

buffer = 50000      # preserving regions of 50 Kb at both extremities of the chromosome
flank = 1000   # minimal distance between two simulated inv


# defining 4 size ranges: 50 bp-1Kb ; 1-5 Kb ; 5Kb - 100 Kb; > 100 Kb
bins = [1000, 5000, 100000]
# Nb of INV in each bins ; trying equal repartition
max_bins = [0]
max_bins.extend([int(numINV / (len(bins) + 1))] * len(bins))

max_bins[0] = numINV - sum(max_bins)

#storing the current number of INV in each bin
n_bins = [0] * (len(bins) + 1)

# print(max_bins)
# print(n_bins)

dictINV = {}   # key = INV position ; value = INV size

while len(dictINV.keys()) < numINV:

    #--------------------------------------------------------
    # Choose size 
    #--------------------------------------------------------
    r_invLEN = random.randint(minINV, maxINV)

    if r_invLEN < bins[0]:
        i_bin = 0
    elif r_invLEN < bins[1]:
        i_bin = 1
    elif r_invLEN < bins[2]:
        i_bin = 2
    else:
        i_bin = 3
    
    while n_bins[i_bin] == max_bins[i_bin]:
    
        r_invLEN = random.randint(minINV, maxINV)

        if r_invLEN < bins[0]:
            i_bin = 0
        elif r_invLEN < bins[1]:
            i_bin = 1
        elif r_invLEN < bins[2]:
            i_bin = 2
        else:
            i_bin = 3

    #--------------------------------------------------------
    # Choose position 
    #--------------------------------------------------------
    r_invPOS = random.randint(buffer, refLEN-buffer-r_invLEN)

    overlap = False

    for p, l in dictINV.items():

        if any([p - flank <= r_invPOS <= p+l-1 + flank, 
                p - flank <= r_invPOS+r_invLEN-1 <= p+l-1 + flank, 
                r_invPOS <= p < p+l-1 <= r_invPOS+r_invLEN-1]):
            
            overlap = True
    
    n_try = 0
    while overlap and n_try < 20:

        r_invPOS = random.randint(buffer, refLEN-buffer)

        overlap = False

        for p, l in dictINV.items():

            if any([p - flank <= r_invPOS <= p+l-1 + flank, 
                    p -flank <= r_invPOS+r_invLEN-1 <= p+l-1 + flank, 
                    r_invPOS <= p < p+l-1 <= r_invPOS+r_invLEN-1]):
                
                overlap = True
        
        n_try += 1

    
    if not overlap:
        dictINV[r_invPOS] = r_invLEN
        n_bins[i_bin] += 1

        # print(r_invPOS, r_invLEN, n_bins)

#===========================================================
# Write as VCF
#===========================================================
def format_vcf_line(ref, pos, typ, leng, tid):

    info = ";".join([
        f"END={str(pos+leng-1)}",
        f"SVTYPE={typ}",
        f"SVSIZE={str(leng)}"
    ])

    vcfLINE = "\t".join([
        ref,
        str(pos),
        f"{typ}{str(tid)}",
        "N",
        "N",
        ".",
        ".",
        info
    ])

    return vcfLINE

#-----------------------------------------------------------
# Prepare header
#-----------------------------------------------------------
headerLINES = [
    '##fileformat=VCFv4.3',
    f'##fileDate={date.today()}',
    '##ALT=<ID=INV,Description="Inversion">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position on reference genome">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV">',
    '##INFO=<ID=SVSIZE,Number=1,Type=Integer,Description="Size in nt of the SV">',
    f'##contig=<ID={refID},length={str(refLEN)}>',
    '#' + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
]


#-----------------------------------------------------------
# Write VCF	 + bed file
#-----------------------------------------------------------

outBED=outVCF+".bed"

with open(outVCF, "w") as vcf, open(outBED, "w") as bed:
     vcf.write("\n".join(headerLINES)+"\n")
     tid = 0

     for pos, leng in sorted(dictINV.items()):

        tid += 1
        vcf.write(format_vcf_line(refID, pos, "INV", leng, tid)+"\n")
        bed.write("\t".join([refIDbed, str(pos), str(pos+leng-1)])+"\n")

print("size repartition 50 bp-1Kb ; 1-5 Kb ; 5Kb - 100 Kb; > 100 Kb")
print(n_bins)
