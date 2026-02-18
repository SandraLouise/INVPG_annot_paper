#! /bin/python3

import sys
import random
import time
from datetime import date

in_fasta, out_fasta, pct, seed = sys.argv[1:]
pct = float(pct)

random.seed(int(seed))

#===========================================================
# 1. Import initial sequence
#===========================================================
iniSeq = ""
iniID = ""

with open(in_fasta, "r") as file:
    for line in file:

        if line[0] == ">":

            iniID = line[1:].rstrip()
        
        else:

            iniSeq = iniSeq + line.rstrip()

#===========================================================
# 2. Add Single Nucleotide Polymorphism
#===========================================================
a = ["A", "C", "T", "G"] 

def define_snp(iniN, a):

    mutations = [c for c in a if c != iniN]
    mutN = random.choice(mutations)

    return mutN

def meth_random_pos(iniSeq, pct, a):

    # Calculate nb of bases to modify
    #-----------------------------------------------------------
    n = int(len(iniSeq) * (pct / 100))

    # Pick SNPs positions
    #-----------------------------------------------------------
    pos = []

    for i in range(n):

        p = int(random.uniform(0, 1) * len(iniSeq))

        # Make sure that one same position is not picked twice
        while p in pos or iniSeq[p] == "N":
            p = random.randrange(len(iniSeq))
        
        pos.append(p)

    # Generate mutated sequence
    #-----------------------------------------------------------
    list_char = list(iniSeq)

    for p in pos:

        iniN = list_char[p]
        mutN = define_snp(iniN, a)

        list_char[p] = mutN

    mutSeq = "".join(list_char)

    return mutSeq, pos

def meth_pos_proba(iniSeq, pct, a):

    list_char = list(iniSeq)
    pos = []

    for i in range(len(iniSeq)):

        r = random.uniform(0, 1)

        if r < (pct/100):
            mutN = define_snp(list_char[i], a)
            list_char[i] = mutN
            pos.append(i)

    mutSeq = "".join(list_char)

    return mutSeq, pos

def meth_hybrid(iniSeq, pct, a):

    n = int(len(iniSeq) * (pct / 100))
    list_char = list(iniSeq)
    ref = ""
    pos = []
    snps = {}

    for i in range(len(iniSeq)):

        if len(pos) == n:
            break

        if list_char[i] == "N":
            continue

        r = random.uniform(0, 1)

        if r < (pct/100):
            mutN = define_snp(list_char[i], a)
            ref =list_char[i]
            list_char[i] = mutN
            snps[i] = [ref, mutN]
            pos.append(i)
    
    while len(pos) < n:

        p = int(random.uniform(0, 1) * len(iniSeq))
        if p not in pos and list_char[p] != "N":

            mutN = define_snp(list_char[p], a)
            ref=list_char[p]
            list_char[p] = mutN
            snps[p] = [ref, mutN]
            pos.append(p)

    mutSeq = "".join(list_char)

    return mutSeq, snps

mutSeq, snps = meth_hybrid(iniSeq, pct, a)

#===========================================================
# 3. Export mutated seq
#===========================================================
mutID = iniID + f"_{str(pct/100)}"

#Write mutated sequence in fasta file
with open(out_fasta, "w") as out:

    out.write(f">{mutID}\n")
    out.write(mutSeq + "\n")

#name to put in VCF
refName = iniID.split("#")[0]

with open(out_fasta+".vcf", "w") as out2:
    headerLINES = [
    '##fileformat=VCFv4.3',
    f'##fileDate={date.today()}',
    f'##contig=<ID={refName},length={len(iniSeq)}>',
    '#' + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
    ]

    out2.write("\n".join(headerLINES)+"\n")


    i = 0
    for position in sorted(snps.keys()):
        i += 1
        out2.write(f"{refName}\t{position}\tsnp{i}\t{snps[position][0]}\t{snps[position][1]}\t.\t.\t.\n")

