#! /bin/bash

###################################################################
## Authors: Sandra Romain, Claire Lemaitre
##
## This script generates a 2-haplotypes dataset for inv-pg-annot experiments
## It takes as input a reference sequence (1rst haplotype), a rate of SNP divergence (parameter DIV), and a VCF with 100 inversion coordinates
## It creates the second haplotype with SNPs and the 100 inv, and formats the corresponding fasta files and input files as needed for the 4 pangenome graph pipelines
##
###################################################################

# Parameters
# ----------
REF=$1	# fasta file of ref genome
DIV=$2   #0.1 for 0,1% div ; 1 for 1%div
INV=$3   # VCF file with the INV
REF_ID=$4   #name of ref individual (eg. chm13)
SIMU_ID=$5   # name of simulated indiv (eg. 2hap) ; for fasta header and for pangenome graph input files
OUTDIR=$6   # dir to put the genomes
SEED=$7    # random seed for generating the snps

#ex : ./create_genome_datasets_2haplotypes.sh chr21.fa 0.1 chr21_sim_100inv.vcf chm13 2hap genomes_2hap 10 


. envpython-3.11.9.sh 

BINDIR="../data/SIMULATED" 


# 1. generate SNPs
# -----------------
prefix=$SIMU_ID
outsnp=$prefix"_div"$DIV"_snp_only.fa"
python3 $BINDIR/02_generate_synthetic_chrom_with_snps.py $REF $outsnp $DIV $SEED

#note: snp list recorded in file $outsnp".vcf"


#2. generate the haplotype with snps + 100 INV
# --------------------------------------------
PCTinv=100
outfinal=$prefix"_div"$DIV
header=$SIMU_ID"#0#chr21"
python3 $BINDIR/03a_simulate_inversions_in_1_synthetic_chrom.py $outsnp $INV $outfinal $PCTinv $header $SEED

# cleaning
rm -f $outsnp

# 2 files: $outfinal.fa $outfinal.vcf


#3. Organize, format the genome files + write the input files for running the different pangenome graph pipelines (MGC, PGGB, minigraph, Cactus)

SIMU=$outfinal".fa"

mkdir -p $OUTDIR
mv $SIMU $OUTDIR/$SIMU


# 3a. MGC pipeline
# ----------------
MGC_file="mgc_input_${SIMU_ID}_div$DIV.txt"
REF_PATH=`realpath $REF`
SIMU_PATH=`realpath $OUTDIR/$SIMU`
echo -e "$REF_ID\t$REF_PATH" > $MGC_file
echo -e "$SIMU_ID\t$SIMU_PATH" >> $MGC_file

# 3b. PGGB and minigraph pipelines
# --------------------------------
PGGB_file="pggb_input_${SIMU_ID}_div$DIV.txt"
echo $REF_PATH > $PGGB_file
echo $SIMU_PATH >> $PGGB_file

MINIGRAPH_file="minigraph_input_${SIMU_ID}_div$DIV.txt"
cp $PGGB_file $MINIGRAPH_file

# the same file can be used for minigraph

# 3c. CACTUS
# ----------
#For cactus: need the FASTA files to be formatted specifically: no '#' in the header, so storing reformatted files in a specific directory

OUTDIR_CACTUS=$OUTDIR/"cactus_genomes"
mkdir -p $OUTDIR_CACTUS
REF_NAME=`basename $REF `
echo ">$REF_ID" > $OUTDIR_CACTUS/$REF_NAME
grep -v "^>" $REF >> $OUTDIR_CACTUS/$REF_NAME

SIMU_NAME=`basename $OUTDIR/$SIMU`
echo ">$SIMU_ID" > $OUTDIR_CACTUS/$SIMU_NAME
grep -v "^>" $OUTDIR/$SIMU >> $OUTDIR_CACTUS/$SIMU_NAME

CACTUS_file="cactus_input_${SIMU_ID}_div$DIV.txt"
CACTUS_REF_PATH=`realpath $OUTDIR_CACTUS/$REF_NAME`
CACTUS_SIMU_PATH=`realpath $OUTDIR_CACTUS/$SIMU_NAME`
echo "($REF_ID:1.0,$SIMU_ID:1.0);"  > $CACTUS_file
echo -e "$REF_ID\t$CACTUS_REF_PATH" >> $CACTUS_file
echo -e "$SIMU_ID\t$CACTUS_SIMU_PATH" >> $CACTUS_file

