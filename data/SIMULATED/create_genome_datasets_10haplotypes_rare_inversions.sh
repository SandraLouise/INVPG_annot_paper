#! /bin/bash

###################################################################
## Authors: Sandra Romain, Claire Lemaitre
##
## This script generates a 10-haplotypes dataset for inv-pg-annot experiments
## It takes as input a reference sequence (1rst haplotype), a rate of SNP divergence (parameter DIV), and a VCF with 100 inversion coordinates
## It creates 9 other haplotypes with each random samples of 50% of the SNPs and each haplotype has 1/9 of the inversion set, so that each initial inversion is present in exactly one haplotype.
## It formats the corresponding fasta files and input files as needed for the 4 pangenome graph pipelines
##
## Note: SNPs and INVs are sampled independently, ie. SNPs alleles inside the inversions are not fixed with the inversions
## Note: the number of haplotypes is not a parameter, it could be changed by changing a variable value inside the script BUT it won't work for the input cactus file that requires a phylogenetic tree of the haplotypes (here it is hard-coded for 10 haplotypes)
##
###################################################################



# Parameters
# ----------
REF=$1	# fasta file of ref genome
DIV=$2   #0.1 for 0,1% div ; 1 for 1%div
INV=$3   # VCF file with the INVs
REF_ID=$4   #name of ref individual (eg. chm13)
SIMU_ID=$5   # prefix name of simulated indiv (eg. tenhaplo) ; for fasta header and for pangenome graph input files
OUTDIR=$6   # dir to put the genomes
SEED=$7    # random seed for generating the snps

#ex : ./create_genome_datasets_10haplotypes.sh chr21.fa 0.1 chr21_sim_100inv.vcf chm13 tenhaplorare genomes_10hap 10 

# Fixed input & env
nHAP=9
snpPCT=50

. envpython-3.11.9.sh 
BINDIR="../data/SIMULATED" 

OUTDIR_VCF="${OUTDIR}_VCF"
mkdir -p $OUTDIR
mkdir -p $OUTDIR_VCF

allINV=$SIMU_ID"_all_inv.txt"
rm -f $allINV


# 0. Prepare Output files = input config files for the pangenome graph pipelines
# ------------------------------------------------------------------------------

#    0.a. MGC pipeline
#    ----------------
MGC_file="mgc_input_${SIMU_ID}_div$DIV.txt"
REF_PATH=`realpath $REF`
echo -e "$REF_ID\t$REF_PATH" > $MGC_file

#    0.b. PGGB and minigraph pipelines
#    --------------------------------
PGGB_file="pggb_input_${SIMU_ID}_div$DIV.txt"
echo $REF_PATH > $PGGB_file

#    0.c. CACTUS
#    ----------
#For cactus: need the FASTA files to be formatted specifically: no '#' in the header, so storing reformatted files in a specific directory

OUTDIR_CACTUS=$OUTDIR/"cactus_genomes"
mkdir -p $OUTDIR_CACTUS
REF_NAME=`basename $REF `
echo ">$REF_ID" > $OUTDIR_CACTUS/$REF_NAME
grep -v "^>" $REF >> $OUTDIR_CACTUS/$REF_NAME

CACTUS_file="cactus_input_${SIMU_ID}_div$DIV.txt"
CACTUS_REF_PATH=`realpath $OUTDIR_CACTUS/$REF_NAME`
echo "($REF_ID:1.0,((((${SIMU_ID}1:1.0,${SIMU_ID}2:1.0):1.0,${SIMU_ID}3:1.0):1.0,(${SIMU_ID}4:1.0,${SIMU_ID}5:1.0):1.0):1.0,((${SIMU_ID}6:1.0,${SIMU_ID}7:1.0):1.0,(${SIMU_ID}8:1.0,${SIMU_ID}9:1.0):1.0)):1.0);"  > $CACTUS_file 
echo -e "$REF_ID\t$CACTUS_REF_PATH" >> $CACTUS_file



# For inversions, we want that each haplotype has a different subset of the full inversion set
# First we need to randomize the order of the INVs
seedFile="toerase_seed.txt"    # we use a seedfile, to be reproducible
echo "forrandom $SEED" > $seedFile
echo "blabla" >> $seedFile  # I dont know why but on the cluster only one line in the seed file does not work...
randomINVFILE=$INV".random"
grep -v "^#" $INV | sort -R --random-source=$seedFile > $randomINVFILE

NBTOT=$(cat $randomINVFILE  | wc -l)
NB_PER_IND=$(($NBTOT / $nHAP)) 

# 1. Generate a SNP list with 2x$DIV    - If DIV != 0
# ----------------------------------

if [ "$DIV" != "0" ]; then
divSNP=$(echo "scale=2; $DIV * 2" | bc)
if [[ $divSNP == .* ]]; then
    divSNP="0$divSNP"
fi
outsnp=$SIMU_ID"_div"$divSNP"_snp_only.fa"
python3 $BINDIR/02_generate_synthetic_chrom_with_snps.py $REF $outsnp $divSNP $SEED

# We can discard the file $outsnp ; but we keep the output VCF file $outsnp".vcf"
rm -f $outsnp
SNP=$outsnp".vcf"
fi


# 2. Simulate 9 haplotypes with 50% SNPs and a specific subset of INVs
# --------------------------------------------------
cumStart=0
for i in $(seq 1 $nHAP); do

# set indiv name + seed
prefix=$SIMU_ID$i
indivSEED=$(($SEED +$i))

#      2.a input 50% of SNPs
#      ---------------------
if [ "$DIV" != "0" ]; then
outsnp=$prefix"_div"$DIV"_snp_only.fa"
python3 $BINDIR/02_generate_synthetic_chrom_with_snps_from_VCF.py $REF $SNP $outsnp $snpPCT $indivSEED

#note: snp list recorded in file $outsnp".vcf"
mv $outsnp".vcf" $OUTDIR_VCF/
fi

#      2.b input 50% of INVs
#      ----------------------
outfinal=$prefix"_div"$DIV"_"$invPCT"INV"
header=$prefix"#0#chr21"

#Creating a vcf file with N inversions ($NB_PER_IND, or N=11 for 100 initial inversions) for this individual : a slice of the shuffled initial inversion file
#Each individual will get a different slice, so that each inversion is present in exactly one haplotype
indivINV="temp_inv.vcf"
tempStart=$(($cumStart + $NB_PER_IND))
if [ $i == $nHAP ]; then
nb=$(($NBTOT-$cumStart))
tail -n $nb $randomINVFILE > $indivINV
else
head -n $tempStart $randomINVFILE | tail -n $NB_PER_IND > $indivINV
fi
cumStart=$tempStart


if [ "$DIV" != "0" ]; then
python3 $BINDIR/03a_simulate_inversions_in_1_synthetic_chrom.py $outsnp $indivINV $outfinal 100 $header $indivSEED
# 2 files: $outfinal.fa $outfinal.vcf

# cleaning
rm -f $outsnp
else
python3 $BINDIR/03a_simulate_inversions_in_1_synthetic_chrom.py $REF $indivINV $outfinal 100 $header $indivSEED
# 2 files: $outfinal.fa $outfinal.vcf
fi
# Concatenate the inversion file of all simulated haplotypes for checking at the end that each inv is present in at least one indiv
grep -v "^#" $outfinal.vcf | cut -f 1-3 >> $allINV

mv $outfinal.vcf $OUTDIR_VCF/
rm -f $indivINV

#      2.c Organize haplotype file + write file name in pg tool input files
#      --------------------------------------------------------------------
SIMU=$outfinal".fa"
mv $SIMU $OUTDIR/$SIMU

SIMU_PATH=`realpath $OUTDIR/$SIMU`
# MGC
echo -e "$prefix\t$SIMU_PATH" >> $MGC_file
# PGGB
echo $SIMU_PATH >> $PGGB_file

# Cactus
SIMU_NAME=`basename $OUTDIR/$SIMU`
echo ">$prefix" > $OUTDIR_CACTUS/$SIMU_NAME
grep -v "^>" $OUTDIR/$SIMU >> $OUTDIR_CACTUS/$SIMU_NAME
CACTUS_SIMU_PATH=`realpath $OUTDIR_CACTUS/$SIMU_NAME`
echo -e "$prefix\t$CACTUS_SIMU_PATH" >> $CACTUS_file

done


# Minigraph, same file as pggb
MINIGRAPH_file="minigraph_input_${SIMU_ID}_div$DIV.txt"
cp $PGGB_file $MINIGRAPH_file


rm -f $randomINVFILE
rm -f $seedFile

# 3. Checking inversion repartition
# ---------------------------------

sort $allINV | uniq -c > $SIMU_ID"_div"$DIV"_inv_distrib.txt"
#rm -f $allINV
cat $SIMU_ID"_div"$DIV"_inv_distrib.txt" | tr -s ' ' | cut -d ' ' -f2 | sort | uniq -c  > $SIMU_ID"_div"$DIV"_inv_histo.txt"

nb=$(wc -l $SIMU_ID"_div"$DIV"_inv_distrib.txt" | cut -d ' ' -f1)
echo "Simulation of $nHAP haplotypes: done"
echo "Genome files are in dir $OUTDIR/"
echo "Input files for PG pipelines are: $MGC_file , $PGGB_file , $CACTUS_file and $MINIGRAPH_file"
echo "n=$nb INV were simulated in at least one haplotype (see the repartition in file ${SIMU_ID}_div${DIV}_inv_distrib.txt and the corresponding histogram in file ${SIMU_ID}_div${DIV}_inv_histo.txt)"



