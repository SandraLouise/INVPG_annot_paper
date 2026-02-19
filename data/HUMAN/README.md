# Human dataset

> [!NOTE]\
> Simply cloning this repo does not create the directory structure. You may need to adapt the command line instructions to your own data to proceed.


## Inversions set

In directory [inversion_sets](inversion_sets)

- Source: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20210917_SSEQplusWHintegrativePhasing_inversionCallset/variants_freeze4inv_sv_inv_hg38_processed_arbigent_filtered_manualDotplot_filtered_PAVgenAdded_withInvCategs.tsv
- Reference: Porubsky, D., Höps, W., Ashraf, H., Hsieh, P., Rodriguez-Martin, B., Yilmaz, F., ... & Korbel, J. O. (2022). Recurrent inversion polymorphisms in humans associate with genetic instability and genomic disorders. Cell, 185(11), 1986-2005.

### Files

- `variants_freeze4inv_sv_inv_hg38_processed_arbigent_filtered_manualDotplot_filtered_PAVgenAdded_withInvCategs.tsv`: original file
- `HGSVC2_inv_hg38.chr7.NA19240_HG00733_HG03486_HG02818.tsv`: processed (filtered) file for chr 7
- `HGSVC2_hg38_chr7_inv.bed`: bed file with filtered coordinates of inversions on chr 7
- `HGSVC2_inv_hg38.chrX.NA19240_HG00733_HG03486_HG02818.tsv`: processed (filtered) file for chr X
- `HGSVC2_hg38_chrX_inv.bed`: bed file with filtered coordinates of inversions on chr X

### Processing

Remove unbalanced inversion calls, extract inversions on one chromosome with genotype in NA19240, HG00733, HG03486, HG02818.

Example for chr 7:

```bash
head -n 1 variants_freeze4inv_sv_inv_hg38_processed_arbigent_filtered_manualDotplot_filtered_PAVgenAdded_withInvCategs.tsv | awk -v OFS='\t' '{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$24,$37,$42,$53,$59,$60,$61,$62,$63 }' > HGSVC2_inv_hg38.chr7.NA19240_HG00733_HG03486_HG02818.tsv
grep -w '^chr7' variants_freeze4inv_sv_inv_hg38_processed_arbigent_filtered_manualDotplot_filtered_PAVgenAdded_withInvCategs.tsv | grep 'pass' | grep -w 'inv' | awk -v OFS='\t' '{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$24,$37,$42,$53,$59,$60,$61,$62,$63 }' >> HGSVC2_inv_hg38.chr7.NA19240_HG00733_HG03486_HG02818.tsv
```

## Genome set

- Source: [HPRC](https://data.humanpangenome.org/assemblies)
- Reference: GRCh38

### Processing

#### Step A - extract multifasta for selected genomes

Script `00_extract_multifasta_from_agc.py`.

Extract from the HPRC `.agc` file (obtainable with `curl -o HPRC-yr1.agc https://zenodo.org/record/5826274/files/HPRC-yr1.agc?download=1`) with added [GRCh38](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC_PLUS/GRCh38/assemblies/).

#### Step B - assign contigs to a selected chromosome 

1. List ref chromosomes to ignore for chrom-contig association (extra contigs)

```bash
grep '^>' ../../haplotypes_full/GRCh38.fa | awk '{ print $1 }' | sed -r 's/>//g' | grep -F '_' > ragtag_exclude_extra.txt
echo 'chrEBV' >> ragtag_exclude_extra.txt
```

2. Run RagTag[^1] for each query haplotype

Script: `01_ragtag_chrom_contig_assoc.sh`
Uses second script `02_run_ragtag_scaffold.sh` to launch each ragtag run.

3. Extract chr contigs for each haplotype

Script: `03_extract_select_contigs.sh`
Lists each haplotypes contigs conrresponding to a given chromosome, then extracts them using samtools faidx.

## References

[^1]: Alonge, M., Lebeigle, L., Kirsche, M., Jenike, K., Ou, S., Aganezov, S., ... & Soyk, S. (2022). Automated assembly scaffolding using RagTag elevates a new tomato system for high-throughput genome editing. Genome biology, 23(1), 258.
