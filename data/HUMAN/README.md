# Human dataset

> [!NOTE]\
> Simply cloning this repo does not create the directory structure. You may need to adapt the command line instructions to your own data to proceed.

- **Data files in:** `../data/human_data`
- **Script files in:** `../data/human_data/scripts`

## Inversions set

- Source: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20210917_SSEQplusWHintegrativePhasing_inversionCallset/variants_freeze4inv_sv_inv_hg38_processed_arbigent_filtered_manualDotplot_filtered_PAVgenAdded_withInvCategs.tsv
- Reference: Porubsky, D., Höps, W., Ashraf, H., Hsieh, P., Rodriguez-Martin, B., Yilmaz, F., ... & Korbel, J. O. (2022). Recurrent inversion polymorphisms in humans associate with genetic instability and genomic disorders. Cell, 185(11), 1986-2005.

### Files

- `variants_freeze4inv_sv_inv_hg38_processed_arbigent_filtered_manualDotplot_filtered_PAVgenAdded_withInvCategs.tsv`: original file
- `HGSVC2_inv_hg38.chr7.NA19240_HG00733_HG03486_HG02818.tsv`: processed (filtered) file

### Processing

Remove unbalanced inversion calls, extract inversions on chromosome 7 with genotype in NA19240, HG00733, HG03486, HG02818.

```bash
head -n 1 variants_freeze4inv_sv_inv_hg38_processed_arbigent_filtered_manualDotplot_filtered_PAVgenAdded_withInvCategs.tsv | awk -v OFS='\t' '{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$24,$37,$42,$53,$59,$60,$61,$62,$63 }' > HGSVC2_inv_hg38.chr7.NA19240_HG00733_HG03486_HG02818.tsv
grep -w '^chr7' variants_freeze4inv_sv_inv_hg38_processed_arbigent_filtered_manualDotplot_filtered_PAVgenAdded_withInvCategs.tsv | grep 'pass' | grep -w 'inv' | awk -v OFS='\t' '{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$24,$37,$42,$53,$59,$60,$61,$62,$63 }' >> HGSVC2_inv_hg38.chr7.NA19240_HG00733_HG03486_HG02818.tsv
```

## Genome set

- Source: HPRC
- Reference:

### Files

- `haplotypes_full/`: contains full diploid genome of GRCh38, NA19240, HG00733, HG03486, HG02818; one file per haplotype.
- `haplotypes_chr7/`: contains chromosome 7 of GRCh38, NA19240, HG00733, HG03486, HG02818; one file per haplotype.
- `haplotypes_chr7/intermediate/`: intermediate files for chromosome-contig association (RagTag output)

### Processing

#### Step A - extract multifasta for selected genomes

Script `00_extract_multifasta_from_agc.py`.
Extract from the HPRC `.agc` file (obtainable with `curl -o HPRC-yr1.agc https://zenodo.org/record/5826274/files/HPRC-yr1.agc?download=1`) with added [GRCh38](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC_PLUS/GRCh38/assemblies/).

#### Step B - assign contigs to chromosome 7

1. List ref chromosomes to ignore for chrom-contig association (extra contigs)

```bash
grep '^>' ../../haplotypes_full/GRCh38.fa | awk '{ print $1 }' | sed -r 's/>//g' | grep -F '_' > ragtag_exclude_extra.txt
echo 'chrEBV' >> ragtag_exclude_extra.txt
```

2. Run RagTag[^1] for each query haplotype

Script: `01_ragtag_chrom_contig_assoc.sh`
Uses second script `02_run_ragtag_scaffold.sh` to launch each ragtag run on the cluster.

3. Extract chr7 contigs for each haplotype

Script: `03_extract_chr7_contigs.sh`
Lists each haplotypes contigs conrresponding to chr7, then extracts them using samtools faidx.

## References

[^1]: Alonge, M., Lebeigle, L., Kirsche, M., Jenike, K., Ou, S., Aganezov, S., ... & Soyk, S. (2022). Automated assembly scaffolding using RagTag elevates a new tomato system for high-throughput genome editing. Genome biology, 23(1), 258.