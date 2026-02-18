# Coenonympha dataset

## Inversions set

### Files

- ``: original file
- ``: processed (filtered) file

### Processing

## Genome set

### Files

- `haplotypes_full/`: contains full diploid genomes of 4 species (1 per species) + reference genome BST1.0.
- `haplotypes_scaffold_*/`: contains contigs associated to chromosome * of BST1 (CARC), LAU8 (CGAR), FTN7 (CDAR), ORR1 (CCEP); one file per haplotype.

### Processing

#### Step A - assign contigs to chromosomes

Done in same fashion as human data, using RagTag then `03_extract_chr_contigs_COENO.sh`.

Need to add reference chromosome (for each <CHR> = {scaffold_6, scaffold_10, scaffold_14, scaffold_21, scaffold_26}) with:

```bash
. envsamtools-1.15.sh
cd ../data/coeno_data/
samtools faidx haplotypes_full/CARC.0.fa <CHR> | sed -r 's/>/>CARC.0./g' > haplotypes_scaffold_6/CARC.0_<CHR>.fa
```

#### Step B - graph construction

1. Build mgc graph

Build the pipeline file for minigraph-cactus (tab-separated .txt file). Uses the script `04_cactus_file.py`.

```bash
ENV_PANCAT="pancat"
conda activate $ENV_PANCAT

# Run this once for each chromosome
python scripts/04_cactus_file.py haplotypes_<CHR> ../../results/pangenomes/coeno_pg/mgc/<CHR>/pipeline/pipeline.txt

conda deactivate
```

Then, constuct the pangenome graph with `05_build_mgc_graph_COENO.sh`:

```bash
sbatch scripts/05_build_mgc_graph_COENO.sh <CHR>
```

2. Build pggb graph

First compute mash distance for each chromosome (using mash: https://github.com/marbl/Mash).

```bash
cd ../data/coeno_data
mkdir distances
mash triangle ../whole_genome_haplotypes/data/chrom_6_8haps.fa.gz > distances/scaffold_6.mash_triangle.txt
```

Get maximum distance with `compute_max_divergence.sh`.

```bash
./scripts/compute_max_divergence.sh
cat distances/coeno.divergence.txt
```

```text
scaffold_10	0.0709514
scaffold_14	0.0579033
scaffold_21	0.0674584
scaffold_26	0.0593303
scaffold_6	0.0626325
```

Choose value for `-p` as <= 100 - max(dist) * 100. Here max(dist) = 0.0709514 so `-p` should be <= 92.90486.
We set `-p 90`.

Build the graph with `06_build_pggb_graph_COENO.sh`.

```bash
cd ../data/coeno_data
cat haplotypes_<CHR>/* | sed -r 's/>CARC.0./>CARC#0#/g' | sed -r 's/>BST1.1./>BST1#1#/g' | sed -r 's/>BST1.2./>BST1#2#/g' | sed -r 's/>LAU8.1./>LAU8#1#/g' | sed -r 's/>LAU8.2./>LAU8#2#/g' | sed -r 's/>FTN7.1./>FTN7#1#/g' | sed -r 's/>FTN7.2./>FTN7#2#/g' | sed -r 's/>ORR1.1./>ORR1#1#/g' | sed -r 's/>ORR1.2./>ORR1#2#/g' > haplotypes_<CHR>_merge/<CHR>.fa

sbatch scripts/06_build_pggb_graph_COENO.sh <CHR>

```