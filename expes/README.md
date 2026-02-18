# Commands and scripts to reproduce the pangenome graphs and the bubble annotation

> [!NOTE]\
> Simply cloning this repo does not create the directory structure. You may need to adapt the command line instructions to your own data to proceed.

For a given set of experiments, we create a directory `expe/` (outside of the directory where we stored or simulated the input haplotype sequences) with the following structure:

```
expeXX/
└── graphs/
└── bubbles/
└── annotations/
└── evaluations/
```


## 1. Graph construction

### Requirement

- the input haplotypes should be in individual fasta files (not necessarily in the same directory, not necessarily in a directory of their own.
- 2 (uncompatible) types of fasta are necessary:
	- for PGGB, the headers of fasta entries must respect the [panSN spec](https://github.com/pangenome/PanSN-spec).  Example: `chm13#0#chr21`. (species name + "#" + haplotype phase + "#" + scaffold/chromosome name)
	- for Cactus, the headers must not contain any "#" characters
	- for MGC and minigraph: we use the PGGB format

### Construction scripts

For each pangenome graph builder (minigraph, mgc, cactus and PGGB), we provide a single bash script in the [graph\_construction\_scripts directory](graph_construction_scripts/) named `build_[tool]_graph.sh` that takes 2 parameters as input:

 - `INPUT_FILE`: the path to a text file containing the absolute paths to the haplotype fasta files and optionally other information necessary for the corresponding tool.
 - `FINAL_GFA`: the output gfa file path.

```
sbatch build_[tool]_graph.sh <INPUT_FILE> <FINAL_GFA>
```

#### Notes: 

- all intermediary files and directories are cleaned after the run, the only output is the final GFA file in correct GFA format (GFA1.0 format with all P-lines).
- works with any number of haplotypes


### Versions of the tools:

| tool | version |
| --- | --- | --- |
| minigraph | `0.21-r606` |
| mc | `cactus_2_9_9` |
| cactus | `6.0.0` (with `--version`)  |
| pggb | `0.7.4` |


## 2. Bubble detection

### VG-deconstruct for Cactus, Minigraph-Cactus, and PGGB graphs

Note: `vg-deconstruct` can detect path-explicit inversions consisting of a single node. It is no longer necessary to use another script for these bubbles.

Env conda : `vg1.65.0`

```bash
vg deconstruct -p REFPATH -a GRAPH.gfa > BUBBLES.vcf
```

- `REFPATH`: Genomes path names in the graph can be checked with the following command : 

   ```bash
   vg paths -Mx GRAPH.gfa
   ```
   
   > Note: that looking at the name in second field of P lines is not sufficient (it can be different) : `grep '^P' GRAPH.gfa | cut -c1-50` 
   
### Minigraph

> gfatools v.0.4-r214-dirty

```bash
sbatch gfatools_pipeline_SBATCH.sh <path_to_output_dir> <path_to_rGFA> <path_to_ref_seq> <path_to_hap1_seq> ... <path_to_hapn_seq>
```

**Utiliser les chemins absolus des fichiers d'input.**


## 3. Bubble annotation

> [INVPG-annot](https://github.com/SandraLouise/INVPG_annot)

For all graph types:

```bash
invpg -v BUBBLES.vcf -g GRAPH.gfa -d DIV_PERCENTAGE
```

- using DIV_PERCENTAGE = 10.