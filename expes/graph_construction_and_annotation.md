# Commands and scripts to reproduce the pangenome graphs and the bubble annotation

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
	- for PGGB, the headers of fasta entries must respect the [panSN spec](https://github.com/pangenome/PanSN-spec).  Example: `CARC#0#chr6` or `chm13#0#chr21`. (species name + "#" + haplotype phase + "#" + scaffold/chromosome name)
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

> May 2025:

| tool | version | env |
| --- | --- | --- |
| minigraph | `0.21-r606` | `minigraph/`| 
| mgc | `cactus_2_8_2` or `7.0.0` (with `--version`) | `cactus_2_8_2/cactus/cactus_env/`| 
| cactus | `6.0.0` (with `--version`)  | `cactus/cactus_env/bin/activate` (and `cactus-bin-v2.6.7/venv_cactus-v2.6.7/` for `hal2vg`) | 
| pggb | `0.6.0` (April 2024) | `pggb`| 





## 2. Bubble detection

### VG-deconstruct for Cactus, Minigraph-Cactus, and PGGB graphs

> vg v1.61.0

Note: `vg-deconstruct` can detect path-explicit inversions consisting of a single node. It is no longer necessary to use another script for these bubbles (for history, the script and examples are still explained at the end of the Readme).

Env conda : `vg1.61.0`

```bash
vg deconstruct -p REFPATH -a GRAPH.gfa > BUBBLES.vcf
```

- `REFPATH`: Genomes path names in the graph can be checked with the following command : 

   ```bash
   vg paths -Mx GRAPH.gfa
   ```
   
   > Note: that looking at the name in second field of P lines is not sufficient (it can be different) : `grep '^P' GRAPH.gfa | cut -c1-50` 
   
  Examples of `REFPATH` for several graphs (coeno simulations div 0.1%):  
  
  ```
  MGC: CARC#0#chr6
  PGGB: CARC#0#chr6#0
  Cactus: CARC#scaffold_6
  ```
    - MGC and PGGB were built with the same `.fa` files (with ref header (panSN): "CARC#0#chr6"), it seems that both re-use this header but PGGB adds "#0" 
    - Cactus: combination of the name put in `cactus_input.txt` ("CARC") and the header of sequence ("scaffold_6") in `.fa`

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

## 4. Comparison of annotated and true inversions, recall

TODO with automated script.