# Scripts to compute inversion recall and statistics for different values of the min-cov parameter

2 scripts :

- `run_annot_eval_mincov.sh`: runs the annotation and evaluation for a given tool, expe and a given min-cov value.
- `expe_robustness_mincov.sh`: calls the previous script for a list of tools and a list of min-cov values

### Main output
 The file `[OUTPUT_DIR]/evaluation/[EXPE_ID].eval`

## Script `run_annot_eval_mincov.sh`

It assumes the graph and bubbles have already been detected, and are stored and organized as if produced by the automation script (see git: [../automation/README.md](../automation/README.md).

#### Input
 - `PATH_INPUT`: absolute path of the main output dir of the automation script that generated the graphs
 - `EXPE_ID`: eg. "2hap_div0.1"
 - `TOOL`: name of the graph pangenome tool 
 - `PATH_TRUTH`: absolute path to the inversion truth file
 - `PATH_OUTPUT`: new dir for these results ; should be different from   
 - `PATH_INPUT`: (otherwise merged stat files will be overwritten)
 - `MINCOV`: value of the in-cov parameter for inv-pg-annot

#### Command

Example of command:

```
run_annot_eval_mincov.sh ../paper_PG_INV/results/chm13_chr21 2hap_div0.1 cactus ../paper_PG_INV/results/chm13_chr21/chr21_sim_100inv.bed ../invpg-annot/chm13_chr21_mincov 0.1
```

## Script `expe_robustness_mincov.sh`

#### Input
 - `PATH_INPUT`: absolute path of the main output dir of the automation script that generated the graphs
 - `EXPE_ID`: eg. "2hap_div0.1"
 - `PATH_TRUTH`: absolute path to the inversion truth file
 - `PATH_OUTPUT`: new dir for these results ; should be different from   
 - `PATH_INPUT`: (otherwise merged stat files will be overwritten)

The list of min-cov values are hard-coded in the script.

#### Command

Genocluster dir: `/scratch/clemaitr/invpg-annot/chm13_chr21_mincov/`

```
expe_robustness_mincov.sh ../paper_PG_INV/results/chm13_chr21 2hap_div0.1 ../paper_PG_INV/results/chm13_chr21/chr21_sim_100inv.bed ../invpg-annot/chm13_chr21_mincov cactus minigraph mgc pggb_v0.7.4
```
