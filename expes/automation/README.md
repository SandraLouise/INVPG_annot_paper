# Automated pipeline run

> [!NOTE]\
> Simply cloning this repo does not create the directory structure. You may need to adapt the command line instructions to your own data to proceed.

The automated pipeline allows to reproduce the results of the paper. It encompasses the graph construction, bubble calling, inversion annotation and annotation evaluation.

## Main scripts

The automated pipeline is run by two main scripts:
- `expe_run.sh`: manages the automation, launches the pipeline (with SBATCH) for a given dataset and with the graph building tools given by the user.
    - git: [expes/automation/expe_run.sh](expes/automation/expe_run.sh)
- `run_pipeline.sh`: called by `expe_run.sh` for each graph building tool, contains the pipeline.
    - git: [expes/automation/run_pipeline.sh](expes/automation/run_pipeline.sh)

## Usage

### Input

- `dir_input_files`: absolute path to the directory containing the genome input files
    - **/!\ the path of the reference genome must be listed first in all input files**
- `expe_id`: identifier of the experience
    - **/!\ must match the filename of the input files** (_e.g._ "2hap_div0.1" for input files named "*tool*_input_2hap_div0.1.txt")
- `path_output`: absolute path to the output directory
- `path_truth`: absolute path to the inversion truth file
    - **/!\ for now best to use a BED formatted file** (see [evaluation_commands.md](https://gitlab.inria.fr/sromain/invpg-annot_publication/-/blob/main/expes/eval_and_stats/evaluation_commands.md?ref_type=heads#requires-truth_file-to-be-in-bed-format))
- `tools`: list of the graph building tools to use (separated by spaces), among {`cactus`, `minigraph`, `mgc`, `pggb`}.

### Command

```bash
expe_run.sh dir_input_files expe_id path_output path_truth tools
```

Examples:
```bash
# Human chm13 chr21
expe_run.sh . 2hap_div0.1 . chr21_sim_100inv.bed cactus minigraph mgc pggb
expe_run.sh . 2hap_div0 . chr21_sim_100inv.bed cactus minigraph mgc pggb

# Coeno carc chr6
expe_run.sh ../paper_PG_INV/data/new_sim_coeno/50_inv/ 2hap_div0.1 ../paper_PG_INV/results/carc_chr6/50_inv/ ../paper_PG_INV/results/carc_chr6/chr6_sim_50inv.bed cactus minigraph mgc pggb
expe_run.sh ../paper_PG_INV/data/new_sim_coeno/100_inv/ 2hap_div0.1 ../paper_PG_INV/results/carc_chr6/100_inv/ ../paper_PG_INV/results/carc_chr6/chr6_sim_100inv.bed cactus minigraph mgc pggb
```

### Output

The output is sorted into 4 subfolders in the given output directory, which are created if not yet existing:
- `graphs`: contains the graphs created with each given graph building tools, named as "`expe_id`_`tool`.gfa"
- `bubbles`: contains the VCFs of the bubbles extracted for each created graphs, named as "`expe_id`_`tool`.vcf"
- `inversions`: contains both the annotation (in BED format) and the annotation statistics for each created graph, named as "`expe_id`\_`tool`.bed" and "`expe_id`_`tool`.stats"
- `evaluation`: contains the graph, annotation and evaluation statistics for each created graph, "`expe_id`_`tool`.eval"

### Prepare merged results for paper figures

#### .eval file

Example with human chm13 chr21 simulation at div 0%:
```bash
cd ../paper_PG_INV/results/chm13_chr21/evaluation
cat 2hap_div0_pggb_v0.7.4.eval | sed 's/pggb/pggb_v0.7.4/g' > merged_2hap_div0.eval
cat $(ls 2hap_div0_*.eval | grep -v 'pggb_v0.7.4') | grep -v '^PG' >> merged_2hap_div0.eval
```

#### .intersect file

Example with human chm13 chr21 simulation:

```bash
cd ../paper_PG_INV/results/chm13_chr21/evaluation
for bed in $(ls ../inversions/2hap_div*_*.bed); do fbname=$(basename "$bed" .bed); arrEXPE=(${fbname//_/ }); pg=${arrEXPE[2]}; div=${arrEXPE[1]}; bedtools intersect -a ../chr21_sim_100inv.bed -b $bed -wao -f 0.5 -r | awk -v OFS='\t' -v pg="$pg" -v div="$div" '{print pg,div,$2,$3,$5,$6,$7}'; done > merged_2hap.intersect
sed -i 's/div//g' merged_2hap.intersect
```

## Secondary scripts used by the pipeline

- [build_cactus_graph.sh](expes/graph_construction_scripts/build_cactus_graph.sh)
- [build_minigraph_graph.sh](expes/graph_construction_scripts/build_minigraph_graph.sh)
- [build_mgc_graph.sh](expes/graph_construction_scripts/build_mgc_graph.sh)
- [build_pggb_graph.sh](expes/graph_construction_scripts/build_pggb_graph.sh)
- [minigraph_call_pipeline.sh](expes/bubble_calling/minigraph_call_pipeline.sh): calls bubbles in minigraph graphs and convert them to VCF format
- [graph_annot_statistics.sh](expes/eval_and_stats/graph_annot_statistics.sh): computes the graph, annotation and evaluation statistics from a GFA file, the INVPG-annot stats output , and the inversion truth file
- [redundancy_stats.py](expes/eval_and_stats/redundancy_stats.py): called by `graph_annot_statistics.sh`, computes the annotation redundancy statistics