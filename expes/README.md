# Commands and scripts to reproduce pangenome graph construction, bubble annotation and evaluation

> [!NOTE]\
> Simply cloning this repo does not create the directory structure. You may need to adapt the command line instructions to your own data to proceed.

The automated pipeline allows to reproduce the results of the paper. It encompasses the graph construction, bubble calling, inversion annotation and annotation evaluation.

For a given set of experiments, we create a directory `expe/` (outside of the directory where we stored or simulated the input haplotype sequences) with the following structure:

```
expeXX/
└── graphs/
└── bubbles/
└── annotations/
└── evaluations/
```


## Main scripts

The automated pipeline is run by two main scripts:
- `expe_run.sh`: manages the automation, launches the pipeline (with SBATCH) for a given dataset and with the graph building tools given by the user.
    - git: [expe_run.sh](expe_run.sh)
- `run_pipeline.sh`: called by `expe_run.sh` for each graph building tool, contains the pipeline.
    - git: [run_pipeline.sh](run_pipeline.sh)

## Usage

### Input

- `dir_input_files`: absolute path to the directory containing the genome input files
    - **/!\ the path of the reference genome must be listed first in all input files**
- `expe_id`: identifier of the experience
    - **/!\ must match the filename of the input files** (_e.g._ "2hap_div0.1" for input files named "*tool*_input_2hap_div0.1.txt")
- `path_output`: absolute path to the output directory
- `path_truth`: absolute path to the inversion truth file
    - /!\ for now best to use a BED formatted file (see [evaluation_commands.md](eval_and_stats/evaluation_commands.md?ref_type=heads#requires-truth_file-to-be-in-bed-format))
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

### Output

The output is sorted into 4 subfolders in the given output directory, which are created if not yet existing:
- `graphs`: contains the graphs created with each given graph building tools, named as "`expe_id`_`tool`.gfa"
- `bubbles`: contains the VCFs of the bubbles extracted for each created graphs, named as "`expe_id`_`tool`.vcf"
- `inversions`: contains both the annotation (in BED format) and the annotation statistics for each created graph, named as "`expe_id`\_`tool`.bed" and "`expe_id`_`tool`.stats"
- `evaluation`: contains the graph, annotation and evaluation statistics for each created graph, "`expe_id`_`tool`.eval"
```


## Secondary scripts used by the pipeline

- [build_cactus_graph.sh](graph_construction_scripts/build_cactus_graph.sh)
- [build_minigraph_graph.sh](graph_construction_scripts/build_minigraph_graph.sh)
- [build_mgc_graph.sh](graph_construction_scripts/build_mgc_graph.sh)
- [build_pggb_graph.sh](graph_construction_scripts/build_pggb_graph.sh)
- [minigraph_call_pipeline.sh](bubble_calling/minigraph_call_pipeline.sh): calls bubbles in minigraph graphs and convert them to VCF format
- [graph_annot_statistics.sh](eval_and_stats/graph_annot_statistics.sh): computes the graph, annotation and evaluation statistics from a GFA file, the INVPG-annot stats output , and the inversion truth file
- [redundancy_stats.py](eval_and_stats/redundancy_stats.py): called by `graph_annot_statistics.sh`, computes the annotation redundancy statistics


## Input data requirements

- the input haplotypes should be in individual fasta files (not necessarily in the same directory, not necessarily in a directory of their own.
- 2 (uncompatible) types of fasta are necessary:
	- for PGGB, the headers of fasta entries must respect the [panSN spec](https://github.com/pangenome/PanSN-spec).  Example: `chm13#0#chr21`. (species name + "#" + haplotype phase + "#" + scaffold/chromosome name)
	- for Cactus, the headers must not contain any "#" characters
	- for MGC and minigraph: we use the PGGB format
