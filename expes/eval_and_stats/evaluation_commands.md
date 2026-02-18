# Commands for statistics and evaluation

## Requires `truth_file` to be in BED format

To convert a standard VCF into BED:

- case of "END=" tag in INFO field, at first position in the INFO field:
    ```bash
    grep -v '^#' truth.vcf | awk -v OFS='\t' '{split($8,info,";"); print $1,$2,info[1]}' | sed 's/END=//g' > truth.bed
    ```
- case of "END=" tag in INFO field, at `i` position (1-based) in the INFO field:
    ```bash
    grep -v '^#' truth.vcf | awk -v OFS='\t' '{split($8,info,";"); print $1,$2,info[i]}' | sed 's/END=//g' > truth.bed
    ```

## Automated script `graph_annot_statistics.sh`

Parameters:
- Graph `GFA`
- True inversions `truth.bed`
- Annot inversions `annot.bed`
- INVPG-annot stats file `annot.stats`
- Expe identifier (for later merge of all eval files) `expeID`

```bash
graph_annot_statistics.sh GFA truth.bed annot.bed annot.stats expeID > expeID.eval
```

### Example

```bash
graph_annot_statistics.sh ../results/pangenomes/sim_pg/minigraph/div01.gfa ../data/sim_data/inversions/simulated_inv.bed ../INVPG_annot/test_dev/test_stats/sim_minigraph_div01.bed ../INVPG_annot/test_dev/test_stats/sim_minigraph_div01.stats minigraph_div01 > ../INVPG_annot/test_dev/test_stats/sim_minigraph_div01.eval
```

```text
ExpeID	minigraph_div01
Graph_size	21425278
Node_number	253
Edge_number	428
Total_bubbles	107
Large_bubbles	80
Inversion_bubbles	74
Path-explicit	69
Alignment-rescued	5
Recall	65
False_positives	0
NonRedundant_bubbles	65
Redundant_bubbles	0
Imprecise_bubbles	9
```

## Individual commands

### Graph Statistics

**Requires:** 
- graph `gfa`
- `vcf`output by `vg deconstruct`/`minigraph pipeline`

|                   |                      |
| :---              | :---                 |
| Graph size        | `grep '^S' <gfa> \| awk '{sum+=length($3)} END{print sum}'` |
| # Nodes   | `grep -c '^S' <gfa>` |
| # Edges   | `grep -c '^L' <gfa>` |
| # Bubbles | `grep -cv '^#' <vcf>` |

### Annotation Statistics

**Requires:**
- INVPG-annot `.stats` file. The stats file also reports the total number of bubbles in the input VCF.

|                                 |                    |
| :---                            | :---               |
| # Large bubbles                 | `annotation.stats` |
| # Inversion bubbles             | `annotation.stats` |
| # Path-explicit annotations     | `annotation.stats` |
| # Alignment-rescued annotation  | `annotation.stats` |

### Annotation Evaluation

**Requires:**
- **bedtools**
- **redundancy_stats.py**
- `annot.bed`:  output annotation 
- `truth.bed`:  true inversions (in bed format)

|                                |                      |
| :---                           | :---                 |
| Recall (50% reciproq. overlap) | `bedtools intersect -a <truth.bed> -b <annot.bed> -f 0.5 -r \| wc -l` <br/> / #TrueINV * 100 (if #TrueINV != 100) |
| Redundancy                     | `bedtools intersect -a <truth.bed> -b <annot.bed> -wao > intersect.tsv` <br/> `python3 redundancy_stats.py intersect.tsv > redundancy.stats` |

#### Redundancy stats
|                              |       |
| :---                         | :---  |
| # Non-redundant & TP bubbles | `redundancy.stats` |
| # Redundant bubbles          | `redundancy.stats` |
| # Imprecise bubbles          | `redundancy.stats` |