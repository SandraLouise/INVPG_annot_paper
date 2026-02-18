# Running the R script

## Dependancies

- R (script written and tested on R 4.5.0)
- ggplot2
- gridExtra
- optparse

## Usage

```bash
Rscript --vanilla draw_figures.R -e merged.eval -i merged.intersect
```

or

```bash
Rscript --vanilla draw_figures.R -e merged.eval -i merged.intersect -o output_dir
```

Help page:
```txt
Rscript --vanilla draw_figures.R --help

Options:
	-e CHARACTER, --evalFile=CHARACTER
		eval file name

	-i CHARACTER, --intersectFile=CHARACTER
		intersect file name

	-o CHARACTER, --out=CHARACTER
		output directory name [default= ./]

	-h, --help
		Show this help message and exit
```
### Input

The results (`.eval` and `.intersect`) for each experiment (PG x dataset) must be merged into single files to be plotted.

For `.intersect` files, a simple `cat *.intersect > merged.intersect` suffices.

For `.eval` files (having a header line):
```bash
cat expe_1.eval > merged.eval
# then for each remaining expe
grep -v '^PG' expe_x.eval >> merged.eval
```

**Note:** the updated pipeline should now generate these merged files (if it ran all expe in one go).

### Output

Draws the figures 3 and 4 of the bioRxiv paper in both pdf and png format.