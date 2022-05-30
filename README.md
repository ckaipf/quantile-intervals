# quantile-intervals

Workflow to visualize the distribution of a continuous variable (`bedgraph`) for each position within genomic intervals (`GFF`).

<p align="center">
<img src="example_plot.png" width="400" height="400">
</p>

## Notes
* Run the script:
```
nextflow main.nf
```
* Workflow takes a list of genomic intervals as a `GFF` file and a genome-wide continuous variable in `bedgraph` format.
* Script was developed especially for the distribution of a variable around a specific position.
  + Intervals get aligned at their centers. Is is not necessary that they are of equal length. However, they should be aligned correctly to enable the detection of distributional differences.
* For input formats consult the example.
  + Example data was taken from Dar and Sorek (2018).
* Different normalization strategies are implemented and can be set within the file `main.nf`.

## Requirements
* `nextflow >= 21.04.3.5560`
* `bedtools >= 2.27.1`
* `python >= 3.7.3`
* `R >= 4.1.2`
* `tidyverse >= 1.3.1`
