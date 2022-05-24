# quantile-intervals

`bedtools` and `tidyverse` workflow to calculate the distribution of a normalized variable (`bedgraph`) for each position within genomic intervals (`GFF`).

<p align="center">
<img src="example_plot.png" width="300" height="300">
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
* `R >= 4.1.2`
* `tidyverse >= 1.3.1`
* `reshape2 >= 1.4.4`
