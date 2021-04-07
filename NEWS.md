# cinaR 0.2.2
* Added:
  - `comparison.scheme` either set to `OVO` (default) which will compare each contrasts to each other. If it's set to `OVA` each contrast will be compared to rest. 

# cinaR 0.2.1
* Added:
  - An extra argument `sv.number` is added if user wants to set the number of surrogate variables.  
  More info [here](https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf).
  - `verbose` argument, default set to `TRUE`. I think, I have implemented this in a really cool way!
  - [`cinaRgenesets`](https://github.com/eonurk/cinaR-genesets) added to vignette.

# cinaR 0.2.0
* Added:
  - `heatmap_differential` function to be able to plot the heatmaps
of a given comparison.
  - `show_comparisons` is added to see the available comparisons 
before you plot differential heatmaps.

* Changed:
  - `heatmap_plot` is now renamed to `heatmap_var_peaks`, which makes more sense.

# cinaR 0.1.1
This version is even cooler, lol:

- Comparisons are now sorted according to `contrast` order
- Added bulk RNA-seq support/documentation
- Removed many dependencies
- Small bug fixes and improvements

# cinaR 0.1.0
- `cinaR` is published in CRAN repository!
