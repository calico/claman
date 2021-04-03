# Claman

*Calico Lipidomics and Metabolomics Analysis*

# Overview

Claman provides methods for analyzing mzroll metabolomics/lipidomics datasets generated using [MAVEN](https://github.com/eugenemel/maven).

Each mzroll (**.mzrollDB**) file is an sqlite database that represents a single experimental metabolomics or lipidomics experiment generated using a shared analytical method (same mode, fragmentation, chromatography, etc). The major experimental data used for interpreting an LC-MS experiment are abundances of peaks across samples for each distinct peakgroup (ions with a characteristic mass and retention time). Claman extracts this data from mzrolls and stores it as an **mzroll_list**; a list of peakgroups, samples, and peaks. This data is further augmented with feature-attributes such as pathway information with sample-attributes which reflect an experiments design. The addition of sample attributes, allows users to integrate mzrolls which include a common set of samples but profiled using different methods (such a positive-, and negative-mode measurements).

mzroll_lists can then be filtered, floored to the limit of detection, and normalized using a suite of built-in methods that can correct for differences in sample loading, compare treatments to a reference condition, and align batches. Next, analytes of interest can be identified using peakgroup-level statistical analysis. Significant findings can then be explored using peakgroup-level visualizations and by identifying pathways that are enriched for significant metabolite changes.

# Initial Setup

The package can be installed from GitHub using:

```r
# install.packages("remotes")
remotes::install_github('calico/claman')

# for vignettes
remotes::install_github('calico/claman', build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
```

# Tutorial

After installing this package, the package functions can be accessed by just loading the package in R:

```r
library(claman)
```

Once calicomics is installed, run package vignette(s):

```r
vignette(package = "claman", topic = "metabolomics_example")
```