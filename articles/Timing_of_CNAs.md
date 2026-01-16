# Timing of CNAs

Understanding when and how genetic resistance to therapy emerges is a
central problem in cancer genomics. Longitudinal sequencing data offer a
powerful way to address this question, but require specialized tools to
reconstruct evolutionary histories and interpret complex genomic events.
This vignette demonstrates how can be used to analyze such data using a
real-world case study from acute myeloid leukemia.

``` r
library(TOSCA)
#> Loading required package: ggplot2
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## Examine the input

In this vignette, we walk through a complete analysis of longitudinal
cancer sequencing data from a single patient with acute myeloid leukemia
who developed resistance to immunological therapy. This real-world
example highlights how can be used to identify and interpret clinically
relevant genomic events.

``` r
data("UPN06")
```

Create the TOSCA object

Fit the appropriate model

Check the quality of the fit

Visualise the result of the fit
