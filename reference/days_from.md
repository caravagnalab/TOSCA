# Computes the days of the event of interest from a date of interest

Computes the days of the event of interest from a date of interest

## Usage

``` r
days_from(x, parameter, from, cna_length, quantiles = c(5, 95))
```

## Arguments

- x:

  TOSCA obj

- parameter:

  the inferred time of interest

- from:

  the date of interest (string, "YY-mm-dd")

- cna_length:

  if the inferred date refers to a CNA, the length of the CNA

- quantiles:

## Value

the distance in days between an arbitrary date (parameter from) and a
chosen inferred time

## Examples

``` r
data("exampleFit")
days_from(exampleFit, parameter = "t_cna", from = "2010-01-04", cna_length = 4.5e+08, quantiles=c(5, 95))
#> # A tibble: 1 × 4
#>   median    mean      q5        q95      
#>   <drtn>    <drtn>    <drtn>    <drtn>   
#> 1 -619 days -621 days -691 days -556 days
```
