# Plot clinical timeline + posterior times with histogram

Plot clinical timeline + posterior times with histogram

## Usage

``` r
plot_timing(x)
```

## Arguments

- x:

  TOSCA obj

## Value

posterior distributions of the inferred timing variables on the clinical
timeline

## Examples

``` r
data("exampleFit")
plot_timing(exampleFit)
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> No id variables; using all as measure variables
```
