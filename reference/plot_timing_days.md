# Plot clinical timeline + posterior times with histogram in days from a chosen date

Plot clinical timeline + posterior times with histogram in days from a
chosen date

## Usage

``` r
plot_timing_days(x, days_from = NULL, event_name = "Custom zero")
```

## Arguments

- x:

  TOSCA obj

## Value

posterior distributions of the inferred timing variables on the clinical
timeline (alternative visualisation)

## Examples

``` r
data("exampleFit")
plot_timing(exampleFit)
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> No id variables; using all as measure variables
```
