# Plot posterior predictive checks

Plot posterior predictive checks

## Usage

``` r
plot_ppc(x)
```

## Arguments

- x:

  TOSCA object

## Value

plot of the posterior predictive distribution versus the real value of
all mutation groups provided in input

## Examples

``` r
data("exampleFit")
plot_ppc(exampleFit)
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Warning: Arguments in `...` must be used.
#> ✖ Problematic argument:
#> • size = 0.05
#> ℹ Did you misspell an argument name?
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the TOSCA package.
#>   Please report the issue to the authors.
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Warning: Arguments in `...` must be used.
#> ✖ Problematic argument:
#> • size = 0.05
#> ℹ Did you misspell an argument name?
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Warning: Arguments in `...` must be used.
#> ✖ Problematic argument:
#> • size = 0.05
#> ℹ Did you misspell an argument name?
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Warning: Arguments in `...` must be used.
#> ✖ Problematic argument:
#> • size = 0.05
#> ℹ Did you misspell an argument name?
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Warning: Arguments in `...` must be used.
#> ✖ Problematic argument:
#> • size = 0.05
#> ℹ Did you misspell an argument name?
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Warning: Arguments in `...` must be used.
#> ✖ Problematic argument:
#> • size = 0.05
#> ℹ Did you misspell an argument name?
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Warning: Arguments in `...` must be used.
#> ✖ Problematic argument:
#> • size = 0.05
#> ℹ Did you misspell an argument name?
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
#> Warning: font family 'Times New Roman' not found in PostScript font database
```
