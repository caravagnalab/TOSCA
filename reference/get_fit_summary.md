# Get mean, mode and q5, q95

Get mean, mode and q5, q95

## Usage

``` r
get_fit_summary(
  x,
  parameter = NULL,
  cna_length = NULL,
  drug_name = NULL,
  quantiles = c(5, 95)
)
```

## Arguments

- x:

  TOSCA obj

- parameter:

  string, optional. Parameter of interest. If nothing is provided the
  function returns the mean, mode and q5, q95, of all parameters.
  Possible values: t_eca, t_mrca, t_cna, t_driver, omega, mu_th_step,
  mu_driver.

- cna_length:

  if parameter == t_cna, length of the CNA of interest

- drug_name:

  if parameter == mu_th_step, name of the drug of interest

- quantiles:

  vector of quantiles, default (5, 95)

## Value

a dataframe summarising the inference results

## Examples

``` r
data("exampleFit")
get_fit_summary(exampleFit, parameter = "t_cna", cna_length = 4.5e+08, quantiles=c(5, 95))
#> # A tibble: 1 × 7
#>   variable median    mean       rhat ess_bulk q5        q95      
#>   <chr>    <chr>     <chr>     <dbl>    <dbl> <chr>     <chr>    
#> 1 t_cna[2] 2008-4-25 2008-4-23  1.00   40515. 2008-2-13 2008-6-27
get_fit_summary(exampleFit, parameter = "mu_th_step", drug_name = "Drug 2", quantiles=c(2, 98))
#> # A tibble: 1 × 7
#>   variable  median       mean              rhat ess_bulk q2            q98      
#>   <chr>     <chr>        <chr>            <dbl>    <dbl> <chr>         <chr>    
#> 1 mu_Drug 2 1.131765e-07 1.2205458505e-07 1.000   14291. 6.6796186e-08 2.428349…
```
