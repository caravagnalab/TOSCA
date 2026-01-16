# Get inference data

Get inference data

## Usage

``` r
get_inference_data(x, model, dormancy = F, max_mrca = NA, reg_dormancy = 0)
```

## Arguments

- x:

  TOSCA obj

- model:

  string, either "CNA" or "Driver"

- dormancy:

  boolean, if model == "CNA", it is possible to include dormancy

- max_mrca:

  date in "YYYY-mm-dd" format, maximum date at which the MRCA could have
  been born (if \< than the date of the second sample)

- reg_dormancy:

  boolean (0 or 1), optional regularisation of the growth imposed in
  case of dormancy

## Value

the list of arguments passed to the stan fit function
