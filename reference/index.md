# Package index

## Available data

Datasets included in the package.

- [`exampleData_CNA_dormancy`](https://caravagnalab.github.io/TOSCA/reference/exampleData_CNA_dormancy.md)
  : Dataset with a toy example for the CNA model with dormancy
- [`exampleData_CNA`](https://caravagnalab.github.io/TOSCA/reference/exampleData_CNA.md)
  : Dataset with a toy example for the CNA model
- [`exampleData_Driver`](https://caravagnalab.github.io/TOSCA/reference/exampleData_Driver.md)
  : Dataset with a toy example for the Driver model
- [`exampleFit`](https://caravagnalab.github.io/TOSCA/reference/exampleFit.md)
  : Example of a TOSCA object with compleated inference
- [`UPN06`](https://caravagnalab.github.io/TOSCA/reference/UPN06.md) :
  Input data for patient UPN06
- [`D9MRCY`](https://caravagnalab.github.io/TOSCA/reference/D9MRCY.md) :
  Input data for patient D9MRCY

## Create TOSCA object

Functions to create the TOSCA object and check the format of input data.

- [`init()`](https://caravagnalab.github.io/TOSCA/reference/init.md) :
  Initialise TOSCA object

- [`convert_date_real()`](https://caravagnalab.github.io/TOSCA/reference/convert_date_real.md)
  : Converts real into date

- [`convert_real_date()`](https://caravagnalab.github.io/TOSCA/reference/convert_real_date.md)
  : Converts date into real

- [`print(`*`<TOSCA>`*`)`](https://caravagnalab.github.io/TOSCA/reference/print.TOSCA.md)
  :

  Print for class `'TOSCA'`.

## Fit model

Function to fit model.

- [`fit()`](https://caravagnalab.github.io/TOSCA/reference/fit.md) :
  Infer the timing of the events of interest fitting the appropriate
  TOSCA model
- [`get_inference_data()`](https://caravagnalab.github.io/TOSCA/reference/get_inference_data.md)
  : Get inference data

## Plot fit

Function to visualise the results of the fit.

- [`check_ppc()`](https://caravagnalab.github.io/TOSCA/reference/check_ppc.md)
  : Posterior Predictive checks
- [`plot_expected_N()`](https://caravagnalab.github.io/TOSCA/reference/plot_expected_N.md)
  : Posterior distribution of the number of cells collected in the first
  and second samples.
- [`plot_ppc()`](https://caravagnalab.github.io/TOSCA/reference/plot_ppc.md)
  : Plot posterior predictive checks
- [`plot_ppc_single_mut()`](https://caravagnalab.github.io/TOSCA/reference/plot_ppc_single_mut.md)
  : Plot the posterior predictive distribution and compares it to the
  real number of mutations
- [`plot_prior_vs_posterior()`](https://caravagnalab.github.io/TOSCA/reference/plot_prior_vs_posterior.md)
  : Produce collective plot of prior vs posterior for all inferred
  parameters
- [`plot_prior_vs_posterior_single_parameter()`](https://caravagnalab.github.io/TOSCA/reference/plot_prior_vs_posterior_single_parameter.md)
  : Prior vs Posterior distribution of inferred parameter
- [`plot_timing()`](https://caravagnalab.github.io/TOSCA/reference/plot_timing.md)
  : Plot clinical timeline + posterior times with histogram
- [`plot_timing_MAP()`](https://caravagnalab.github.io/TOSCA/reference/plot_timing_MAP.md)
  : Plot clinical timeline + posterior times with MAP
- [`days_from()`](https://caravagnalab.github.io/TOSCA/reference/days_from.md)
  : Computes the days of the event of interest from a date of interest

## Examine fit

Get the results from the inference.

- [`get_fit_summary()`](https://caravagnalab.github.io/TOSCA/reference/get_fit_summary.md)
  : Get mean, mode and q5, q95
