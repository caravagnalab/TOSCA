# Package index

## Available data

Datasets included in the package.

- [`exampleData_CNA_dormancy`](https://caravagnalab.github.io/TOSCA/reference/exampleData_CNA_dormancy.md)
  : Dataset with a toy example for the CNA model with dormancy
- [`exampleData_CNA`](https://caravagnalab.github.io/TOSCA/reference/exampleData_CNA.md)
  : Dataset with a toy example for the CNA model
- [`exampleData_Driver`](https://caravagnalab.github.io/TOSCA/reference/exampleData_Driver.md)
  : Dataset with a toy example for the Driver model

## Create TOSCA object

Functions to create the TOSCA object and check the format of input data.

- [`init()`](https://caravagnalab.github.io/TOSCA/reference/init.md) :
  Initialise TOSCA object
- [`check_clinical_input()`](https://caravagnalab.github.io/TOSCA/reference/check_clinical_input.md)
  : Check clinical input
- [`check_genomic_input()`](https://caravagnalab.github.io/TOSCA/reference/check_genomic_input.md)
  : Check that the input mutation dataframe is valid
- [`check_parameters_input()`](https://caravagnalab.github.io/TOSCA/reference/check_parameters_input.md)
  : Checks that the parameters dataframe is valid
- [`check_required_cols()`](https://caravagnalab.github.io/TOSCA/reference/check_required_cols.md)
  : Check that for the required columns for each dataframe are present
- [`check_valid_date()`](https://caravagnalab.github.io/TOSCA/reference/check_valid_date.md)
  : Check that the input date is valid
- [`check_valid_karyotype()`](https://caravagnalab.github.io/TOSCA/reference/check_valid_karyotype.md)
  : Check that the input karytypes in the mutations dataframe are valid
- [`convert_date_real()`](https://caravagnalab.github.io/TOSCA/reference/convert_date_real.md)
  : Converts date into real number
- [`convert_real_date()`](https://caravagnalab.github.io/TOSCA/reference/convert_real_date.md)
  : Converts real number into date

## Fit model

Function to fit model.

- [`fit()`](https://caravagnalab.github.io/TOSCA/reference/fit.md) :
  Infer the timing of the events of interest fitting the appropriate
  TOSCA model

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
- [`get_name_mu_step()`](https://caravagnalab.github.io/TOSCA/reference/get_name_mu_step.md)
  : Title
- [`plot_prior_vs_posterior()`](https://caravagnalab.github.io/TOSCA/reference/plot_prior_vs_posterior.md)
  : Produce collective plot of prior vs posterior for all inferred
  parameters
- [`plot_prior_vs_posterior_single_parameter()`](https://caravagnalab.github.io/TOSCA/reference/plot_prior_vs_posterior_single_parameter.md)
  : Prior vs Posterior distribution of inferred parameter
- [`plot_timing()`](https://caravagnalab.github.io/TOSCA/reference/plot_timing.md)
  : Plot clinical timeline + posterior times with histogram
- [`plot_timing_MAP()`](https://caravagnalab.github.io/TOSCA/reference/plot_timing_MAP.md)
  : Plot clinical timeline + posterior times with MAP
- [`plot_timing_density()`](https://caravagnalab.github.io/TOSCA/reference/plot_timing_density.md)
  : Plot clinical timeline + posterior times with density
