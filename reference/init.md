# Initialise TOSCA object

Initialise TOSCA object

## Usage

``` r
init(mutations, parameters, samples, therapies, sample_name = "Example")
```

## Arguments

- mutations:

  data.frame with the information on mutation falling on the timed
  branches. 1. Name (str), name of choosing for the cluster of
  mutations 2. Length (int), length of the genomic segment they sit
  on 3. Karyptype (str), karyotype of the genomic segment they sit as
  "Major:minor" (ex "1:1") 4. Type (str), one between the following:
  "clock-like primary", "clock-like relapse", "alpha", "beta", "driver",
  the name of the mutagenic drug to which they are associated as
  reported in the therapies dataframe 5. Value (int), number of
  mutations

- parameters:

  data.frame with the parameters required for inference. 1. Name
  (str), 2. Value (real), 3. Index (int)

- samples:

  data.frame with the name and collection date of the longitudinal
  samples + date of birth (optional). 1. Name (str), name of the
  sample 2. Date (date) date of collection in the format "YYYY-mm-dd"

- therapies:

  data.frame with the name and exposure period to therapies affecting
  the mutation accumulation process. 1. Name (str) 2. Class (str), one
  of the following : "Mutagenic step", "Mutagenic cauchy", "Driver
  responsive", "Chemotherapy inducing dormancy" 3. Start (date), start
  of the cicle of therapy in the format "YYYY-mm-dd" 4. End (date), end
  of the cicle of therapy in the format "YYYY-mm-dd"

- sample_name:

  string with the name of the sample

## Value

TOSCA object

## Examples

``` r
library(TOSCA)
library(dplyr)
library(ggplot2)
data("exampleData_CNA")
mutations = exampleData_CNA$Mutations
parameters = exampleData_CNA$Parameters
samples = exampleData_CNA$Samples
therapies = exampleData_CNA$Therapies

x = init(mutations=mutations, samples=samples, therapies=therapies, parameters=parameters)
#> Error in init(mutations = mutations, samples = samples, therapies = therapies,     parameters = parameters): Input must be a data.frame.
```
