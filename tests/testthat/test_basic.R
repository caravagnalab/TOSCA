library(testthat)
library(TOSCA)

mutations = exampleData_CNA$Mutations
parameters = exampleData_CNA$Parameters
samples = exampleData_CNA$Samples
therapies = exampleData_CNA$Therapies

test_that("init() returns a TOSCA object", {
  tosca_obj <- init(mutations, parameters, samples, therapies)
  expect_s3_class(tosca_obj, "TOSCA")
})
