library(TOSCA)
library(dplyr)
library(ggplot2)
data("exampleData_CNA")
set.seed(123)

mutations=exampleData_CNA$Mutations
parameters=exampleData_CNA$Parameters
samples=exampleData_CNA$Samples
therapies=exampleData_CNA$Therapies

x=init(mutations=mutations,samples=samples,therapies=therapies,parameters=parameters)
exampleFit=TOSCA::fit(x,model_name='CNA',
n_iterations=10000,
n_chains=4,
warm_up=5000)
usethis::use_data(exampleFit, overwrite = TRUE)
