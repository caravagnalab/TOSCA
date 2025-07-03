# Timeline
t_eca = .2
Sample_1 = .25

step_1_start_cycle_1 = .3
step_1_end_cycle_1 = .35

step_1_start_cycle_2 = .5
t_driver = .55
step_1_end_cycle_2 = .65

t_mrca = 1
Sample_2 = 5

# Parameters
omega = 5
omega_alpha = 50
omega_beta = 10
hist(rgamma(10000,omega_alpha,omega_beta))

l_diploid = 3e9
mu_clock = 1e-9
mu_clock_driver = 5e-9

mu_driver = 5e-8
mu_driver_alpha = 10
mu_driver_beta = 2e8
hist(rgamma(10000,mu_driver_alpha, mu_driver_beta))

expected_mutation_rate_th = function(muts_th, omega, delta_therapy, length_genome=3e9){
  mu_th_LB = muts_th / (2*length_genome*omega*delta_therapy)
  mu_th_UB = muts_th / (2*length_genome)
  mu_th = mean(mu_th_LB, mu_th_UB)
  mu_th
}

expected_mutation_rate_th(mu_driver, 5,(.35-.3)+(.65-.5))

get_hyperparameters_mutation_rate = function(mu_th, delta_therapy, var=.05){
  alpha_mu_th = mu_th*delta_therapy*var
  beta_mu_th = delta_therapy*var
  return(c('alpha'=alpha_mu_th, 'beta'=beta_mu_th))
}

# Driver legato a evento esogeno
m_clock = 2*omega*l_diploid*(mu_clock*(t_driver-t_eca) + mu_clock_driver*(t_mrca-t_driver))
m_driver = 2*omega*l_diploid*mu_driver*(step_1_end_cycle_2-t_driver)
N = exp(omega*(Sample_2-t_mrca))
N_min = N*.9
N_max = N*1.1

exampleData = list(
  'Clinical Timepoints' =
    data.frame(
      'Clinical.name'= c('Sample','Sample','Therapy driver','Therapy driver'),
      'Clinical.type'= c('1','2', '2','2'), # Driver == 1 endogeno, ==2 esogeno
      'Clinical.value.start'= c(Sample_1, Sample_2, step_1_start_cycle_1, step_1_start_cycle_2),
      'Clinical.value.end'= c(NA, NA, step_1_end_cycle_1, step_1_end_cycle_2),
      'Clinical.index' = c(NA, NA, '1','2')
    ),
  'Mutations' =
    data.frame(
      'Mutation.name'= c('m_clock','m_driver'),
      'Mutation.type'= c('clock','2'),
      'Mutation.index'= c(NA, NA),
      'Mutation.value'= c(m_clock, m_driver)
    ),
  'Parameters' =
    data.frame('Parameter.name' = c('l_diploid','mu_clock','mu_clock_driver','mu_driver_alpha','mu_driver_beta','omega_alpha','omega_beta','k_step','N_min','N_max','exponential_growth'),
               'Parameter.value' = c(l_diploid, mu_clock, mu_clock_driver, mu_driver_alpha, mu_driver_beta, omega_alpha, omega_beta, 10, N_min, N_max, 1),
               'Parameter.index' = c(rep(NA, 11))
               )
)

#saveRDS(exampleData, file = "exampleData.rds")

usethis::use_data(exampleData, overwrite = TRUE)





