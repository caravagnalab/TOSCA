# Timeline
t_eca = .2
t_mrca_sample_1 = .23
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

# Driver legato a evento esogeno
m_clock_primary = 2*omega*l_diploid*mu_clock*(t_mrca_sample_1-t_eca)
m_clock = 2*omega*l_diploid*(mu_clock*(t_driver-t_eca) + mu_clock_driver*(t_mrca-t_driver))
m_driver = 2*omega*l_diploid*mu_driver*(step_1_end_cycle_2-t_driver)
N2 = exp(omega*(Sample_2-t_mrca))
N1 = exp(omega*(Sample_1-t_mrca_sample_1))
N_min = N1*.9
N_max = N2*1.1

exampleData_Driver = list(
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
      'Mutation.name'= c('m_clock', 'm_clock','m_driver'),
      'Mutation.type'= c('primary', 'relapse','2'),
      'Mutation.index'= c(NA, NA, NA),
      'Mutation.value'= c(m_clock_primary, m_clock, m_driver),
      'Mutation.coeff'= rep(NA,3),
      'Mutation.source'= c('clock','clock','driver')
    ),
  'Parameters' =
    data.frame('Parameter.name' = c('l_diploid','mu_clock','mu_clock_driver','mu_driver_alpha','mu_driver_beta','omega_alpha','omega_beta','k_step','N_min','N_max','exponential_growth', 'mrca_alpha','mrca_beta',
                                    "mu_driver_clock", "mu_driver", "phi_clock","phi_driver"),
               'Parameter.value' = c(l_diploid, mu_clock, mu_clock_driver, mu_driver_alpha, mu_driver_beta, omega_alpha, omega_beta, 10, N_min, N_max, 1,1,1,
                                     mu_clock_driver, mu_driver,.1,.1),
               'Parameter.index' = c(rep(NA, 17))
               )
)

#saveRDS(exampleData, file = "exampleData.rds")

usethis::use_data(exampleData_Driver, overwrite = TRUE)





