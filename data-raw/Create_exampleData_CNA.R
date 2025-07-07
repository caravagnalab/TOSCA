# Timeline
t_eca = .2
Sample_1 = .25

step_1_start_cycle_1 = .3
step_1_end_cycle_1 = .35

step_1_start_cycle_2 = .5
t_CNA_1 = .55
step_1_end_cycle_2 = .65

step_2_start_cycle_1 = .8
step_2_end_cycle_1 = .9

t_CNA_2 = 1

t_mrca = 1.2
Sample_2 = 5

# Parameters
omega = 5
omega_alpha = 50
omega_beta = 10
hist(rgamma(10000,omega_alpha,omega_beta))

l_diploid = 3e9
l_CNA1 = 4.5e7
l_CNA2 = 2e7
mu_clock = 1e-8
mu_th_step_1 = 5e-8
mu_th_step_2 = 1e-7

mu_th_step1_alpha = 10
mu_th_step1_beta = 2e8
hist(rgamma(10000,mu_th_step1_alpha, mu_th_step1_beta))

mu_th_step2_alpha = 3
mu_th_step2_beta = 2e7
hist(rgamma(10000,mu_th_step2_alpha, mu_th_step2_beta))

# Driver legato a evento esogeno
m_clock = 2*omega*l_diploid*mu_clock*(t_mrca-t_eca)

coeff_cna_1 = 2
m_alpha_cna_1 = omega*l_CNA1*(mu_clock*( (step_1_start_cycle_1-t_eca) + (step_1_start_cycle_2 - step_1_end_cycle_1) ) +
                                mu_th_step_1*(step_1_end_cycle_1-step_1_start_cycle_1))
m_beta_cna_1 = coeff_cna_1*omega*l_CNA1*(
  mu_clock * ( (step_1_end_cycle_2-step_1_end_cycle_2)  + (t_mrca-step_2_end_cycle_1)) +
  mu_th_step_1 * ( (step_1_end_cycle_2-t_CNA_1) )  +
  mu_th_step_2 * (step_2_end_cycle_1-step_2_start_cycle_1)
)

coeff_cna_2 = 4
m_alpha_cna_2 = omega*l_CNA2*(
  mu_clock * ( (step_1_start_cycle_1-t_eca) + (step_1_start_cycle_2-step_1_end_cycle_1) + (step_2_start_cycle_1-step_1_end_cycle_2) + (t_CNA_2-step_2_end_cycle_1))+
  mu_th_step_1 * ((step_1_end_cycle_1-step_1_start_cycle_1) + (step_1_end_cycle_2-step_1_start_cycle_2)) +
  mu_th_step_2 *(step_2_end_cycle_1-step_2_start_cycle_1)
)
m_beta_cna_2 = coeff_cna_2*omega*l_CNA2*(
  mu_clock*(t_mrca-t_CNA_2)
)

m_step_1 = 2*omega*l_diploid*mu_th_step_1*((step_1_end_cycle_1-step_1_start_cycle_1) + (step_1_end_cycle_2-step_1_start_cycle_2))
m_step_2 = 2*omega*l_diploid*mu_th_step_2*(step_2_end_cycle_1-step_2_start_cycle_1)

N = exp(omega*(Sample_2-t_mrca))
N_min = N*.9
N_max = N*1.1

exampleData_CNA = list(
  'Clinical Timepoints' =
    data.frame(
      'Clinical.name'= c('Sample','Sample','Therapy step','Therapy step','Therapy step'),
      'Clinical.type'= c('1','2', '1','1','2'),
      'Clinical.value.start'= c(Sample_1, Sample_2, step_1_start_cycle_1, step_1_start_cycle_2, step_2_start_cycle_1),
      'Clinical.value.end'= c(NA, NA, step_1_end_cycle_1, step_1_end_cycle_2, step_2_end_cycle_1),
      'Clinical.index' = c(NA, NA, '1','2','1')
    ),
  'Mutations' =
    data.frame(
      'Mutation.name'= c('m_clock','m_cna','m_cna','m_cna','m_cna', 'm_th', 'm_th'),
      'Mutation.type'= c('clock','alpha','beta','alpha','beta', 'step','step'),
      'Mutation.index'= c(NA, '1','1','2','2','1','2'),
      'Mutation.value'= c(m_clock, m_alpha_cna_1, m_beta_cna_1, m_alpha_cna_2, m_beta_cna_2, m_step_1, m_step_2)
    ),
  'Parameters' =
    data.frame('Parameter.name' = c('l_diploid',
                                    'l_CNA','l_CNA',
                                    'mu_th_step_alpha','mu_th_step_beta',
                                    'mu_th_step_alpha','mu_th_step_beta',
                                    'coeff_CNA','coeff_CNA',
                                    'mu_clock','omega_alpha','omega_beta','k_step','N_min','N_max',
                                    'exponential_growth', 'mrca_alpha','mrca_beta'),
               'Parameter.value' = c(l_diploid,
                                     l_CNA1,l_CNA2,
                                     mu_th_step1_alpha,mu_th_step1_beta,
                                     mu_th_step2_alpha,mu_th_step2_beta,
                                     coeff_cna_1, coeff_cna_2,
                                     mu_clock, omega_alpha, omega_beta, 10, N_min, N_max, 1,1,1),
               'Parameter.index' = c(NA,'1','2','1','1','2','2', '1', '2', rep(NA,9))
    )
)

#saveRDS(exampleData, file = "exampleData.rds")

usethis::use_data(exampleData_CNA, overwrite = TRUE)





