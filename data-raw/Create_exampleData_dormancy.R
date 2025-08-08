t_eca = 1
t_mrca_sample_1 = 2
t_sampling_1 = 3

t_th_step_1_cycle_1_start = 4
t_th_step_1_cycle_1_end = 5

t_th_step_1_cycle_2_start = 6
t_th_step_1_cycle_2_end = 7

t_th_step_2_cycle_1_start = 8
t_wgd = 8.5
t_th_step_2_cycle_1_end = 9

t_mrca = 10
t_sampling_2 = 11

omega = 5
omega_alpha = 50
omega_beta = 10

mu_clock = 1e-8
mu_th_1 = 1e-7
mu_th_2 = 5e-8

mu_th_step1_alpha = 3 # 10
mu_th_step1_beta = 2e7 #2e8
hist(rgamma(10000,mu_th_step1_alpha, mu_th_step1_beta))

mu_th_step2_alpha = 10
mu_th_step2_beta = 2e8
hist(rgamma(10000,mu_th_step2_alpha, mu_th_step2_beta))

N = exp(omega*(t_sampling_2-t_mrca))
N_min = N*.9
N_max = N*1.1

N_sample_1 = exp(omega*(t_sampling_1-t_mrca_sample_1))
N_min_sample_1 = N_sample_1*.9
N_max_sample_1 = N_sample_1*1.1

l_cnloh = 1e9
l_tetraploid = 2e9

alpha_tetraploid_clock = 2*l_tetraploid*omega*mu_clock*(t_wgd-t_eca)
beta_tetraploid_clock = 4*l_tetraploid*omega*mu_clock*(t_mrca-t_wgd)
alpha_cnloh_clock = l_cnloh*omega*mu_clock*(t_wgd-t_eca)
beta_cnloh_clock = 2*l_cnloh*omega*mu_clock*(t_mrca-t_wgd)

alpha_tetraploid_step1 = 2*l_tetraploid*omega*mu_th_1*((t_th_step_1_cycle_1_end-t_th_step_1_cycle_1_start) + (t_th_step_1_cycle_2_end-t_th_step_1_cycle_2_start))
beta_tetraploid_step1 = 0 #4*l_tetraploid*omega*mu_th_1*()
alpha_cnloh_step1 = l_cnloh*omega*mu_th_1*((t_th_step_1_cycle_1_end-t_th_step_1_cycle_1_start) + (t_th_step_1_cycle_2_end-t_th_step_1_cycle_2_start))
beta_cnloh_step1 = 0 #2*l_cnloh*omega*mu_th_1*()

alpha_tetraploid_step2 = 2*l_tetraploid*omega*mu_th_2*(t_wgd-t_th_step_2_cycle_1_start)
beta_tetraploid_step2 = 4*l_tetraploid*omega*mu_th_2*(t_th_step_2_cycle_1_end-t_wgd)
alpha_cnloh_step2 = l_cnloh*omega*mu_th_2*(t_wgd-t_th_step_2_cycle_1_start)
beta_cnloh_step2  = 2*l_cnloh*omega*mu_th_2*(t_th_step_2_cycle_1_end-t_wgd)

l_diploid = 3e9
m_clock_primary = 2*omega*l_diploid*mu_clock*(t_mrca_sample_1-t_eca)

exampleData_WGD = list(
  'Clinical Timepoints' =
    data.frame(
      'Clinical.name'= c('Sample','Sample','Therapy step','Therapy step','Therapy step'),
      'Clinical.type'= c('1','2', '1','1','2'),
      'Clinical.value.start'= c(t_sampling_1, t_sampling_2, t_th_step_1_cycle_1_start, t_th_step_1_cycle_2_start, t_th_step_2_cycle_1_start),
      'Clinical.value.end'= c(NA, NA, t_th_step_1_cycle_1_end, t_th_step_1_cycle_2_end, t_th_step_2_cycle_1_end),
      'Clinical.index' = c(NA, NA, '1','2','1')
    ),
  'Mutations' =
    data.frame(
      'Mutation.name'= c('m_clock',rep('m_wgd', 12)),
      'Mutation.type'= c('primary',rep(c('alpha','beta'),6)),
      'Mutation.index'= c(NA,NA,NA,NA,NA,'1','1', '1','1', '2','2','2','2'),
      'Mutation.value'= c(
        m_clock_primary,
        alpha_tetraploid_clock, beta_tetraploid_clock,
        alpha_cnloh_clock, beta_cnloh_clock,
        alpha_tetraploid_step1, beta_tetraploid_step1,
        alpha_cnloh_step1, beta_cnloh_step1,
        alpha_tetraploid_step2, beta_tetraploid_step2,
        alpha_cnloh_step2, beta_cnloh_step2),
      'Mutation.coeff' = c(NA, rep(c('4','4','2','2'),3)),
      'Mutation.source' = c(rep('clock',5), rep('step',8))
    ),
  'Parameters' =
    data.frame('Parameter.name' = c('l_tetraploid', 'l_cnloh','l_diploid',
                                    'mu_th_step_alpha','mu_th_step_beta',
                                    'mu_th_step_alpha','mu_th_step_beta',
                                    'mu_clock','omega_alpha','omega_beta','k_step','N_min','N_max','N_min','N_max',
                                    'exponential_growth', 'mrca_alpha','mrca_beta'),
               'Parameter.value' = c(l_tetraploid, l_cnloh, l_diploid,
                                     mu_th_step1_alpha,mu_th_step1_beta,
                                     mu_th_step2_alpha,mu_th_step2_beta,
                                     mu_clock, omega_alpha, omega_beta, 10, N_min, N_max, N_min_sample_1, N_max_sample_1, 1,1,1),
               'Parameter.index' = c(NA, NA, NA,'1','1','2','2',rep(NA,11))
    )
)

#saveRDS(exampleData, file = "exampleData.rds")

usethis::use_data(exampleData_WGD, overwrite = TRUE)

