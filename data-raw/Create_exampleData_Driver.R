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
N_min1 = N1*.9
N_max1 = N1*1.1
N_min2 = N2*0.9
N_max2 = N2*1.1

exampleData_Driver = list(
  "Samples" = data.frame(
    "Name"=c("Diagnosis", "Relapse"),
    "Date"=c(TOSCA:::convert_date_real(Sample_1),
             TOSCA:::convert_date_real(Sample_2))),
  "Therapies" = data.frame(
    "Name"=c('Timolozomide','Timolozomide'),
    "Class"= c("Driver responsive","Driver responsive"),
    "Start"=c(TOSCA:::convert_date_real(step_1_start_cycle_1),
              TOSCA:::convert_date_real(step_1_start_cycle_2)),
    "End"=c(TOSCA:::convert_date_real(step_1_end_cycle_1),
            TOSCA:::convert_date_real(step_1_end_cycle_2))),

  "Mutations" = data.frame(
    "Name" = c("sbs1_primary", "sbs1_relapse", "sbs11"),
    "Length" = c(l_diploid, l_diploid, l_diploid),
    "Karyptype" = c("1:1","1:1", "1:1"),
    "Type"=c("clock-like primary", "clock-like relapse", "driver"),
    "Value"=c(m_clock_primary, m_clock, m_driver)
  ),

  'Parameters' =
    data.frame('Name' = c('mu_clock',
                          'omega_alpha','omega_beta','k_step','N_min','N_max', 'N_min','N_max' ,'exponential_growth', 'mrca_alpha','mrca_beta',
                          "mu_driver_clock", "mu_driver", "phi_clock","phi_driver"),
               'Value' = c(mu_clock, omega_alpha, omega_beta, 10, N_min1, N_max1, N_min2, N_max2, 1,1,1, mu_clock_driver, mu_driver,.1,.1),
               'Index' = c(rep(NA, 4), "Diagnosis", "Diagnosis", "Relapse", "Relapse", rep(NA, 7))
    ))




usethis::use_data(exampleData_Driver, overwrite = TRUE)





