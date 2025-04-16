### Example data
library(dplyr)

## Timeline
t_eca = 3 # time of ECA (unknown)
t_obs_1 = 4 # first sampling time (known)
t_th_1_start = 5 # start of the therapy 1 (known)
t_cna = 14 # time of the CNA (unknown)
t_th_1_end = 16 # end of therapy 1 (known)
t_mrca = 20 # time of MRCA (unknown)
t_obs_2 = 22 # second sampling time (known)

## Rates
mu = 6e-8 # clock-like mut rate (known)
mu_th_1 = 50*mu # mutation rate under therapy 1 (unknown)
omega = 10 # growth rate (unknown)

# Major and minor allele
nA = 2
nB = 0 # other possibilities: 1,2
n= 2

N_min= 1e8 # exp(omega*(t_obs_2 - t_mrca)) - exp(omega*(t_obs_2 - t_mrca))*.5
N_max= 1e9 # exp(omega*(t_obs_2 - t_mrca)) + exp(omega*(t_obs_2 - t_mrca))*.5
k= 10000
diploid_length = sum(CNAqc:::get_reference('hg38')$length)
CNA_length = 2e7

m_clock = omega*mu*diploid_length*(t_mrca-t_eca)
m_alpha = omega*CNA_length*(mu*(t_cna - t_eca) +
                   mu_th_1*(t_cna - t_th_1_start))
m_beta = n*omega*CNA_length*(mu*(t_mrca-t_th_1_end) +
                    mu_th_1*(t_th_1_end-t_cna))
m_th_1 = omega*mu_th_1*diploid_length*(t_th_1_end - t_th_1_start)

alpha_mu_th_1 = mu_th_1*(t_th_1_end - t_th_1_start) * .05
beta_mu_th_1 = (t_th_1_end - t_th_1_start) * .05

# lb : omega giving rise to N_min if t_mrca = t_th_1_end
lb_omega= log(N_min) / (t_obs_2 - t_th_1_end)
# ub : omega giving rise to N_max if t_mrca = t_obs_2 - epsilon
epsilon = 1
ub_omega= log(N_max) / (t_obs_2 - (t_obs_2 - epsilon))

E = (lb_omega + ub_omega) / 2
V = 10
beta_omega = V / E
alpha_omega = E*beta_omega

convert_date_real = function(x) {
  y = 2000 + x %>% floor()
  py = (x - (x %>% floor())) * 365
  m = (py / 30) %>% floor
  d = (py - (m * 30)) %>% floor
  date_string = paste(y, m +1, d + 1, sep = '-')

  date_string = ifelse (d>30, paste(y, m +1, 30, sep = '-'), date_string)
  date_string = ifelse ((m==1 & d>27), paste(y, m +1, 28, sep = '-'), date_string)
  date_string = ifelse (m>=12, paste(y+1, 1, 1, sep = '-'), date_string)

  return(date_string)
}

convert_real_date = function(date = NULL, ref_year = 2000) {
  ref_month = 1
  ref_day = 1

  year = as.integer(strsplit(date, '-')[[1]][1])
  month = as.integer(strsplit(date, '-')[[1]][2])
  day = as.integer(strsplit(date, '-')[[1]][3])

  return((year - ref_year) + (month / 12 - ref_month / 12) + (day / 365 - ref_day / 365))
}

exampleData = list(
  'Clinical Timepoints' =
    data.frame('Timepoint' = c('Sample 1', 'Sample 2', 'Therapy 1'),
               'Start' = c(convert_date_real(t_obs_1), convert_date_real(t_obs_2), convert_date_real(t_th_1_start)),
               'End'= c(NA, NA, convert_date_real(t_th_1_end))
    ),
  'Mutations' =
    data.frame('Mutation type' = c('m_clock', 'm_alpha', 'm_beta', 'm_th_1'),
               'Number of mutations' = c(m_clock, m_alpha, m_beta, m_th_1)),
  'Parameters' =
    data.frame('Param name' = c('N_max','N_min','mu','omega_alpha','omega_beta',
                                'alpha_mu_th_1','beta_mu_th_1', 'k', 'diploid_length', 'CNA_length'),
               'Value' = c(N_max,N_min,mu,alpha_omega,beta_omega,alpha_mu_th_1,beta_mu_th_1, k,
                           diploid_length, CNA_length))
)

#saveRDS(exampleData, file = "exampleData.rds")

usethis::use_data(exampleData, overwrite = TRUE)
