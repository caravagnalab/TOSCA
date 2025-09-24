# Timeline
birth = 0
t_eca = 4
t_mrca_primary = 4.5
SAMPLE_1 = 6
TH1_START = 6.1
TH1_END = 6.5
CHEMO_START = 6.6
t_dormancy_start = CHEMO_START + 0.01
CHEMO_END = 6.8
TH2_START = 10
t_dormancy_end = 10.2
t_cna = 10.5
TH2_END = 11.2
t_mrca = 11.3
SAMPLE_2 = 13


# Parameters
omega = 5
omega_alpha = 50
omega_beta = 10
# hist(rgamma(10000,omega_alpha,omega_beta))

l_diploid = 3e9
l_CNA1 = 4.5e8
mu_clock = 1e-8
mu_th_step_1 = 5e-8
mu_th_step_2 = 1e-7

mu_th_step1_alpha = 10
mu_th_step1_beta = 2e8
hist(rgamma(10000,mu_th_step1_alpha, mu_th_step1_beta))

mu_th_step2_alpha = 3
mu_th_step2_beta = 2e7
hist(rgamma(10000,mu_th_step2_alpha, mu_th_step2_beta))


m_clock = 2*omega*l_diploid*mu_clock*( (t_dormancy_start - t_eca) + (t_mrca - t_dormancy_end))
m_clock_primary = 2*omega*l_diploid*mu_clock*(t_mrca_primary-t_eca)

coeff_cna_1 = 2
m_alpha_cna_1 = omega*l_CNA1*(
  mu_clock * ( (TH1_START - t_eca) + (t_dormancy_start-TH1_END)) +
  mu_th_step_1 * (TH1_END-TH1_START) +
  mu_th_step_2 * (t_cna-t_dormancy_end)
)

m_beta_cna_1 = coeff_cna_1*omega*l_CNA1*(
  mu_clock * (t_mrca - t_cna) +
  mu_th_step_2 * (TH2_END-t_cna)
)


m_step_1 = 2*omega*l_diploid*mu_th_step_1*(TH1_END-TH1_START)
m_step_2 = 2*omega*l_diploid*mu_th_step_2*(TH2_END-TH2_START)

N = exp(omega*(SAMPLE_2-t_mrca))
N_min = N*.9
N_max = N*1.1

N_sample_1 = exp(omega*(SAMPLE_1-t_mrca_primary))
N_min_sample_1 = N_sample_1*.9
N_max_sample_1 = N_sample_1*1.1


real_to_date = function(x, date = NULL, ref_year = 2000) {

  y = ref_year + x %>% floor()
  py = (x - (x %>% floor())) * 365
  m = (py / 30) %>% floor
  d = (py - (m * 30)) %>% floor
  date_string = paste(y, m +1, d + 1, sep = '-')

  date_string = ifelse (d>30, paste(y, m +1, 30, sep = '-'), date_string)
  date_string = ifelse ((m==1 & d>27), paste(y, m +1, 28, sep = '-'), date_string)
  date_string = ifelse (m>=12, paste(y+1, 1, 1, sep = '-'), date_string)

  return(date_string)
}

exampleData_CNA_dormancy = list(
  "Samples" = data.frame(
    "Name"=c("Birth","Diagnosis", "Relapse"),
    "Date"=c(
      real_to_date(birth, ref_year = 2000),
      real_to_date(SAMPLE_1),
      real_to_date(SAMPLE_2))),
  "Therapies" = data.frame(
    "Name"=c("Drug 1", "Chemo", "Drug 2"),
    "Class"= c("Mutagenic","Chemotherapy inducing dormancy","Mutagenic"),
    "Start"=c(real_to_date(TH1_START),
                real_to_date(CHEMO_START),
              real_to_date(TH2_START)),
    "End"=c(real_to_date(TH1_END),
            real_to_date(CHEMO_END),
            real_to_date(TH2_END))),
  "Mutations" = data.frame(
    "Name" = c("m_primary", "m_relapse", "m_pre", "m_post", "m_drug1", "m_drug2"),
    "Length" = c(l_diploid, l_diploid, l_CNA1, l_CNA1, l_diploid, l_diploid),
    "Karyotype" = c("1:1","1:1", "2:0", "2:0", "1:1", "1:1"),
    "Type"=c("clock-like primary", "clock-like relapse", "alpha", "beta", "Drug 1", "Drug 2"),
    "Value"=c(m_clock_primary, m_clock, m_alpha_cna_1, m_beta_cna_1, m_step_1, m_step_2)
  ),
  'Parameters' =
    data.frame('Name' = c('mu_clock','omega_alpha','omega_beta','k_step',
                          'N_min','N_max','N_min','N_max',
                          'exponential_growth', 'exponential_growth', 'mrca_alpha','mrca_beta',
                          "alpha_th_step","beta_th_step", "alpha_th_step","beta_th_step",
                          "phi_clock",rep("phi_th_step",2),
                          "phi_cna_alpha", "phi_cna_beta") # rep("phi_cna",2)
               ,
    'Value' = c(mu_clock, omega_alpha, omega_beta, 10,
                N_min, N_max, N_min_sample_1, N_max_sample_1,
                1,1, 1,1,
                mu_th_step1_alpha, mu_th_step1_beta, mu_th_step2_alpha, mu_th_step2_beta,
                rep(0.05, 5)),
    'Index' = c(rep(NA,4), "Relapse", "Relapse", "Diagnosis","Diagnosis", "Diagnosis", "Relapse",NA, NA,
                "Drug 1","Drug 1", "Drug 2","Drug 2",NA, "Drug 1", "Drug 2",
                l_CNA1, l_CNA1)
    )
)


usethis::use_data(exampleData_CNA_dormancy, overwrite = TRUE)





