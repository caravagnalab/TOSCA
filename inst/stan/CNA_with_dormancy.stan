functions {

  real lambda_therapy(real ti, real tf,real t_therapy_i,real t_therapy_f, real k){

    vector[2] a1;
    vector[2] a2;
    vector[2] a3;
    vector[2] a4;

    real f1;
    real f2;
    real f3;
    real f4;

    a1[1] = k*(tf-t_therapy_i);
    a1[2] = 1;
    a2[1] = k*(tf-t_therapy_f);
    a2[2] = 1;
    a3[1] = k*(ti-t_therapy_i);
    a3[2] = 1;
    a4[1] = k*(ti-t_therapy_f);
    a4[2] = 1;

    f1 = 1/k*(log(exp(a1[1]-max(a1)) + exp(-max(a1))) + max(a1));
    f2 = 1/k*(log(exp(a2[1]-max(a2)) + exp(-max(a2))) + max(a2));
    f3 = 1/k*(log(exp(a3[1]-max(a3)) + exp(-max(a3))) + max(a3));
    f4 = 1/k*(log(exp(a4[1]-max(a4)) + exp(-max(a4))) + max(a4));
    return f1 - f2 - f3 + f4;

  }

  real delta(real t,real dorm_s,real dorm_e,real k){
    return  dorm_e + (t-dorm_e)/(1+exp(k*(t-dorm_e))) - dorm_s;
    }

  real traslation(real t, real dorm_start, real dorm_end ) {
    real t_cn =t;
    if (t> dorm_start){
      t_cn= t_cn + dorm_end - dorm_start;
      }
      return t_cn;
	}
}

data{

 int <lower=0> m_clock_primary;
 int <lower=0> m_clock;
 real <lower=0> l_diploid;
 real <lower=0> mu_clock;

 int <lower=0> m_alpha;
 int <lower=0> m_beta;
 real <lower=0> l_CNA;
 int <lower=0> coeff_alpha; // 1 for 2:0, 2 for 2:2
 int <lower=0> coeff_beta; // 2 for 2:0, 4 for 2:2

 int <lower=0, upper=1> extra_therapy;
 real <lower=0> start_th_step;
 real <lower=0> end_th_step;
 real <lower=0> mu_th_step;
 real <lower=0> m_th_step;

	real <lower=0> omega_alpha;
  real <lower=0> omega_beta;
  real <lower=0> k_step;

  real <lower=0> Sample_1;
  real <lower=0> Sample_2;
  real <lower=0> chemo_start;
  real <lower=0> chemo_end;

  int <lower=0, upper=1> exponential_growth;
  array[2] real<lower=0> N_min;
  array[2] real<lower=0> N_max;

  real <lower=0> alpha_mrca;
  real <lower=0> beta_mrca;

  // Overdispersion parameters for all mutation types
  real <lower=0> phi_clock;
  real <lower=0> phi_th_step;
  real <lower=0> phi_cna;
}

parameters{

  real <lower=0, upper= Sample_1> t_eca;
  real <lower=t_eca, upper= Sample_1> t_mrca_primary;
  real <lower= chemo_start, upper= chemo_end> t_dormancy_start;
  real <lower= chemo_end, upper= Sample_2> t_dormancy_end;
  real <lower= t_dormancy_end, upper= Sample_2> t_mrca;
  real <lower= t_eca, upper= t_mrca - (t_dormancy_end - t_dormancy_start)> t_cna_tr;

  real <lower= 0> omega;
}


model{
  t_eca ~ uniform(0, Sample_1);
  t_mrca_primary ~ uniform(t_eca, Sample_1);
  t_dormancy_start~ uniform(chemo_start, chemo_end);
  t_dormancy_end~ uniform(chemo_end, Sample_2);
  t_mrca ~ uniform(t_dormancy, Sample_2);
  t_cna_tr ~ uniform(t_eca, t_mrca - (t_dormancy_end - t_dormancy_start) );

  omega ~ gamma(omega_alpha,omega_beta);

  real shape_clock = 1 / phi_clock;
  real shape_step = 1 / phi_th_step;
  real shape_cna = 1 / phi_cna;

  m_clock ~  neg_binomial_2_rng(
    2*mu_clock*omega*l_diploid*((t_dormancy_start - t_eca) + (t_mrca-t_dormancy_end)) +.1,
    shape_clock
  );

  if (extra_therapy){
    m_th_step  ~  neg_binomial_2_rng(
      2*mu_th_step*omega*l_diploid*lambda_therapy(t_dormancy_end,t_mrca, start_th_step, end_th_step,k_step) +.1,
      shape_step
      );

    m_alpha ~  neg_binomial_2_rng(coeff_alpha*omega*l_CNA*(mu_clock*(t_cna_tr - t_eca)
                                                         + mu_th_step*lambda_therapy(t_eca, t_cna_tr, start_th_step - delta(start_th_step, t_dormancy_start, t_dormancy_end, k_step),
                                                                                                      end_th_step - delta(end_th_step, t_dormancy_start, t_dormancy_end, k_step), k_step)),
                                  shape_cna);
    m_beta ~ neg_binomial_2_rng(coeff_beta*omega*l_CNA*(mu_clock*( t_mrca - (t_dormancy_end - t_dormancy_start) - t_cna_tr )
                                                      + mu_th_step*lambda_therapy(t_cna_tr,t_mrca-(t_dormancy_end - t_dormancy_start), start_th_step - delta(start_th_step,t_dormancy_start,t_dormancy_end,k_step),
                                                                                                                                       end_th_step - delta(end_th_step,t_dormancy_start,t_dormancy_end,k_step),k_step)),
                                shape_cna);

  }else{
    m_alpha ~  neg_binomial_2_rng(coeff_alpha*omega*l_CNA*(mu_clock*(t_cna_tr - t_eca),
                                  shape_cna);
    m_beta ~ neg_binomial_2_rng(coeff_beta*omega*l_CNA*(mu_clock*( t_mrca - (t_dormancy_end - t_dormancy_start) - t_cna_tr ),
                                shape_cna);

  }

  if (exponential_growth==1){
    target += -N_min[1]*exp(-omega*(Sample_1 - t_mrca_primary)) + log(1-exp(-(N_max[1]-N_min[1])*exp(-omega*(Sample_1 - t_mrca_primary))));
    target += -N_min[2]*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max[2]-N_min[2])*exp(-omega*(Sample_2 - t_mrca))));
  }
}

generated quantities{

    real t_copynumber;

    int m_clock =  neg_binomial_2_rng(
    2*mu_clock*omega*l_diploid*((t_dormancy_start - t_eca) + (t_mrca-t_dormancy_end)) +.1,
    shape_clock
  );

  if (extra_therapy){
    int m_th_step_rep =  neg_binomial_2_rng(
      2*mu_th_step*omega*l_diploid*lambda_therapy(t_dormancy_end,t_mrca, start_th_step, end_th_step,k_step) +.1,
      shape_step
      );

    int m_alpha_rep =  neg_binomial_2_rng(coeff_alpha*omega*l_CNA*(mu_clock*(t_cna_tr - t_eca)
                                                         + mu_th_step*lambda_therapy(t_eca, t_cna_tr, start_th_step - delta(start_th_step, t_dormancy_start, t_dormancy_end, k_step),
                                                                                                      end_th_step - delta(end_th_step, t_dormancy_start, t_dormancy_end, k_step), k_step)),
                                  shape_cna);
    int m_beta_rep =  neg_binomial_2_rng(coeff_beta*omega*l_CNA*(mu_clock*( t_mrca - (t_dormancy_end - t_dormancy_start) - t_cna_tr )
                                                      + mu_th_step*lambda_therapy(t_cna_tr,t_mrca-(t_dormancy_end - t_dormancy_start), start_th_step - delta(start_th_step,t_dormancy_start,t_dormancy_end,k_step),
                                                                                                                                       end_th_step - delta(end_th_step,t_dormancy_start,t_dormancy_end,k_step),k_step)),
                                shape_cna);

  }else{
    int m_alpha_rep =  neg_binomial_2_rng(coeff_alpha*omega*l_CNA*(mu_clock*(t_cna_tr - t_eca),
                                  shape_cna);
    int m_beta_rep =  neg_binomial_2_rng(coeff_beta*omega*l_CNA*(mu_clock*( t_mrca - (t_dormancy_end - t_dormancy_start) - t_cna_tr ),
                                shape_cna);

  }

    t_cna = traslation(t_cna_tr,  t_dormancy_start,  t_dormancy_end);

    if (exponential_growth==1){
    real N_sample_2 = exp(omega*(Sample_2-t_mrca));
    real N_sample_1 = exp(omega*(Sample_1-t_eca));
  }
}
