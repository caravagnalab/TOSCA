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

	real couchy(real location, real scale, real a,real b){
    real d= ((1/pi()) * atan((b-location)/scale) + .5) - ((1/pi()) * atan((a-location)/scale) + .5);
  return(d);
	}
}

data{

  // Clock-like mutations
  int <lower=0> m_clock_primary;
  // int <lower=0> m_clock;
  real <lower=0> l_diploid;
  real <lower=0> mu_clock;

	// mutations on CNA
  // int <lower=0> n_cna;
  // array[n_cna] int<lower=0> m_alpha;
  // array[n_cna] int<lower=0> m_beta;
  // array[n_cna] real<lower=0> l_wgd;
  // array[n_cna] int<lower=0> coeff_alpha; // 1 for 2:0, 2 for 2:2
  // array[n_cna] int<lower=0> coeff_beta; // 2 for 2:0, 4 for 2:2
  //
  // real<lower=0> start_th_step;
  // real<lower=0> end_th_step;
  // int <lower=0> m_th_step;
  // int <lower=0> mu_th_step;

  int <lower=0> m_alpha;
  int <lower=0> m_beta;
  int <lower=0> l_wgd;

	// other parameters
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
}

parameters{

  real <lower=0, upper= Sample_1> t_eca;
  real <lower=t_eca, upper= Sample_1> t_mrca_primary;
  real <lower= chemo_start, upper= chemo_end> t_dormancy_start;
  real <lower= t_dormancy_start, upper= Sample_2> t_dormancy_end;
  real <lower= t_dormancy_end, upper= Sample_2> t_mrca;
  real <lower=t_eca, upper=t_mrca - (t_dormancy_end - t_dormancy_start)> t_wgd_tr;

  real <lower= 0> omega;

}

model{

  t_eca ~ uniform(0, Sample_1);
  t_mrca_primary ~ uniform(t_eca, Sample_1);
  t_dormancy_start ~ uniform(chemo_start, chemo_end);
  t_dormancy_end ~ uniform(t_dormancy_start, Sample_2);
  t_mrca ~ uniform(t_dormancy_end, Sample_2);
  t_wgd_tr ~ uniform(t_eca, t_mrca - (t_dormancy_end - t_dormancy_start));

  omega ~ gamma(omega_alpha,omega_beta);

  m_clock_primary ~ poisson(2*mu_clock*omega*l_diploid*(t_mrca_primary-t_eca));

  m_alpha ~  poisson(2*omega*l_wgd*mu_clock*(t_wgd_tr-t_eca));
  m_beta ~  poisson(4*omega*l_wgd*mu_clock*( (t_mrca - (t_dormancy_end - t_dormancy_start) ) - t_wgd_tr));

  if (exponential_growth==1){
    target += -N_min[1]*exp(-omega*(Sample_1 - t_mrca_primary)) + log(1-exp(-(N_max[1]-N_min[1])*exp(-omega*(Sample_1 - t_mrca_primary))));
    target += -N_min[2]*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max[2]-N_min[2])*exp(-omega*(Sample_2 - t_mrca))));
  }
}

generated quantities{

  int m_clock_primary_rep = poisson_rng(2*mu_clock*omega*l_diploid*(t_mrca_primary-t_eca));
  int m_alpha_rep =  poisson_rng(2*omega*l_wgd*mu_clock*(t_wgd_tr-t_eca));
  int m_beta_rep =  poisson_rng(4*omega*l_wgd*mu_clock*( (t_mrca - (t_dormancy_end - t_dormancy_start) ) - t_wgd_tr));

  real t_wgd = traslation(t_wgd_tr,  t_dormancy_start,  t_dormancy_end );

  if (exponential_growth==1){
    real N_sample_2 = exp(omega*(Sample_2-t_mrca));
    real N_sample_1 = exp(omega*(Sample_1-t_eca));
  }

}
