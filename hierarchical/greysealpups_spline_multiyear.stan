// Grey seal pup production model: multi-colony, multi-year hierarchical gam
// 
// Year-specific deviations using mgcv factor smooth: s(Year_f, ClineDist, bs="fs")
// Colony random effects for birth date
// Estimated leave process (same across all colonies and years)
// Fixed survival and carcass persistence parameters (same across all colonies and years)
// Penalized spline priors using precision matrices from jagam

functions {
  vector skewnormBDFunction(real mu, real sigma, real alpha, int nsteps, vector days){
    vector[nsteps] xvec;
    real omega_skn = sigma/sqrt(1-(2/pi())*( square(alpha)/(1+square(alpha)) ) );
    real xi_skn = mu - omega_skn * (alpha/(sqrt(1+square(alpha)))) * sqrt(2/pi());
    for (i in 1:nsteps) xvec[i] = exp(skew_normal_lpdf(days[i] | xi_skn, omega_skn, alpha));
    return xvec / sum(xvec);
  }
  
  vector hazFunction(real mu, real std, int nsteps, int dstep) {
    real prob;
    vector[nsteps] hazprob;
    hazprob[1] = normal_cdf(1 * dstep | mu, std);
    for (s in 2:nsteps) {
      prob = normal_cdf(s * dstep | mu, std) - normal_cdf((s - 1) * dstep | mu, std);
      hazprob[s] = (1 - normal_cdf((s - 1) * dstep | mu, std)) < 1e-6
                     ? 1.0
                     : prob / (1 - normal_cdf((s - 1) * dstep | mu, std));
    }
    return hazprob;
  }
  
  matrix transMatrix(real xphi, real xgamma, real xzeta, real xchi, int stages, int dstep) {
    matrix[stages, stages] matarray = rep_matrix(0, stages, stages);
    matarray[1, 1] = pow(xphi, dstep) * (1 - xgamma) * (1 - xzeta);
    matarray[2, 1] = pow(xphi, dstep) * xgamma * (1 - xzeta);
    matarray[4, 1] = 1 - pow(xphi, dstep);
    matarray[2, 2] = pow(xphi, dstep) * (1 - xzeta);
    matarray[3, 2] = pow(xphi, dstep) * xzeta;
    matarray[4, 2] = 1 - pow(xphi, dstep);
    matarray[4, 4] = pow(xchi, dstep);
    matarray[5, 4] = 1 - pow(xchi, dstep);
    matarray[3, 3] = 1;
    matarray[5, 5] = 1;
    return matarray;
  }
}

data {
  int<lower=1> nsurveys;
  int<lower=1> nstates;
  int<lower=1> ndays;
  int<lower=1> nobsstates;
  int<lower=1> nvisits;
  int<lower=1> ncols;
  int<lower=1> nyears;
  
  int<lower=1> K;                     // spline basis dimension (overall + factor)
  matrix[nsurveys, K] X_spline;       // full design matrix from jagam
  matrix[K, K] S1_block;              // penalty matrix from jagam
  matrix[K, K] S2_block;             
  
  array[nsurveys] int<lower=1, upper=ncols> colid_svy;
  array[nsurveys] int<lower=1, upper=nyears> year_svy;
  array[nsurveys] int<lower=0> nvisits_svy;
  
  array[nvisits] int obs_vector_w;
  array[nvisits] int obs_vector_m;
  array[nvisits] int obs_vector_d;
  array[nvisits] int<lower=0> obsdays_svy;
  
  int<lower=0, upper=1> est_pobs_d;
  int<lower=0, upper=1> est_pobs_w;
  int<lower=0, upper=1> est_pobs_m;
  int<lower=0, upper=1> est_pcm;
}

transformed data {
  int daysPerBin = 5;
  int nIntervals_s = ndays %/% daysPerBin;
  array[nvisits] int obsdays_bin;
  for (v in 1:nvisits)
    obsdays_bin[v] = min(obsdays_svy[v] %/% daysPerBin + 1, nIntervals_s);
  
  int nIntervals_c = min(20, nIntervals_s);
  vector[nIntervals_c] cohort_days;
  for (z in 1:nIntervals_c)
    cohort_days[z] = z * daysPerBin - 0.5 * daysPerBin;
  
  real xchi = inv_logit(2.5);
  real muMoult = 23;
  real sdMoult = 5;
  vector[nIntervals_s] xgamma = hazFunction(muMoult, sdMoult, nIntervals_s, daysPerBin);
}

parameters {
  // spline coefficients (overall + year-specific)
  vector[K] beta_spline_raw;
  vector<lower=1e-8>[2] lambda_spline;  // smoothing penalties
  
  vector[nsurveys] nborn_col_raw;
  
  real muBday_global;
  real<lower=0> sigma_col;
  vector[ncols] z_col;
  
  // Birth distribution parameters
  real<lower=0> sdBday_raw;
  real<lower=0> alphaBday_raw;
  
  // Leave distribution
  real muLeave_raw;
  real<lower=0> sdLeave_raw;
  
  // Observation probabilities
  real<lower=0.7, upper=1> pobs_d_est;
  real<lower=0.7, upper=1> pobs_w_est;
  real<lower=0.7, upper=1> pobs_m_est;
  real<lower=0.7, upper=1> pcm_est;
}

transformed parameters {
  vector[nsurveys] nborn_col = exp(5.0 + 2.0 * nborn_col_raw);
  
  // Spline effect for each survey 
  vector[nsurveys] beta_spline_mu_s = X_spline * beta_spline_raw;

  // Colony random effects
  vector[ncols] col_eff_raw = sigma_col * z_col;
  vector[nsurveys] col_eff;
  for (s in 1:nsurveys)
    col_eff[s] = col_eff_raw[colid_svy[s]];
 
  real sdBday = sdBday_raw * 6;
  real alphaBday = alphaBday_raw * 6;
  
  real muLeave = exp(3.5 + 0.1 * muLeave_raw);
  real sdLeave = 1 + 5.0 * sdLeave_raw;
  vector[nIntervals_s] xzeta = hazFunction(muLeave, sdLeave, nIntervals_s, daysPerBin);
  
  // Fixed survival (same for all colonies and years)
  vector[nsurveys] survInt_svy = rep_vector(4.0, nsurveys);
  real survSlope = 0.4;
  
  real pobs_d = est_pobs_d == 1 ? pobs_d_est : 0.85;
  real pobs_w = est_pobs_w == 1 ? pobs_w_est : 0.98;
  real pobs_m = est_pobs_m == 1 ? pobs_m_est : 0.91;
  real pcm    = est_pcm == 1 ? pcm_est : 0.84;
}

model {
  muBday_global ~ normal(45, 10);

  // Smoothing penalties
  lambda_spline[1] ~ cauchy(0, 1);  
  lambda_spline[2] ~ cauchy(0, 1);
  
  // Build precision matrix from jagam penalties
  matrix[K, K] K1;
  {
    matrix[K, K] S_total = lambda_spline[1] * S1_block + lambda_spline[2] * S2_block;
    K1 = S_total + diag_matrix(rep_vector(0.01, K));  // ridge 
  }

  // Penalized spline prior
  beta_spline_raw ~ multi_normal_prec(rep_vector(0, K), K1);

  nborn_col_raw ~ std_normal();
  z_col ~ std_normal();
  sigma_col ~ normal(0, 1);

  sdBday_raw ~ gamma(0.1, 0.1);
  alphaBday_raw ~ std_normal();
  
  muLeave_raw ~ std_normal();
  sdLeave_raw ~ gamma(1, 1);
  
  
  // Daily survival rate
  array[nsurveys] vector[ndays] xphi;
  for (surv in 1:nsurveys) {
    for (a in 1:ndays) {
      xphi[surv][a] = inv_logit(survInt_svy[surv] + survSlope * (a - 1));
    }
  }

  // Observation probabilities
  if (est_pobs_d == 1) pobs_d_est ~ beta(35, 5);
  if (est_pobs_w == 1) pobs_w_est ~ beta(45, 5);
  if (est_pobs_m == 1) pobs_m_est ~ beta(45, 5);
  if (est_pcm == 1)    pcm_est ~ beta(35, 5);
  
  // Observation matrix
  matrix[nstates, nobsstates + 1] obsMatrix = rep_matrix(0, nstates, nobsstates + 1);
  obsMatrix[1, 1] = pobs_w;
  obsMatrix[1, 4] = 1 - pobs_w;
  obsMatrix[2, 1] = pobs_m * (1 - pcm);
  obsMatrix[2, 2] = pobs_m * pcm;
  obsMatrix[2, 4] = 1 - pobs_m;
  obsMatrix[3, 4] = 1;
  obsMatrix[4, 3] = pobs_d;
  obsMatrix[4, 4] = 1 - pobs_d;
  obsMatrix[5, 4] = 1;
  
  // Likelihood
  int pos = 1;
  for (s in 1:nsurveys) {  # each survey is unique colony-year combination
    // Birth date mean for this survey (global + colony effect + year-specific spline)
    vector[nIntervals_c] Pb = skewnormBDFunction(
      muBday_global + col_eff[s] + beta_spline_mu_s[s],
      sdBday, alphaBday, nIntervals_c, cohort_days);
    
    vector[nIntervals_c] cohortSizeBin = nborn_col[s] * Pb;
    
    array[nIntervals_s] matrix[nstates, nstates] matarray;
    for (i in 1:nIntervals_s) {
      matarray[i] = transMatrix(
        xphi[s][(i - 1) * daysPerBin + 1], xgamma[i], xzeta[i], xchi, nstates, daysPerBin);
    }
    
    array[nIntervals_s] vector[nstates] stageN;
    for (a in 1:nIntervals_s) stageN[a] = rep_vector(0, nstates);
    
    for (h in 1:nIntervals_c) {
      vector[nstates] tempCohortVec = rep_vector(0, nstates);
      tempCohortVec[1] = cohortSizeBin[h];
      stageN[h] += tempCohortVec;
      for (a in 1:(nIntervals_s - h)) {
        tempCohortVec = matarray[a] * tempCohortVec;
        stageN[h + a] += tempCohortVec;
      }
    }
    
    for (v in 1:nvisits_svy[s]) {
      row_vector[nobsstates + 1] lambda_obs_vector = to_row_vector(stageN[obsdays_bin[pos]]) * obsMatrix;
      obs_vector_w[pos] ~ poisson(lambda_obs_vector[1]);
      obs_vector_m[pos] ~ poisson(lambda_obs_vector[2]);
      obs_vector_d[pos] ~ poisson(lambda_obs_vector[3]);
      pos += 1;
    }
  }
}
