functions{
//This is the jth BP basis of degree M; A transformation to the argument t to (t+1)/2 is done so that the support is on [-1,1]; 0.5*(M+1) is normalizing constant  
vector Btilde(int M, int j, int N, vector t){
  vector[N] a;
  for(jj in 1:N){
    a[jj] = (0.5/(M+1))*exp(beta_lpdf((t[jj]+1)/2 |j+1, M-j+1));
  }
  return(a);
}

// Berstein (of degree M) basis functions evaluated at vector t; return matrix whose columns evaluate each basis (BP of M+1 basis) at the input vector t
matrix Balpha(int M, int N, vector t) {
  matrix[N,M+1] bb; 
    for(j in 1:(M+1)){
      bb[,j] = Btilde(M,j-1,N,t);
    } 
  return(bb);
}
}

data {
  int<lower=0> nsubs;                         // total no. of individuals
  // int<lower=0> ntooth;                        // maximum no. of teeth for an individual
  int<lower=0> nfixef;                        // total no. of fixed effects
  int<lower=0> ngroup;                        // total no. of clusters
  int<lower=0> nobs;                          // total no. of observations (uncensored)
  matrix[nobs, nfixef] xmat;                  // fixed effects design matrix
  array[nobs] int group;                            // cluster ids
  array[nobs] int sub;                              // subject ids
  // int tooth[nobs];                            // tooth ids
  vector[nobs] yvec;                          // observations
  
  int<lower=0> ncens;                         // Number of censored data
  matrix[ncens,nfixef] xmat_cens;             // censored data's design matrix
  array[ncens] int group_cens;                      // censored data's group index
  array[ncens] int sub_cens;                        // subject ids
  // int tooth_cens[ncens];                      // tooth ids
  vector[ncens] yvec_censor;                  // censored observations

  real<lower=0> scale_icept;                  // prior std for the intercept

  int<lower=0> M;                            //degree of BP basis
  matrix[M+1,M+1] A;                         //matrix that converts first order differences (as well as the initial element) to the original vector
}

transformed data {
  // both fix eff. and rand. eff. are apriori centered at 0; the GPD distribution for standard deviation of fix eff. and rand. eff are apriori centered at 0
  vector[nfixef] meanGroup;
  real mean_GPD_group;                       // mean par for GPD

  // u[1] = -1;
  // for(ss in 2:(L+1))
  //   u[ss] = u[ss-1] + 2.0/L;
  meanGroup = rep_vector(0.0, nfixef);  
  mean_GPD_group = 0;
}

parameters {
  real intercept;                           // overall intercept (outside the single index function)
  matrix[nfixef,ngroup] fixef_group;       // fix. eff. group deviation
  vector<lower=0>[nfixef] group_sd_local;  // local sd in GLT shrinkage for fix. eff. group 
  // real intercept[T];                    // time specific intercept
  real<lower=0> scale_GPD_group;                // scale par for GPD
  real<lower=0.5> shape_GPD_group;          // shape par for GPD
  real mu_xi;                               // mean for shape par for GPD
  real<lower=0> sigma_xi;                   // sd for shape par for GPD

  real shape_GEV;                          // shape par for GEV
  // real scale_GEV_raw;             // (log) scale par for GEV; constrained in transformed parameter section 
  real scale_GEV;
  
  vector[nsubs] ranef;                  // subject level rand. eff.
  real<lower=0> stdErrRanef;           // std err. in subject level. rand. eff.

  vector[M] alpha;                 // first order differences of BP basis coefficients; 

}

transformed parameters {
  vector[ncens] yHat_cens;
  vector[nobs] yHat;
  // vector[ntooth] intercept_0;                           // tooth level intercept (zero centered)
  matrix[nfixef,ngroup] sc_fix;                            // scaled fixef. plus. group deviation
  matrix[nobs,ngroup] si_linear;        // single index matrix (with different group coefficients) for the failure time data
  matrix[ncens,ngroup] si_linear_cens; // single index matrix for the censored time data
  vector[nobs] si_linear_select;       // single index w.r.t. selected group for the failure time data
  vector[ncens] si_linear_select_cens; // single index w.r.t. selected group for the censored time data

  for(gg in 1:ngroup){
       sc_fix[,gg] =  fixef_group[,gg]/norm2(fixef_group[,gg]);
       si_linear[,gg] = xmat*sc_fix[,gg];
       si_linear_cens[,gg] = xmat_cens*sc_fix[,gg];
  }

  for(ii in 1:nobs){
    si_linear_select[ii] = si_linear[ii,group[ii]];
    yHat[ii] = ranef[sub[ii]];
  }

  for(i in 1:ncens){
    si_linear_select_cens[i] = si_linear_cens[i,group_cens[i]];
    yHat_cens[i] = ranef[sub_cens[i]];
  }

  yHat += Balpha(M,nobs,si_linear_select)*(A*append_row(0,abs(alpha))) + intercept;     //the first coefficient is constraint as 0 because it is essentially intercept
  yHat_cens += Balpha(M,ncens,si_linear_select_cens)*(A*append_row(0,abs(alpha))) + intercept;

}

model {
  // vector[ncens] yHat_cens;
  // vector[nobs] yHat;
  // // vector[ntooth] intercept_0;                           // tooth level intercept (zero centered)
  // matrix[nfixef,ngroup] sc_fix;                            // scaled fixef. plus. group deviation
  // matrix[nobs,ngroup] si_linear;        // single index matrix (with different group coefficients) for the failure time data
  // matrix[ncens,ngroup] si_linear_cens; // single index matrix for the censored time data
  // vector[nobs] si_linear_select;       // single index w.r.t. selected group for the failure time data
  // vector[ncens] si_linear_select_cens; // single index w.r.t. selected group for the censored time data
  
  // sample overall intercept (not used for DP model)
  intercept ~ normal(0.0, scale_icept);

  // sample rand. eff.  
  ranef ~ normal(0.0,stdErrRanef);   //ranef is a vector
    
  // sample mSIM coefficients
  alpha ~ normal(0.0,10.0);

  // sample the shape parameter of GPD for group local sd 
  sigma_xi ~ normal(0.0,10.0);
  
  // sample fix. eff. and GLT parameters
  scale_GPD_group ~ inv_gamma(nfixef/shape_GPD_group+1.0, 1.0);
  shape_GPD_group ~ lognormal(mu_xi, sigma_xi); // T[0.5,];
  
  for(gg in 1:ngroup)
    fixef_group[,gg] ~ normal(meanGroup, group_sd_local);
  target += -nfixef*log(scale_GPD_group) - (shape_GPD_group+1)/shape_GPD_group*sum(log1p(shape_GPD_group*(group_sd_local-mean_GPD_group)*pow(scale_GPD_group,-1)));

  // for(gg in 1:ngroup){
  //      sc_fix[,gg] =  fixef_group[,gg]/norm2(fixef_group[,gg]);
  //      si_linear[,gg] = xmat*sc_fix[,gg];
  //      si_linear_cens[,gg] = xmat_cens*sc_fix[,gg];
  // }
  // 
  // for(ii in 1:nobs){
  //   si_linear_select[ii] = si_linear[ii,group[ii]];    
  //   yHat[ii] = ranef[sub[ii]];
  // }
  // 
  // for(i in 1:ncens){
  //   si_linear_select_cens[i] = si_linear_cens[i,group_cens[i]];
  //   yHat_cens[i] = ranef[sub_cens[i]];
  // }
  // 
  // yHat += Balpha(M,nobs,si_linear_select)*(A*append_row(0,abs(alpha))) + intercept;     //the first coefficient is constraint as 0 because it is essentially intercept
  // yHat_cens += Balpha(M,ncens,si_linear_select_cens)*(A*append_row(0,abs(alpha))) + intercept;
  
  // sample shape par for GEV 
  shape_GEV ~ uniform(-0.5, 0.5);
  if(max(append_row(yHat-yvec, yHat_cens-yvec_censor)*shape_GEV) > 0)
    scale_GEV ~ normal(0.0, 10.0) T[log(max(append_row(yHat-yvec, yHat_cens-yvec_censor)*shape_GEV)),];
  else     
    scale_GEV ~ normal(0.0, 10.0);

  target += -nobs*scale_GEV - (shape_GEV+1)/shape_GEV*sum(log1p(shape_GEV*(yvec-yHat)*exp(-scale_GEV)));
  for(ii in 1:nobs)
     target += - pow(1.0 + shape_GEV*(yvec[ii]-yHat[ii])*exp(-scale_GEV),-inv(shape_GEV));

  for(ii in 1:ncens)
    target += log1m_exp(-pow(1.0+shape_GEV*(yvec_censor[ii]-yHat_cens[ii])*exp(-scale_GEV), -inv(shape_GEV)));  //survival function for GEV distribution
}

generated quantities{
  vector[nsubs] log_lik;
  // vector[ncens] yHat_cens;
  // vector[nobs] yHat;
  // matrix[nfixef,ngroup] sc_fix;                            // scaled fixef. plus. group deviation
  // matrix[nobs,ngroup] si_linear;        // single index matrix (with different group coefficients) for the failure time data
  // matrix[ncens,ngroup] si_linear_cens; // single index matrix for the censored time data
  // vector[nobs] si_linear_select;       // single index w.r.t. selected group for the failure time data
  // vector[ncens] si_linear_select_cens; // single index w.r.t. selected group for the censored time data


  // for(gg in 1:ngroup){
  //      sc_fix[,gg] =  fixef_group[,gg]/norm2(fixef_group[,gg]);
  //      si_linear[,gg] = xmat*sc_fix[,gg];
  //      si_linear_cens[,gg] = xmat_cens*sc_fix[,gg];
  // }
  // 
  // for(ii in 1:nobs){
  //   si_linear_select[ii] = si_linear[ii,group[ii]];    
  //   yHat[ii] = ranef[sub[ii]];
  // }
  // 
  // for(i in 1:ncens){
  //   si_linear_select_cens[i] = si_linear_cens[i,group_cens[i]];
  //   yHat_cens[i] = ranef[sub_cens[i]];
  // }
  // 
  // yHat += Balpha(M,nobs,si_linear_select)*(A*append_row(0,abs(alpha))) + intercept;     //the first coefficient is constraint as 0 because it is essentially intercept
  // yHat_cens += Balpha(M,ncens,si_linear_select_cens)*(A*append_row(0,abs(alpha))) + intercept;

  log_lik = rep_vector(0.0, nsubs);
  for(ii in 1:nobs){
    log_lik[sub[ii]] += - pow(1.0 + shape_GEV*(yvec[ii]-yHat[ii])*exp(-scale_GEV),-inv(shape_GEV));
  }
  for(ii in 1:ncens){
    log_lik[sub_cens[ii]] += log1m_exp(-pow(1.0+shape_GEV*(yvec_censor[ii]-yHat_cens[ii])*exp(-scale_GEV), -inv(shape_GEV)));
  }

}
