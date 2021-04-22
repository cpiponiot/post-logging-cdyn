
// defining data
data {
  // indices
  int<lower=0> N;                 // nb of observations with small trees 
  int<lower=0> P;                 // nb of plots
  int<lower=0> L;                 // number of logged plots
  int<lower=0> t[N];              // time since logging
  int<lower=1,upper=P> np[N];     // plot each observation belongs to  
  int<lower=1> logged[L];         // logged plots ID
  
  // response variables
  real<lower=-1> cGrowth[N];          // cumulative volume gain
  real<lower=-1> cMort[N];          // cumulative volume mortality (negative values)
  real<lower=0> Vol[N];             // total volume
  real deltaV[L];                 // logging-induced volume loss
  real maxV[P];
  
  // parameters bounds
  real lowBound[9]; 
  real upBound[9]; 
  real supSigma[5];
}

// defining model parameters
parameters {
  
  real<lower=lowBound[2], upper=upBound[2]> aM;                    // asymptotic volume mortality
  real<lower=lowBound[5],upper=upBound[5]> theta;                 // volume respiration rate
  real<lower=lowBound[8],upper=upBound[8]> vmax;
  real<lower=lowBound[3],upper=upBound[3]> bG;                // volume productivity increase rate
  real<lower=lowBound[4],upper=bG*(vmax*theta/aM+1)> bM;                  // volume mortality increase rate
  real<lower=lowBound[6],upper=upBound[6]> t0;                   // pre-logging maturity
  real<lower=lowBound[7],upper=t0> tlog[L];                // post-logging maturity
  // standard deviations
  real<lower=0,upper=supSigma[1]> sigma_G;              // cumulative volume growth
  real<lower=0,upper=supSigma[2]> sigma_M;              // cumulative volume mortality
  real<lower=0,upper=supSigma[3]> sigma_V;              // total volume 
  real<lower=0,upper=supSigma[4]> sigma_deltaV;         // logging-induced volume loss
  
}

transformed parameters{
  real<lower=lowBound[1], upper=upBound[1]> aG[P];
  real t1[P];
  
  for (i in 1:P) { 
    aG[i] = aM + theta*vmax;
  }
    
 // initial maturity: t0 for control plots, tlog for logged plots
  for (i in 1:P) {t1[i] = t0;} 
  for (l in 1:L) {t1[logged[l]] = tlog[l];} 
  
}

// the model
model{
  real muG[N];
  real muM[N];
  real muV[N];
  real mudeltaV[L];

  // calculate volume dynamics predictions knowing the parameters values
  for (i in 1:N)
  { 
    muG[i] = aG[np[i]]*bG/(theta-bG) * ( (exp(-bG*t1[np[i]])-exp(-bG*(t[i]+t1[np[i]])))/bG - (exp(-theta*t1[np[i]])-exp(-theta*(t[i]+t1[np[i]])))/theta) + aM*(t[i] - (theta/bM*(exp(-bM*t1[np[i]])-exp(-bM*(t[i]+t1[np[i]]))) -  bM/theta*(exp(-theta*t1[np[i]])-exp(-theta*(t[i]+t1[np[i]]))))/(theta-bM)) ;
    
    muM[i] = aM * ( t[i] - (exp(-bM*t1[np[i]])-exp(-bM*(t[i]+t1[np[i]])))/bM ) ;
    
    muV[i] = aG[np[i]]/theta * (1 - (theta*exp(-bG*(t[i]+t1[np[i]])) - bG*exp(-theta*(t[i]+t1[np[i]])))/(theta-bG) ) - aM/theta * (1 - (theta*exp(-bM*(t[i]+t1[np[i]])) - bM*exp(-theta*(t[i]+t1[np[i]])))/(theta-bM) ) ; 
    
    // define model likelihoods
    target += normal_lpdf(log(cGrowth[i]+1) | log(muG[i]+1), sigma_G);
    target += normal_lpdf(log(cMort[i]+1) | log(muM[i]+1), sigma_M);
    target += normal_lpdf(log(Vol[i]+1)   | log(muV[i]+1), sigma_V);
  }
  
  // logging-induced volume loss predictions
  for (i in 1:L)
  {
     mudeltaV[i] = aG[logged[i]]/((theta-bG)) * ((exp(-bG*(tlog[i]))- exp(-bG*(t0))) - bG/theta*(exp(-theta*(tlog[i])) - exp(-theta*(t0))) ) - aM/((theta-bM)) * ((exp(-bM*(tlog[i]))- exp(-bM*(t0))) - bM/theta*(exp(-theta*(tlog[i])) - exp(-theta*(t0)) ))   ;

    // likelihood
    target += normal_lpdf(log(deltaV[i]+1)  | log(mudeltaV[i]+1), sigma_deltaV)*3;
  }
  
  // priors
  t0 ~ normal(200,50);
  sigma_deltaV ~ normal(0,1);
  sigma_V ~ normal(0,1); 
  sigma_G ~ normal(0,1);
  sigma_M ~ normal(0,1);
}

