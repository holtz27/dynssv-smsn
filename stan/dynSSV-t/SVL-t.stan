data{
  int<lower=0> T;
  vector[T] y;
}
parameters{
  real mu;
  real<lower=0, upper=1> phiT_h;
  real<lower=0> s_h;
  vector[T] h_std;
  vector<lower=0>[T] U;
  real<lower=0, upper=1> rhoT;
  real<lower=2> v;
}
transformed parameters{
  real<lower=-1,upper=1> phi_h= (2*phiT_h - 1);
  real<lower=1> k2 = v/(v-2);
  real<lower=0> s2_h= pow(s_h, 2);
  real<lower=-1, upper=1> rho=(2*rhoT-1);
  vector[T] h=h_std*s_h*sqrt(1-rho*rho);

  h[1] /= sqrt(1-phi_h*phi_h);
  h+=mu;
  for(t in 2:T){
    h[t]+=phi_h*(h[t-1]-mu)+s_h*rho*exp(-0.5*h[t-1])*sqrt(k2*U[t-1])*y[t-1];
  }
  vector<lower=0>[T] sigma_t=exp(0.5*h)./sqrt(k2*U);
}
model{
  // Prioris h
  mu ~ normal(0, sqrt(10));
  phiT_h ~ beta(20, 1.5);
  s2_h ~ inv_gamma(2.5, 0.025);
  v ~ gamma(2.0, 0.1);
  //rhoT ~ beta(3.0,6.0);
  rho ~ uniform(-1.0, 1.0);
  
  // model
  U ~ gamma(0.5*v, 0.5*v);
  h_std ~ std_normal();
  y ~ normal(0, sigma_t);
}
