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
  real<lower=2> v;
}
transformed parameters{
  real<lower=-1,upper=1> phi_h;
  real<lower=0> s2_h;
  vector[T] h = h_std * s_h;
  phi_h = (2*phiT_h - 1);
  s2_h = pow(s_h, 2);
  
  h[1] /= sqrt(1 - phi_h * phi_h);
  h += mu;
  for (t in 2:T){
    h[t] += phi_h*(h[t-1]-mu);
  }

  vector[T] sigma_t = sqrt((v-2)/v)*exp(0.5*h)./sqrt(U);
}
model{
  // Prioris h
  mu ~ normal(0, sqrt(10));
  phiT_h ~ beta(20, 1.5);
  s2_h ~ inv_gamma(2.5, 0.025);
  
  //tails
  v ~ gamma(2.0, 0.1);
  
  // model
  U ~ gamma(0.5*v, 0.5*v);
  h_std ~ std_normal();
  y ~ normal(0, sigma_t);
}
