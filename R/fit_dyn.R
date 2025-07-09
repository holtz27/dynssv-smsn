source('https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R')
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/Topicos/source/criterion.R')

# compiling stan models
path='~/topicos/st/models/pcp.stan'
model_stan=rstan::stan_model(file=path)

# loading data
source('~/topicos/aplication/insample/bitcoin.R')
log.ret=log.ret-mean(log.ret)

warmup=5e1
iter=2e1

###################
## fitting model ##
###################
draws=rstan::sampling(model_stan, 
                      data=list(T=length(log.ret), 
                                y=as.numeric(log.ret),
                                lambda=-log(0.5)/0.5),
                      chains=1,
                      warmup=warmup,
                      iter=warmup+iter,
                      cores=1)
pars=c('mu','phi_h','s_h','a1','h','a','s_a', 'v','mu_t','sigma_t', 'U')
x=rstan::extract(draws, pars=pars)
theta=rbind(x$mu, x$phi_h, x$s_h, x$a1, x$s_a, x$v)
names=c('mu','phi','sh','a1','sa','v')

summary=num_analisys(draws=theta, names=names, digits=4, hdp=TRUE)
summary

pdf('dyn.pdf',width=20, height=10)
trace_plots(theta, names=names)
dev.off()

# Info Criterios
info.crit=c(dic=dic(data=log.ret, mu_t=x$mu_t, sigma_t=x$sigma_t),
            waic=waic(data=log.ret, mu_t=x$mu_t, sigma_t=x$sigma_t)$estimates['waic',1],
            loo=loo(data=log.ret, mu_t=x$mu_t, sigma_t=x$sigma_t)$estimates['looic',1])
info.crit

a=x$a
h=x$h
u=x$U
save(summary, info.crit, h, a, u, file='dyn_st.RData')