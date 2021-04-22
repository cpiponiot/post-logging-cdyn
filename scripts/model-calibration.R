#### Fit carbon model to Paracou data ####

library(data.table)
library(ggplot2)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

sapply(list.files("R", full.names = TRUE), source)

load("data/cdata.rda")

## cumulative biomass growth and loss
setorder(cdata, plot, year)
cdataPre = cdata[year<1990]
cdata = cdata[year>=1990]
cdata[, `:=`(cgain = cumsum(gain*dT), 
             closs = cumsum(loss*dT)), 
      .(plot)]

data_list <- list(N = nrow(cdata), 
                  P = length(unique(cdata$plot)), 
                  L = length(unique(cdata$plot[cdata$treat!=0])),
                  t = cdata$year - 1987, 
                  np = cdata$plot, 
                  logged = unique(cdata$plot[cdata$treat!=0]), 
                  cGrowth = cdata$cgain, 
                  cMort = cdata$closs, 
                  Vol = cdata$agb, 
                  deltaV = cdataPre[treat>0, .(agb[year==1986]-min(agb)), .(plot)]$V1, 
                  maxV = tapply(cdata$agb, cdata$plot, max),
                  lowBound = rep(0, 9), 
                  upBound = c(50, 10, 0.05, 0.1, 0.05, 300, 300, 3*max(cdata$agb), 500), 
                  supSigma = c(5,5,5,10,40))

if (!file.exists("stan/stan_cDyn.rda")) {
  stan_cDyn <- stan("stan/vdyn_model.stan", data = data_list, iter = 500, chains = 3)
  save(stan_cDyn, file = "stan/stan_cDyn.rda")
  traceplot(stan_cDyn, pars = c("lp__","theta","aM","bG","bM","t0"))
  traceplot(stan_cDyn, pars = c("mu_vmax"))
  traceplot(stan_cDyn, pars = c("vmax"))
} else load("stan/stan_cDyn.rda")

if (!file.exists("stan/stan_cDyn_1vmax.rda")) {
  stan_cDyn_1vmax <- stan("stan/vdyn_model_1vmax.stan", data = data_list, 
                          iter = 2000, chains = 3)
  save(stan_cDyn_1vmax, file = "stan/stan_cDyn_1vmax.rda")
  traceplot(stan_cDyn_1vmax, pars = c("lp__","theta","aM","bG","bM","t0"))
  traceplot(stan_cDyn_1vmax, pars = c("vmax"))
} else load("stan/stan_cDyn_1vmax.rda")

### max likelihood predictions
dfpred = data.table(expand.grid(mat = 1:200, plot = 1:12))
pars = rstan::extract(stan_cDyn_1vmax)
maxL = which.max(pars$lp__)
t0mL = pars$t0[maxL]
pars = data.table(plot = 1:12, 
                  aG = pars$aG[maxL,], 
                  aM = pars$aM[maxL], 
                  theta = pars$theta[maxL], 
                  bG = pars$bG[maxL], 
                  bM = pars$bM[maxL], 
                  t1 = pars$t1[maxL,], 
                  t0 = pars$t0[maxL])
dfpred = merge(dfpred, pars, by = "plot")

dfpred[, agb := aG / theta * (1 - (theta * exp(-bG * (mat)) - bG *
                                     exp(-theta * (mat))) / (theta - bG)) - aM / theta * (1 - (theta *
                                                                                                 exp(-bM * (mat)) - bM * exp(-theta * (mat))) / (theta - bM))]
dfpred[, awp := aG*(1-exp(-bG*(mat))) - theta*agb]
dfpred[, awm := aM*(1-exp(-bM*(mat)))]

dfpred = merge(dfpred, unique(cdata[, c("plot", "treat")]))

cdata = merge(cdata, pars[, c("plot", "t1")])
cdata$mat <- cdata$t1 + (cdata$year - 1987)

ggplot(dfpred, aes(x = mat, y = agb, color = as.factor(treat))) + 
  geom_vline(xintercept = t0mL, lty = 2)+ 
  geom_line() + 
  geom_point(data = cdata) 
ggplot(dfpred, aes(x = mat, y = awp, color = as.factor(treat))) + 
  geom_vline(xintercept = t0mL, lty = 2)+ 
  geom_line() + 
  geom_point(data = cdata) 
ggplot(dfpred, aes(x = mat, y = awm, color = as.factor(treat))) + 
  geom_vline(xintercept = t0mL, lty = 2)+ 
  geom_line() + 
  geom_point(data = cdata)  

