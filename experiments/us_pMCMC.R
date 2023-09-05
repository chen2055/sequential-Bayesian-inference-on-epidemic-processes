##this script runs a particle Metroplis Hassting on US ILI data

##load the process model, measurement model, particle filter algorithm
#setwd("/Users/yanchen/Documents/git/stat656_project")
source("funs/process_funs.R")
source("funs/measure_funs.R")
source("funs/particle_filter.R")

####read in US ili data
us.ili.full=read.csv("data/us_ili_full.csv",stringsAsFactors = F)
#remove 2000-2001, 1999-2000, 2019-
us.ili.select = us.ili.full[-c(158:209,628:679,1149:nrow(us.ili.full)),c("Year","Week","AH3")]
rownames(us.ili.select) = 1:nrow(us.ili.select)
#y = us.ili.select[301:508,"AH3"]
y = us.ili.select[301:352,"AH3"]

all_times = 1:length(y)
reset_n_week = 52 #the number of weeks to restart the initiation of a stochastic process

pop.size <- 20000 #total population size
t0 = 0

params0 <- list(r0=1.5,gamma=0.15,mu = 0.0002,delta.t = 1/7, N = pop.size, rho = 0.8,noiseSD = 20, size = 20,I0_ratio = 5/pop.size,a = 1.1, b = 0.6, c = 0.35,a1 = 0.15,a2 = 0.05,b1 = 0.05,b2 = 0.05,sigma_b = 0.5,seas_option = 1,dom.freqs = c(0.8695, 1.9129),measure_option = "dpois")

#use Poisson Likelihood
my_pf=pfPos(1000,simx0,SIR.weekly.model,dpoisLik,y,all_times)

#use Negative Binomial Likelihood
#my_pf=pfPos(1000,simx0,SIR.weekly.model,dnbLik,y,all_times)
#params0$measure_option = "dnb"

est = c("r0","gamma",'rho',"I0_ratio")
pest = c()
qest = c("r0","gamma",'rho',"I0_ratio")

mc_par = params0[est]
mc_par$r0 = mc_par$r0/7
mc_par[qest] = qlogis(as.numeric(mc_par[qest]))

iters=500
p=length(est)

thmatq=matrix(0,nrow=iters+1,ncol=p)
thmatp=matrix(0,nrow=iters+1,ncol=p)
colnames(thmatq)=est
colnames(thmatp)=est
th = as.numeric(mc_par)
thmatq[1,]=th
mh_allpars = params0
mh_allpars[est] <- sapply(th,plogis)
mh_allpars$r0 = mh_allpars$r0*7
init_result=my_pf(mh_allpars)
pos = init_result$pos
#cov_appr = diag(c(0.01,0.01,0.04,1))
cov_appr = diag(c(0.01,0.015,0.04,0.5))
#cov_mh = cov_appr* (2.4^2 / p)
cov_mh = cov_appr
library(mvtnorm)
# Main pMCMC loop
for (i in 1:iters) {
  message(paste(i,""),appendLF=FALSE)
  thprop= rmvnorm(n=1, mean=thmatq[i,], sigma=cov_mh, method="chol")
  mh_allpars[est] <- sapply(thprop,plogis)
  mh_allpars$r0 = mh_allpars$r0*7
  temp_result=my_pf(mh_allpars)
  pos_prop = temp_result$pos
  print(pos_prop)
  if (log(runif(1)) < pos_prop - pos) {
    print("accept proposal")
    print(paste("pos improved",round(pos_prop - pos,2)))
    th=thprop
    pos=pos_prop
  }
  thmatq[i+1,]=th
}
message("Done!")
# Compute and plot some basic summaries
library(smfsb)
pdf("figures/thmat_summary_us_seasonal.pdf")
burnin = 100
thmatp= plogis(thmatq[(burnin+1):(iters + 1),])
thmatp[,"r0"] = thmatp[,"r0"]*7
mcmcSummary(thmatp)
dev.off()
#saveRDS(thmatp,"us_seas_pMH_thmatp.RDS")