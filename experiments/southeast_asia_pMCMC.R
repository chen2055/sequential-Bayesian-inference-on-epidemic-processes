#this script runs a particle Metropolis Hasting on the Southeast Asia dataset

##load the process model, measurement model, particle filter algorithm
#setwd("/Users/yanchen/Documents/git/stat656_project")
source("funs/process_funs.R")
source("funs/measure_funs.R")
source("funs/particle_filter.R")

#obtain dom.freqs
source("funs/spec.R")
#dom.freqs = c(0.8695, 1.9129)

#read in southeast asia data
ili.sea.agg = read.csv("data/sea_ili_full.csv",stringsAsFactors = F)
ili.region.b = ili.sea.agg[366:nrow(ili.sea.agg),]
rownames(ili.region.b)=1:nrow(ili.region.b)
## y is the measurements we are going to use and is usually in increments of 52 weeks
#y = round(ili.region.b[54:(54+52*4-1),"AH3"])
y = round(ili.region.b[(54+104):(54+52*4-1),"AH3"])
weeks = 1:52 #specifies that the granularity of our model is in week
all_times = 1:length(y) #the index for the time horizon we are studying
#reset_n_week = 26
reset_n_week = 20 #the number of weeks to restart the initiation of a stochastic process


t0 = 0
pop.size = 2000
params0 <- list(r0=1.4,gamma=0.15,mu = 0.0002,delta.t = 1/7, N = pop.size, rho = 0.6,noiseSD = 20, size = 20,I0_ratio = 5/pop.size,a = 1.1, b = 0.6, c = 0.3,a1 = 0.05,a2 = 0.05,b1 = 0.25,b2 = -0.2,sigma_b = 0.5,seas_option = 2,dom.freqs = c(0.8695, 1.9129),measure_option = "dpois")


#use the Poisson measurement model
my_pf=pfPos(1000,simx0,SIR.weekly.model,dpoisLik,y,all_times,h = 0)

#==== pMCMC
est = c("c","a1",'a2',"b1",'b2','r0','gamma','rho','I0_ratio')
pest = c("c","a1",'a2',"b1",'b2')
qest = c("r0","gamma",'rho',"I0_ratio")

#qlogis conversion
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
ls_th = as.list(th)
names(ls_th) = est
thmatq[1,]=th
mh_allpars = params0
#plogis conversion
mh_allpars[pest] <- ls_th[pest]
mh_allpars[qest] <- sapply(ls_th[qest],plogis)
mh_allpars$r0 = mh_allpars$r0*7
init_result=my_pf(mh_allpars)
pos = init_result$pos#get initial posterior
#cov_appr = diag(c(0.01,0.01,0.04,1))
#cov_appr = diag(c(rep(0.05,5),c(0.01,0.015,0.04,0.5)))
cov_appr = diag(c(rep(0.01,5),c(0.01,0.015,0.04,0.1)))
#cov_mh = cov_appr* (2.4^2 / p)
cov_mh = cov_appr
library(mvtnorm)
# Main pMCMC loop
for (i in 1:iters) {
  message(paste(i,""),appendLF=FALSE)
  thprop= rmvnorm(n=1, mean=thmatq[i,], sigma=cov_mh, method="chol")
  ls_thprop = as.list(thprop)
  names(ls_thprop) = est
  mh_allpars[pest] <- ls_thprop[pest]
  mh_allpars[qest] <- sapply(ls_thprop[qest],plogis)
  mh_allpars$r0 = mh_allpars$r0*7
  #print(as.numeric(mh_allpars[est]))
  temp_result=my_pf(mh_allpars)
  pos_prop = temp_result$pos
  print(pos_prop)
  if (log(runif(1)) < pos_prop - pos) {
    print("accept proposal")
    print(paste("pos improved",round(pos_prop - pos,2)))
    th=thprop
    pos=pos_prop
  }
  #}
  thmatq[i+1,]=th
  #thmatp[i+1,]=plogis(th)
  #thmatp[i+1,"r0"]=thmatp[i+1,"r0"]*7
}
message("Done!")
library(smfsb)# plot some basic summaries
#setwd("/Users/yanchen/Documents/git/stat656_project/figures")
pdf("figures/thmat_summary_southeastasia_smallp.pdf")
burnin = 100
thmatp= thmatq[(burnin+1):(iters + 1),]
thmatp[,qest]= plogis(thmatq[(burnin+1):(iters + 1),qest])
thmatp[,"r0"] = thmatp[,"r0"]*7
#function from smfsb to draw the trace, ACF and histogram plots for draws from MCMC
mcmcSummary(thmatq[(burnin+1):(iters + 1),])
dev.off()
#View(mcmcSummary)

library(coda)
cat("The effective sample sizes of the draws for the pMH draws are:")
round(effectiveSize(thmatp))

#setwd("/Users/yanchen/Documents/git/stat656_project/data")
#saveRDS(thmatp,"southeastasia_seas_pMH_thmatp.RDS")


