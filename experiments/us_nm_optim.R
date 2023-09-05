##this script runs a Nelder-Mead Optimization on US ILI data

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

params0 <- list(r0=1.34,gamma=0.15,mu = 0.0002,delta.t = 1/7, N = pop.size, rho = 0.64,noiseSD = 20, size = 20,I0_ratio = 5/pop.size,a = 1.1, b = 0.6, c = 0.1,a1 = 0.05,a2 = 0.05,b1 = 0.05,b2 = 0.05,sigma_b = 0.5,seas_option = 0,dom.freqs = c(0.8695, 1.9129),measure_option = "dpois")

#use Poisson Likelihood
my_pf=pfPos(1000,simx0,SIR.weekly.model,dpoisLik,y,all_times)

#use Negative Binomial Likelihood
#my_pf=pfPos(1000,simx0,SIR.weekly.model,dnbLik,y,all_times)
#params0$measure_option = "dnb"

neg.pos <- function (par, pest = c("c","a1"),qest = c("r0","gamma",'rho',"I0_ratio")) {
  allpars <- params0
  allpars[pest] <- par[pest]
  allpars[qest] <- sapply(par[qest],plogis)
  allpars$r0 = allpars$r0*7
  result = my_pf(allpars)
  -result$pos
}

#no seasonal parameters
#est = c("r0","gamma",'rho',"I0_ratio")
#pest = c()
#qest = c("r0","gamma",'rho',"I0_ratio")

#seasonal parameters
est = c("c","a1","r0","gamma",'rho','I0_ratio')
pest = c("c","a1")
qest = c("r0","gamma",'rho','I0_ratio')

temp_par = params0[est]
temp_par$r0 = temp_par$r0/7
temp_par[qest] = qlogis(as.numeric(temp_par[qest]))
#temp_pos = neg.pos(temp_par)

nm_fit <- optim(
  par=temp_par,
  pest = c("c","a1"),
  qest = c("r0","gamma",'rho','I0_ratio'),
  fn=neg.pos,
  method="Nelder-Mead",
  control=list(maxit=50,trace = 1)
)

#convert the parameters to between (0,1) after Nelder-Mead optimization
convert_nm_optim_par<-function(nm_fit,pest = c("c","a1"),qest = c("r0","gamma",'rho','I0_ratio')){
  nm_optim_par = nm_fit$par
  nm_optim_allpars = params0
  nm_optim_allpars[pest] <- nm_optim_par[pest]
  nm_optim_allpars[qest] <- sapply(nm_optim_par[qest],plogis)
  nm_optim_allpars$r0 = nm_optim_allpars$r0*7
  nm_optim_allpars
  }

nm_optim_allpars_seas = convert_nm_optim_par(nm_fit)
optim_result_seas = my_pf(nm_optim_allpars_seas)

par(mfrow=c(1, 1))
plot(y,type = "l",main = paste("dnb",",logPos=",round(optim_result_seas$pos),",r0=",round(nm_optim_allpars_seas$r0,2),",rho=",round(nm_optim_allpars_seas$rho,2)), axes=FALSE,ylim=c(0,max(max(y),max(optim_result_seas$xmean*nm_optim_allpars_seas$rho))),xlab = "week",cex.main = 0.8)
box()
axis(2, cex.axis=0.8, las=2)
lines(optim_result_seas$xmean[2:(length(all_times)+1)]*nm_optim_allpars_seas$rho,col = "red")
axis(side = 1,cex.axis=0.8)

#saveRDS(nm_optim_allpars_seas,"data/us_seas_nm_optim_allpars_dnb_logbeta.RDS")
