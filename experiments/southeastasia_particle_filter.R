#this script runs a particle filter for the Southeast Asia dataset

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
my_params = params0
result = my_pf(my_params)

#par(mfrow=c(1, 1))
plot(y,type = "l",main = paste("dpois",",logPos=",round(result$pos),",r0=",round(my_params$r0,2),",rho=",round(my_params$rho,2),",reset_week=",reset_n_week), axes=FALSE,ylim=c(0,max(max(y),max(result$xmean*my_params$rho))),xlab = "week",cex.main = 0.8)
box()
axis(2, cex.axis=0.8, las=2)
lines(result$xmean[2:(length(all_times)+1)]*my_params$rho,col = "red")
axis(side = 1,cex.axis=0.8)