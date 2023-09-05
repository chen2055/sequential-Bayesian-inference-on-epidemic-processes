##this script runs a particle filter on US ILI data with or without seasonality

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

my_pf=pfPos(1000,simx0,SIR.weekly.model,dpoisLik,y,all_times)
#my_pf=pfPos(1000,simx0,SIR.weekly.model,dnbLik,y,all_times,h = 0)
my_params = params0
#if adding seasonality, uncomment the next line
#my_params[c("seas_option")] = 1
result = my_pf(my_params)

#par(mfrow=c(1, 1))
plot(y,type = "l",main = paste("dpois",",logPos=",round(result$pos),",r0=",round(my_params$r0,2),",rho=",round(my_params$rho,2)), axes=FALSE,ylim=c(0,max(max(y),max(result$xmean*my_params$rho))),xlab = "week",cex.main = 0.8,lwd = 2)
box()
axis(2, cex.axis=0.8, las=2)
lines(result$xmean[2:(length(all_times)+1)]*my_params$rho,col = "red",lwd = 2)
axis(side = 1,cex.axis=0.8)

result$v_pos[is.infinite(result$v_pos)]=-1e10
plot_pos = -log(abs(result$v_pos))
par(new = T)
plot(plot_pos,col="blue"
     ,type='l',lty = "dotted"
     ,axes = F,xlab = '',ylab = ''
     ,ylim = range(plot_pos,na.rm = T),lwd = 2)
axis(side = 4,cex.axis=0.8)
abline(v = 11,col = "blue",lwd = 2)
