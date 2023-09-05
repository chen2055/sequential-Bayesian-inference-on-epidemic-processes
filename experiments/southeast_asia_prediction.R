#this script makes predictions based on posterior draws for the Southeast Asia dataset

##load the process model, measurement model, particle filter algorithm
setwd("/Users/yanchen/Documents/git/stat656_project")
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

#===== predictions
y = round(ili.region.b[(54+104):(54+140-1),"AH3"])
all_times = 1:length(y)
n_pred = 3
my_pf_pred=pfPos(1000,simx0,SIR.weekly.model,dpoisLik,y,all_times,h = n_pred)

pred_result = my_pf_pred(params0)

#par(mfrow=c(1, 1))
plot(y,type = "l",main = paste("dpois",",logPos=",round(pred_result$pos),",PPC"), axes=FALSE,ylim=c(0,max(max(y),max(pred_result$xmean*params0$rho))),xlab = "week",cex.main = 0.8)
box()
axis(2, cex.axis=0.8, las=2)
lines(pred_result$xmean[2:(length(all_times)+1-n_pred)]*params0$rho,col = "red")
lines((length(all_times)-n_pred):length(all_times),pred_result$xmean[(length(all_times)+1-n_pred):(length(all_times)+1)]*params0$rho,col = "blue",lwd = 2)
axis(side = 1,cex.axis=0.8)


#setwd("/Users/yanchen/Documents/git/stat656_project/data")
thmatp = readRDS("data/southeastasia_seas_pMH_thmatp.RDS")

ls_xmean = list()#store xmean
ls_pos = list()#store posterior
pp_pos = matrix(NA,nrow = round(nrow(thmatp)/20))#store total post predictive posterior
v_rho = matrix(NA,nrow = round(nrow(thmatp)/20))#store under-reporting ratio
est = c("c","a1",'a2',"b1",'b2','r0','gamma','rho','I0_ratio')#all parameters in PMH
pest = c("c","a1",'a2',"b1",'b2') #unconstrained parameters
qest = c("r0","gamma",'rho',"I0_ratio")#parameters between (0,1)
count = 0
for(i in 2:nrow(thmatp)){
  if(i%%20 == 2){
    count = count + 1
    v_rho[count] = thmatp[i,"rho"]
    th = thmatp[i,]
    ls_th = as.list(th)
    params0[est] <- ls_th[est]
    temp_result=my_pf_pred(params0)
    ls_xmean[[count]]=temp_result$xmean
    ls_pos[[count]]=temp_result$v_pos
    pp_pos[count] = temp_result$pos
  }
}


#par(mfrow=c(1, 1))
plot(y,type = "l",main = "Post Predictive Checks", axes=FALSE,xlab = "week",ylim = c(0,max(sapply(ls_xmean,max))),cex.main = 0.8, lwd = 3)
box()
axis(2, cex.axis=0.8, las=2)
axis(side = 1,cex.axis=0.8)
count = 0
for(i in 2:nrow(thmatp)){
  if(i%%20 == 2){
    count = count + 1
    lines(ls_xmean[[count]][2:(length(all_times)+1-n_pred)]*v_rho[count],col = "grey")
    lines((length(all_times)-n_pred):length(all_times),ls_xmean[[count]][(length(all_times)+1-n_pred):(length(all_times)+1)]*v_rho[count],col = count + 1,lwd = 2)
  }
}


count = 0
ls_pred = list()
pp_error = matrix(NA,nrow = round(nrow(thmatp)/20))#store the error between actual and post predictions
for(i in 2:nrow(thmatp)){
  if(i%%20 == 2){
    count = count + 1
    ls_pred[[count]] = ls_xmean[[count]][(length(all_times)+1-n_pred):(length(all_times)+1)]*v_rho[count]
    pp_error[count] = mean(abs(ls_pred[[count]]  - y[(length(y)-n_pred):length(y)]))
  }
}

#plot(pp_error)
