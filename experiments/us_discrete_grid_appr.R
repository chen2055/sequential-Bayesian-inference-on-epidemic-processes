##this script runs a discrete grid approximation for the US ILI dataset

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

#Compute the log likelihood function over grid values of mu and sigma_sq
r0_sequence_coarse = seq(1.1,2.0,0.1)
rho_sequence_coarse = seq(0.1,1,0.1)

log_pos_values = matrix(NA, nrow=length(r0_sequence), ncol=length(rho_sequence))

for(i in 1:length(r0_sequence)){
  for(j in 1:length(rho_sequence)){
    temp_params = params0
    temp_params[c("r0","rho")] = c(r0_sequence[i], rho_sequence[j])
    result = my_pf(temp_params)
    log_pos_values[i,j] = result$pos
  }
}

#visualize log posterior function
contour(r0_sequence_coarse, rho_sequence_coarse[3:10], log_pos_values[,3:10], main=expression(paste("Joint log posterior of "*r0*" and "*rho*" under dnb")), xlab=expression(r0), ylab=expression(rho), col="blue")

#finer discrete grid so that we could sample from the approximated posterior
r0_sequence_fine = seq(1.25,1.5,0.01)
rho_sequence_fine = seq(0.5,0.7,0.01)

grid_log_posterior_values = matrix(NA, nrow=length(r0_sequence_fine), ncol=length(rho_sequence_fine))

for(i in 1:length(r0_sequence_fine)){
  for(j in 1:length(rho_sequence_fine)){
    temp_params = params0
    temp_params[c("r0","rho")] = c(r0_sequence_fine[i], rho_sequence_fine[j])
    result = my_pf(temp_params)
    grid_log_posterior_values[i,j] = result$pos
  }
}

par(mfrow=c(1, 2))

contour(r0_sequence_coarse, rho_sequence_coarse[3:10], log_pos_values[,3:10], xlab=expression(r0), ylab=expression(rho), col="blue")

contour(r0_sequence_fine, rho_sequence_fine, grid_log_posterior_values, xlab=expression(r0), ylab=expression(rho), col="brown")



#convert log posterior to posterior
grid_posterior_values = exp(grid_log_posterior_values - min(grid_log_posterior_values))
grid_posterior_values = grid_posterior_values/(sum(grid_posterior_values))

contour(r0_sequence_fine, rho_sequence_fine, grid_posterior_values, xlab=expression(r0), ylab=expression(rho), col="black")

c_grid_posterior_values = c(grid_posterior_values)

grid_appr_df = expand.grid(r0_sequence,rho_sequence)

# sample from the discretized grid approximated posterior
set.seed(84735)
post_sample_inds <- sample(1:length(c_grid_posterior_values), size = 10000, prob = c_grid_posterior_values, replace = TRUE)
draws_discrete = grid_appr_df[post_sample_inds,]

par(mfrow=c(1, 2))
hist(draws_discrete[,1], main=expression("Discrete for "*r0*""), probability=TRUE, col="gray", border="white",xlab = expression(r0))
d1 <- density(na.omit(draws_discrete[,1]))
lines(d1, col="blue")

hist(draws_discrete[,2], main=expression("Discrete for "*rho*""), probability=TRUE, col="gray", border="white",xlab = expression(rho))
d1 <- density(na.omit(draws_discrete[,2]))
lines(d1, col="blue")