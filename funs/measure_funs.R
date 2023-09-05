#this script includes funtions related to measurement model
library(truncnorm)
library(extraDistr)

#prior distribution
dprior <- function(params = params0){
  #seas_option: 0: no seasonality;
  #             1: only one frequency
  #             2: two frequencies
  #measure_option: 
  #             dpois: use dpois
  #             dnb: use dnbnoim
  with(
    as.list(params),
    {
      base_prior_value = log(dtruncnorm(r0, a=1, b=Inf, mean = 1.25, sd = 0.15))+dunif(rho,log = T)+dbeta(gamma,2,5,log = T)+log(dtruncnorm(I0_ratio, a=0.1/N, b=0.3, mean = 10/N, sd = 10/N))+log(dtruncnorm(mu, a=0, b=0.01, mean = 0.00035, sd = 0.00015))+log(dtruncnorm(noiseSD, a=1, b=Inf, mean = 20, sd = 10))
      
      measure_prior_value = 0
      if(measure_option=="dnb"){
        #if use negative binomial, there is an extra over-dispersion parameter
        measure_prior_value = log(dtruncnorm(size, a=1, b=Inf, mean = 20, sd = 10))
      }
      
      season_prior_value = 0
      if(seas_option==1){
        #season_prior_value = sum(dnorm(c(a,b),log = T))
        season_prior_value = dnorm(c,log = T)+sum(dnorm(c(a1),mean = 0,sd = sigma_b,log = T))+log(dtruncnorm(sigma_b, a=0, b=Inf, mean = 0.5, sd = 0.1))
      }
      if(seas_option == 2){
        season_prior_value = dnorm(c,log = T)+sum(dnorm(c(a1,a2,b1,b2),mean = 0,sd = sigma_b,log = T))+log(dtruncnorm(sigma_b, a=0, b=Inf, mean = 0.5, sd = 0.1))
      }
      prior_value = base_prior_value + measure_prior_value + season_prior_value
      prior_value
    }
  )
}

#dnorm likelihood function for y
dnormLik <- function(x,y,log=TRUE,params = params0){
  with(
    as.list(params),
    {
      ll=dnorm(y,x*rho,noiseSD,log=TRUE)
      if (log)
        ll
      else
        exp(ll)
    }
  )
}

#dpois likelihood function for y
dpoisLik <- function(x,y,log=TRUE,params = params0){
  with(
    as.list(params),
    {
      if(x!=0){
        pois_mu = x*rho
      }else{
        pois_mu = 0.1
      }
      ll=dpois(y,pois_mu,log=TRUE)
      if (log)
        ll
      else
        exp(ll)
    }
  )
}

#negative binomial likelihood function for y
dnbLik <- function(x,y,log=TRUE,params = params0){
  with(
    as.list(params),
    {
      if(x!=0){
        db_mu = x*rho
      }else{
        db_mu = 0.1
      }
      ll=dnbinom(y, size = size, mu = x*rho, log =T)
      if (log)
        ll
      else
        exp(ll)
    }
  )
}
