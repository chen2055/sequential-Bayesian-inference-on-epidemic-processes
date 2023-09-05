#this script includes function related to the process model

# define a sampler for the prior on the initial state
simx0 <- function(n,t0 = 0,params = params0){
  #n: number of particles
  #t0: starting week
  #params: parameters
  with(
    as.list(params),
    {
      I0 = max(I0_ratio*N,0)
      S0 = max(N - I0,0)
      R0 = N - I0 - S0
      inci0 = 0
      week = t0
      SI_mat=cbind(rpois(n,S0),rpois(n,I0))
      SIR_mat = cbind(SI_mat,pmax(0,N - (SI_mat[,1]) - SI_mat[,2]))
      SIR_mat = cbind(SIR_mat,rep(inci0,nrow(SIR_mat)))
      SIR_mat = cbind(SIR_mat,rep(week,nrow(SIR_mat)))
      colnames(SIR_mat)=c("S","I","R","inci","week")
      rownames(SIR_mat) = 1:nrow(SIR_mat)
      SIR_mat
    }
  )
}

#step function
step_sir <- function(x,params = params0){
  #x:current state
  #params:parameters
  #seas_option: 0: no seasonality;
  #             1: only one frequency
  #             2: two frequencies
  #dom.freqs: two dominant frequencies of southeast Asia
  S <- x[1]
  I <- x[2]
  R <- x[3]
  inci <- x[4]
  week<-x[5]
  with( #use with as in deterministic model to simplify code
    as.list(params),
    {
      bS <- rbinom(n = 1, size = round(N), prob = 1-exp(-mu*delta.t)) #birth of Susceptible
      w_gamma = gamma * 7
      #default seasonality option: no seasonality
      beta = (w_gamma + mu)*r0
      if(seas_option ==1){
        cos1 = cos(2*pi*(week+40)/52.17)
        #beta = a + b*cos1
        beta = exp(c + a1*cos1)
      }
      
      if(seas_option ==2){
        cos1 = cos(2*pi*dom.freqs[1]*week/52.17)
        cos2 = cos(2*pi*dom.freqs[2]*week/52.17)
        sin1 = sin(2*pi*dom.freqs[1]*week/52.17)
        sin2 = sin(2*pi*dom.freqs[2]*week/52.17)
        beta = exp(c + a1*cos1 + a2*cos2+ b1*sin1 + b2*sin2)
      }
      
      if(S>0 && I>0){
        dN_SI <- rbinom(n=1,size=round(S),prob=1-exp(-beta*I/N*delta.t))
        dS <- rbinom(n = 1, size = round(S), prob = 1-exp(-mu*delta.t)) #death of Susceptible
      }else{
        dN_SI<-0
        dS<-0
      }
      if(I>0){
        dN_IR <- rbinom(n=1,size=round(I),prob=1-exp(-w_gamma*delta.t))
        dI <- rbinom(n = 1, size = round(I), prob = 1-exp(-mu*delta.t)) #death of infected
      }else{
        dN_IR<-0
        dI<-0
      }
      
      if(R>0){
        dR <- rbinom(n = 1, size = round(R), prob = 1-exp(-mu*delta.t)) #death of recovered
      }else{
        dR <- 0
      }
      
      S <- S - dN_SI + bS - dS
      I <- I + dN_SI - dN_IR - dI
      R <- R + dN_IR - dR
      inci <- max(dN_SI - bS + dS,0)
      week <- week + delta.t
      c(S,I,R,inci,week)
    }
  )
}

#weekly step function
SIR.weekly.model <- function (x, params= params0){ #function to simulate stochastic SIR
  dim_x = length(x)
  output <- array(dim=c(8,dim_x)) #set up array to store results
  colnames(output) <- c("S","I","R","inci","week") #name variables
  output[1,] <- x #first record of output is initial condition
  for (k in 1:7) { #iterate for nstep steps
    output[k+1,] <- x <- step_sir(x,params)
  }
  #update incidence
  output[8,4] = max(0,output[1,1] -  output[8,1])
  output[8,] #return output
}
