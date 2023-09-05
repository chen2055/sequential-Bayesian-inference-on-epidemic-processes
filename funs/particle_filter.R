#function to implement particle filter
pfPos <- function (n, simx0, stepFun, dataLik,y,all_times=1:52,h = 1) 
  #n: number of particles
  #simx0: intial distribution of x
  #stepFun: simulate next x
  #dataLik: likelihood of y
  #y:data
  #all_times:time/step sequences
  #h:predictive steps
{
  return(function(...) {
    # if(measure_option == "dpois"){
    #   dataLik = dpoisLik
    # }else{
    #   dataLik = dnbLik
    # }
    nstep = length(all_times)
    xmean<-matrix(NA,nrow = nstep+1)
    xmat = simx0(n, t0, ...)
    xmean[1]= mean(xmat[,4])
    
    v_pos = rep(NA,nstep)
    prior_value = dprior(...)
    pos = prior_value
    for (i in 1:(nstep-h)) {
      #print(paste("step=======",i))
      if(i%%reset_n_week==1 && i >1){
        xmat = simx0(n, t0, ...)
      }else{
        xmat = t(apply(xmat, 1, stepFun, ...))
      }
      
      if(max(xmat[,4])!=0){
        #if max of sim x >0, use dpois
        w = apply(as.matrix(xmat[,4]), 1, dataLik, y = y[i], log = F,...)
      }else{
        #else, use dnorm
        w = apply(as.matrix(xmat[,4]), 1, dnormLik, y = y[i], log = F,...)
      }
      
      #print(paste("mean of sim:",mean(xmat[,4])))
      #print(paste("data",y[i]))
      # if(i==11){
      #   print(summary(w))
      #   hist(xmat[,4])
      #   test_xmat = xmat[,4]
      # }
      
      #if (max(w) < 1e-20) {
      if (sum(w) == 0) {
        #print("Particle filter failed")
        if(max(xmat[,4])!=0){
          #if max of sim x >0, use dpois
          w = apply(as.matrix(xmat[,4]), 1, dataLik, y = y[i], log = T,...) #
        }else{
          #else, use dnorm
          w = apply(as.matrix(xmat[,4]), 1, dnormLik, y = y[i], log = T,...)
        }
        v_pos[i] = max(w)
        stored_w = w
        max_ind = which.max(w)
      }else{
        v_pos[i] = log(mean(w))
        w = w/sum(w) #normalize weights
        rows = sample(1:n, n, replace = TRUE, prob = w)
        xmat = xmat[rows, ]
        stored_w = w[rows]
      }
      pos = pos + v_pos[i]
      xmean[i+1]= mean(xmat[,4])
    }
    ##predictions
    if(h>0){
      for (i in (nstep-h+1):nstep) {
        xmat = t(apply(xmat, 1, stepFun, ...))
        if(max(xmat[,4])!=0){
          #if max of sim x >0, use dpois
          w = apply(as.matrix(xmat[,4]), 1, dataLik, y = y[i], log = F,...)
        }else{
          #else, use dnorm
          w = apply(as.matrix(xmat[,4]), 1, dnormLik, y = y[i], log = F,...)
        }
        v_pos[i] = log(mean(w))
        xmean[i+1]= mean(xmat[,4])
        #use the last available weights to make predictions
        # w = w/sum(w) #normalize weights
        # if (sum(stored_w) != 0) {
        #   stored_w = stored_w/sum(stored_w) #normalize weights
        #   xmean[i+1]= sum(xmat[,4]*stored_w)
        # }else{
        #   xmean[i+1]= mean(xmat[,4])
        # }
      }
    }
    
    #list(pos = pos, v_pos=v_pos,xmean = xmean,test_xmat=test_xmat)
    list(pos = pos, v_pos=v_pos,xmean = xmean)
  })
}
