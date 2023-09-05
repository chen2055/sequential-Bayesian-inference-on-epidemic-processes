#library(ggplot2)
library(quantmod) #function findPeaks

ili.sea.agg = read.csv("data/sea_ili_full.csv",stringsAsFactors = F)
ili.region.b = ili.sea.agg[366:nrow(ili.sea.agg),]
rownames(ili.region.b)=1:nrow(ili.region.b)

df.to.fit = ili.region.b[,c("Year","Week","AH3")]
nsamp = nrow(df.to.fit)
start.row = 1
start.year0 = 2004
v.start.row = c()
v.reach.row = c()
col_name = "AH3"
df.to.fit$rescale = df.to.fit[,col_name]
week.unit = 52
#scale the magnitude of cases within each year
for(t in 1:ceiling(nsamp/week.unit)){
  reach.row = min(start.row + week.unit, nsamp)
  v.start.row = c(v.start.row,start.row)
  v.reach.row = c(v.reach.row,reach.row)
  h3n2.select = df.to.fit[start.row:reach.row,col_name]
  h3n2.rescale = scale(h3n2.select,center = min(h3n2.select,na.rm = T),scale =max(h3n2.select,na.rm = T) - min(h3n2.select,na.rm = T) )
  df.to.fit[start.row:reach.row,"rescale"]=h3n2.rescale
  if(reach.row < nsamp){start.row = reach.row + 1}
}
df.to.fit =df.to.fit[1:nsamp,]
df.to.fit$rescale[is.na(df.to.fit$rescale)]=0
x = df.to.fit$rescale
len_x = length(x)

#====
## convoluted daniel kernel with window of 3
k3vd <- kernel("daniell", c(3,3))
xk3vd <- c(x[1:6],kernapply(x, k3vd),x[(len_x-5):len_x])

#====
## spetral analysis
x.spec <- spectrum(x,log="no",spans=c(3, 3),plot=FALSE)
del<-1/52.17 # sampling interval
spx <- x.spec$freq/del
spy <- 2*x.spec$spec
#plot(spy~spx,xlab="frequency",ylab="spectral density",type="l",lwd = 2)

#====
## After getting the spetral density, the next task is to choose dominant frequencies, the best number of dominant frequencies may be more than 1
dom.freqs = spx[intersect(findPeaks(spy,0.05)-1,which(spy>1.2))] #frequencies that have density greater than 1.2 and are local peaks