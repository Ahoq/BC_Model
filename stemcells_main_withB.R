#Main file for full non-dimensionalized system
#install.packages("deSolve")
library(deSolve)
source('~/Dropbox/Cancer Modeling/TGF-B R files/stemcells_parameters_withB.R')
source('~/Dropbox/Cancer Modeling/TGF-B R files/stemcells_ode_withB.R')


yini=c(S_c=10^6*1/79, E_a=10^6*20/79, Tumor=10^6, B = 10^4, R = 10^6*20/79)
allout=c()

dt = 1e-1

  Tvals=seq(timestart,timeend, by=dt)
  out=lsoda(y = yini, times = Tvals, func = ODE, parms = NULL, jactype="fullint")
  LL=length(out[,1])
  allout=rbind(allout, out[1:(LL-1),])
  #yini=out[length(out[,1]),2:3]


# Rescalings Putting Dimensions back in. 
times=allout[,1]#*tau
S_c = allout[,2]#*Delta
E_a = allout[,3]#*epsilon
Tumor = allout[,4]#*epsilon
B = allout[,5]#*epsilon
R = allout[,6]

all.plots <- par(mfrow=c(2, 3))
pl1 = plot(times,S_c,col='purple',type='l',lwd=2,ylab='Stem Cells',xlab='Time in Days')
pl2 = plot(times,E_a,col='blue',type='l',lwd=2,ylab='Activated T-Cells',xlab='Time in Days')
pl3 = plot(times,Tumor,col='green',type='l',lwd=2,ylab='Tumor Cells',xlab='Time in Days')
pl4 = plot(times,B,col='yellow',type='l',lwd=2,ylab='TGF-B Cells',xlab='Time in Days')
pl5 = plot(times,R,col='orange',type='l',lwd=2,ylab='Tregs',xlab='Time in Days')
par(all.plots)

