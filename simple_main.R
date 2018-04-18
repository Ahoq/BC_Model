 #Main file for full non-dimensionalized system
library(deSolve)
source('~/Dropbox/REUF_CancerVaccines/RCode/new_fullParam.R')
source('~/Dropbox/REUF_CancerVaccines/RCode/simple_Param.R')
#source('~/Dropbox/REUF_CancerVaccines/RCode/vblood.R')
#source('~/Dropbox/REUF_CancerVaccines/RCode/vtumor.R')
source('~/Dropbox/REUF_CancerVaccines/RCode/simple_fullDE.R')


yini=c(D_blood=0, Ea_blood=0, Em_blood=0, D_spleen=0, W=0, P=0, Ea_spleen=0, Em_spleen=0, Ea_tumor=0, Tumor=0, D_tumor=0)
yini=yini*c(1/Delta,1/epsilon,1/epsilon,1/Delta,1/Delta,1/epsilon,1/epsilon,1/epsilon,1/epsilon,1/sigma,1/Delta)
allout=c()

#vaccine is currently defined in new_fullParam

dt = 1e-7
for (i in 1:(length(vaccine[1,])-1)){
Tvals=seq(vaccine[1,i],vaccine[1,i+1], by=dt)
vb=vaccine[2,i]
out=lsoda(y = yini, times = Tvals, func = simple_fullDE, parms = vb, jactype="fullint")
LL=length(out[,1])
allout=rbind(allout, out[1:(LL-1),])
yini=out[length(out[,1]),2:12]
}

# Rescalings
times=allout[,1]*tau
D_blood = allout[,2]*Delta
Ea_blood = allout[,3]*epsilon
Em_blood = allout[,4]*epsilon
D_spleen = allout[,5]*Delta
W = allout[,6]*Delta
P = allout[,7]*epsilon
Ea_spleen = allout[,8]*epsilon
Em_spleen = allout[,9]*epsilon
Ea_tumor = allout[,10]*epsilon
Tumor = allout[,11]*sigma
D_tumor = allout[,12]*Delta


plot(times,D_blood,col='purple',type='l',lwd=2,ylab='D_blood',xlab='Time in Days')
plot(times,Ea_blood,col='blue',type='l',lwd=2,ylab='Ea_blood',xlab='Time in Days')
plot(times,Em_blood,col='green',type='l',lwd=2,ylab='Em_blood',xlab='Time in Days')
plot(times,D_spleen,col='purple',type='l',lwd=2,ylab='DC_spleen',xlab='Time in Days')
plot(times,W,col='red',type='l',lwd=2,ylab='W',xlab='Time in Days')
plot(times,P,col='orange',type='l',lwd=2,ylab='P',xlab='Time in Days')
plot(times,Ea_spleen,col='blue',type='l',lwd=2,ylab='Ea_spleen',xlab='Time in Days')
plot(times,Em_spleen,col='green',type='l',lwd=2,ylab='Em_spleen',xlab='Time in Days')
#plot(times,Ea_tumor,col='blue',type='l',lwd=2,ylab='Ea_tumor',xlab='Time in Days')
#plot(times,Tumor,col='black',type='l',lwd=2,ylab='Tumor',xlab='Time in Days')
#plot(times,D_tumor,col='purple',type='l',lwd=2,ylab='D_tumor',xlab='Time in Days')

#plot(out[,1],out[,3],col='red',type='l',lwd=2,
#    xlab='Time',ylab='W', main='Reproduction State')
#plot(out[,1],out[,3], col='blue', type='l', lwd=2)
