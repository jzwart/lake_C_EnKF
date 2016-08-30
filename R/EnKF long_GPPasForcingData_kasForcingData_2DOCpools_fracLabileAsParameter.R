# setting up Long Data for EnKF; EnKF runs at end of code
# JAZ; 2016-03-28

# rm(list=ls())
# # setwd('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/MCMC/Chris_model/')
# source('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/R Code/loadLongData.R')
# library(MASS)
# library(data.table)
# library(R.utils)
# sourceDirectory('/Users/Jake/Desktop/R functions/',modifiedOnly = F)
# 
# data$Iota<-data$Iota/32/1000*1000/(60*60*24) # converting iota from mg O2 L-1 day-1 (mmol photons m-2 sec-1)-1 to mol C m-3 day-1 (mmol photons m-2 day-1)-1 
# data$sdIota<-data$sdIota/32/1000*1000/(60*60*24)
# data$GPP<-data$GPP/32/1000*1000*data$epiVol # converting GPP from mg O2 L-1 day-1 to mol C day-1 in epi 
# data$sdGPP<-data$sdGPP/32/1000*1000*data$epiVol # converting sd GPP from mg O2 L-1 day-1 to mol C m-3 day-1 
# data$R<-data$R/32/1000*1000*data$epiVol # converting R from mg O2 L-1 day-1 to mol C day-1 in epi 
# data$sdR<-data$sdR/32/1000*1000*data$epiVol # converting R from mg O2 L-1 day-1 to mol C day-1 in epi 
# data$doc<-data$doc/12/1000*1000*data$epiVol # converting from mg L-1 to mol C in epi
# data$dic<-data$dic*data$epiVol #converting from mol C m-3 to mol C in epi
# 
# data<-data[min(which(!is.na(data$docIn))):max(which(!is.na(data$docIn))),]
# data<-data[min(which(!is.na(data$GPP))):max(which(!is.na(data$GPP))),] # earliest to latest GPP observations 
# data<-data[!is.na(data$GPP),] # get rid of days without GPP obs 
# # quick fix for missing outlet discharge data 
# data$QoutInt<-data$Qout
# data$QoutInt[1]<-data$Qout[min(which(!is.na(data$Qout)))]
# data$QoutInt[length(data$datetime)]<-data$Qout[max(which(!is.na(data$Qout)))] 
# data$timeStep<-seq(1:length(data$datetime))
# data$QoutInt<-approx(data$timeStep,data$QoutInt,data$timeStep)$y
# data<-data[!is.na(data$docIn),]
# 
# # QA/QC dic data 
# data$dic[data$dic<0]<-NA # changin negative DIC to NA 
# data$dic[data$dic>mean(data$dic,na.rm=T)+3*sd(data$dic,na.rm=T)]<-NA # if dic is 3 standard deviations beyond mean, turn to NA 
# 
# ##### In lake processes
# photoOx=0.0042   # photooxidation rate constant, [mol c m-2 day-1]; Graneli et al. 1996, L&O, 41: 698-706
# exude=0.03	# fraction of GPP released as refractory DOC, [day-1]; Biddanda & Benner 1997 via Hanson et al. 2004
# exudeLabile=0.07 # fraction of GPP released as labile DOC, [day-1]; Hanson etal 2004 
# exudeTotal=exude+exudeLabile
# leafLeach=0.12/14	# fraction of leaf load released as DOC, [day-1]; France et al. 1997; saw 6-18% loss over two weeks  do this on an annual leaf load * leach rate OR make this seasonal; this also has to be input to sediment/hypo OM
# floc=0.005		# fraction of DOC that floculates, [day-1]; von Wachenfeldt & Tranvik 2008 via Jones et al. 2012
# GPPrespired=0.85	# fraction of GPP that is respired quickly, [day-1]; Quay et al. 1996, Cole et al. 2002 via Hanson et al. 2004; Hanson uses .8; cole et al. ranges from .85-.90
# DOCrespired=0.005 # fraction of DOC that is respired - background R, [day-1]; Houser et al. 2001 via Hanson et al. 2004
# 
# # temperature dependence for r; r units = day-1 
# rFunc<-function(r20,temp){
#   r<-r20*1.047^(temp-20) # scaled R based on Holtgrieve et al. 2010
#   return(r)
# }
# # rFunc<-function(r20,temp){
# #   r<-r20*1.09^(temp-20) # scaled R based on Holtgrieve et al. 2010
# #   return(r)
# # }
# 
# # wind dependent k; k units = m day-1 ; returns kCO2
# kFunc<-function(c2,wnd,temp){ # c2 is the intercept, c1 the slope of the Cole Caraco 1998 wind-based k600 model 
#   k600=c2+(0.215*(wnd^1.7)) #units in cm h-1
#   k600=k600*24/100 # units in m day-1 
#   kCO2=k600.2.kGAS.base(k600,temp,gas='CO2') # k.CO2 in m day-1
#   return(kCO2)
# }

###################
#Run Kalman filter
rm(list=ls())
load('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/R Data/EnKF_LongData.RData')

# spin up, just repeating the first X days at the begining of the timeseries; autocorrelation effects?? 
spinUpLength<-0
spinUp<-data[1:spinUpLength,]
data2<-rbind(spinUp,data)

fracLabile0<-0.08 # fraction of initial DOC pool that is labile
fracLabile<-0.1 # fraction loaded t-DOC that is labile; estimate from Berggren et al. 2010 Eco Letts & ISME

# need some initial conditions for observations (draw from normal distribution?)
data2$dic[1]<-data2$dic[min(which(!is.na(data2$dic)))]
data2$doc[1]<-data2$doc[min(which(!is.na(data2$doc)))]

## *********************** EnKF ***************************## Gao et al. 2011 is a useful reference 
nEn<-100 # number of ensembles 
nStep<-length(data2$datetime)

# draws from priors to create ensemble 
parGuess <- c(0.0015,0.08,0.2) #r20; units: day-1; fraction labile of loaded DOC 
min<-c(0.0001,0.06,0.15)
max<-c(0.003,0.11,0.25)

rPDF<-abs(rnorm(n=nEn,mean = parGuess[1],sd = (max[1]-min[1])/5)) # max-min / 5 is rule of thumb sd; forcing positive for negative draws 
hist(rPDF)
rPDF_fast<-abs(rnorm(n=nEn,mean = parGuess[2],sd = (max[2]-min[2])/5))
hist(rPDF_fast)
fracPDF<-abs(rnorm(n=nEn,mean=parGuess[3],sd=(max[3]-min[3])/5))
hist(fracPDF)

# setting up initial parameter values for all ensemble members 
rVec<-matrix(rPDF) # each row is an ensemble member 
rVec_fast<-matrix(rPDF_fast)
fracVec<-matrix(fracPDF)

# initial B transition matrix for each ensemble  
B<-array(NA,dim=c(4,4,nStep,nEn)) # array of transition matrix B[a,b,c,d]; 
                                      # where [a,b] is transition matrix DIC & DOCr & DOCl, c=timeStep, and d=ensemble member 
# initial parameters of B at timestep 1 
for(i in 1:nEn){
  B[,,1,i]<-matrix(c(1-data2$QoutInt[1]/data2$epiVol[1]-data2$kCO2[1]/data2$thermo.depth[1],
                rFunc(rVec[i],data2$wtr[1]),rFunc(rVec_fast[i],data2$wtr[1]),0,
                0,1-rFunc(rVec[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1],0,0,
                0,0,1-rFunc(rVec_fast[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1],0,
                0,1-rFunc(rVec[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1],1-rFunc(rVec_fast[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1],0),
                nrow=4,byrow=T)
}

# setting up state values for all ensemble members (all the same initial pool size )
# updated on 2016-07-16 to draw from a normal distribution using first observation and SD
dicVec<-data2$dic
docVec<-data2$doc
if(is.na(dicVec[1])){
  dicVec[1]<-dicVec[min(which(!is.na(dicVec)))]
  docVec[1]<-docVec[min(which(!is.na(docVec)))]
}

# estimate of observation error covariance matrix
# error in DOC / DIC pool esitmate comes from concentration estimates and volume of epilimnion
# DOC / DIC concentration error from replication day in WL on 2014-07-30. CV for DOC is 0.1325858, DIC is 0.1367684
# Epi: area sd = 1000m2; depth sd = 0.5 m;
# iota sd from metab estimates 
docSD<-(data2$doc/data2$epiVol)*0.1325858 # DOC concentration sd in mol C
docSD<-docSD # DOC concentration sd in mol C; modifying for sensativity analysis
dicSD<-(data2$dic/data2$epiVol)*0.1367684 # DIC concentration sd in mol C 
dicSD<-dicSD*4
areaSD<-0 # constant SD in m 
depthSD<-0 # constant SD in m 

H<-array(0,dim=c(2,2,nStep))
# propogation of error for multiplication 
docPoolSD<-data2$doc*sqrt((docSD/(data2$doc/data2$epiVol))^2+(areaSD/data2$A0)^2+(depthSD/data2$thermo.depth)^2)
dicPoolSD<-data2$dic*sqrt((dicSD/(data2$dic/data2$epiVol))^2+(areaSD/data2$A0)^2+(depthSD/data2$thermo.depth)^2)
H[1,1,]<-dicPoolSD^2 #variance of DIC
H[2,2,]<-docPoolSD^2 #variance of DOC 
H[1,1,]<-ifelse(is.na(H[1,1,]),mean(H[1,1,],na.rm=T),H[1,1,]) # taking care of na's in DIC; setting to mean of variance if NA 
H[2,2,]<-ifelse(is.na(H[2,2,]),mean(H[2,2,],na.rm=T),H[2,2,]) # taking care of na's in DOC; setting to mean of variance if NA 

y=array(rbind(dicVec,docVec),dim=c(2,1,nStep))
y=array(rep(y,nEn),dim=c(2,1,nStep,nEn)) # array of observations y[a,b,c,d]; where a=dic/doc, b=column, c=timeStep, and d=ensemble member 

X<-array(NA,dim =c(4,1,nStep,nEn)) # model estimate matrices X[a,b,c,d]; where a=dic/doc_r/doc_l, b=column, c=timestep, and d=ensemble member 

#initializing estimate of state, X, with first observation 
# drawing from normal distribution for inital pool sizes 
X[1,1,1,]<-rnorm(n=nEn,y[1,1,1,],sd=dicPoolSD[1])
X[2,1,1,]<-rnorm(n=nEn,y[2,1,1,]*(1-fracLabile0),sd=docPoolSD[1]*(1-fracLabile0)) # labile pool is 90% of initial DOC pool
X[3,1,1,]<-rnorm(n=nEn,y[2,1,1,]*fracLabile0,sd=docPoolSD[1]*fracLabile0) # recalcitrant pool is 10% of initial DOC pool
X[4,1,1,]<-X[2,1,1,]+X[3,1,1,] # recalcitrant pool + labile pool 

# operator matrix saying 1 if there is observation data available, 0 otherwise 
h<-array(0,dim=c(2,7,nStep))
for(i in 1:nStep){
  h[1,3,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic 
  h[2,7,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) #doc total (we only have data on total DOC pool)
}

P <- array(0,dim=c(3,3,nStep))

#Define matrix C, parameters of covariates [2x6]
C<-array(NA,dim=c(4,6,nStep,nEn)) # array of transition matrix C[a,b,c,d]; 
#intializing first time step 
for(i in 1:nEn){
  C[,,1,i]<-matrix(c(data2$kCO2[1],1,-data2$GPP[1],0,0,0,0,0,0,data2$GPP[1],0,(1-fracVec[i]),
                     0,0,0,0,data2$GPP[1],fracVec[i],
                     0,0,0,data2$GPP[1],data2$GPP[1],1),nrow=4,byrow=T)
}

ut<-array(NA,dim=c(6,1,nStep,nEn)) # array of transition matrix C[a,b,c,d]; 
#intializing first time step 
for(i in 1:nEn){
  ut[,,1,i]<-matrix(c(data2$DICeq[1]*data2$epiVol[1]/data2$thermo.depth[1],
                      data2$dicIn[1],(1-GPPrespired),
                      exude,exudeLabile,data2$docIn[1]),ncol=1)
}

pars<-array(rep(NA,nEn),dim=c(3,1,nStep,nEn)) # parameters: r20
pars[1,1,1,]<-rVec
pars[2,1,1,]<-rVec_fast
pars[3,1,1,]<-fracVec

# set up a list for all matrices 
z=list(B=B,y=y,X=X,C=C,ut=ut,pars=pars)

# i is ensemble member; t is timestep 
i=1
t=2

# set up Y vector for which we concatonate parameters, states, and observed data 
Y<-array(NA,c(7,1,nStep,nEn))
Y[1,1,1,]<-rVec # r20 parameter 
Y[2,1,1,]<-rVec_fast # r20 labile parameter 
Y[3,1,1,]<-fracVec # fraction labile of loaded DOC
Y[4,1,1,]<-z$X[1,1,1,] # DIC state  
Y[5,1,1,]<-z$X[2,1,1,] # DOC recalcitrant state  
Y[6,1,1,]<-z$X[3,1,1,] # DOC labile state  
Y[7,1,1,]<-z$X[4,1,1,] # DOC total state  


#Iterate through time
for(t in 2:nStep){
  for(i in 1:nEn){
    # Forecasting; need to update parameters, 
    z$pars[1:3,1,t,i]<-Y[1:3,1,t-1,i] # r20
    z$X[1:4,1,t-1,i]<-Y[4:7,1,t-1,i] # updating state variables from Y 
    
    #Predictions
    z$X[,,t,i]<-z$B[,,t-1,i]%*%z$X[,,t-1,i] + z$C[,,t-1,i]%*%z$ut[,,t-1,i] # forecasting state variable predictions 
    z$B[,,t,i]<-matrix(c(1-data2$QoutInt[t]/data2$epiVol[t]-data2$kCO2[t]/data2$thermo.depth[t], # zmix is fraction of pool available for exchange with atm; parameters estimates are the same as previous time step
             rFunc(z$pars[1,1,t,i],data2$wtr[t]),rFunc(z$pars[2,1,t,i],data2$wtr[t]),0,
             0,1-rFunc(z$pars[1,1,t,i],data2$wtr[t])-data2$QoutInt[t]/data2$epiVol[t],0,0,
             0,0,1-rFunc(z$pars[2,1,t,i],data2$wtr[t])-data2$QoutInt[t]/data2$epiVol[t],0,
             0,1-rFunc(z$pars[1,1,t,i],data2$wtr[t])-data2$QoutInt[t]/data2$epiVol[t],1-rFunc(z$pars[2,1,t,i],data2$wtr[t])-data2$QoutInt[t]/data2$epiVol[t],0),
           nrow=4,byrow=T)
    z$y[,,t,i]<-z$y[,,t,i] # observation of states stay the same 
    # parameters are of previous timestep for matrix C, parameters of covariates [2x6]
    z$C[,,t,i]<-matrix(c(data2$kCO2[t],1,-data2$GPP[t],0,0,0,0,0,0,data2$GPP[t],0,(1-z$pars[3,1,t,i]),
             0,0,0,0,data2$GPP[t],z$pars[3,1,t,i],
             0,0,0,data2$GPP[t],data2$GPP[t],1),nrow=4,byrow=T)
    z$ut[,,t,i] <- matrix(c(data2$DICeq[t]*data2$epiVol[t]/data2$thermo.depth[t],
             data2$dicIn[t],(1-GPPrespired),
             exude,exudeLabile,data2$docIn[t]),ncol=1)
    
    # forecast Y vector 
    Y[1:3,1,t,i]<-Y[1:3,1,t-1,i] #r20, r20_fast, fraction labile loaded DOC parameters same as previous timestep
    Y[4:7,1,t,i]<-z$X[1:4,1,t,i] #forecasted states 

    } # end forecast for each ensemble for timestep t
    
    #begin data assimilation if there are any observations 
    if(any(!is.na(z$y[,,t,i]))==TRUE){ # update vector as long as there is one observation of state 
      
      #mean of vector Y for all ensembles at time step t 
      YMean<-matrix(apply(Y[,,t,],MARGIN = 1,FUN=mean),nrow=length(Y[,1,1,1]))
      delta_Y<-Y[,,t,]-matrix(rep(YMean,nEn),nrow=length(Y[,1,1,1]))# difference in ensemble state and mean of all ensemble states 

      K<-((1/(nEn-1))*delta_Y%*%t(delta_Y)%*%t(h[,,t]))%*%
        qr.solve(((1/(nEn-1))*h[,,t]%*%delta_Y%*%t(delta_Y)%*%t(h[,,t])+H[,,t]))   # Kalman gain 7x2 matrix; 7x3 if iota is included 

      # yObs as a 2x1 matrix  ; 3x1 if iota included 
      yObs<-array(NA,c(2,1,nEn))
      yObs[1,1,]<-z$y[1,1,t,1] #DIC obs
      yObs[2,1,]<-z$y[2,1,t,1] #DOC obs
      yObs[,,]<-ifelse(is.na(yObs[,,]),0,yObs[,,])
      
      # update Y vector 
      for(i in 1:nEn){
        Y[,,t,i]<-Y[,,t,i]+K%*%(yObs[,,i]-h[,,t]%*%Y[,,t,i])
      }
      
      # checking for parameter convergence 
      parsMean<-matrix(apply(Y[1:length(z$pars[,1,1,1]),,t,],MARGIN = 1,FUN=mean),nrow=length(z$pars[,1,1,1]))
      pars_matrix<-array(NA,dim=c(length(z$pars[,1,1,1]),length(z$pars[,1,1,1]),nEn))
      for(i in 1:nEn){
        delta_pars<-(Y[1:length(z$pars[,1,1,1]),1,t,i]-parsMean)
        pars_matrix[,,i]<-delta_pars%*%t(delta_pars)
      }
      
      Ptemp<-matrix(0,ncol=length(z$pars[,1,1,1]),nrow=length(z$pars[,1,1,1]))
      for(i in 1:nEn){
        Ptemp<-Ptemp+pars_matrix[,,i]
      }
      
      P[,,t]<-(1/(nEn-1))*Ptemp
    }

} # End iteration

# 
# # plotting ***************************************************
windows()
DOCout<-apply(Y[7,1,,],MARGIN = 1,FUN=mean)
ylim=range(Y[7,1,,]/data2$epiVol*12)
plot(DOCout/data2$epiVol*12,type='l',ylim=ylim,ylab='DOC mg/L')
for(i in 1:nEn){
  lines(Y[7,1,,i]/data2$epiVol*12,col='gray',ylab='')
}
lines(DOCout/data2$epiVol*12,ylab='')
par(new=T)
plot(z$y[2,1,,1]/data2$epiVol*12,ylim=ylim,col='red',pch=16,ylab='')
arrows(seq(1:nStep),(z$y[2,1,,1]/data2$epiVol*12)-docSD*12,seq(1:nStep),(z$y[2,1,,1]/data2$epiVol*12)+docSD*12,code=3,length=0.1,angle=90,col='red')

windows()
DICout<-apply(z$X[1,1,,],MARGIN = 1,FUN=mean)
ylim=range(z$X[1,1,,]/data2$epiVol*12)
plot(DICout/data2$epiVol*12,type='l',ylim=ylim,ylab='DIC mg/L')
for(i in 1:nEn){
  lines(z$X[1,1,,i]/data2$epiVol*12,col='gray',ylab='')
}
lines(DICout/data2$epiVol*12,ylab='')
par(new=T)
plot(z$y[1,1,,1]/data2$epiVol*12,ylim=ylim,col='red',pch=16,ylab='')
arrows(seq(1:nStep),(z$y[1,1,,1]/data2$epiVol*12)-dicSD*12,seq(1:nStep),(z$y[1,1,,1]/data2$epiVol*12)+dicSD*12,code=3,length=0.1,angle=90,col='red')

windows()
rOut<-apply(rFunc(Y[1,1,,],data2$wtr),MARGIN = 1,FUN=mean)
ylim=range(rFunc(Y[1,1,,],data2$wtr))
plot(rOut,type='l',ylim=ylim,ylab='r day-1')
for(i in 1:nEn){
  lines(rFunc(Y[1,1,,i],data2$wtr),col='gray',ylab='')
}
lines(rOut,ylab='')
lines(apply(Y[1,1,,],MARGIN=1,FUN=mean),col='red') # r20 is in red

windows()
rOut<-apply(rFunc(Y[2,1,,],data2$wtr),MARGIN = 1,FUN=mean)
ylim=range(rFunc(Y[2,1,,],data2$wtr))
plot(rOut,type='l',ylim=ylim,ylab='r fast day-1')
for(i in 1:nEn){
  lines(rFunc(Y[2,1,,i],data2$wtr),col='gray',ylab='')
}
lines(rOut,ylab='')
lines(apply(Y[2,1,,],MARGIN=1,FUN=mean),col='red') # r20 is in red

windows()
rOut<-apply(Y[3,1,,],MARGIN = 1,FUN=mean)
ylim=range(Y[3,1,,])
plot(rOut,type='l',ylim=ylim,ylab='frac Labile')
for(i in 1:nEn){
  lines(Y[3,1,,i],col='gray',ylab='')
}
lines(rOut,ylab='')

# windows()
# rOut<-apply(Y[1,1,,],MARGIN = 1,FUN=mean)
# ylim=range(Y[1,1,,])
# plot(rOut,type='l',ylim=ylim,ylab='r day-1')
# for(i in 1:nEn){
#   lines(Y[1,1,,i],col='gray',ylab='')
# }
# lines(rOut,ylab='')


# # parameter convergence 
windows()
plot(P[1,1,]/rOut,pch=16,main='r20 convergence')

# plotting pools not concentration 
windows()
DOCout<-apply(Y[7,1,,],MARGIN = 1,FUN=mean)
ylim=range(c(Y[7,1,,],z$y[2,1,,1]),na.rm=T)
plot(DOCout,type='l',ylim=ylim,ylab='DOC mg/L')
for(i in 1:nEn){
  lines(Y[7,1,,i],col='gray',ylab='')
}
lines(DOCout,ylab='')
par(new=T)
plot(z$y[2,1,,1],ylim=ylim,col='red',pch=16,ylab='')
arrows(seq(1:nStep),z$y[2,1,,1]-docPoolSD,seq(1:nStep),z$y[2,1,,1]+docPoolSD,code=3,length=0.1,angle=90,col='red')

windows()
DICout<-apply(z$X[1,1,,],MARGIN = 1,FUN=mean)
ylim=range(c(z$X[1,1,,],z$y[1,1,,1]),na.rm=T)
plot(DICout,type='l',ylim=ylim,ylab='DIC mg/L')
for(i in 1:nEn){
  lines(z$X[1,1,,i],col='gray',ylab='')
}
lines(DICout,ylab='')
par(new=T)
plot(z$y[1,1,,1],ylim=ylim,col='red',pch=16,ylab='')
arrows(seq(1:nStep),z$y[1,1,,1]-dicPoolSD,seq(1:nStep),z$y[1,1,,1]+dicPoolSD,code=3,length=0.1,angle=90,col='red')

# DOC labile pool 
windows()
DOCout<-apply(Y[6,1,,],MARGIN = 1,FUN=mean)
ylim=range(Y[6,1,,],na.rm=T)
plot(DOCout,type='l',ylim=ylim,ylab='DOC mol')
for(i in 1:nEn){
  lines(Y[6,1,,i],col='gray',ylab='')
}
lines(DOCout,ylab='')
# DOC Recalcitrant pool 
windows()
DOCout<-apply(Y[5,1,,],MARGIN = 1,FUN=mean)
ylim=range(Y[5,1,,],na.rm=T)
plot(DOCout,type='l',ylim=ylim,ylab='DOC mol')
for(i in 1:nEn){
  lines(Y[5,1,,i],col='gray',ylab='')
}
lines(DOCout,ylab='')
# fraction labile pool 
windows()
DOCoutL<-apply(Y[6,1,,],MARGIN = 1,FUN=mean)
DOCoutR<-apply(Y[5,1,,],MARGIN = 1,FUN=mean)
fracLout<-DOCoutL/DOCoutR
ylim=range(fracLout,na.rm=T)
plot(fracLout,type='l',ylim=ylim,ylab='fraction labile pool')
for(i in 1:nEn){
  lines(Y[6,1,,i]/Y[5,1,,i],col='gray',ylab='')
}
lines(fracLout,ylab='')

# emergent r20 (r20 slow + r20 fast )
windows()
r20L<-apply(Y[2,1,,],MARGIN = 1,FUN=mean)
r20R<-apply(Y[1,1,,],MARGIN = 1,FUN=mean)
DOCoutL<-apply(Y[6,1,,],MARGIN = 1,FUN=mean)
DOCoutR<-apply(Y[5,1,,],MARGIN = 1,FUN=mean)
fracLout<-DOCoutL/DOCoutR
r20All<-r20L*fracLout+r20R*(1-fracLout)
ylim=range(r20All,na.rm=T)
ylim[1]<-0
ylim[2]<-ylim[2]*1.3
plot(r20All,type='l',ylim=ylim,ylab='r20 (emergent)')
for(i in 1:nEn){
  curFracLout<-Y[6,1,,i]/Y[5,1,,i]
  curr20All<-Y[2,1,,i]*curFracLout+Y[1,1,,i]*(1-curFracLout)
  lines(curr20All,col='gray',ylab='')
}
lines(r20All,ylab='')



#plotting r vs. SWin
# 
windows()
rOut<-apply(rFunc(Y[1,1,,],data2$wtr),MARGIN = 1,FUN=mean)
ylim=range(rFunc(Y[1,1,,],data2$wtr))
plot(rOut,type='l',ylim=ylim,ylab='r day-1')
for(i in 1:nEn){
  lines(rFunc(Y[1,1,,i],data2$wtr),col='gray',ylab='')
}
lines(rOut,ylab='')
lines(apply(Y[1,1,,],MARGIN=1,FUN=mean),col='red') # r20 is in red
par(new=T)
plot(data2$docIn,type='l',col='brown',ylab='',yaxt='n')
# parameter considered to be converged when sd is 1/4 of initial sd 
dSD0<-sd(Y[1,1,1,]) #initial sd in d20 
dSD<-apply(Y[1,1,,],MARGIN = 1,FUN=sd)
dConverged<-dSD/dSD0 # fraction of initial sd 
# 
# windows()
# plot(ma(data2$NEP[spinUpLength:length(data2$NEP)],n = 7)~apply(Y[1,1,spinUpLength:length(data2$NEP),],MARGIN = 1,FUN=mean),pch=16,cex=1.5)
# abline(lm(ma(data2$NEP[spinUpLength:length(data2$NEP)],n=7)~apply(Y[1,1,spinUpLength:length(data2$NEP),],MARGIN = 1,FUN=mean)),lwd=2,lty=2)
# summary(lm(ma(data2$NEP[spinUpLength:length(data2$NEP)],n=7)~apply(Y[1,1,spinUpLength:length(data2$NEP),],MARGIN = 1,FUN=mean)))
# 
# windows()
# plot(data2$NEP[spinUpLength:length(data2$NEP)]~apply(rFunc(Y[1,1,spinUpLength:length(data2$NEP),],data2$wtr[spinUpLength:length(data2$NEP)]),MARGIN = 1,FUN=mean))
# abline(lm(data2$NEP[spinUpLength:length(data2$NEP)]~apply(rFunc(Y[1,1,spinUpLength:length(data2$NEP),],data2$wtr[spinUpLength:length(data2$NEP)]),MARGIN = 1,FUN=mean)))

# source('/Users/Jake/Desktop/R functions/backgroundR.R')
# windows()
# plot(backgroundR(data2$R,data2$GPP,n=10)$slope[spinUpLength:nStep],type='o')
# par(new=T)
# plot(apply(Y[1,1,spinUpLength:length(data2$NEP),],MARGIN = 1,FUN=mean,type='o',col='red'))
windows()
plot(data2$docIn,type='l',ylab='')
par(new=T)
plot(apply(Y[1,1,,],MARGIN = 1,FUN=mean),type='o',col='red',ylab='',yaxt='n')
axis(4)

# RMSE for states 
windows()
plot(sqrt((apply(Y[3,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12-z$y[1,1,,1]/data2$epiVol*12)^2),
     pch=16,ylab='DIC RMSE',type='o',xlab='time')
windows()
plot(sqrt((apply(Y[6,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12-z$y[2,1,,1]/data2$epiVol*12)^2),
     pch=16,ylab='DOC RMSE',type='o',xlab='time')

# mean RMSE for states 
sqrt(mean((apply(Y[4,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12-z$y[1,1,,1]/data2$epiVol*12)^2,na.rm=T))
sqrt(mean((apply(Y[7,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12-z$y[2,1,,1]/data2$epiVol*12)^2,na.rm=T))

windows()
plot(sqrt((apply(Y[3,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12-z$y[1,1,,1]/data2$epiVol*12)^2)~
  sqrt((apply(Y[6,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12-z$y[2,1,,1]/data2$epiVol*12)^2),
  pch=16,cex=1.5)
# 
# abline(lm(sqrt((apply(Y[2,1,spinUpLength:nStep,],MARGIN = 1,FUN=mean)/data2$epiVol[spinUpLength:nStep]*12-z$y[1,1,spinUpLength:nStep,1]/data2$epiVol[spinUpLength:nStep]*12)^2)~
#   sqrt((apply(Y[3,1,spinUpLength:nStep,],MARGIN = 1,FUN=mean)/data2$epiVol[spinUpLength:nStep]*12-z$y[2,1,spinUpLength:nStep,1]/data2$epiVol[spinUpLength:nStep]*12)^2)))
# 
# summary(lm(sqrt((apply(Y[2,1,spinUpLength:nStep,],MARGIN = 1,FUN=mean)/data2$epiVol[spinUpLength:nStep]*12-z$y[1,1,spinUpLength:nStep,1]/data2$epiVol[spinUpLength:nStep]*12)^2)~
#      sqrt((apply(Y[3,1,spinUpLength:nStep,],MARGIN = 1,FUN=mean)/data2$epiVol[spinUpLength:nStep]*12-z$y[2,1,spinUpLength:nStep,1]/data2$epiVol[spinUpLength:nStep]*12)^2)))
# 
out<-apply(Y[3,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12
out<-cbind(out,(z$y[1,1,,1]/data2$epiVol*12))
resCO2<-out[,1]-out[,2]

out<-apply(Y[6,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12
out<-cbind(out,(z$y[2,1,,1]/data2$epiVol*12))
resDOC<-out[,1]-out[,2]
windows()
plot(resDOC~resCO2,pch=16,cex=2,ylab='residuals DOC (modeled - obs)', xlab='residuals CO2 (modeled - obs)',cex.lab=1.5,cex.axis=2)
summary(lm(resDOC~resCO2),lwd=2,lty=2)

par(mar=c(5,5,3,6))
plot(data2$totalAlloC,type='l',ylab='DOC load',cex.axis=2,cex.lab=1.5,xlab='time step')
par(new=T)
plot(resDOC,pch=16,xlab='',xaxt='n',ylab='',yaxt='n',cex=2)
axis(4,ylab='residuals DOC (modeled - obs)',cex.axis=2)
par(new=T)
plot(resCO2,pch=16,xlab = '',ylab='',xaxt='n',yaxt='n',cex=2,col='red')

plot(data2$totalAlloC~resDOC,cex=2,pch=16)

########################## Plotting for MS ###################################

# Figure 1 
png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig1_DOC_CO2.png', 
    res=300, width=14, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=2
ylim=c(0.5,3)
par(mar=c(5,6,4,2),mfrow=c(1,2))
DOCout<-apply(Y[7,1,,],MARGIN = 1,FUN=mean)
ylim=range(c(Y[7,1,,],z$y[2,1,,1]+docPoolSD,z$y[2,1,,1]-docPoolSD),na.rm=T)
plot(DOCout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),type='l',ylim=ylim,lwd=lwd,cex.axis=1.5,
     ylab=expression(DOC~(mol~C)),xlab='',cex.lab=cex)
for(i in 1:nEn){
  lines(Y[7,1,(spinUpLength+1):nStep,i]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
        col='gray',ylab='',lwd=lwd)
}
lines(DOCout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),ylab='',lwd=lwd)
par(new=T)
plot(z$y[2,1,(spinUpLength+1):nStep,1]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
     ylim=ylim,col='red',pch=16,ylab='',cex=cex,xaxt='n',yaxt='n',xlab='')
arrows(as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[2,1,(spinUpLength+1):nStep,1]-docPoolSD[(spinUpLength+1):nStep],
       as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[2,1,(spinUpLength+1):nStep,1]+docPoolSD[(spinUpLength+1):nStep],
       code=3,length=0.1,angle=90,col='red',lwd=lwd)
legend("topleft", legend=c("Estimated State","Ensemble Mean",'Observed State'),
       col=c('gray','black','red'),pt.bg=c('gray','black','red'), 
       ncol=1,lwd=c(4,4,4),bty='n',lty=c(1,1,1))

DICout<-apply(Y[4,1,,],MARGIN = 1,FUN=mean)
ylim=range(c(Y[4,1,,],z$y[1,1,,1]+dicPoolSD,z$y[1,1,,1]-dicPoolSD),na.rm=T)
plot(DICout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),type='l',ylim=ylim,lwd=lwd,cex.axis=1.5,
     ylab=expression(CO[2]~(mol~C)),xlab='',cex.lab=cex)
for(i in 1:nEn){
  lines(Y[4,1,(spinUpLength+1):nStep,i]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
        col='gray',ylab='',lwd=lwd)
}
lines(DICout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),ylab='',lwd=lwd)
par(new=T)
plot(z$y[1,1,(spinUpLength+1):nStep,1]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
     ylim=ylim,col='red',pch=16,ylab='',cex=cex,xaxt='n',yaxt='n',xlab='')
arrows(as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[1,1,(spinUpLength+1):nStep,1]-dicPoolSD[(spinUpLength+1):nStep],
       as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[1,1,(spinUpLength+1):nStep,1]+dicPoolSD[(spinUpLength+1):nStep],
       code=3,length=0.1,angle=90,col='red',lwd=lwd)
dev.off()

# Figure 1 plotting concentration 
png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig1_DOC_CO2_concentration.png', 
    res=300, width=14, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=2
ylim=c(0.5,3)
par(mar=c(5,6,4,2),mfrow=c(1,2))
DOCout<-apply(Y[7,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
ylim=range(c(Y[7,1,,]/data2$epiVol*12,z$y[2,1,,1]/data2$epiVol*12+docPoolSD/data2$epiVol*12,z$y[2,1,,1]/data2$epiVol*12-docPoolSD/data2$epiVol*12),na.rm=T)
plot(DOCout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),type='l',ylim=ylim,lwd=lwd,cex.axis=1.5,
     ylab=expression(DOC~(mg~C~L^-1)),xlab='',cex.lab=cex)
for(i in 1:nEn){
  lines(Y[7,1,(spinUpLength+1):nStep,i]/data2$epiVol*12~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
        col='gray',ylab='',lwd=lwd)
}
lines(DOCout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),ylab='',lwd=lwd)
par(new=T)
plot(z$y[2,1,(spinUpLength+1):nStep,1]/data2$epiVol*12~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
     ylim=ylim,col='red',pch=16,ylab='',cex=cex,xaxt='n',yaxt='n',xlab='')
arrows(as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[2,1,(spinUpLength+1):nStep,1]/data2$epiVol*12-docPoolSD[(spinUpLength+1):nStep]/data2$epiVol*12,
       as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[2,1,(spinUpLength+1):nStep,1]/data2$epiVol*12+docPoolSD[(spinUpLength+1):nStep]/data2$epiVol*12,
       code=3,length=0.1,angle=90,col='red',lwd=lwd)
legend("topright", legend=c("Estimated State","Ensemble Mean",'Observed State'),
       col=c('gray','black','red'),pt.bg=c('gray','black','red'), 
       ncol=1,lwd=c(4,4,4),bty='n',lty=c(1,1,1))

DICout<-apply(Y[4,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
ylim=range(c(Y[4,1,,]/data2$epiVol*12,z$y[1,1,,1]/data2$epiVol*12+dicPoolSD/data2$epiVol*12,z$y[1,1,,1]/data2$epiVol*12-dicPoolSD/data2$epiVol*12),na.rm=T)
plot(DICout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),type='l',ylim=ylim,lwd=lwd,cex.axis=1.5,
     ylab=expression(CO[2]~(mg~C~L^-1)),xlab='',cex.lab=cex)
for(i in 1:nEn){
  lines(Y[4,1,(spinUpLength+1):nStep,i]/data2$epiVol*12~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
        col='gray',ylab='',lwd=lwd)
}
lines(DICout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),ylab='',lwd=lwd)
par(new=T)
plot(z$y[1,1,(spinUpLength+1):nStep,1]/data2$epiVol*12~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
     ylim=ylim,col='red',pch=16,ylab='',cex=cex,xaxt='n',yaxt='n',xlab='')
arrows(as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[1,1,(spinUpLength+1):nStep,1]/data2$epiVol*12-dicPoolSD[(spinUpLength+1):nStep]/data2$epiVol*12,
       as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[1,1,(spinUpLength+1):nStep,1]/data2$epiVol*12+dicPoolSD[(spinUpLength+1):nStep]/data2$epiVol*12,
       code=3,length=0.1,angle=90,col='red',lwd=lwd)
dev.off()

# Figure 1 plotting concentration observed vs. predicted with error bars 
png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig1_DOC_CO2_concentration_obsVpred.png', 
    res=300, width=14, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=2
par(mar=c(5,6,4,2),mfrow=c(1,2))
DOCout<-apply(Y[7,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
DOCpredictSD<-apply(Y[7,1,,]/data2$epiVol*12,MARGIN = 1,FUN=sd)
DOCout<-cbind(DOCout,z$y[2,1,(spinUpLength+1):nStep,1]/data2$epiVol*12)
ylim=range(c(na.omit(DOCout+DOCpredictSD),na.omit(DOCout-DOCpredictSD),z$y[2,1,,1]/data2$epiVol*12+docPoolSD/data2$epiVol*12,
             z$y[2,1,,1]/data2$epiVol*12-docPoolSD/data2$epiVol*12),na.rm=T)
xlim=ylim

plot(DOCout[spinUpLength+1:nStep,1]~DOCout[spinUpLength+1:nStep,2],ylim=ylim,cex=2,cex.axis=1.5,xlim=xlim,pch=16,
     ylab=expression(Predicted~DOC~(mg~C~L^-1)),cex.lab=cex,xlab=expression(Observed~DOC~(mg~C~L^-1)))
arrows(DOCout[,2]-docPoolSD[(spinUpLength+1):nStep]/data2$epiVol*12,DOCout[,1], # error bars for observed DOC concentration 
       DOCout[,2]+docPoolSD[(spinUpLength+1):nStep]/data2$epiVol*12,DOCout[,1],
       code=3,length=0.1,angle=90,col='gray60',lwd=lwd)
arrows(DOCout[,2],DOCout[,1]-DOCpredictSD, # error bars for predicted DOC concentration 
       DOCout[,2],DOCout[,1]+DOCpredictSD,
       code=3,length=0.1,angle=90,col='gray60',lwd=lwd)
points(DOCout[spinUpLength+1:nStep,1]~DOCout[spinUpLength+1:nStep,2],cex=2,pch=16)
abline(0,1,lty=2,lwd=2,col='gray')

DICout<-apply(Y[4,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
DICpredictSD<-apply(Y[4,1,,]/data2$epiVol*12,MARGIN = 1,FUN=sd)
DICout<-cbind(DICout,z$y[1,1,(spinUpLength+1):nStep,1]/data2$epiVol*12)
ylim=range(c(na.omit(DICout+DICpredictSD),na.omit(DICout-DICpredictSD),z$y[1,1,,1]/data2$epiVol*12+dicPoolSD/data2$epiVol*12,
             z$y[1,1,,1]/data2$epiVol*12-dicPoolSD/data2$epiVol*12),na.rm=T)
xlim=ylim
plot(DICout[spinUpLength+1:nStep,1]~DICout[spinUpLength+1:nStep,2],ylim=ylim,cex=2,cex.axis=1.5,xlim=xlim,pch=16,
     ylab=expression(Predicted~pCO[2]~(mg~C~L^-1)),cex.lab=cex,xlab=expression(Observed~pCO[2]~(mg~C~L^-1)))
arrows(DICout[,2]-dicPoolSD[(spinUpLength+1):nStep]/data2$epiVol*12,DICout[,1], # error bars for observed DOC concentration 
       DICout[,2]+dicPoolSD[(spinUpLength+1):nStep]/data2$epiVol*12,DICout[,1],
       code=3,length=0.1,angle=90,col='gray60',lwd=lwd)
arrows(DICout[,2],DICout[,1]-DICpredictSD, # error bars for predicted DOC concentration 
       DICout[,2],DICout[,1]+DICpredictSD,
       code=3,length=0.1,angle=90,col='gray60',lwd=lwd)
points(DICout[spinUpLength+1:nStep,1]~DICout[spinUpLength+1:nStep,2],cex=2,pch=16)
abline(0,1,lty=2,lwd=2,col='gray')
dev.off()

# Figure 1 plotting pools of observed vs. predicted with error bars 
png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig1_DOC_CO2_pools_obsVpred.png', 
    res=300, width=14, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=2
par(mar=c(5,6,4,2),mfrow=c(1,2))
DOCout<-apply(Y[7,1,,],MARGIN = 1,FUN=mean)
DOCpredictSD<-apply(Y[7,1,,],MARGIN = 1,FUN=sd)
DOCout<-cbind(DOCout,z$y[2,1,(spinUpLength+1):nStep,1])
ylim=range(c(na.omit(DOCout+DOCpredictSD),na.omit(DOCout-DOCpredictSD),z$y[2,1,,1]+docPoolSD,
             z$y[2,1,,1]-docPoolSD),na.rm=T)
xlim=ylim

plot(DOCout[spinUpLength+1:nStep,1]~DOCout[spinUpLength+1:nStep,2],ylim=ylim,cex=2,cex.axis=1.5,xlim=xlim,pch=16,
     ylab=expression(Predicted~DOC~(mol~C)),cex.lab=cex,xlab=expression(Observed~DOC~(mol~C)))
arrows(DOCout[,2]-docPoolSD[(spinUpLength+1):nStep],DOCout[,1], # error bars for observed DOC concentration 
       DOCout[,2]+docPoolSD[(spinUpLength+1):nStep],DOCout[,1],
       code=3,length=0.1,angle=90,col='gray60',lwd=lwd)
arrows(DOCout[,2],DOCout[,1]-DOCpredictSD, # error bars for predicted DOC concentration 
       DOCout[,2],DOCout[,1]+DOCpredictSD,
       code=3,length=0.1,angle=90,col='gray60',lwd=lwd)
points(DOCout[spinUpLength+1:nStep,1]~DOCout[spinUpLength+1:nStep,2],cex=2,pch=16)
abline(0,1,lty=2,lwd=2,col='gray')

DICout<-apply(Y[4,1,,],MARGIN = 1,FUN=mean)
DICpredictSD<-apply(Y[4,1,,],MARGIN = 1,FUN=sd)
DICout<-cbind(DICout,z$y[1,1,(spinUpLength+1):nStep,1])
ylim=range(c(na.omit(DICout+DICpredictSD),na.omit(DICout-DICpredictSD),z$y[1,1,,1]+dicPoolSD,
             z$y[1,1,,1]-dicPoolSD),na.rm=T)
xlim=ylim
plot(DICout[spinUpLength+1:nStep,1]~DICout[spinUpLength+1:nStep,2],ylim=ylim,cex=2,cex.axis=1.5,xlim=xlim,pch=16,
     ylab=expression(Predicted~pCO[2]~(mol~C)),cex.lab=cex,xlab=expression(Observed~pCO[2]~(mol~C)))
arrows(DICout[,2]-dicPoolSD[(spinUpLength+1):nStep],DICout[,1], # error bars for observed DOC concentration 
       DICout[,2]+dicPoolSD[(spinUpLength+1):nStep],DICout[,1],
       code=3,length=0.1,angle=90,col='gray60',lwd=lwd)
arrows(DICout[,2],DICout[,1]-DICpredictSD, # error bars for predicted DOC concentration 
       DICout[,2],DICout[,1]+DICpredictSD,
       code=3,length=0.1,angle=90,col='gray60',lwd=lwd)
points(DICout[spinUpLength+1:nStep,1]~DICout[spinUpLength+1:nStep,2],cex=2,pch=16)
abline(0,1,lty=2,lwd=2,col='gray')
dev.off()



# Figure 2 
png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig2_d.png', 
    res=300, width=7, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=2
ylim=c(0.5,3)
par(mar=c(5,6,4,2))
rOut<-apply(rFunc(Y[1,1,,],data2$wtr),MARGIN = 1,FUN=mean)
ylim=range(rFunc(Y[1,1,(spinUpLength+1):nStep,],data2$wtr[(spinUpLength+1):nStep]))
plot(rOut[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
     type='l',ylim=ylim,ylab=expression(d~(day^-1)),lwd=lwd,xlab='',cex.lab=cex,cex.axis=1.5)
for(i in 1:nEn){
  lines(rFunc(Y[1,1,(spinUpLength+1):nStep,i],data2$wtr[(spinUpLength+1):nStep])~
          as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),col='gray',ylab='',lwd=lwd,xlab='')
}
lines(rOut[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),ylab='',lwd=lwd,
      xlab='')
# plotting red for when parameter is converged 
lines(rOut[dConverged<0.25]~as.POSIXct(data2$datetime[dConverged<0.25]),ylab='',lwd=lwd,col='red',
      xlab='')
dev.off()

# Figure 3 
png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig3_d20_carbonLoad.png', 
    res=300, width=10, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=2
ylim=c(0.5,3)
par(mar=c(5,6,4,6))
rOut<-apply(Y[1,1,,],MARGIN = 1,FUN=mean)
ylim=range(rFunc(Y[1,1,(spinUpLength+1):nStep,],data2$wtr[(spinUpLength+1):nStep]))
ylim=c(0.00,0.015)
plot(rOut[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
     type='l',ylim=ylim,ylab=expression(d[20]~(day^-1)),lwd=lwd,xlab='',cex.lab=cex,cex.axis=cex)
for(i in 1:nEn){
  lines(Y[1,1,(spinUpLength+1):nStep,i]~
          as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),col='gray',ylab='',lwd=lwd,xlab='')
}
lines(rOut[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),ylab='',lwd=lwd,
      xlab='')
par(new=T)
plot(data2$docIn[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),type='l',lwd=4,
     col='blue',yaxt='n',xaxt='n',ylab='',xlab='')
axis(4,cex.axis=cex)
mtext(expression(loadDOC~(mol~C~day^-1)),4,line = 4,cex=cex)
legend("top", legend=c("d20 Estimate",'d20 Ensemble Mean','DOC loading'),
       col=c('gray','black','blue'),pt.bg=c('gray','black','blue'), 
       ncol=1,lwd=c(4,4,4),bty='n',lty=c(1,1,1))

dev.off()


# Figure 4 
# moving average weight of GPP and R to estimate NEP weighted 
GPPweight<-ma.weighted(data2$GPP[(spinUpLength+1):nStep],
                       cv = data2$sdGPP[(spinUpLength+1):nStep]/data2$GPP[(spinUpLength+1):nStep],n=7)
Rweight<-ma.weighted(data2$R[(spinUpLength+1):nStep],
                       cv = data2$sdR[(spinUpLength+1):nStep]/data2$R[(spinUpLength+1):nStep],n=7)
NEPweight<-GPPweight-Rweight

png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig4_d20_NEP.png', 
    res=300, width=7, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=2
ylim=c(0.5,3)
par(mar=c(5,6,4,6))
plot(NEPweight~apply(Y[1,1,(spinUpLength+1):nStep,],MARGIN = 1,FUN=mean),
     pch=16,cex=cex,cex.axis=cex,cex.lab=cex,ylab=expression(NEP~(mol~C~day^-1)),
     xlab=expression(d[20]~(day^-1)))
abline(lm(NEPweight~apply(Y[1,1,(spinUpLength+1):nStep,],MARGIN = 1,FUN=mean)),lwd=lwd,lty=2)
dev.off()

summary(lm(NEPweight~apply(Y[1,1,(spinUpLength+1):nStep,],MARGIN = 1,FUN=mean)))

# Fig 5
png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig5_d20_convergence.png', 
    res=300, width=7, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=2
ylim=c(0.5,3)
par(mar=c(5,6,4,6))
plot(P[1,1,][P[1,1,]>0]~seq(1:nStep)[P[1,1,]>0],
     pch=16,cex=cex,cex.axis=cex,cex.lab=cex,ylab=expression(d[20]~Convergence~(P[t])),
     xlab=expression(Model~Timestep))
dev.off()


# Figure 6 plotting model uncertainty versus observed uncertainty  
png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig6_DOC_CO2_uncertainty_obsVpred.png', 
    res=300, width=7, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=4
par(mar=c(5,6,4,2))
DOCout<-apply(Y[4,1,,],MARGIN = 1,FUN=mean)
DOCpredictSD<-apply(Y[4,1,,],MARGIN = 1,FUN=sd)
DOCcv<-DOCpredictSD/DOCout
DOCcv<-cbind(DOCcv,docPoolSD/z$y[2,1,,1])
DICout<-apply(Y[2,1,,],MARGIN = 1,FUN=mean)
DICpredictSD<-apply(Y[2,1,,],MARGIN = 1,FUN=sd)
DICcv<-DICpredictSD/DICout
DICcv<-cbind(DICcv,dicPoolSD/z$y[1,1,,1])
ylim=range(c(na.omit(DOCcv),na.omit(DICcv)),na.rm=T)

plot(DOCcv[,1]~as.POSIXct(data2$datetime),ylim=ylim,cex=2,cex.axis=1.5,pch=16,type='l',lwd=lwd,
     ylab=expression(CV~DOC~and~pCO[2]~(mol~C)),cex.lab=cex,xlab='')
lines(DICcv[,1]~as.POSIXct(data2$datetime),ylim=ylim,cex=2,cex.axis=1.5,pch=16,type='l',lwd=lwd,col='gray60')
abline(h=mean(DOCcv[,2],na.rm=T),lwd=lwd,lty=2)
abline(h=mean(DICcv[,2],na.rm=T),lwd=lwd,lty=2,col='gray60')
legend("right", legend=c("DOC Pred CV",'pCO2 Pred CV'),
       col=c('black','gray60'),pt.bg=c('black','gray60'), 
       ncol=1,lwd=c(4,4),bty='n',lty=c(1,1))
dev.off()

# does d times DOC equal NEP? 
rOut<-apply(rFunc(Y[1,1,,],data2$wtr),MARGIN = 1,FUN=mean)
predNEP<--DOCout*rOut-((1-GPPrespired)*data2$GPP)
predNEP<-cbind(predNEP,NEPweight)
plot(predNEP[,2])
points(predNEP[,1])
plot(predNEP[,1]~predNEP[,2],pch=16)
abline(lm(predNEP[,1]~predNEP[,2]))



png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig4_PredNEP_NEP.png', 
    res=300, width=7, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=2
par(mar=c(5,6,4,6))
ylim=range(predNEP,na.rm=T)
xlim=ylim
plot(predNEP[,2]~predNEP[,1],
     pch=16,cex=cex,cex.axis=cex,cex.lab=cex,ylab=expression(Observed~NEP~(mol~C~day^-1)),
     xlab=expression(Predicted~NEP~(mol~C~day^-1)))
# abline(lm(predNEP[,2]~predNEP[,1]),lwd=lwd,lty=2)
dev.off()

# background R 
backR<-backgroundR(data2$R,data2$GPP,n=60)
plot(backR$intercept,type='l')
par(new=T)
plot(predNEP[,1],type='l',col='red')


