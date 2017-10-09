# creating fake data for MDF 
# JAZ; 2016-11-22 


rm(list=ls())

load('/afs/crc.nd.edu/user/j/jzwart/ModelDataFusion/R_Data/EnKF_LongData_20170223.RData')

splitFunc<-function(epiDens,streamDens,fracIn){ # function that tells how much load goes into epi 
  fracInEpi=exp(-fracIn*(streamDens-epiDens))
  return(fracInEpi)
}
splitFunc<-function(epiDens,streamDens,fracIn){ # function that tells how much load goes into epi 
  fracInEpi=exp(-fracIn*(streamDens-epiDens))
  return(fracIn)
}

require(deSolve)
require(LakeMetabolizer)
require(snow)
require(sp)
require(rgeos)
require(parallel)

# spin up, just repeating the first X days at the begining of the timeseries; autocorrelation effects?? 
spinUpLength<-0
spinUp<-data[1:spinUpLength,]
data2<-rbind(spinUp,data)
# should we cut down to only when high frequency discharge out was known? occurs on 2014-06-02
data2<-data2[min(which(!is.na(data2$highFreqWaterHeight))):max(which(!is.na(data2$highFreqWaterHeight))),]
data2<-data2[!is.na(data2$ma_gpp),]
data2<-data2[data2$datetime<as.POSIXct('2015-01-01'),] # only keeping 2014
data2<-data2[min(which(!is.na(data2$doc))):nrow(data2),]

# adding in trash pump discharge out starting on Aug. 22, 2015 until the end of time series in 2015 
data2$QoutInt[data2$datetime>=as.POSIXct('2015-08-22')]<-data2$QoutInt[data2$datetime>=as.POSIXct('2015-08-22')]+400

fracLabile0<-0.01 # fraction of initial DOC pool that is labile
fracLabile<-0.08 # fraction loaded t-DOC that is labile; estimate from Berggren et al. 2010 Eco Letts & ISME

# need some initial conditions for observations (draw from normal distribution?)
data2$dic[1]<-data2$dic[min(which(!is.na(data2$dic)))]
data2$doc[1]<-data2$doc[min(which(!is.na(data2$doc)))]

## *********************** EnKF ***************************## Gao et al. 2011 is a useful reference 
nEn<-100 # number of ensembles 
nStep<-length(data2$datetime)

# draws from priors to create ensemble 
parGuess <- c(0.004,0.3,0.3,0.1) #r20; units: day-1; fraction labile of loaded DOC; fraction inlet that goes into epi; turnover rate parameters set constant frac labile estimated 
min<-c(0.004,0.31,0.5,0.9)
max<-c(0.004,0.31,0.5,0.9)

rPDF<-abs(rnorm(n=nEn,mean = parGuess[1],sd = (max[1]-min[1])/5)) # max-min / 5 is rule of thumb sd; forcing positive for negative draws 
rPDF_fast<-abs(rnorm(n=nEn,mean = parGuess[2],sd = (max[2]-min[2])/5))
fracPDF<-abs(rnorm(n=nEn,mean=parGuess[3],sd=(max[3]-min[3])/5))
fracInPDF<-abs(rnorm(n=nEn,mean=parGuess[4],sd=(max[4]-min[4])/5))

# setting up initial parameter values for all ensemble members 
rVec<-matrix(rPDF) # each row is an ensemble member 
rVec_fast<-matrix(rPDF_fast)
fracVec<-matrix(fracPDF)
fracInVec<-matrix(fracInPDF)
halfSat<-4 # half saturation constant for turnover rate of doc 

mm<-function(r,doc,halfSat,epiVol){
  doc<-doc/epiVol*12
  rout<-r*doc/(halfSat+doc)
  return(rout)
}

# setting up state values for all ensemble members (all the same initial pool size )
# updated on 2016-07-16 to draw from a normal distribution using first observation and SD
dicVec<-data2$dic
docVec<-data2$doc
if(is.na(dicVec[1])){
  dicVec[1]<-dicVec[min(which(!is.na(dicVec)))]
  docVec[1]<-docVec[min(which(!is.na(docVec)))]
}

data2$entrainHypo<-as.numeric(data2$entrainVol<0)
data2$entrainEpi<-as.numeric(data2$entrainVol>0)

# initial B transition matrix for each ensemble  
B<-array(NA,dim=c(4,4,nStep,nEn)) # array of transition matrix B[a,b,c,d]; 
# where [a,b] is transition matrix DIC & DOCr & DOCl, c=timeStep, and d=ensemble member 
# initial parameters of B at timestep 1 
for(i in 1:nEn){
  B[,,1,i]<-matrix(c(1-data2$QoutInt[1]/data2$epiVol[1]-data2$kCO2[1]/data2$thermo.depth[1]+
                       (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                     mm(rFunc(rVec[i],data2$wtr[1]),docVec[1],halfSat,data2$epiVol[1]),mm(rFunc(rVec_fast[i],data2$wtr[1]),docVec[1],halfSat,data2$epiVol[1]),0,
                     0,1-mm(rFunc(rVec[i],data2$wtr[1]),docVec[1],halfSat,data2$epiVol[1])-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],0,0,
                     0,0,1-mm(rFunc(rVec_fast[i],data2$wtr[1]),docVec[1],halfSat,data2$epiVol[1])-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],0,
                     0,1-mm(rFunc(rVec[i],data2$wtr[1]),docVec[1],halfSat,data2$epiVol[1])-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                     1-mm(rFunc(rVec_fast[i],data2$wtr[1]),docVec[1],halfSat,data2$epiVol[1])-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],0),
                   nrow=4,byrow=T)
}

# estimate of observation error covariance matrix
# error in DOC / DIC pool esitmate comes from concentration estimates and volume of epilimnion
# DOC / DIC concentration error from replication day in WL on 2014-07-30. CV for DOC is 0.1325858, DIC is 0.1367684
# Epi: area sd = 1000m2; depth sd = 0.5 m;
# iota sd from metab estimates 
docSDadjust<-1
dicSDadjust<-1
docSD<-(data2$doc/data2$epiVol)*0.1325858 # DOC concentration sd in mol C
docSD<-ifelse(is.na(docSD),docSD,docSD[data2$datetime=='2014-07-30']) # making SD the same for all observations; not based on concentration 2016-11-22 
docSD<-docSD*docSDadjust # DOC concentration sd in mol C; modifying for sensativity analysis
dicSD<-(data2$dic/data2$epiVol)*0.1367684 # DIC concentration sd in mol C 
dicSD<-ifelse(is.na(dicSD),dicSD,dicSD[data2$datetime=='2014-07-30']) # making SD the same for all observations; not based on concentration 2016-11-22 
dicSD<-dicSD*dicSDadjust
areaSD<-4000 # constant SD in m 
depthSD<-0.25 # constant SD in m 

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
h<-array(0,dim=c(2,8,nStep))
for(i in 1:nStep){
  h[1,4,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic 
  h[2,8,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) #doc total (we only have data on total DOC pool)
}

P <- array(0,dim=c(4,4,nStep))

#Define matrix C, parameters of covariates [2x6]
C<-array(NA,dim=c(4,7,nStep,nEn)) # array of transition matrix C[a,b,c,d]; 
#intializing first time step 
for(i in 1:nEn){
  C[,,1,i]<-matrix(c(data2$kCO2[1],1,-data2$ma_gpp[1],0,0,0,data2$hypo_dicInt[1],0,0,0,data2$ma_gpp[1],0,(1-fracVec[i]),data2$hypo_docInt[1]*(1-fracLabile0),
                     0,0,0,0,data2$ma_gpp[1],fracVec[i],data2$hypo_docInt[1]*fracLabile0,
                     0,0,0,data2$ma_gpp[1],data2$ma_gpp[1],1,data2$hypo_docInt[1]),nrow=4,byrow=T)
}

ut<-array(NA,dim=c(7,1,nStep,nEn)) # array of transition matrix C[a,b,c,d]; 
#intializing first time step 
for(i in 1:nEn){
  ut[,,1,i]<-matrix(c(data2$DICeq[1]*data2$epiVol[1]/data2$thermo.depth[1],
                      data2$dicIn[1]-data2$streamWaterdisch[1]*data2$Inlet_dic[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])),(1-GPPrespired),
                      exude,exudeLabile,data2$docIn[1]-data2$streamDOCdisch[1]/12/1000*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])),
                      (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainEpi[1]),ncol=1)
}

pars<-array(rep(NA,nEn),dim=c(4,1,nStep,nEn)) # parameters: r20
pars[1,1,1,]<-rVec
pars[2,1,1,]<-rVec_fast
pars[3,1,1,]<-fracVec
pars[4,1,1,]<-fracInVec

# set up a list for all matrices 
z=list(B=B,y=y,X=X,C=C,ut=ut,pars=pars)

# i is ensemble member; t is timestep 
i=1
t=2

# set up Y vector for which we concatonate parameters, states, and observed data 
Y<-array(NA,c(8,1,nStep,nEn))
Y[1,1,1,]<-rVec # r20 parameter 
Y[2,1,1,]<-rVec_fast # r20 labile parameter 
Y[3,1,1,]<-fracVec # fraction labile of loaded DOC
Y[4,1,1,]<-fracInVec # fraction of inlet into epi
Y[5,1,1,]<-z$X[1,1,1,] # DIC state  
Y[6,1,1,]<-z$X[2,1,1,] # DOC recalcitrant state  
Y[7,1,1,]<-z$X[3,1,1,] # DOC labile state  
Y[8,1,1,]<-z$X[4,1,1,] # DOC total state  


#Iterate through time
for(t in 2:nStep){
  # Forecasting; need to update parameters, 
  z$pars[1:4,1,t,i]<-Y[1:4,1,t-1,i] # r20
  z$X[1:4,1,t-1,i]<-Y[5:8,1,t-1,i] # updating state variables from Y 
  
  #Predictions
  z$X[,,t,i]<-z$B[,,t-1,i]%*%z$X[,,t-1,i] + z$C[,,t-1,i]%*%z$ut[,,t-1,i] # forecasting state variable predictions 
  z$B[,,t,i]<-matrix(c(1-data2$QoutInt[t]/data2$epiVol[t]-data2$kCO2[t]/data2$thermo.depth[t]+
                         (data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t], # zmix is fraction of pool available for exchange with atm; parameters estimates are the same as previous time step
                       mm(rFunc(z$pars[1,1,t,i],data2$wtr[t]),z$X[2,1,t,i],halfSat,data2$epiVol[t]),mm(rFunc(z$pars[2,1,t,i],data2$wtr[t]),z$X[2,1,t,i],halfSat,data2$epiVol[t]),0,
                       0,1-mm(rFunc(z$pars[1,1,t,i],data2$wtr[t]),z$X[2,1,t,i],halfSat,data2$epiVol[t])-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],0,0,
                       0,0,1-mm(rFunc(z$pars[2,1,t,i],data2$wtr[t]),z$X[2,1,t,i],halfSat,data2$epiVol[t])-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],0,
                       0,1-mm(rFunc(z$pars[1,1,t,i],data2$wtr[t]),z$X[2,1,t,i],halfSat,data2$epiVol[t])-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],
                       1-mm(rFunc(z$pars[2,1,t,i],data2$wtr[t]),z$X[2,1,t,i],halfSat,data2$epiVol[t])-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],0),
                     nrow=4,byrow=T)
  z$y[,,t,i]<-z$y[,,t,i] # observation of states stay the same 
  # parameters are of previous timestep for matrix C, parameters of covariates [2x6]
  z$C[,,t,i]<-matrix(c(data2$kCO2[t],1,-data2$ma_gpp[t],0,0,0,data2$hypo_dicInt[t],0,0,0,data2$ma_gpp[t],0,(1-z$pars[3,1,t,i]),data2$hypo_docInt[t]*(1-fracLabile0),
                       0,0,0,0,data2$ma_gpp[t],z$pars[3,1,t,i],data2$hypo_docInt[t]*fracLabile0,
                       0,0,0,data2$ma_gpp[t],data2$ma_gpp[t],1,data2$hypo_docInt[t]),nrow=4,byrow=T)
  z$ut[,,t,i] <- matrix(c(data2$DICeq[t]*data2$epiVol[t]/data2$thermo.depth[t],
                          data2$dicIn[t]-data2$streamWaterdisch[t]*data2$Inlet_dic[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])),(1-GPPrespired),
                          exude,exudeLabile,data2$docIn[t]-data2$streamDOCdisch[t]/12/1000*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])),
                          (data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainEpi[t]),ncol=1)
  
  # forecast Y vector 
  Y[1:4,1,t,i]<-Y[1:4,1,t-1,i] #r20, r20_fast, fraction labile loaded DOC parameters same as previous timestep
  Y[5:8,1,t,i]<-z$X[1:4,1,t,i] #forecasted states 
} # End iteration

# True state is in Y vector with decay rate of DOC set to 0.005 day-1 
# plot(Y[1,1,,1])
# plot(Y[2,1,,1])
# plot(Y[3,1,,1])
# plot(Y[4,1,,1])# DIC state 
# plot(Y[5,1,,1])
# plot(Y[6,1,,1])
# plot(Y[7,1,,1])# DOC total state  

true<-Y

obs<-1 # adjustment factor for observation error 
freq<-c(1,7,14,21,28,35,63) # frequency of observations 
reps<-c(1,2,3,4,5,6,10) # number of replicate sampling 
n<-100 ## number of times to run MDF to account for random sampling 

source('/afs/crc.nd.edu/user/j/jzwart/ModelDataFusion/R_Code/wrongModels_20171009.R')

starting<-startingValues()
models<-c('oneDOCpool_noMM_noTemp','oneDOCpool_noMM_Temp','oneDOCpool_MM_noTemp','oneDOCpool_MM_Temp','twoDOCpools_noMM_noTemp','twoDOCpools_noMM_Temp',
          'twoDOCpools_MM_noTemp','twoDOCpools_MM_Temp')

runs<-expand.grid(obs,freq,reps,models) # data frame with all combinations of factors for model runs
runs$unique<-paste(runs[,4],runs[,1],runs[,2],runs[,3],sep='-')

sim<-function(run){
  cur<-strsplit(run,'-',fixed=T)
  curMod<-cur[[1]][1]
  curObs<-as.numeric(cur[[1]][2])
  curFreq<-as.numeric(cur[[1]][3])
  curReps<-as.numeric(cur[[1]][4])
  if(length(grep('one',curMod))>0){
    out<-do.call(curMod,args = list(obs=curObs,freq=curFreq,reps=curReps,rVec=starting$rVec,fracInVec=starting$fracInVec,
                                    covar_inflat_vec=starting$covar_inflat_vec,n=n))
    curOut<-data.frame(model=curMod,obs=curObs,freq=curFreq,reps=curReps,DOCrmse=out$DOCrmse,DOCbias=out$DOCbias,DOCprecision=out$DOCprecision,
                       DICrmse=out$DICrmse,DICbias=out$DICbias,DICprecision=out$DICprecision,
                       obsDOCrmse=out$obsDOCrmse,obsDOCbias=out$obsDOCbias,obsDOCprecision=out$obsDOCprecision,
                       obsDICrmse=out$obsDICrmse,obsDICbias=out$obsDICbias,obsDICprecision=out$obsDICprecision,
                       mdfobsDOCrmse=out$mdfobsDOCrmse,mdfobsDOCbias=out$mdfobsDOCbias,mdfobsDOCprecision=out$mdfobsDOCprecision,
                       mdfobsDICrmse=out$mdfobsDICrmse,mdfobsDICbias=out$mdfobsDICbias,mdfobsDICprecision=out$mdfobsDICprecision,
                       mdfDOCcv=out$mdfDOCcv,mdfDICcv=out$mdfDICcv)
    file<-paste(curMod,curObs,curFreq,curReps,sep='_')
    save(curOut,file=paste('/afs/crc.nd.edu/user/j/jzwart/ModelDataFusion/Results/',file,'.Rda',sep = ''))
  }else if(length(grep('two',curMod))>0){
    out<-do.call(curMod,args = list(obs=curObs,freq=curFreq,reps=curReps,rVec_slow = starting$rVec_slow,
                                    rVec_fast = starting$rVec_fast,fracVec = starting$fracVec,fracInVec=starting$fracInVec,
                                    covar_inflat_vec=starting$covar_inflat_vec,n=n))
    curOut<-data.frame(model=curMod,obs=curObs,freq=curFreq,reps=curReps,DOCrmse=out$DOCrmse,DOCbias=out$DOCbias,DOCprecision=out$DOCprecision,
                       DICrmse=out$DICrmse,DICbias=out$DICbias,DICprecision=out$DICprecision,
                       obsDOCrmse=out$obsDOCrmse,obsDOCbias=out$obsDOCbias,obsDOCprecision=out$obsDOCprecision,
                       obsDICrmse=out$obsDICrmse,obsDICbias=out$obsDICbias,obsDICprecision=out$obsDICprecision,
                       mdfobsDOCrmse=out$mdfobsDOCrmse,mdfobsDOCbias=out$mdfobsDOCbias,mdfobsDOCprecision=out$mdfobsDOCprecision,
                       mdfobsDICrmse=out$mdfobsDICrmse,mdfobsDICbias=out$mdfobsDICbias,mdfobsDICprecision=out$mdfobsDICprecision,
                       mdfDOCcv=out$mdfDOCcv,mdfDICcv=out$mdfDICcv)
    file<-paste(curMod,curObs,curFreq,curReps,sep='_')
    save(curOut,file=paste('/afs/crc.nd.edu/user/j/jzwart/ModelDataFusion/Results/',file,'.Rda',sep = ''))
  }
}

runs=runs$unique

mclapply(runs,sim,mc.cores=48)
