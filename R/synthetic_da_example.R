# creating fake data for MDF
# JAZ; 2016-11-22

load('Data/EnKF_LongData_20170223.RData')

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

set.seed(42)

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

reps=2
freq=7
obs=1

# draws from priors to create ensemble
parGuess <- c(0.007,0.1,1) #r20; units: day-1
min<-c(0.001,0.3,-3)
max<-c(0.02,0.5,3)

rPDF<-abs(rnorm(n=nEn,mean = parGuess[1],sd = (max[1]-min[1])/5)) # max-min / 5 is rule of thumb sd; forcing positive for negative draws
fracInPDF<-abs(rnorm(n=nEn,mean=parGuess[2],sd=(max[2]-min[2])/5))
covar_inflat_PDF <- abs(rnorm(n=nEn, mean=parGuess[3], sd = (max[3]-min[3])/5))
# fracInPDF<-rbeta(n=nEn,shape1 = 2,shape2 = .5)

# setting up initial parameter values for all ensemble members
rVec<-matrix(rPDF) # each row is an ensemble member
fracInVec<-matrix(fracInPDF)
covar_inflat_vec<-matrix(covar_inflat_PDF)
hist(fracInVec)

# initial B transition matrix for each ensemble
B<-array(t(c(rep(NA,nEn),rep(NA,nEn),rep(NA,nEn),rep(NA,nEn))),dim=c(2,2,nStep,nEn)) # array of transition matrix B[a,b,c,d];
# where [a,b] is transition matrix DIC & DOC, c=timeStep, and d=ensemble member

# initial parameters of B at timestep 1
for(i in 1:nEn){
  B[,,1,i]<-matrix(c(1-data2$QoutInt[1]/data2$epiVol[1]-data2$kCO2[1]/data2$thermo.depth[1]+
                       (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                     rVec[i],0,1-rVec[i]-data2$QoutInt[1]/data2$epiVol[1]+
                       (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1]),
                   nrow=2,byrow=T)
}
# observations come from model generated data
docCV<-0.133 # observation error to add to true data
dicCV<-0.137 # observation error to add to true data
reps<-reps # number of replicates for sampling from true data distribution
dicVec<-true[5,1,,1]
docVec<-true[8,1,,1]
dicSD_samp<-rep(NA,nStep)
docSD_samp<-rep(NA,nStep)
docSD<-(data2$doc/data2$epiVol)*0.1325858 # DOC concentration sd in mol C
docSD<-ifelse(is.na(docSD),docSD,docSD[data2$datetime=='2014-07-30']) # making SD the same for all observations; not based on concentration 2016-11-22
dicSD<-(data2$dic/data2$epiVol)*0.1367684 # DIC concentration sd in mol C
dicSD<-ifelse(is.na(dicSD),dicSD,dicSD[data2$datetime=='2014-07-30'])
for(t in 1:nStep){
  curdic<-rnorm(n=reps,mean=true[5,1,t,1]/data2$epiVol[t],sd = dicSD[data2$datetime=='2014-07-30']*obs)
  curdoc<-rnorm(n=reps,mean=true[8,1,t,1]/data2$epiVol[t],sd = docSD[data2$datetime=='2014-07-30']*obs)
  dicVec[t]<-mean(curdic)*data2$epiVol[t]
  docVec[t]<-mean(curdoc)*data2$epiVol[t]
  dicSD_samp[t]<-sd(curdic)*data2$epiVol[t]
  docSD_samp[t]<-sd(curdoc)*data2$epiVol[t]
}
# how many observations to have
freq<-freq # sampling frequency in days
dicVec[2:nStep]<-ifelse(data2$timeStep[2:nStep]%%freq==0,dicVec[2:nStep],NA)
docVec[2:nStep]<-ifelse(data2$timeStep[2:nStep]%%freq==0,docVec[2:nStep],NA)
if(is.na(dicVec[1])){
  dicVec[1]<-dicVec[min(which(!is.na(dicVec)))]
  docVec[1]<-docVec[min(which(!is.na(docVec)))]
}


docSD<-(data2$doc/data2$epiVol)*0.1325858 # DOC concentration sd in mol C
docSD<-ifelse(is.na(docVec),NA,docSD[data2$datetime=='2014-07-30']) # making SD the same for all observations; not based on concentration 2016-11-22
docSD<-docSD*obs # DOC concentration sd in mol C; modifying for sensativity analysis
dicSD<-(data2$dic/data2$epiVol)*0.1367684 # DIC concentration sd in mol C
dicSD<-ifelse(is.na(dicVec),NA,dicSD[data2$datetime=='2014-07-30']) # making SD the same for all observations; not based on concentration 2016-11-22
dicSD<-dicSD*obs
areaSD<-4000 # constant SD in m
depthSD<-0.25 # constant SD in m

H<-array(0,dim=c(2,2,nStep))
# propogation of error for multiplication
docPoolSD<-docVec*sqrt((docSD/(docVec/data2$epiVol))^2+(areaSD/data2$A0)^2+(depthSD/data2$thermo.depth)^2)
dicPoolSD<-dicVec*sqrt((dicSD/(dicVec/data2$epiVol))^2+(areaSD/data2$A0)^2+(depthSD/data2$thermo.depth)^2)
H[1,1,]<-dicPoolSD^2 #variance of DIC
H[2,2,]<-docPoolSD^2 #variance of DOC
H[1,1,]<-ifelse(is.na(H[1,1,]),mean(H[1,1,],na.rm=T),H[1,1,]) # taking care of na's in DIC; setting to mean of variance if NA
H[2,2,]<-ifelse(is.na(H[2,2,]),mean(H[2,2,],na.rm=T),H[2,2,]) # taking care of na's in DOC; setting to mean of variance if NA

y=array(rbind(dicVec,docVec),dim=c(2,1,nStep))
y=array(rep(y,nEn),dim=c(2,1,nStep,nEn)) # array of observations y[a,b,c,d]; where a=dic/doc, b=column, c=timeStep, and d=ensemble member

X<-array(NA,dim =c(2,1,nStep,nEn)) # model estimate matrices X[a,b,c,d]; where a=dic/doc, b=column, c=timestep, and d=ensemble member

#initializing estimate of state, X, with first observation
# drawing from normal distribution for inital pool sizes
X[1,1,1,]<-rnorm(n=nEn,y[1,1,1,],sd=dicPoolSD[1])
X[2,1,1,]<-rnorm(n=nEn,y[2,1,1,],sd=docPoolSD[1])

# operator matrix saying 1 if there is observation data available, 0 otherwise
h<-array(0,dim=c(2,5,nStep))
for(i in 1:nStep){
  h[1,4,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic
  h[2,5,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) #doc
}


P <- array(0,dim=c(3,3,nStep))

#Define matrix C, parameters of covariates [2x6]
C<-array(NA,dim=c(2,6,nStep,nEn)) # array of transition matrix C[a,b,c,d];
#intializing first time step
for(i in 1:nEn){
  C[,,1,i]<-matrix(c(data2$kCO2[1],1,-data2$ma_gpp[1],0,0,data2$hypo_dicInt[1],0,0,0,data2$ma_gpp[1],1,data2$hypo_docInt[1]),nrow=2,byrow=T)
}

ut<-array(NA,dim=c(6,1,nStep,nEn)) # array of transition matrix C[a,b,c,d];
#intializing first time step
for(i in 1:nEn){
  ut[,,1,i]<-matrix(c(data2$DICeq[1]*data2$epiVol[1]/data2$thermo.depth[1],
                      data2$dicIn[1]-data2$streamWaterdisch[1]*data2$Inlet_dic[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])),(1-GPPrespired),
                      exudeTotal,data2$docIn[1]-data2$streamDOCdisch[1]/12/1000*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])),
                      (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainEpi[1]),ncol=1)
}

pars<-array(rep(NA,nEn),dim=c(3,1,nStep,nEn)) # parameters: r20
pars[1,1,1,]<-rVec
pars[2,1,1,]<-fracInVec
pars[3,1,1,]<-covar_inflat_vec

# set up a list for all matrices
z=list(B=B,y=y,X=X,C=C,ut=ut,pars=pars)

# i is ensemble member; t is timestep
i=1
t=2

# set up Y vector for which we concatonate parameters, states, and observed data
Y<-array(NA,c(5,1,nStep,nEn))
Y[1,1,1,]<-rVec # r20 parameter
Y[2,1,1,]<-fracInVec # fraction that goes into epi
Y[3,1,1,]<-covar_inflat_vec
Y[4,1,1,]<-z$X[1,1,1,] # DIC state
Y[5,1,1,]<-z$X[2,1,1,] # DOC state

#Iterate through time
for(t in 2:nStep){
  for(i in 1:nEn){
    # Forecasting; need to update parameters,
    z$pars[1:3,1,t,i]<-Y[1:3,1,t-1,i] # r20
    z$X[1:2,1,t-1,i]<-Y[4:5,1,t-1,i] # updating state variables from Y

    #Predictions
    z$X[,,t,i]<-z$B[,,t-1,i]%*%z$X[,,t-1,i] + z$C[,,t-1,i]%*%z$ut[,,t-1,i] # forecasting state variable predictions
    z$B[,,t,i]<-matrix(c(1-data2$QoutInt[t]/data2$epiVol[t]-data2$kCO2[t]/data2$thermo.depth[t]+
                           (data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[2,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],
                         z$pars[1,1,t,i],0,1-z$pars[1,1,t,i]-data2$QoutInt[t]/data2$epiVol[t]+
                           (data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[2,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t]),
                       nrow=2,byrow=T)
    z$y[,,t,i]<-z$y[,,t,i] # observation of states stay the same
    # parameters are of previous timestep for matrix C, parameters of covariates [2x6]
    z$C[,,t,i]<-matrix(c(data2$kCO2[t],1,-data2$ma_gpp[t],0,0,data2$hypo_dicInt[t],0,0,0,data2$ma_gpp[t],1,data2$hypo_docInt[t]),nrow=2,byrow=T)
    z$ut[,,t,i]<-matrix(c(data2$DICeq[t]*data2$epiVol[t]/data2$thermo.depth[t],
                          data2$dicIn[t]-data2$streamWaterdisch[t]*data2$Inlet_dic[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[2,1,t,i])),(1-GPPrespired),
                          exudeTotal,data2$docIn[t]-data2$streamDOCdisch[t]/12/1000*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[2,1,t,i])),
                          (data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[2,1,t,i])))*data2$entrainEpi[t]),ncol=1)

    # forecast Y vector
    Y[1:3,1,t,i]<-Y[1:3,1,t-1,i] #r20,same as previous timestep
    Y[4:5,1,t,i]<-z$X[1:2,1,t,i] #forecasted states

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

    covar_inflat <- mean(Y[3,1,t,])

    K<-((1/(nEn-1))*covar_inflat*delta_Y%*%t(delta_Y)%*%t(h[,,t]))%*%
      qr.solve(((1/(nEn-1))*covar_inflat*h[,,t]%*%delta_Y%*%t(delta_Y)%*%t(h[,,t])+H[,,t]))   # Kalman gain 7x2 matrix; 7x3 if iota is included

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

DOCout<-apply(Y[5,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
DOCpredictSD<-apply(Y[5,1,,]/data2$epiVol*12,MARGIN = 1,FUN=sd)
DOCout<-cbind(DOCout,z$y[2,1,(spinUpLength+1):nStep,1]/data2$epiVol*12)
DOCout<-cbind(DOCout,(true[8,1,,1]/data2$epiVol*12))

DICout<-apply(Y[4,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
DICpredictSD<-apply(Y[4,1,,]/data2$epiVol*12,MARGIN = 1,FUN=sd)
DICout<-cbind(DICout,z$y[1,1,(spinUpLength+1):nStep,1]/data2$epiVol*12)
DICout<-cbind(DICout,(true[5,1,,1]/data2$epiVol*12))

sqrt(mean((DOCout[,1]-DOCout[,2])^2,na.rm=T)) #RMSE mod-obs
sqrt(mean((DOCout[,1]-DOCout[,3])^2,na.rm=T)) # RMSE mod-true
sqrt(mean((DOCout[,2]-DOCout[,3])^2,na.rm=T)) #RMSE obs-true

sqrt(mean((DICout[,1]-DICout[,2])^2,na.rm=T)) #RMSE mod-obs
sqrt(mean((DICout[,1]-DICout[,3])^2,na.rm=T)) # RMSE mod-true
sqrt(mean((DICout[,2]-DICout[,3])^2,na.rm=T)) #RMSE obs-true

mean(DOCout[,1]-DOCout[,2],na.rm=T)^2 #Bias mod-obs
mean(DOCout[,1]-DOCout[,3],na.rm=T)^2 #Bias mod-true
mean(DOCout[,2]-DOCout[,3],na.rm=T)^2 #Bias obs-true

mean(DICout[,1]-DICout[,2],na.rm=T)^2 #Bias mod-obs
mean(DICout[,1]-DICout[,3],na.rm=T)^2 #Bias mod-true
mean(DICout[,2]-DICout[,3],na.rm=T)^2 #Bias obs-true

summary(lm(DOCout[,1]~DOCout[,2]))$r.squared #Precision mod-obs
summary(lm(DOCout[,1]~DOCout[,3]))$r.squared #Precision mod-true
summary(lm(DOCout[,2]~DOCout[,3]))$r.squared #Precision obs-true

summary(lm(DICout[,1]~DICout[,2]))$r.squared #Precision mod-obs
summary(lm(DICout[,1]~DICout[,3]))$r.squared #Precision mod-true
summary(lm(DICout[,2]~DICout[,3]))$r.squared #Precision obs-true


# DOC pCO2 throught time Fig. 3
##############
png('Figures/Fig3_all.png',
res=300, width=14, height=14, units = 'in')
# windows()
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels
leg=0.05

cex=3
cex.lab=3
cex.axis=3
cex.leg=1.5
lwd=3
ylim=c(0.5,3)
xlim=range(as.POSIXct(data2$datetime[spinUpLength+1:nStep]))
par(mar=c(5,7,4,2),mfrow=c(2,2))
DOCout<-apply(Y[5,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
ylim=range(c(Y[5,1,,]/data2$epiVol*12,z$y[2,1,,1]/data2$epiVol*12+docPoolSD/data2$epiVol*12,z$y[2,1,,1]/data2$epiVol*12-docPoolSD/data2$epiVol*12),na.rm=T)
plot(DOCout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),type='l',ylim=ylim,lwd=lwd,cex.axis=cex.axis,
     ylab=expression(DOC~(mg~C~L^-1)),xlab='',cex.lab=cex.lab)
for(i in 1:nEn){
  lines(Y[5,1,(spinUpLength+1):nStep,i]/data2$epiVol*12~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
        col='gray',ylab='',lwd=lwd)
}
lines(DOCout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),ylab='',lwd=3,col='gray30')
lines(true[8,1,,1]/data2$epiVol*12~as.POSIXct(data2$datetime),ylim = ylim,col='black',lwd=3,lty=6)
par(new=T)
plot(z$y[2,1,(spinUpLength+1):nStep,1]/data2$epiVol*12~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),cex.axis=cex.axis,ylab=expression(DOC~(mg~C~L^-1)),
     ylim=ylim,col='black',pch=16,cex=cex,xlab='',cex.lab=cex.lab)
arrows(as.POSIXct(data2$datetime),(z$y[2,1,,1]/data2$epiVol*12)-docPoolSD/data2$epiVol*12,as.POSIXct(data2$datetime),
       (z$y[2,1,,1]/data2$epiVol*12)+docPoolSD/data2$epiVol*12,code=3,length=0.1,angle=90,col='black',lwd=lwd)
legend(x = as.POSIXct('2014-06-10'), y = 22.5, legend=c("Estimated State","Ensemble Mean",'True State','Observed State'),
       col=c('gray','gray30','black','black'),pt.bg=c('gray','gray30','black','black'), cex=cex.leg,
       ncol=1,lwd=c(4,4,4,0),bty='n',lty=c(1,1,3,0),pt.cex = c(0,0,0,2),pch=c(0,0,0,16))
text(x=xlim[2]-leg*(xlim[2]-xlim[1]),y = ylim[1]+leg*(ylim[2]-ylim[1]),labels = 'A', cex = cex.lab)

DICout<-apply(Y[4,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
ylim=range(c(Y[4,1,,]/data2$epiVol*12,z$y[1,1,,1]/data2$epiVol*12+dicPoolSD/data2$epiVol*12,z$y[1,1,,1]/data2$epiVol*12-dicPoolSD/data2$epiVol*12),na.rm=T)
plot(DICout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),type='l',ylim=ylim,lwd=lwd,cex.axis=cex.axis,
     ylab=expression(CO[2]~(mg~C~L^-1)),xlab='',cex.lab=cex.lab)
for(i in 1:nEn){
  lines(Y[4,1,(spinUpLength+1):nStep,i]/data2$epiVol*12~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
        col='gray',ylab='',lwd=lwd)
}
lines(DICout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),ylab='',lwd=lwd,col='gray30')
lines(true[5,1,,1]/data2$epiVol*12~as.POSIXct(data2$datetime),ylim = ylim,col='black',lwd=lwd,lty=6)
par(new=T)
plot(z$y[1,1,(spinUpLength+1):nStep,1]/data2$epiVol*12~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),cex.axis=cex.axis,
     ylim=ylim,col='black',pch=16,cex=cex,ylab=expression(CO[2]~(mg~C~L^-1)),xlab='',cex.lab=cex.lab)
arrows(as.POSIXct(data2$datetime),(z$y[1,1,,1]/data2$epiVol*12)-dicPoolSD/data2$epiVol*12,as.POSIXct(data2$datetime),
       (z$y[1,1,,1]/data2$epiVol*12)+dicPoolSD/data2$epiVol*12,code=3,length=0.1,angle=90,col='black',lwd=lwd)
legend("topleft", legend=c("Estimated State","Ensemble Mean",'True State','Observed State'),
       col=c('gray','gray30','black','black'),pt.bg=c('gray','gray30','black','black'), cex=cex.leg,
       ncol=1,lwd=c(4,4,4,0),bty='n',lty=c(1,1,3,0),pt.cex = c(0,0,0,2),pch=c(0,0,0,16))
text(x=xlim[2]-leg*(xlim[2]-xlim[1]),y = ylim[1]+leg*(ylim[2]-ylim[1]),labels = 'B', cex = cex.lab)


#d20
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

dSD0<-sd(Y[1,1,1,]) #initial sd in d20
dSD<-apply(Y[1,1,,],MARGIN = 1,FUN=sd)
dConverged<-dSD/dSD0 # fraction of initial sd
ylim=c(0.5,3)
par(mar=c(5,7,4,2))
rOut<-apply(Y[1,1,,],MARGIN = 1,FUN=mean)
ylim=range(rFunc(Y[1,1,(spinUpLength+1):nStep,],data2$wtr[(spinUpLength+1):nStep]))
ylim=c(0.00,0.015)
plot(rOut[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
     type='l',ylim=ylim,ylab=expression(d[20]~(day^-1)),lwd=lwd,xlab='',cex.lab=cex,cex.axis=cex.axis)
for(i in 1:nEn){
  lines(Y[1,1,(spinUpLength+1):nStep,i]~
          as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),col='gray',ylab='',lwd=lwd,xlab='')
}
lines(rOut[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),ylab='',lwd=3,col='gray30',
      xlab='')

r20L<-true[2,1,,1]
r20R<-true[1,1,,1]
DOCoutL<-true[7,1,,1]
DOCoutR<-true[6,1,,1]
DOCoutT<-true[8,1,,1]
fracLout<-DOCoutL/DOCoutT
fracRout<-DOCoutR/DOCoutT
r20All<-r20L*fracLout+r20R*fracRout
lines(r20All~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),col='black',lty=6,
      type='l',ylim=ylim,ylab=expression(d[20]~(day^-1)),lwd=lwd,xlab='',cex.lab=cex,cex.axis=cex)
legend("top", legend=c("Estimated d20","Ensemble Mean",'True d20'),
       col=c('gray','gray30','black'),pt.bg=c('gray','gray30','black'), cex=cex.leg,
       ncol=1,lwd=c(4,4,4),bty='n',lty=c(1,1,3))
text(x=xlim[2]-leg*(xlim[2]-xlim[1]),y = ylim[1]+leg*(ylim[2]-ylim[1]),labels = 'C', cex = cex.lab)

# CV
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

par(mar=c(5,7,4,2))
DOCout<-apply(Y[5,1,,],MARGIN = 1,FUN=mean)
DOCpredictSD<-apply(Y[5,1,,],MARGIN = 1,FUN=sd)
DOCcv<-DOCpredictSD/DOCout
DOCcv<-cbind(DOCcv,docPoolSD/z$y[2,1,,1])
DICout<-apply(Y[4,1,,],MARGIN = 1,FUN=mean)
DICpredictSD<-apply(Y[4,1,,],MARGIN = 1,FUN=sd)
DICcv<-DICpredictSD/DICout
DICcv<-cbind(DICcv,dicPoolSD/dicSDadjust/z$y[1,1,,1])
ylim=range(DICcv,DOCcv,na.rm=T)

plot(DOCcv[,1]~as.POSIXct(data2$datetime),ylim=ylim,cex=2,cex.axis=1.5,pch=16,type='l',lwd=lwd,
     ylab=expression(CV~DOC~and~CO[2]~(mol~C)),cex.lab=cex,cex.axis = cex.axis, xlab='')
lines(DICcv[,1]~as.POSIXct(data2$datetime),ylim=ylim,cex=2,cex.axis=1.5,pch=16,type='l',lwd=lwd,col='gray60')
points(DOCcv[,2]~as.POSIXct(data2$datetime),lwd=lwd,lty=2,pch=16,cex=cex)
points(DICcv[,2]~as.POSIXct(data2$datetime),lwd=lwd,lty=2,pch=16,cex=cex,col='grey60')
legend(x = as.POSIXct('2014-08-01'), y = 0.15, legend=c("DA DOC CV",'DA CO2 CV','Obs DOC CV','Obs CO2 CV'),
       col=c('black','gray60','black','gray60'),pt.bg=c('black','gray60','black','gray60'), cex=cex.leg,
       ncol=1,lwd=c(4,4,0,0),bty='n',lty=c(1,1,0,0),pch = c(0,0,16,16),pt.cex=c(0,0,2,2))
text(x=xlim[2]-.01*(xlim[2]-xlim[1]),y = ylim[1]+.01*(ylim[2]-ylim[1]),labels = 'D', cex = cex.lab)

dev.off()
###########


