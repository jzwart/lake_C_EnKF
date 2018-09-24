# # # # setting up Long Data for EnKF; EnKF runs at end of code
# # # # JAZ; 2016-03-28
# validating with every other observation

# this model explicitly accounts for changes in pH and carbonate equilibrium in model

###################
#Run Kalman filter
library(dplyr)
load('Data/EL_20170408.RData') # in this Rdata version, dic = dic and co2 = co2; in the previous one EnKF_LongData_20170223.RData, dic = co2 and there was no total dic data
ph = readRDS('Data/EL_ph.rds') %>%
  mutate(date= as.character(date))
source('R/EnKF_dic_co2.R')
library(AquaEnv)

splitFunc<-function(epiDens,streamDens,fracIn){ # not based on density difference between, just fraction split
  fracInEpi=exp(-fracIn*(streamDens-epiDens))
  return(fracIn)
}

frac_co2 <- function(dic, ph, temp, vol_epi, sal = 0){ # function that takes temperature, DIC pool, and pH to calculate CO2 fraction
  temp = temp # deg C
  DIC_pool = dic # mol C
  Vepi = vol_epi # m^3
  pH = ph
  sal = sal # salinity

  aquaOut = AquaEnv::aquaenv(S = sal, # salinity
                             t = temp, # water temp
                             d = 0, # depth
                             lat = 46, # latitude
                             SumCO2 = DIC_pool/Vepi/rLakeAnalyzer::water.density(temp,sal),
                             pH = pH)

  CO2=aquaOut$CO2[1]*water.density(temp,sal)*Vepi
  fracCO2=CO2/DIC_pool
  return(fracCO2)
}

set.seed(42)

data2 <- data
data2 <- left_join(data2, ph, by = c('datetime' = 'date'))
data2$hypo_dicInt<-data2$hypo_dicInt*0.25 # entrained CO2 is much less than where we measure CO2

# cut down to only when high frequency discharge out was known. occurs on 2014-06-02
data2<-data2[min(which(!is.na(data2$highFreqWaterHeight))):max(which(!is.na(data2$highFreqWaterHeight))),]
data2<-data2[!is.na(data2$ma_gpp),]
data2<-data2[data2$datetime<as.POSIXct('2015-01-01'),] # only keeping 2014
data2<-data2[min(which(!is.na(data2$doc))):nrow(data2),]

# adding in trash pump discharge out starting on Aug. 22, 2015 until the end of time series in 2015
# data2$QoutInt[data2$datetime>=as.POSIXct('2015-08-22')]<-data2$QoutInt[data2$datetime>=as.POSIXct('2015-08-22')]+400

fracLabile0<-0.01 # fraction of initial DOC pool that is labile
fracLabile<-0.08 # fraction loaded t-DOC that is labile; estimate from Berggren et al. 2010 Eco Letts & ISME

# need some initial conditions for observations
data2$dic[1]<-data2$dic[min(which(!is.na(data2$dic)))]
data2$doc[1]<-data2$doc[min(which(!is.na(data2$doc)))]
data2$co2[1]<-data2$co2[min(which(!is.na(data2$co2)))]

## *********************** EnKF ***************************## Gao et al. 2011 is a useful reference
nEn<-100 # number of ensembles
nStep<-length(data2$datetime)

# draws from priors to create ensemble
parGuess <- c(0.004,0.3,0.30,0.1,1) #r20; units: day-1; fraction labile of loaded DOC; fraction inlet that goes into epi; turnover rate parameters set constant frac labile estimated
min<-c(0.004,0.31,0.005,0.3,-3)
max<-c(0.004,0.31,0.5,0.5,3)

rPDF<-abs(rnorm(n=nEn,mean = parGuess[1],sd = (max[1]-min[1])/5)) # max-min / 5 is rule of thumb sd;
rPDF_fast<-abs(rnorm(n=nEn,mean = parGuess[2],sd = (max[2]-min[2])/5))
fracPDF<-abs(rnorm(n=nEn,mean=parGuess[3],sd=(max[3]-min[3])/5))
fracInPDF<-abs(rnorm(n=nEn,mean=parGuess[4],sd=(max[4]-min[4])/5))
covar_inflat_PDF <- rnorm(n=nEn, mean=parGuess[5], sd=(max[5]-min[5])/5)

# setting up initial parameter values for all ensemble members
rVec<-matrix(rPDF) # each row is an ensemble member
rVec_fast<-matrix(rPDF_fast)
fracVec<-matrix(fracPDF)
fracInVec<-matrix(fracInPDF)
covar_inflat_vec <- matrix(covar_inflat_PDF)

data2$entrainHypo<-as.numeric(data2$entrainVol<0)
data2$entrainEpi<-as.numeric(data2$entrainVol>0)

# setting up state values for all ensemble members (all the same initial pool size )
# updated on 2016-07-16 to draw from a normal distribution using first observation and SD
dicVec<-data2$dic
docVec<-data2$doc
co2Vec<-data2$co2
if(is.na(dicVec[1])){
  dicVec[1]<-dicVec[min(which(!is.na(dicVec)))]
  docVec[1]<-docVec[min(which(!is.na(docVec)))]
  co2Vec[1]<-co2Vec[min(which(!is.na(co2Vec)))]
}

# estimate of observation error covariance matrix
# error in DOC / DIC pool esitmate comes from concentration estimates and volume of epilimnion
# DOC / DIC concentration error from replication day in WL on 2014-07-30. CV for DOC is 0.1325858, DIC is 0.1367684
# Epi: area sd = 4000m2; depth sd = 0.25 m;
# iota sd from metab estimates
docSDadjust<-1
dicSDadjust<-1
co2SDadjust<-1
docSD<-(data2$doc/data2$epiVol)*0.1325858 # DOC concentration sd in mol C
docSD<-ifelse(is.na(docSD),docSD,docSD[data2$datetime=='2014-07-30']) # making SD the same for all observations; not based on concentration 2016-11-22
docSD<-docSD*docSDadjust # DOC concentration sd in mol C; modifying for sensativity analysis
dicSD<-(data2$dic/data2$epiVol)*0.1367684 # DIC concentration sd in mol C
dicSD<-ifelse(is.na(dicSD),dicSD,dicSD[data2$datetime=='2014-07-30']) # making SD the same for all observations; not based on concentration 2016-11-22
dicSD<-dicSD*dicSDadjust
co2SD<-(data2$co2/data2$epiVol)*0.1367684 # DIC concentration sd in mol C
co2SD<-ifelse(is.na(co2SD),co2SD,co2SD[data2$datetime=='2014-07-30']) # making SD the same for all observations; not based on concentration 2016-11-22
co2SD<-co2SD*dicSDadjust
areaSD<-4000 # constant SD in m
depthSD<-0.25 # constant SD in m

H<-array(0,dim=c(3,3,nStep))
# propogation of error for multiplication
docPoolSD<-data2$doc*sqrt((docSD/(data2$doc/data2$epiVol))^2+(areaSD/data2$A0)^2+(depthSD/data2$thermo.depth)^2)
dicPoolSD<-data2$dic*sqrt((dicSD/(data2$dic/data2$epiVol))^2+(areaSD/data2$A0)^2+(depthSD/data2$thermo.depth)^2)
co2PoolSD<-data2$co2*sqrt((co2SD/(data2$co2/data2$epiVol))^2+(areaSD/data2$A0)^2+(depthSD/data2$thermo.depth)^2)
H[1,1,]<-dicPoolSD^2 #variance of DIC
H[2,2,]<-co2PoolSD^2 # variance of CO2
H[3,3,]<-docPoolSD^2 #variance of DOC

H[1,1,]<-ifelse(is.na(H[1,1,]),mean(H[1,1,],na.rm=T),H[1,1,]) # taking care of na's in DIC; setting to mean of variance if NA
H[2,2,]<-ifelse(is.na(H[2,2,]),mean(H[2,2,],na.rm=T),H[2,2,]) # taking care of na's in CO2; setting to mean of variance if NA
H[3,3,]<-ifelse(is.na(H[3,3,]),mean(H[3,3,],na.rm=T),H[3,3,]) # taking care of na's in DOC; setting to mean of variance if NA

y=array(rbind(dicVec,co2Vec,docVec),dim=c(3,1,nStep))
y=array(rep(y,nEn),dim=c(3,1,nStep,nEn)) # array of observations y[a,b,c,d]; where a=dic/doc/co2, b=column, c=timeStep, and d=ensemble member

X<-array(NA,dim =c(5,1,nStep,nEn)) # model estimate matrices X[a,b,c,d]; where a=dic/co2/doc_r/doc_l/doc, b=column, c=timestep, and d=ensemble member

#initializing estimate of state, X, with first observation
# drawing from normal distribution for inital pool sizes
X[1,1,1,]<-rnorm(n=nEn,y[1,1,1,],sd=dicPoolSD[1])
X[2,1,1,]<-rnorm(n=nEn,y[2,1,1,],sd=co2PoolSD[1])
X[3,1,1,]<-rnorm(n=nEn,y[3,1,1,]*(1-fracLabile0),sd=docPoolSD[1]*(1-fracLabile0)) # recalcitrant pool is 90% of initial DOC pool
X[4,1,1,]<-rnorm(n=nEn,y[3,1,1,]*fracLabile0,sd=docPoolSD[1]*fracLabile0) # labile pool is 10% of initial DOC pool
X[5,1,1,]<-X[3,1,1,]+X[4,1,1,] # recalcitrant pool + labile pool

# initial B transition matrix for each ensemble
B<-array(NA,dim=c(5,5,nStep,nEn)) # array of transition matrix B[a,b,c,d];
# where [a,b] is transition matrix DIC & DOCr & DOCl, c=timeStep, and d=ensemble member
# initial parameters of B at timestep 1
for(i in 1:nEn){
  B[,,1,i]<-matrix(c(1-data2$QoutInt[1]/data2$epiVol[1]-
                       frac_co2(dic = X[1,,1,i], ph = data2$pH_sensor_corr[1], temp = data2$wtr[1], vol_epi = data2$epiVol[1])*
                       data2$kCO2[1]/data2$thermo.depth[1]+
                       (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*
                       data2$entrainHypo[1]/data2$epiVol[1],0,rFunc(rVec[i],data2$wtr[1]),rFunc(rVec_fast[i],data2$wtr[1]),0, # end of row 1
                     1*frac_co2(dic = X[1,,1,i], ph = data2$pH_sensor_corr[1], temp = data2$wtr[1], vol_epi = data2$epiVol[1]),
                       0,0,0,0, # end of row 2
                     0,0,1-rFunc(rVec[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1]+
                       (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*
                       data2$entrainHypo[1]/data2$epiVol[1],0,0, # end of row 3
                     0,0,0,1-rFunc(rVec_fast[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1]+
                       (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*
                       data2$entrainHypo[1]/data2$epiVol[1],0, # end of row 4
                     0,0,1-rFunc(rVec[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1]+
                       (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*
                       data2$entrainHypo[1]/data2$epiVol[1],1-rFunc(rVec_fast[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1]+
                       (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*
                       data2$entrainHypo[1]/data2$epiVol[1],0), # end of row 5
                   nrow=5,byrow=T)
}

# operator matrix saying 1 if there is observation data available, 0 otherwise
h<-array(0,dim=c(3,10,nStep))
for(i in 1:nStep){
  h[1,6,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic
  h[2,7,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) # co2
  h[3,10,i]<-ifelse(!is.na(y[3,1,i,1]),1,0) #doc total (we only have data on total DOC pool)
}

P <- array(0,dim=c(5,5,nStep))
S <- array(0,dim=c(5,5,nStep))
PS <- array(0, dim=c(10,10,nStep))

#Define matrix C, parameters of covariates [2x6]
C<-array(NA,dim=c(5,7,nStep,nEn)) # array of transition matrix C[a,b,c,d];
#intializing first time step
for(i in 1:nEn){
  C[,,1,i]<-matrix(c(data2$kCO2[1],1,-data2$ma_gpp[1],0,0,0,data2$hypo_dicInt[1], # end of row 1
                     0,0,0,0,0,0,0, # end of row 2
                     0,0,0,data2$ma_gpp[1],0,(1-fracVec[i]),data2$hypo_docInt[1]*(1-fracLabile0), # end of row 3
                     0,0,0,0,data2$ma_gpp[1],fracVec[i],data2$hypo_docInt[1]*fracLabile0, # end of row 4
                     0,0,0,data2$ma_gpp[1],data2$ma_gpp[1],1,data2$hypo_docInt[1]), # end of row 5
                   nrow = 5, byrow =T)
}

ut<-array(NA,dim=c(7,1,nStep,nEn)) # array of transition matrix C[a,b,c,d];
#intializing first time step
for(i in 1:nEn){
  ut[,,1,i]<-matrix(c(data2$DICeq[1]*data2$epiVol[1]/data2$thermo.depth[1],
                      data2$dicIn[1]-data2$streamDICdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])),
                      (1-GPPrespired),
                      exude,
                      exudeLabile,
                      data2$docIn[1]-data2$streamDOCdisch[1]/12/1000*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])),
                      (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainEpi[1]),
                    ncol=1)
}

pars<-array(rep(NA,nEn),dim=c(5,1,nStep,nEn)) # parameters: r20
pars[1,1,1,]<-rVec
pars[2,1,1,]<-rVec_fast
pars[3,1,1,]<-fracVec
pars[4,1,1,]<-fracInVec
pars[5,1,1,]<-covar_inflat_vec

# set up a list for all matrices
z=list(B=B,y=y,X=X,C=C,ut=ut,pars=pars)

# i is ensemble member; t is timestep
i=1
t=2

# set up Y vector for which we concatonate parameters, states, and observed data
Y<-array(NA,c(10,1,nStep,nEn))
Y[1,1,1,]<-rVec # r20 parameter
Y[2,1,1,]<-rVec_fast # r20 labile parameter
Y[3,1,1,]<-fracVec # fraction labile of loaded DOC
Y[4,1,1,]<-fracInVec # fraction of Inlet stream that goes into epi
Y[5,1,1,]<-pars[5,1,1,]# inflation factor
Y[6,1,1,]<-z$X[1,1,1,] # DIC state
Y[7,1,1,]<-z$X[2,1,1,] # CO2 state
Y[8,1,1,]<-z$X[3,1,1,] # DOC recalcitrant state
Y[9,1,1,]<-z$X[4,1,1,] # DOC labile state
Y[10,1,1,]<-z$X[5,1,1,] # DOC total state

assim_number <- 1 # counter for checking if observation should be assimilated or not (assimilating everyother obs, validating on left out obs); assimilate odd numbers
assim_obs <- rep(0,nStep)

out = EnKF_2pools_dic_co2(Y, z, i, t) # run EnKF dic / co2 version
Y = out$Y
assim_obs = out$assim_obs

# # mean RMSE for states
print(sqrt(mean((apply(Y[6,1,assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[assim_obs==0]*12-z$y[1,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12)^2,na.rm=T)))
print(sqrt(mean((apply(Y[7,1,assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[assim_obs==0]*12-z$y[2,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12)^2,na.rm=T)))
print(sqrt(mean((apply(Y[10,1,assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[assim_obs==0]*12-z$y[3,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12)^2,na.rm=T)))

cor(apply(Y[6,1,!is.na(z$y[1,1,,1])&assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[!is.na(z$y[1,1,,1])&assim_obs==0]*12,z$y[1,1,!is.na(z$y[1,1,,1])&assim_obs==0,1]/data2$epiVol[!is.na(z$y[1,1,,1])&assim_obs==0]*12)^2
cor(apply(Y[7,1,!is.na(z$y[2,1,,1])&assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[!is.na(z$y[2,1,,1])&assim_obs==0]*12,z$y[2,1,!is.na(z$y[2,1,,1])&assim_obs==0,1]/data2$epiVol[!is.na(z$y[2,1,,1])&assim_obs==0]*12)^2
cor(apply(Y[10,1,!is.na(z$y[3,1,,1])&assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[!is.na(z$y[3,1,,1])&assim_obs==0]*12,z$y[3,1,!is.na(z$y[3,1,,1])&assim_obs==0,1]/data2$epiVol[!is.na(z$y[3,1,,1])&assim_obs==0]*12)^2

#bias
mean(apply(Y[7,1,assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[assim_obs==0]*12-z$y[2,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12,na.rm = T)^2
mean(apply(Y[10,1,assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[assim_obs==0]*12-z$y[3,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12,na.rm=T)^2

# # mean RMSE for states; mol C
print(sqrt(mean((apply(Y[7,1,assim_obs==0,],MARGIN = 1,FUN=mean)-z$y[2,1,assim_obs==0,1])^2,na.rm=T)))
print(sqrt(mean((apply(Y[10,1,assim_obs==0,],MARGIN = 1,FUN=mean)-z$y[3,1,assim_obs==0,1])^2,na.rm=T)))

# r2; mol C
cor(apply(Y[7,1,!is.na(z$y[2,1,,1])&assim_obs==0,],MARGIN = 1,FUN=mean),z$y[2,1,!is.na(z$y[2,1,,1])&assim_obs==0,1])^2
cor(apply(Y[10,1,!is.na(z$y[3,1,,1])&assim_obs==0,],MARGIN = 1,FUN=mean),z$y[3,1,!is.na(z$y[3,1,,1])&assim_obs==0,1])^2

# AIC CO2
fit = apply(Y[7,1,assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[assim_obs==0]*12
obs = z$y[2,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12
AIC(logLik(lm(obs~fit)))

# AIC DOC
fit = apply(Y[10,1,assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[assim_obs==0]*12
obs = z$y[3,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12
AIC(logLik(lm(obs~fit)))

# checking if variance in ensembles captures the obs at 95% confidence
# #
# doc_05 = (apply(Y[9,1,assim_obs==0,],MARGIN = 1,FUN=mean) - 1.96 * apply(Y[9,1,assim_obs==0,],MARGIN = 1,FUN=sd)) / data2$epiVol[assim_obs==0]*12
# doc_95 = (apply(Y[9,1,assim_obs==0,],MARGIN = 1,FUN=mean) + 1.96 * apply(Y[9,1,assim_obs==0,],MARGIN = 1,FUN=sd)) / data2$epiVol[assim_obs==0]*12
# doc_obs = z$y[2,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12
# perc_doc_obs_in_ci = sum(doc_obs < doc_95 & doc_obs > doc_05, na.rm = T) / sum(!is.na(doc_obs)) * 100
#
# co2_05 = (apply(Y[6,1,assim_obs==0,],MARGIN = 1,FUN=mean) - 1.96 * apply(Y[6,1,assim_obs==0,],MARGIN = 1,FUN=sd)) / data2$epiVol[assim_obs==0]*12
# co2_95 = (apply(Y[6,1,assim_obs==0,],MARGIN = 1,FUN=mean) + 1.96 * apply(Y[6,1,assim_obs==0,],MARGIN = 1,FUN=sd)) / data2$epiVol[assim_obs==0]*12
# co2_obs = z$y[1,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12
# perc_co2_obs_in_ci = sum(co2_obs < co2_95 & co2_obs > co2_05, na.rm = T) / sum(!is.na(co2_obs)) * 100
#
# # all obs
# doc_05 = (apply(Y[9,1,,],MARGIN = 1,FUN=mean) - 1.96 * apply(Y[9,1,,],MARGIN = 1,FUN=sd)) / data2$epiVol[]*12
# doc_95 = (apply(Y[9,1,,],MARGIN = 1,FUN=mean) + 1.96 * apply(Y[9,1,,],MARGIN = 1,FUN=sd)) / data2$epiVol[]*12
# doc_obs = z$y[2,1,,1]/data2$epiVol[]*12
# perc_doc_obs_in_ci = sum(doc_obs < doc_95 & doc_obs > doc_05, na.rm = T) / sum(!is.na(doc_obs)) * 100
#
# co2_05 = (apply(Y[6,1,,],MARGIN = 1,FUN=mean) - 1.96 * apply(Y[6,1,,],MARGIN = 1,FUN=sd)) / data2$epiVol[]*12
# co2_95 = (apply(Y[6,1,,],MARGIN = 1,FUN=mean) + 1.96 * apply(Y[6,1,,],MARGIN = 1,FUN=sd)) / data2$epiVol[]*12
# co2_obs = z$y[1,1,,1]/data2$epiVol[]*12
# perc_co2_obs_in_ci = sum(co2_obs < co2_95 & co2_obs > co2_05, na.rm = T) / sum(!is.na(co2_obs)) * 100
#
# # 95% CI of obs
# doc_obs_05 = z$y[2,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12- 1.96 * docPoolSD[assim_obs==0]/data2$epiVol[assim_obs==0]*12
# doc_obs_95 = z$y[2,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12+ 1.96 * docPoolSD[assim_obs==0]/data2$epiVol[assim_obs==0]*12
# perc_doc_ci_in_ci = sum(doc_obs_05 < doc_95 & doc_obs_95 > doc_05, na.rm = T) / sum(!is.na(doc_obs_05)) * 100
#
# co2_obs_05 = z$y[1,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12- 1.96 * dicPoolSD[assim_obs==0]/data2$epiVol[assim_obs==0]*12
# co2_obs_95 = z$y[1,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12+ 1.96 * dicPoolSD[assim_obs==0]/data2$epiVol[assim_obs==0]*12
# perc_co2_ci_in_ci = sum(co2_obs_05 < co2_95 & co2_obs_95 > co2_05, na.rm = T) / sum(!is.na(co2_obs_05)) * 100


#this is kinda cool
windows()
for(t in 1:length(data2$datetime)){
  if(t == 1){
    plot(density(Y[10,1,t,]), ylim = c(0,1e-3))
  }else{
    lines(density(Y[10,1,t,]))
  }
}

for(t in 1:length(data2$datetime)){
  if(t == 1){
    plot(density(Y[5,1,t,]), ylim = c(0,2))
    abline(v = 1)
  }else{
    lines(density(Y[5,1,t,]))
  }
}

# # Figure 2 plotting concentration
###############
png('Figures/Fig2_all_rev2.png',
    res=300, width=14, height=21, units = 'in')
# windows()
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels
leg = 0.05 # distance from corner for panel label

cex=4
cex.lab=3
cex.axis=2
lwd=3
ylim=c(0.5,3)
xlim=range(as.POSIXct(data2$datetime))
par(mar=c(5,7,4,2),mfrow=c(3,2))
DOCout<-apply(Y[10,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
ylim=range(c(Y[10,1,,]/data2$epiVol*12,z$y[3,1,,1]/data2$epiVol*12+docPoolSD/data2$epiVol*12,z$y[3,1,,1]/data2$epiVol*12-docPoolSD/data2$epiVol*12),na.rm=T)
plot(DOCout~as.POSIXct(data2$datetime),type='l',ylim=ylim,lwd=lwd,cex.axis=cex.axis,
     ylab=expression(DOC~(mg~C~L^-1)),xlab='',cex.lab=cex.lab)
for(i in 1:nEn){
  lines(Y[10,1,,i]/data2$epiVol*12~as.POSIXct(data2$datetime),
        col='gray',ylab='',lwd=lwd)
}
lines(DOCout~as.POSIXct(data2$datetime),ylab='',lwd=3,col='gray30')
par(new=T)
plot(z$y[3,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12~as.POSIXct(data2$datetime[assim_obs==0]),cex=cex,
     ylim=ylim,col='black',pch=21,ylab=expression(DOC~(mg~C~L^-1)),xlab='',cex.lab=cex.lab,cex.axis=cex.axis,xlim=xlim)
arrows(as.POSIXct(data2$datetime[assim_obs==0]),z$y[3,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12-docPoolSD[assim_obs==0]/data2$epiVol[assim_obs==0]*12,
       as.POSIXct(data2$datetime[assim_obs==0]),z$y[3,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12+docPoolSD[assim_obs==0]/data2$epiVol[assim_obs==0]*12,
       code=3,length=0.1,angle=90,col='black',lwd=3)
par(new=T)
plot(z$y[3,1,assim_obs==1,1]/data2$epiVol[assim_obs==1]*12~as.POSIXct(data2$datetime[assim_obs==1]),cex=cex,
     ylim=ylim,col='black',pch=16,ylab=expression(DOC~(mg~C~L^-1)),xlab='',cex.lab=cex.lab,cex.axis=cex.axis,xlim=xlim)
arrows(as.POSIXct(data2$datetime[assim_obs==1]),z$y[3,1,assim_obs==1,1]/data2$epiVol[assim_obs==1]*12-docPoolSD[assim_obs==1]/data2$epiVol[assim_obs==1]*12,
       as.POSIXct(data2$datetime[assim_obs==1]),z$y[3,1,assim_obs==1,1]/data2$epiVol[assim_obs==1]*12+docPoolSD[assim_obs==1]/data2$epiVol[assim_obs==1]*12,
       code=3,length=0.1,angle=90,col='black',lwd=3)
legend("topleft", legend=c("Estimated State","Ensemble Mean",'Observed State'),
       col=c('gray','gray30','black'),pt.bg=c('gray','gray30','black'),cex = cex.axis,
       ncol=1,lwd=c(4,4,0),bty='n',lty=c(1,1,0),pt.cex = c(0,0,2),pch = c(0,0,16))
text(x=xlim[2]-leg*(xlim[2]-xlim[1]),y = ylim[1]+leg*(ylim[2]-ylim[1]),labels = 'A', cex = cex.lab)

CO2out<-apply(Y[7,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
ylim=range(c(Y[7,1,,]/data2$epiVol*12,z$y[2,1,,1]/data2$epiVol*12+co2PoolSD/data2$epiVol*12,z$y[2,1,,1]/data2$epiVol*12-co2PoolSD/data2$epiVol*12),na.rm=T)
plot(CO2out~as.POSIXct(data2$datetime),type='l',ylim=ylim,lwd=lwd,cex.axis=cex.axis,
     ylab=expression(CO[2]~(mg~C~L^-1)),xlab='',cex.lab=cex.lab)
for(i in 1:nEn){
  lines(Y[7,1,,i]/data2$epiVol*12~as.POSIXct(data2$datetime),
        col='gray',ylab='',lwd=lwd)
}
lines(CO2out~as.POSIXct(data2$datetime),ylab='',lwd=3,col='gray30')
par(new=T)
plot(z$y[2,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12~as.POSIXct(data2$datetime[assim_obs==0]),cex=cex,
     ylim=ylim,col='black',pch=21,ylab = '',xlab='',cex.lab=cex.lab,cex.axis=cex.axis,xlim=xlim)
arrows(as.POSIXct(data2$datetime[assim_obs==0]),z$y[2,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12-co2PoolSD[assim_obs==0]/data2$epiVol[assim_obs==0]*12,
       as.POSIXct(data2$datetime[assim_obs==0]),z$y[2,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12+co2PoolSD[assim_obs==0]/data2$epiVol[assim_obs==0]*12,
       code=3,length=0.1,angle=90,col='black',lwd=3)
par(new=T)
plot(z$y[2,1,assim_obs==1,1]/data2$epiVol[assim_obs==1]*12~as.POSIXct(data2$datetime[assim_obs==1]),cex=cex,
     ylim=ylim,col='black',pch=16,ylab='',xlab='',cex.lab=cex.lab,cex.axis=cex.axis,xlim=xlim)
arrows(as.POSIXct(data2$datetime[assim_obs==1]),z$y[2,1,assim_obs==1,1]/data2$epiVol[assim_obs==1]*12-co2PoolSD[assim_obs==1]/data2$epiVol[assim_obs==1]*12,
       as.POSIXct(data2$datetime[assim_obs==1]),z$y[2,1,assim_obs==1,1]/data2$epiVol[assim_obs==1]*12+co2PoolSD[assim_obs==1]/data2$epiVol[assim_obs==1]*12,
       code=3,length=0.1,angle=90,col='black',lwd=3)
legend("topleft", legend=c("Estimated State","Ensemble Mean",'Observed State'),
       col=c('gray','gray30','black'),pt.bg=c('gray','gray30','black'),cex = cex.axis,
       ncol=1,lwd=c(4,4,0),bty='n',lty=c(1,1,0),pt.cex = c(0,0,2),pch = c(0,0,16))
text(x=xlim[2]-leg*(xlim[2]-xlim[1]),y = ylim[1]+leg*(ylim[2]-ylim[1]),labels = 'B', cex = cex.lab)

# turnover rate
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

r20L<-apply(rFunc(Y[2,1,,],data2$epiTemp),MARGIN = 1,FUN=mean)
r20R<-apply(rFunc(Y[1,1,,],data2$epiTemp),MARGIN = 1,FUN=mean)
DOCoutL<-apply(Y[9,1,,],MARGIN = 1,FUN=mean)
DOCoutR<-apply(Y[8,1,,],MARGIN = 1,FUN=mean)
DOCoutT<-apply(Y[10,1,,],MARGIN = 1,FUN=mean)
fracLout<-DOCoutL/DOCoutT
fracRout<-DOCoutR/DOCoutT
r20All<-r20L*fracLout+r20R*fracRout
ylim=range(rFunc(Y[2,1,,],data2$epiTemp)*(Y[9,1,,]/Y[10,1,,])+rFunc(Y[1,1,,],data2$epiTemp)*(Y[8,1,,]/Y[10,1,,]))
ylim[1]=0.001
par(mar=c(5,7,4,2))
plot(r20All~as.POSIXct(data2$datetime),
     type='l',ylim=ylim,ylab=expression(d~(day^-1)),lwd=lwd,xlab='',cex.lab=cex,cex.axis=cex.axis)
for(i in 1:nEn){
  curFracLout<-Y[9,1,,i]/Y[10,1,,i]
  curFracRout<-Y[8,1,,i]/Y[10,1,,i]
  curr20All<-rFunc(Y[2,1,,i],data2$epiTemp)*curFracLout+rFunc(Y[1,1,,i],data2$epiTemp)*curFracRout
  lines(curr20All,col='gray',ylab='')
  lines(curr20All~
          as.POSIXct(data2$datetime),col='gray',ylab='',lwd=lwd,xlab='')
}
lines(r20All~as.POSIXct(data2$datetime),ylab='',lwd=lwd,col='gray30',
      xlab='')
legend("topright", legend=c("d Estimate","Ensemble Mean"),
       col=c('gray','gray30'),pt.bg=c('gray','gray30'),cex = cex.axis,
       ncol=1,lwd=c(4,4),bty='n',lty=c(1,1),pt.cex = c(0,0),pch = c(0,0))
text(x=xlim[2]-leg*(xlim[2]-xlim[1]),y = ylim[1]+leg*(ylim[2]-ylim[1]),labels = 'C', cex = cex.lab)

# d 20
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

r20L<-apply(Y[2,1,,],MARGIN = 1,FUN=mean)
r20R<-apply(Y[1,1,,],MARGIN = 1,FUN=mean)
DOCoutL<-apply(Y[9,1,,],MARGIN = 1,FUN=mean)
DOCoutR<-apply(Y[8,1,,],MARGIN = 1,FUN=mean)
DOCoutT<-apply(Y[10,1,,],MARGIN = 1,FUN=mean)
fracLout<-DOCoutL/DOCoutT
fracRout<-DOCoutR/DOCoutT
r20All<-r20L*fracLout+r20R*fracRout
ylim=range(rFunc(Y[2,1,,],data2$epiTemp)*(Y[9,1,,]/Y[10,1,,])+rFunc(Y[1,1,,],data2$epiTemp)*(Y[8,1,,]/Y[10,1,,]))
ylim[1]=0.001
par(mar=c(5,7,4,7))
rOut<-apply(Y[1,1,,],MARGIN = 1,FUN=mean)
plot(r20All~as.POSIXct(data2$datetime),
     type='l',ylim=ylim,ylab=expression(d[20]~(day^-1)),lwd=lwd,xlab='',cex.lab=cex,cex.axis=cex.axis)
for(i in 1:nEn){
  curFracLout<-Y[9,1,,i]/Y[10,1,,i]
  curFracRout<-Y[8,1,,i]/Y[10,1,,i]
  curr20All<-Y[2,1,,i]*curFracLout+Y[1,1,,i]*curFracRout
  lines(curr20All,col='gray',ylab='')
  lines(curr20All~
          as.POSIXct(data2$datetime),col='gray',ylab='',lwd=lwd,xlab='')
}
lines(r20All~as.POSIXct(data2$datetime),ylab='',lwd=lwd,col='gray30',
      xlab='')
par(new=T)
plot(data2$docIn~as.POSIXct(data2$datetime),type='l',lwd=4,
     col='black',yaxt='n',xaxt='n',ylab='',xlab='',lty=6)
axis(4,cex.axis=cex.axis)
mtext(expression(loadDOC~(mol~C~day^-1)),4,line = 5,cex=2)
legend("topleft", legend=c("d20 Estimate",'d20 Ensemble Mean','DOC loading'),
       col=c('gray','gray30','black'),pt.bg=c('gray','gray30','black'),cex = cex.axis,
       ncol=1,lwd=c(4,4,4),bty='n',lty=c(1,1,3))
text(x=xlim[2]-.01*(xlim[2]-xlim[1]),y = ylim[1]+leg*(ylim[2]-ylim[1]),labels = 'D', cex = cex.lab)

l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

dSD0<-sd(Y[3,1,1,]) #initial sd in d20
dSD<-apply(Y[3,1,,],MARGIN = 1,FUN=sd)
dConverged<-dSD/dSD0 # fraction of initial sd
par(mar=c(5,6,4,2))
rOut<-apply(Y[3,1,,],MARGIN = 1,FUN=mean)
ylim=range(Y[3,1,,])
plot(rOut~as.POSIXct(data2$datetime),
     type='l',ylim=ylim,ylab=expression(fracFast),lwd=lwd,xlab='',cex.lab=cex,cex.axis=cex.axis)
for(i in 1:nEn){
  lines(Y[3,1,,i]~as.POSIXct(data2$datetime),col='gray',ylab='',lwd=lwd,xlab='')
}
lines(rOut~as.POSIXct(data2$datetime),ylab='',lwd=lwd,col='gray30',xlab='')
legend("topright", legend=c("fracFast Estimate","Ensemble Mean"),
       col=c('gray','gray30'),pt.bg=c('gray','gray30'),cex = cex.axis,
       ncol=1,lwd=c(4,4),bty='n',lty=c(1,1),pt.cex = c(0,0),pch = c(0,0))
text(x=xlim[2]-leg*(xlim[2]-xlim[1]),y = ylim[1]+leg*(ylim[2]-ylim[1]),labels = 'E', cex = cex.lab)


l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

par(mar=c(5,7,4,2))
DOCout<-apply(Y[10,1,,],MARGIN = 1,FUN=mean)
DOCpredictSD<-apply(Y[10,1,,],MARGIN = 1,FUN=sd)
DOCcv<-DOCpredictSD/DOCout
DOCcv<-cbind(DOCcv,docPoolSD/z$y[3,1,,1])
CO2out<-apply(Y[7,1,,],MARGIN = 1,FUN=mean)
CO2predictSD<-apply(Y[7,1,,],MARGIN = 1,FUN=sd)
CO2cv<-CO2predictSD/CO2out
CO2cv<-cbind(CO2cv,co2PoolSD/co2SDadjust/z$y[2,1,,1])
ylim=range(c(na.omit(DOCcv),na.omit(CO2cv)),na.rm=T)

plot(DOCcv[,1]~as.POSIXct(data2$datetime),ylim=ylim,cex=cex,cex.axis=cex.axis,pch=16,type='l',lwd=lwd,
     ylab=expression(CV~DOC~and~CO[2]~(mol~C)),cex.lab=cex,xlab='')
lines(CO2cv[,1]~as.POSIXct(data2$datetime),ylim=ylim,cex=cex,cex.axis=cex.axis,pch=16,type='l',lwd=lwd,col='gray60')
points(DOCcv[,2]~as.POSIXct(data2$datetime),pch=16,cex=cex)
points(CO2cv[,2]~as.POSIXct(data2$datetime),pch=16,cex=cex,col='gray60')
legend("topright", legend=c("DA DOC CV",'DA CO2 CV','Obs DOC CV','Obs CO2 CV'),
       col=c('black','gray60','black','gray60'),pt.bg=c('black','gray60','black','gray60'),cex = cex.axis,
       ncol=1,lwd=c(4,4,0,0),bty='n',lty=c(1,1,0,0),pt.cex = c(0,0,2,2),pch=c(0,0,16,16))
text(x=xlim[2]-.01*(xlim[2]-xlim[1]),y = ylim[1]+.01*(ylim[2]-ylim[1]),labels = 'F', cex = cex.lab)

dev.off()




DOCout<-apply(Y[9,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
DICout<-apply(Y[6,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)

# how much CO2 missed due to linear interpolation? between obs time point 6 & 7
linInt=z
# linInt$y[1,1,36:43,1]<-approx(36:43,z$y[1,1,36:43,1]/data2$epiVol[36:43]*12,xout = 36:43)$y
linInt$y[1,1,,1]<-approx(1:nStep,z$y[1,1,,1]/data2$epiVol[]*12,xout = 1:nStep)$y
linInt$y[2,1,,1]<-approx(1:nStep,z$y[2,1,,1]/data2$epiVol[]*12,xout = 1:nStep)$y

range(DICout[]/linInt$y[1,1,,1],na.rm = T) # all diff
range(DICout[!as.numeric(h[1,6,])]/linInt$y[1,1,!as.numeric(h[1,6,]),1],na.rm = T) # diff during unmonitored
# max decrease
max((DICout[!as.numeric(h[1,6,])]-linInt$y[1,1,!as.numeric(h[1,6,]),1])/DICout[!as.numeric(h[1,6,])]*100,na.rm = T)
#max increase
max((linInt$y[1,1,!as.numeric(h[1,6,]),1]-DICout[!as.numeric(h[1,6,])])/DICout[!as.numeric(h[1,6,])]*100,na.rm = T)
range(DICout[as.logical(h[1,6,])]/linInt$y[1,1,as.logical(h[1,6,]),1],na.rm = T) # diff during monitored

range(DOCout[]/linInt$y[2,1,,1],na.rm = T) # all diff
range(DOCout[!as.numeric(h[2,9,])]/linInt$y[2,1,!as.numeric(h[2,9,]),1],na.rm = T) # diff during unmonitored
# max decrease
max((DOCout[!as.numeric(h[2,9,])]-linInt$y[2,1,!as.numeric(h[2,9,]),1])/DOCout[!as.numeric(h[2,9,])]*100,na.rm = T)
#max increase
max((linInt$y[2,1,!as.numeric(h[2,9,]),1]-DOCout[!as.numeric(h[2,9,])])/DOCout[!as.numeric(h[2,9,])]*100,na.rm = T)
range(DOCout[as.logical(h[2,9,])]/linInt$y[2,1,as.logical(h[2,9,]),1],na.rm = T) # diff during monitored

#what if we only had bi-weekly obs
linInt$y[1,1,,1]<-approx(1:nStep,z$y[1,1,,1]/data2$epiVol[]*12,xout = 1:nStep)$y
linInt$y[2,1,,1]<-approx(1:nStep,z$y[2,1,,1]/data2$epiVol[]*12,xout = 1:nStep)$y

range(DICout[]/linInt$y[1,1,,1],na.rm = T)
range(DICout[!as.numeric(h[1,6,])]/linInt$y[1,1,!as.numeric(h[1,6,]),1],na.rm = T)
range(DICout[as.logical(h[1,6,])]/linInt$y[1,1,as.logical(h[1,6,]),1],na.rm = T)

range(DOCout[]/linInt$y[2,1,,1],na.rm = T)
range(DOCout[!as.numeric(h[2,9,])]/linInt$y[2,1,!as.numeric(h[2,9,]),1],na.rm = T)
range(DOCout[as.logical(h[2,9,])]/linInt$y[2,1,as.logical(h[2,9,]),1],na.rm = T)

