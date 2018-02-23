# setting up Long Data for EnKF; EnKF runs at end of code
# JAZ; 2016-03-28

#Run Kalman filter
rm(list=ls())
# load('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/R Data/EnKF_LongData_20161218.RData')
load('/Users/jzwart/LakeCarbonEnKF/Data/EnKF_LongData_20170223.RData')

splitFunc<-function(epiDens,streamDens,fracIn){ # function that tells how much load goes into epi 
  fracInEpi=exp(-fracIn*(streamDens-epiDens))
  return(fracInEpi)
}

splitFunc<-function(epiDens,streamDens,fracIn){ # function that tells how much load goes into epi 
  fracInEpi=exp(-fracIn*(streamDens-epiDens))
  return(fracIn)
}

# spin up, just repeating the first X days at the begining of the timeseries; autocorrelation effects?? 
spinUpLength<-0
spinUp<-data[1:spinUpLength,]
data2<-rbind(spinUp,data)
# should we cut down to only when high frequency discharge out was known? occurs on 2014-06-02
data2<-data2[min(which(!is.na(data2$highFreqWaterHeight))):max(which(!is.na(data2$highFreqWaterHeight))),]
data2<-data2[!is.na(data2$ma_gpp),]
data2<-data2[data2$datetime<as.POSIXct('2015-01-01'),] # only keeping 2014
data2<-data2[min(which(!is.na(data2$doc))):nrow(data2),]

data2$hypo_dicInt<-data2$hypo_dicInt*0.25 # entrained CO2 is much less than where we measure CO2 
# data2<-data2[data2$datetime<as.POSIXct('2015-08-25'),] # missing big entrainment event 

# adding in trash pump discharge out starting on Aug. 22, 2015 until the end of time series in 2015 
data2$QoutInt[data2$datetime>=as.POSIXct('2015-08-22')]<-data2$QoutInt[data2$datetime>=as.POSIXct('2015-08-22')]+400

# need some initial conditions for observations (draw from normal distribution?)
data2$dic[1]<-data2$dic[min(which(!is.na(data2$dic)))]
data2$doc[1]<-data2$doc[min(which(!is.na(data2$doc)))]

## *********************** EnKF ***************************## Gao et al. 2011 is a useful reference 
nEn<-100 # number of ensembles 
nStep<-length(data2$datetime)

# draws from priors to create ensemble 
parGuess <- c(0.007,0.1,1) #r20; units: day-1
min<-c(0.001,0.3,-30)
max<-c(0.02,0.5,30)

rPDF<-abs(rnorm(n=nEn,mean = parGuess[1],sd = (max[1]-min[1])/5)) # max-min / 5 is rule of thumb sd; forcing positive for negative draws 
fracInPDF<-abs(rnorm(n=nEn,mean=parGuess[2],sd=(max[2]-min[2])/5))
covar_inflat_PDF <- rnorm(n=nEn, mean=parGuess[3], sd=(max[3]-min[3])/5)
# fracInPDF<-rbeta(n = nEn,shape1 = 4,shape2 = 1)

# setting up initial parameter values for all ensemble members 
rVec<-matrix(rPDF) # each row is an ensemble member 
fracInVec<-matrix(fracInPDF)
covar_inflat_vec<-matrix(covar_inflat_PDF)
hist(fracInVec)
  
  data2$entrainHypo<-as.numeric(data2$entrainVol<0)
  data2$entrainEpi<-as.numeric(data2$entrainVol>0)
  
  # initial B transition matrix for each ensemble  
  B<-array(NA,dim=c(2,2,nStep,nEn)) # array of transition matrix B[a,b,c,d]; 
  # where [a,b] is transition matrix DIC & DOCr & DOCl, c=timeStep, and d=ensemble member 
  # initial parameters of B at timestep 1 
  for(i in 1:nEn){
    B[,,1,i]<-matrix(c(1-data2$QoutInt[1]/data2$epiVol[1]-data2$kCO2[1]/data2$thermo.depth[1]+
                         (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                       rFunc(rVec[i],data2$wtr[1]),0,1-rFunc(rVec[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1]+
                         (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1]),
                     nrow=2,byrow=T)
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
  docSDadjust<-1
  dicSDadjust<-1
  docSD<-(data2$doc/data2$epiVol)*0.1325858 # DOC concentration sd in mol C
  docSD<-ifelse(is.na(docSD),docSD,docSD[data2$datetime=='2014-07-30']) # making SD the same for all observations; not based on concentration 2016-11-22
  # docSD<-ifelse(is.na(docSD),docSD,0.11)
  docSD<-docSD*docSDadjust # DOC concentration sd in mol C; modifying for sensativity analysis
  dicSD<-(data2$dic/data2$epiVol)*0.1367684 # DIC concentration sd in mol C 
  dicSD<-ifelse(is.na(dicSD),dicSD,dicSD[data2$datetime=='2014-07-30']) # making SD the same for all observations; not based on concentration 2016-11-22
  # dicSD<-ifelse(is.na(dicSD),dicSD,0.005)
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
  
  X<-array(NA,dim =c(2,1,nStep,nEn)) # model estimate matrices X[a,b,c,d]; where a=dic/doc_r/doc_l, b=column, c=timestep, and d=ensemble member 
  
  #initializing estimate of state, X, with first observation 
  # drawing from normal distribution for inital pool sizes 
  X[1,1,1,]<-rnorm(n=nEn,y[1,1,1,],sd=dicPoolSD[1])
  X[2,1,1,]<-rnorm(n=nEn,y[2,1,1,],sd=docPoolSD[1]) 
  
  # operator matrix saying 1 if there is observation data available, 0 otherwise 
  h<-array(0,dim=c(2,5,nStep))
  for(i in 1:nStep){
    h[1,4,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic 
    h[2,5,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) #doc total (we only have data on total DOC pool)
  }
  
  P <- array(0,dim=c(3,3,nStep))
  S <- array(0,dim=c(2,2,nStep))
  
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
  Y[2,1,1,]<-fracInVec # frac inlet into epi
  Y[3,1,1,]<-covar_inflat_vec # covariance inlation 
  Y[4,1,1,]<-z$X[1,1,1,] # DIC state  
  Y[5,1,1,]<-z$X[2,1,1,] # DOC total state  
  
  assim_number <- 0 # counter for checking if observation should be assimilated or not (assimilating everyother obs, validating on left out obs); assimilate odd numbers 
  assim_obs <- rep(0,nStep)
  
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
               rFunc(z$pars[1,1,t,i],data2$wtr[t]),0,1-rFunc(z$pars[1,1,t,i],data2$wtr[t])-data2$QoutInt[t]/data2$epiVol[t]+
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
    
    #begin data assimilation if there are any observations and obs is an assimilation obs 
    if(any(!is.na(z$y[,,t,i]))==TRUE & (assim_number %% 2) == 0){ # update vector as long as there is one observation of state 
      assim_number <- assim_number+1
      assim_obs[t] <- 1
      
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
      
      # 
      statesMean<-matrix(apply(Y[(length(z$pars[,1,1,1])+1):length(Y[,1,1,1]),,t,],MARGIN = 1,FUN=mean),nrow=nrow(Y)-nrow(z$pars))
      states_matrix<-array(NA,dim=c(nrow(Y)-nrow(z$pars),nrow(Y)-nrow(z$pars),nEn))
      for(i in 1:nEn){
        delta_states<-(Y[(nrow(z$pars)+1):nrow(Y),1,t,i]-statesMean)
        states_matrix[,,i]<-delta_states%*%t(delta_states)
      }
      statesTemp<-matrix(0,ncol=nrow(Y)-nrow(z$pars),nrow=nrow(Y)-nrow(z$pars))
      for(i in 1:nEn){
        statesTemp<-statesTemp+states_matrix[,,i]
      }
      S[,,t]<-(1/(nEn-1))*statesTemp
      
    }else if(any(!is.na(z$y[,,t,i]))==TRUE & (assim_number %% 2) != 0){
      assim_number <- assim_number+1
    }
    
  } # End iteration
  
  Y14=Y # saving Y parameter for plotting later 
  data14=data2
  z14=z
  assim_obs14=assim_obs
  # take parameter output from 2014 fit and validate with 2015 data
  nStep14=nStep
  
  load('/Users/jzwart/LakeCarbonEnKF/Data/EnKF_LongData_20170223.RData')
  
  # spin up, just repeating the first X days at the begining of the timeseries; autocorrelation effects?? 
  spinUpLength<-0
  spinUp<-data[1:spinUpLength,]
  data2<-rbind(spinUp,data)
  # should we cut down to only when high frequency discharge out was known? occurs on 2014-06-02
  data2<-data2[min(which(!is.na(data2$highFreqWaterHeight))):max(which(!is.na(data2$highFreqWaterHeight))),]
  data2<-data2[!is.na(data2$ma_gpp),]
  data2<-data2[data2$datetime>as.POSIXct('2015-01-01'),] # only keeping 2014
  data2<-data2[min(which(!is.na(data2$doc))):nrow(data2),]
  
  data2$hypo_dicInt<-data2$hypo_dicInt*0.25 # entrained CO2 is much less than where we measure CO2 
  # data2<-data2[data2$datetime<as.POSIXct('2015-08-25'),] # missing big entrainment event 
  
  # adding in trash pump discharge out starting on Aug. 22, 2015 until the end of time series in 2015 
  data2$QoutInt[data2$datetime>=as.POSIXct('2015-08-22')]<-data2$QoutInt[data2$datetime>=as.POSIXct('2015-08-22')]+400
  
  # need some initial conditions for observations (draw from normal distribution?)
  data2$dic[1]<-data2$dic[min(which(!is.na(data2$dic)))]
  data2$doc[1]<-data2$doc[min(which(!is.na(data2$doc)))]
  
  ## *********************** EnKF ***************************## Gao et al. 2011 is a useful reference 
  nEn<-100 # number of ensembles 
  nStep<-length(data2$datetime)
  
  # initial ensembles from 2014 fits 
  rVec<-matrix(Y[1,1,length(Y[1,1,,1]),]) # each row is an ensemble member 
  fracInVec<-matrix(Y[2,1,length(Y[1,1,,1]),])
  covar_inflat_vec <- matrix(Y[3,1,length(Y[1,1,,1]),])
  
  data2$entrainHypo<-as.numeric(data2$entrainVol<0)
  data2$entrainEpi<-as.numeric(data2$entrainVol>0)
  
  # initial B transition matrix for each ensemble  
  B<-array(NA,dim=c(2,2,nStep,nEn)) # array of transition matrix B[a,b,c,d]; 
  # where [a,b] is transition matrix DIC & DOCr & DOCl, c=timeStep, and d=ensemble member 
  # initial parameters of B at timestep 1 
  for(i in 1:nEn){
    B[,,1,i]<-matrix(c(1-data2$QoutInt[1]/data2$epiVol[1]-data2$kCO2[1]/data2$thermo.depth[1]+
                         (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                       rFunc(rVec[i],data2$wtr[1]),0,1-rFunc(rVec[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1]+
                         (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1]),
                     nrow=2,byrow=T)
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
  docSDadjust<-1
  dicSDadjust<-1
  docSD<-(data2$doc/data2$epiVol)*0.1325858 # DOC concentration sd in mol C
  # docSD<-ifelse(is.na(docSD),docSD,docSD[data2$datetime=='2014-07-30']) # making SD the same for all observations; not based on concentration 2016-11-22
  docSD<-ifelse(is.na(docSD),docSD,0.11)
  docSD<-docSD*docSDadjust # DOC concentration sd in mol C; modifying for sensativity analysis
  dicSD<-(data2$dic/data2$epiVol)*0.1367684 # DIC concentration sd in mol C 
  # dicSD<-ifelse(is.na(dicSD),dicSD,dicSD[data2$datetime=='2014-07-30']) # making SD the same for all observations; not based on concentration 2016-11-22
  dicSD<-ifelse(is.na(dicSD),dicSD,0.005)
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
  
  X<-array(NA,dim =c(2,1,nStep,nEn)) # model estimate matrices X[a,b,c,d]; where a=dic/doc_r/doc_l, b=column, c=timestep, and d=ensemble member 
  
  #initializing estimate of state, X, with first observation 
  # drawing from normal distribution for inital pool sizes 
  X[1,1,1,]<-rnorm(n=nEn,y[1,1,1,],sd=dicPoolSD[1])
  X[2,1,1,]<-rnorm(n=nEn,y[2,1,1,],sd=docPoolSD[1]) 
  
  # operator matrix saying 1 if there is observation data available, 0 otherwise 
  h<-array(0,dim=c(2,5,nStep))
  for(i in 1:nStep){
    h[1,4,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic 
    h[2,5,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) #doc total (we only have data on total DOC pool)
  }
  
  P <- array(0,dim=c(3,3,nStep))
  S <- array(0,dim=c(2,2,nStep))
  
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
  Y[2,1,1,]<-fracInVec # frac inlet into epi
  Y[3,1,1,]<-covar_inflat_vec # covariance inlation 
  Y[4,1,1,]<-z$X[1,1,1,] # DIC state  
  Y[5,1,1,]<-z$X[2,1,1,] # DOC total state  
  
  assim_number <- 0 # counter for checking if observation should be assimilated or not (assimilating everyother obs, validating on left out obs); assimilate odd numbers 
  assim_obs <- rep(0,nStep)
  
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
                           rFunc(z$pars[1,1,t,i],data2$wtr[t]),0,1-rFunc(z$pars[1,1,t,i],data2$wtr[t])-data2$QoutInt[t]/data2$epiVol[t]+
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
    
    # begin data assimilation if there are any observations
    #begin data assimilation if there are any observations and obs is an assimilation obs 
    if(any(!is.na(z$y[,,t,i]))==TRUE & (assim_number %% 2) == 0){ # update vector as long as there is one observation of state 
      assim_number <- assim_number+1
      assim_obs[t] <- 1
      
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

      #
      statesMean<-matrix(apply(Y[(length(z$pars[,1,1,1])+1):length(Y[,1,1,1]),,t,],MARGIN = 1,FUN=mean),nrow=nrow(Y)-nrow(z$pars))
      states_matrix<-array(NA,dim=c(nrow(Y)-nrow(z$pars),nrow(Y)-nrow(z$pars),nEn))
      for(i in 1:nEn){
        delta_states<-(Y[(nrow(z$pars)+1):nrow(Y),1,t,i]-statesMean)
        states_matrix[,,i]<-delta_states%*%t(delta_states)
      }
      statesTemp<-matrix(0,ncol=nrow(Y)-nrow(z$pars),nrow=nrow(Y)-nrow(z$pars))
      for(i in 1:nEn){
        statesTemp<-statesTemp+states_matrix[,,i]
      }
      S[,,t]<-(1/(nEn-1))*statesTemp

    }else if(any(!is.na(z$y[,,t,i]))==TRUE & (assim_number %% 2) != 0){
      assim_number <- assim_number+1
    }

  } # End iteration
  
  nStep15=nStep
  
  # 
  data2 <- rbind(data14,data2)
  zz <- z
  YY <- Y
  nStep <- length(Y14[1,1,,1])+nStep
  Y <- array(NA,c(5,1,nStep,nEn))
  for(i in 1:length(Y[,1,1,1])){
    for(j in 1:nEn){
      Y[i,1,,j] <- c(Y14[i,1,,j],YY[i,1,,j])
    }
  }
  
  z$y <- array(NA,c(2,1,nStep,nEn))
  z$X <- array(NA,c(2,1,nStep,nEn))
  for(i in 1:length(z$y[,1,1,1])){
    for(j in 1:nEn){
      z$y[i,1,,j] <- c(z14$y[i,1,,j],zz$y[i,1,,j])
      z$X[i,1,,j] <- c(z14$X[i,1,,j],zz$X[i,1,,j])
    }
  }
  
  docSD=rep(0.11,nStep)
  dicSD=rep(0.005,nStep)
  docPoolSD<-data2$doc*sqrt((docSD/(data2$doc/data2$epiVol))^2+(areaSD/data2$A0)^2+(depthSD/data2$thermo.depth)^2)
  dicPoolSD<-data2$dic*sqrt((dicSD/(data2$dic/data2$epiVol))^2+(areaSD/data2$A0)^2+(depthSD/data2$thermo.depth)^2)
  
  assim_obs <- c(assim_obs,assim_obs14)
  
  # 
  # # plotting ***************************************************
  windows()
  DOCout<-apply(Y[5,1,,],MARGIN = 1,FUN=mean)
  ylim=range(Y[5,1,,]/data2$epiVol*12)
  plot(DOCout/data2$epiVol*12,type='l',ylim=ylim,ylab='DOC mg/L')
  for(i in 1:nEn){
    lines(Y[5,1,,i]/data2$epiVol*12,col='gray',ylab='')
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
  rOut<-apply(splitFunc(data2$epiDens,data2$streamDens,Y[2,1,,]),MARGIN = 1,FUN=mean)
  ylim=range(splitFunc(data2$epiDens,data2$streamDens,Y[2,1,,]))
  plot(rOut,type='l',ylim=ylim,ylab='r day-1')
  for(i in 1:nEn){
    lines(splitFunc(data2$epiDens,data2$streamDens,Y[2,1,,i]),col='gray',ylab='')
  }
  lines(rOut,ylab='')
  
  windows()
  rOut<-apply(splitFunc(data2$epiDens,data2$streamDens,Y[3,1,,]),MARGIN = 1,FUN=mean)
  ylim=range(splitFunc(data2$epiDens,data2$streamDens,Y[3,1,,]))
  plot(rOut,type='l',ylim=ylim,ylab='r day-1')
  for(i in 1:nEn){
    lines(splitFunc(data2$epiDens,data2$streamDens,Y[3,1,,i]),col='gray',ylab='')
  }
  lines(rOut,ylab='')
  
  windows()
  rOut<-apply(Y[1,1,,],MARGIN = 1,FUN=mean)
  ylim=range(Y[1,1,,])
  plot(rOut,type='l',ylim=ylim,ylab='r day-1')
  for(i in 1:nEn){
    lines(Y[1,1,,i],col='gray',ylab='')
  }
  lines(rOut,ylab='')
  
  
  # # parameter convergence 
  windows()
  plot(P[1,1,]/rOut,pch=16,main='r20 convergence')
  
  plot(S[2,1,])
  
  # plotting pools not concentration 
  windows()
  DOCout<-apply(Y[5,1,,],MARGIN = 1,FUN=mean)
  ylim=range(c(Y[5,1,,],z$y[2,1,,1]),na.rm=T)
  plot(DOCout,type='l',ylim=ylim,ylab='DOC mg/L')
  for(i in 1:nEn){
    lines(Y[5,1,,i],col='gray',ylab='')
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
  
   
  # #plotting r vs. SWin
  # # 
  # windows()
  # rOut<-apply(rFunc(Y[1,1,,],data2$wtr),MARGIN = 1,FUN=mean)
  # ylim=range(rFunc(Y[1,1,,],data2$wtr))
  # plot(rOut,type='l',ylim=ylim,ylab='r day-1')
  # for(i in 1:nEn){
  #   lines(rFunc(Y[1,1,,i],data2$wtr),col='gray',ylab='')
  # }
  # lines(rOut,ylab='')
  # lines(apply(Y[1,1,,],MARGIN=1,FUN=mean),col='red') # r20 is in red
  # par(new=T)
  # plot(data2$docIn,type='l',col='brown',ylab='',yaxt='n')
  # # parameter considered to be converged when sd is 1/4 of initial sd 
  # dSD0<-sd(Y[1,1,1,]) #initial sd in d20 
  # dSD<-apply(Y[1,1,,],MARGIN = 1,FUN=sd)
  # dConverged<-dSD/dSD0 # fraction of initial sd 
  # # 
  # # windows()
  # # plot(ma(data2$NEP[spinUpLength:length(data2$NEP)],n = 7)~apply(Y[1,1,spinUpLength:length(data2$NEP),],MARGIN = 1,FUN=mean),pch=16,cex=1.5)
  # # abline(lm(ma(data2$NEP[spinUpLength:length(data2$NEP)],n=7)~apply(Y[1,1,spinUpLength:length(data2$NEP),],MARGIN = 1,FUN=mean)),lwd=2,lty=2)
  # # summary(lm(ma(data2$NEP[spinUpLength:length(data2$NEP)],n=7)~apply(Y[1,1,spinUpLength:length(data2$NEP),],MARGIN = 1,FUN=mean)))
  # # 
  # # windows()
  # # plot(data2$NEP[spinUpLength:length(data2$NEP)]~apply(rFunc(Y[1,1,spinUpLength:length(data2$NEP),],data2$wtr[spinUpLength:length(data2$NEP)]),MARGIN = 1,FUN=mean))
  # # abline(lm(data2$NEP[spinUpLength:length(data2$NEP)]~apply(rFunc(Y[1,1,spinUpLength:length(data2$NEP),],data2$wtr[spinUpLength:length(data2$NEP)]),MARGIN = 1,FUN=mean)))
  # 
  # # source('/Users/Jake/Desktop/R functions/backgroundR.R')
  # # windows()
  # # plot(backgroundR(data2$R,data2$GPP,n=10)$slope[spinUpLength:nStep],type='o')
  # # par(new=T)
  # # plot(apply(Y[1,1,spinUpLength:length(data2$NEP),],MARGIN = 1,FUN=mean,type='o',col='red'))
  # windows()
  # plot(data2$docIn,type='l',ylab='')
  # par(new=T)
  # plot(apply(Y[1,1,,],MARGIN = 1,FUN=mean),type='o',col='red',ylab='',yaxt='n')
  # axis(4)
  # 
  # # RMSE for states ; only 2015 data (evaluation years)
  windows()
  plot(sqrt((apply(Y[4,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12-z$y[1,1,,1]/data2$epiVol*12)^2),
       pch=16,ylab='DIC RMSE',type='o',xlab='time')
  windows()
  plot(sqrt((apply(Y[5,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12-z$y[2,1,,1]/data2$epiVol*12)^2),
       pch=16,ylab='DOC RMSE',type='o',xlab='time')
  # 
  data2$timeStep<-seq(1,nStep)
  # # mean RMSE for states 
  print(sqrt(mean((apply(Y[4,1,assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[assim_obs==0]*12-z$y[1,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12)^2,na.rm=T)))
  print(sqrt(mean((apply(Y[5,1,assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[assim_obs==0]*12-z$y[2,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12)^2,na.rm=T)))
  
  cor(apply(Y[4,1,!is.na(z$y[1,1,,1])&assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[!is.na(z$y[1,1,,1])&assim_obs==0]*12,z$y[1,1,!is.na(z$y[1,1,,1])&assim_obs==0,1]/data2$epiVol[!is.na(z$y[1,1,,1])&assim_obs==0]*12)^2
  cor(apply(Y[5,1,!is.na(z$y[2,1,,1])&assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[!is.na(z$y[2,1,,1])&assim_obs==0]*12,z$y[2,1,!is.na(z$y[2,1,,1])&assim_obs==0,1]/data2$epiVol[!is.na(z$y[2,1,,1])&assim_obs==0]*12)^2
  
  #bias 
  mean(apply(Y[4,1,assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[assim_obs==0]*12-z$y[1,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12,na.rm = T)^2
  mean(apply(Y[5,1,assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[assim_obs==0]*12-z$y[2,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12,na.rm=T)^2
  
  # # mean RMSE for states; mol C 
  print(sqrt(mean((apply(Y[4,1,assim_obs==0,],MARGIN = 1,FUN=mean)-z$y[1,1,assim_obs==0,1])^2,na.rm=T)))
  print(sqrt(mean((apply(Y[5,1,assim_obs==0,],MARGIN = 1,FUN=mean)-z$y[2,1,assim_obs==0,1])^2,na.rm=T)))
  
  # r2; mol C 
  cor(apply(Y[4,1,!is.na(z$y[1,1,,1])&assim_obs==0,],MARGIN = 1,FUN=mean),z$y[1,1,!is.na(z$y[1,1,,1])&assim_obs==0,1])^2
  cor(apply(Y[5,1,!is.na(z$y[2,1,,1])&assim_obs==0,],MARGIN = 1,FUN=mean),z$y[2,1,!is.na(z$y[2,1,,1])&assim_obs==0,1])^2
  
  # AIC CO2
  fit = apply(Y[4,1,assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[assim_obs==0]*12
  obs = z$y[1,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12
  AIC(logLik(lm(obs~fit)))
  
  # AIC DOC 
  fit = apply(Y[5,1,assim_obs==0,],MARGIN = 1,FUN=mean)/data2$epiVol[assim_obs==0]*12
  obs = z$y[2,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12
  AIC(logLik(lm(obs~fit)))
  
  
# range of mean d across run 
range(apply(Y[1,1,,],MARGIN = 1,FUN=mean))

ccf(data2$totalAlloC[1:100],apply(Y[1,1,,1:100],MARGIN = 1,FUN=mean),lag.max = 30)

oneDOC_d<-Y

windows()
plot(sqrt((apply(Y[2,1,spinUpLength:nStep,],MARGIN = 1,FUN=mean)/data2$epiVol[spinUpLength:nStep]*12-z$y[1,1,spinUpLength:nStep,1]/data2$epiVol[spinUpLength:nStep]*12)^2)~
  sqrt((apply(Y[4,1,spinUpLength:nStep,],MARGIN = 1,FUN=mean)/data2$epiVol[spinUpLength:nStep]*12-z$y[2,1,spinUpLength:nStep,1]/data2$epiVol[spinUpLength:nStep]*12)^2),
  pch=16,cex=1.5)

abline(lm(sqrt((apply(Y[2,1,spinUpLength:nStep,],MARGIN = 1,FUN=mean)/data2$epiVol[spinUpLength:nStep]*12-z$y[1,1,spinUpLength:nStep,1]/data2$epiVol[spinUpLength:nStep]*12)^2)~
  sqrt((apply(Y[4,1,spinUpLength:nStep,],MARGIN = 1,FUN=mean)/data2$epiVol[spinUpLength:nStep]*12-z$y[2,1,spinUpLength:nStep,1]/data2$epiVol[spinUpLength:nStep]*12)^2)))

summary(lm(sqrt((apply(Y[2,1,spinUpLength:nStep,],MARGIN = 1,FUN=mean)/data2$epiVol[spinUpLength:nStep]*12-z$y[1,1,spinUpLength:nStep,1]/data2$epiVol[spinUpLength:nStep]*12)^2)~
     sqrt((apply(Y[4,1,spinUpLength:nStep,],MARGIN = 1,FUN=mean)/data2$epiVol[spinUpLength:nStep]*12-z$y[2,1,spinUpLength:nStep,1]/data2$epiVol[spinUpLength:nStep]*12)^2)))

out<-apply(Y[2,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12
out<-cbind(out,(z$y[1,1,,1]/data2$epiVol*12))
resCO2<-out[,1]-out[,2]

out<-apply(Y[4,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12
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
plot(resCO2,pch=16,xlab = '',ylab='',xaxt='n',yaxt='n',cex=2)

plot(data2$totalAlloC~resDOC,cex=2,pch=16)

########################## Plotting for MS ###################################

# Figure 1 
# png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig1_DOC_CO2.png', 
#     res=300, width=14, height=7, units = 'in')
windows()
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=2
ylim=c(0.5,3)
par(mar=c(5,6,4,2),mfrow=c(1,2))
DOCout<-apply(Y[5,1,,],MARGIN = 1,FUN=mean)
ylim=range(c(Y[5,1,,],z$y[2,1,,1]+docPoolSD,z$y[2,1,,1]-docPoolSD),na.rm=T)
plot(DOCout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),type='l',ylim=ylim,lwd=lwd,cex.axis=1.5,
     ylab=expression(DOC~(mol~C)),xlab='',cex.lab=cex)
for(i in 1:nEn){
  lines(Y[5,1,(spinUpLength+1):nStep,i]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
        col='gray',ylab='',lwd=lwd)
}
lines(DOCout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),ylab='',lwd=lwd)
par(new=T)
plot(z$y[2,1,(spinUpLength+1):nStep,1]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),ylab=expression(DOC~(mol~C)),cex.lab=cex,cex.axis=1.5,
     ylim=ylim,col='red',pch=16,cex=cex,xlab='')
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
plot(z$y[1,1,(spinUpLength+1):nStep,1]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),ylab=expression(CO[2]~(mol~C)),cex.lab=cex,cex.axis=1.5,
     ylim=ylim,col='red',pch=16,cex=cex,xlab='')
arrows(as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[1,1,(spinUpLength+1):nStep,1]-dicPoolSD[(spinUpLength+1):nStep],
       as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[1,1,(spinUpLength+1):nStep,1]+dicPoolSD[(spinUpLength+1):nStep],
       code=3,length=0.1,angle=90,col='red',lwd=lwd)
# dev.off()

# Figure 1 plotting concentration 
################
png('/Users/jzwart/LakeCarbonEnKF/Figures/Fig1_all.png',
    res=300, width=14, height=14, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels
leg = 0.05 # distance from corner for panel label 

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
lines(DOCout[spinUpLength+1:nStep]~as.POSIXct(data2$datetime[spinUpLength+1:nStep]),ylab='',lwd=lwd,col='gray30')
par(new=T)
plot(z$y[2,1,(spinUpLength+1):nStep,1]/data2$epiVol*12~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
     ylim=ylim,col='black',pch=16,ylab='',cex=cex,xaxt='n',yaxt='n',xlab='')
arrows(as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[2,1,(spinUpLength+1):nStep,1]/data2$epiVol*12-docPoolSD[(spinUpLength+1):nStep]/data2$epiVol*12,
       as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[2,1,(spinUpLength+1):nStep,1]/data2$epiVol*12+docPoolSD[(spinUpLength+1):nStep]/data2$epiVol*12,
       code=3,length=0.1,angle=90,col='black',lwd=lwd)
legend("topleft", legend=c("Estimated State","Ensemble Mean",'Observed State'),
       col=c('gray','gray30','black'),pt.bg=c('gray','gray30','black'), cex = cex.leg,
       ncol=1,lwd=c(4,4,0),bty='n',lty=c(1,1,0),pt.cex = c(0,0,2),pch=c(0,0,16))
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
par(new=T)
plot(z$y[1,1,(spinUpLength+1):nStep,1]/data2$epiVol*12~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
     ylim=ylim,col='black',pch=16,ylab='',cex=cex,xaxt='n',yaxt='n',xlab='')
arrows(as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[1,1,(spinUpLength+1):nStep,1]/data2$epiVol*12-dicPoolSD[(spinUpLength+1):nStep]/data2$epiVol*12,
       as.POSIXct(data2$datetime[spinUpLength+1:nStep]),z$y[1,1,(spinUpLength+1):nStep,1]/data2$epiVol*12+dicPoolSD[(spinUpLength+1):nStep]/data2$epiVol*12,
       code=3,length=0.1,angle=90,col='black',lwd=lwd)
legend("topleft", legend=c("Estimated State","Ensemble Mean",'Observed State'),
       col=c('gray','gray30','black'),pt.bg=c('gray','gray30','black'), cex = cex.leg,
       ncol=1,lwd=c(4,4,0),bty='n',lty=c(1,1,0),pt.cex = c(0,0,2),pch=c(0,0,16))
text(x=xlim[2]-leg*(xlim[2]-xlim[1]),y = ylim[1]+leg*(ylim[2]-ylim[1]),labels = 'B', cex = cex.lab)

# d 
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

# parameter considered to be converged when sd is 1/4 of initial sd
dSD0<-sd(Y[1,1,1,]) #initial sd in d20
dSD<-apply(Y[1,1,,],MARGIN = 1,FUN=sd)
dConverged<-dSD/dSD0 # fraction of initial sd
ylim=c(0.5,3)
par(mar=c(5,6,4,2))
rOut<-apply(rFunc(Y[1,1,,],data2$wtr),MARGIN = 1,FUN=mean)
ylim=range(rFunc(Y[1,1,(spinUpLength+1):nStep,],data2$wtr[(spinUpLength+1):nStep]))
plot(rOut[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
     type='l',ylim=ylim,ylab=expression(d~(day^-1)),lwd=lwd,xlab='',cex.lab=cex,cex.axis=cex.axis)
for(i in 1:nEn){
  lines(rFunc(Y[1,1,(spinUpLength+1):nStep,i],data2$wtr[(spinUpLength+1):nStep])~
          as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),col='gray',ylab='',lwd=lwd,xlab='')
}
lines(rOut[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),ylab='',lwd=3,col='gray30',
      xlab='')
legend("topleft", legend=c("Estimated State","Ensemble Mean"),
       col=c('gray','gray30'),pt.bg=c('gray','gray30'), cex = cex.leg,
       ncol=1,lwd=c(4,4),bty='n',lty=c(1,1),pt.cex = c(0,0),pch=c(0,0))
text(x=xlim[2]-leg*(xlim[2]-xlim[1]),y = ylim[1]+leg*(ylim[2]-ylim[1]),labels = 'C', cex = cex.lab)

# d20 
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels
dSD0<-sd(Y[1,1,1,]) #initial sd in d20
dSD<-apply(Y[1,1,,],MARGIN = 1,FUN=sd)
dConverged<-dSD/dSD0 # fraction of initial sd
# ylim=c(0.5,3)
par(mar=c(5,7,4,7))
rOut<-apply(Y[1,1,,],MARGIN = 1,FUN=mean)
# ylim=range(rFunc(Y[1,1,(spinUpLength+1):nStep,],data2$wtr[(spinUpLength+1):nStep]))
# ylim=c(0.00,0.012)
plot(rOut[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
     type='l',ylim=ylim,ylab=expression(d[20]~(day^-1)),lwd=lwd,xlab='',cex.lab=cex,cex.axis=cex.axis)
for(i in 1:nEn){
  lines(Y[1,1,(spinUpLength+1):nStep,i]~
          as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),col='gray',ylab='',lwd=lwd,xlab='')
}
lines(rOut[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),ylab='',lwd=3,col='gray30',
      xlab='')
# lines(rOut[dConverged<0.25]~as.POSIXct(data2$datetime[dConverged<0.25]),ylab='',lwd=lwd,col='red',
#       xlab='')
par(new=T)
plot(data2$docIn[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),type='l',lwd=4,
     col='black',yaxt='n',xaxt='n',ylab='',xlab='',lty=6)
axis(4,cex.axis=cex.axis)
mtext(expression(loadDOC~(mol~C~day^-1)),4,line = 5,cex=2.5)
legend("top", legend=c("d20 Estimate",'d20 Ensemble Mean','DOC loading'),
       col=c('gray','gray30','black'),pt.bg=c('gray','gray30','black'), cex = cex.leg,
       ncol=1,lwd=c(4,4,4),bty='n',lty=c(1,1,3))
text(x=xlim[2]-.01*(xlim[2]-xlim[1]),y = ylim[1]+.01*(ylim[2]-ylim[1]),labels = 'D', cex = cex.lab)

dev.off()

################


# Figure 1 plotting concentration observed vs. predicted with error bars 
# png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig1_DOC_CO2_concentration_obsVpred.png', 
#     res=300, width=14, height=7, units = 'in')
windows()
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=2
par(mar=c(5,6,4,2),mfrow=c(1,2))
DOCout<-apply(Y[3,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
DOCpredictSD<-apply(Y[3,1,,]/data2$epiVol*12,MARGIN = 1,FUN=sd)
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

DICout<-apply(Y[2,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
DICpredictSD<-apply(Y[2,1,,]/data2$epiVol*12,MARGIN = 1,FUN=sd)
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
# dev.off()

# Figure 1 plotting pools of observed vs. predicted with error bars 
# png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig1_DOC_CO2_pools_obsVpred.png', 
#     res=300, width=14, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=2
par(mar=c(5,6,4,2),mfrow=c(1,2))
DOCout<-apply(Y[3,1,,],MARGIN = 1,FUN=mean)
DOCpredictSD<-apply(Y[3,1,,],MARGIN = 1,FUN=sd)
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

DICout<-apply(Y[2,1,,],MARGIN = 1,FUN=mean)
DICpredictSD<-apply(Y[2,1,,],MARGIN = 1,FUN=sd)
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
# dev.off()



# Figure 2 
# range of mean d across run 
range(apply(Y[1,1,,],MARGIN = 1,FUN=mean))
sd(apply(Y[1,1,,],MARGIN = 1,FUN=mean))/mean(apply(Y[1,1,,],MARGIN = 1,FUN=mean))
sd(Y[1,1,,])/mean(Y[1,1,,])

range(apply(rFunc(Y[1,1,,],data2$epiTemp),MARGIN = 1,FUN=mean))


mean(apply(rFunc(Y[1,1,,],data2$wtr),MARGIN = 1,FUN=mean))
range(apply(rFunc(Y[1,1,,],data2$wtr),MARGIN = 1,FUN=mean))
png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig1_d.png',
    res=300, width=7, height=7, units = 'in')
# windows()
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

# parameter considered to be converged when sd is 1/4 of initial sd
dSD0<-sd(Y[1,1,1,]) #initial sd in d20
dSD<-apply(Y[1,1,,],MARGIN = 1,FUN=sd)
dConverged<-dSD/dSD0 # fraction of initial sd
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
lines(rOut[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),ylab='',lwd=3,col='gray30',
      xlab='')
legend("topleft", legend=c("Estimated State","Ensemble Mean"),
       col=c('gray','gray30'),pt.bg=c('gray','gray30'), 
       ncol=1,lwd=c(4,4),bty='n',lty=c(1,1),pt.cex = c(0,0),pch=c(0,0))
# plotting red for when parameter is converged 
# lines(rOut[dConverged<0.25]~as.POSIXct(data2$datetime[dConverged<0.25]),ylab='',lwd=lwd,col='red',
#       xlab='')
dev.off()

# Figure 3 
range(apply(Y[1,1,,],MARGIN = 1,FUN=mean))
sd(apply(Y[1,1,,],MARGIN = 1,FUN=mean))/mean(apply(Y[1,1,,],MARGIN = 1,FUN=mean))
png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig1_d20_carbonLoad.png',
    res=300, width=7, height=7, units = 'in')
# windows()
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

dSD0<-sd(Y[1,1,1,]) #initial sd in d20
dSD<-apply(Y[1,1,,],MARGIN = 1,FUN=sd)
dConverged<-dSD/dSD0 # fraction of initial sd
cex=2
lwd=2
# ylim=c(0.5,3)
par(mar=c(5,6,4,6))
rOut<-apply(Y[1,1,,],MARGIN = 1,FUN=mean)
# ylim=range(rFunc(Y[1,1,(spinUpLength+1):nStep,],data2$wtr[(spinUpLength+1):nStep]))
# ylim=c(0.00,0.012)
plot(rOut[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),
     type='l',ylim=ylim,ylab=expression(d[20]~(day^-1)),lwd=lwd,xlab='',cex.lab=cex,cex.axis=1.5)
for(i in 1:nEn){
  lines(Y[1,1,(spinUpLength+1):nStep,i]~
          as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),col='gray',ylab='',lwd=lwd,xlab='')
}
lines(rOut[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),ylab='',lwd=3,col='gray30',
      xlab='')
# lines(rOut[dConverged<0.25]~as.POSIXct(data2$datetime[dConverged<0.25]),ylab='',lwd=lwd,col='red',
#       xlab='')
par(new=T)
plot(data2$docIn[(spinUpLength+1):nStep]~as.POSIXct(data2$datetime[(spinUpLength+1):nStep]),type='l',lwd=4,
     col='black',yaxt='n',xaxt='n',ylab='',xlab='',lty=6)
axis(4,cex.axis=1.5)
mtext(expression(loadDOC~(mol~C~day^-1)),4,line = 4,cex=cex)
legend("top", legend=c("d20 Estimate",'d20 Ensemble Mean','DOC loading'),
       col=c('gray','gray30','black'),pt.bg=c('gray','gray30','black'), 
       ncol=1,lwd=c(4,4,4),bty='n',lty=c(1,1,3))

dev.off()


# Figure 4 
# moving average weight of GPP and R to estimate NEP weighted 
GPPweight<-ma.weighted(data2$GPP[(spinUpLength+1):nStep],
                       cv = data2$sdGPP[(spinUpLength+1):nStep]/data2$GPP[(spinUpLength+1):nStep],n=7)
Rweight<-ma.weighted(data2$R[(spinUpLength+1):nStep],
                       cv = data2$sdR[(spinUpLength+1):nStep]/data2$R[(spinUpLength+1):nStep],n=7)
NEPweight<-GPPweight-Rweight

# png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig4_d20_NEP.png', 
#     res=300, width=7, height=7, units = 'in')
windows()
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
# dev.off()

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
DOCout<-apply(Y[3,1,,],MARGIN = 1,FUN=mean)
DOCpredictSD<-apply(Y[3,1,,],MARGIN = 1,FUN=sd)
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
plot(predNEP[,2],predNEP[,1],
     pch=16,cex=cex,cex.axis=cex,cex.lab=cex,xlab=expression(Observed~NEP~(mol~C~day^-1)),
     ylab=expression(Predicted~NEP~(mol~C~day^-1)))
abline(lm(predNEP[,1]~predNEP[,2]),lwd=lwd,lty=2)
dev.off()

# background R 
backR<-backgroundR(data2$R,data2$GPP,n=30)
plot(backR$intercept,type='l')
par(new=T)
plot(predNEP[,1],type='l',col='red')


# CO2 flux from lake modeled vs. linearly interpolated observations 
windows()
co2FluxSD=dicPoolSD*data2$kCO2/data2$thermo.depth
co2Flux=data2$dic*data2$kCO2/data2$thermo.depth
plot(co2Flux-data2$DICeq*data2$kCO2/data2$thermo.depth,pch=16,ylim=c(100,1100))
arrows(1:length(data2$datetime),co2Flux-co2FluxSD,
       1:length(data2$datetime),co2Flux+co2FluxSD,
       code=3,length=0.1,angle=90,col='red',lwd=2)
DICout<-apply(z$X[1,1,,],MARGIN = 1,FUN=mean)
lines(DICout*data2$kCO2/data2$thermo.depth-data2$DICeq*data2$kCO2/data2$thermo.depth,pch=16,lwd=2)
for(i in 1:nEn){
  lines(z$X[1,1,,i]*data2$kCO2/data2$thermo.depth-data2$DICeq*data2$kCO2/data2$thermo.depth,col='gray',ylab='')
}
lines(DICout*data2$kCO2/data2$thermo.depth-data2$DICeq*data2$kCO2/data2$thermo.depth,pch=16,lwd=2)

points(co2Flux-data2$DICeq*data2$kCO2/data2$thermo.depth,col='red',pch=16,ylim=c(100,1100))
arrows(1:length(data2$datetime),co2Flux-co2FluxSD,
       1:length(data2$datetime),co2Flux+co2FluxSD,
       code=3,length=0.1,angle=90,col='red',lwd=2)






