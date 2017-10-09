twoDOCpools_noMM_Temp<-function(obs,freq,reps,rVec_slow,rVec_fast,fracVec,fracInVec,n){
  docOut<-c()
  docOutBias<-c()
  docOutPrecision<-c()
  dicOut<-c()
  dicOutBias<-c()
  dicOutPrecision<-c()
  obsDICout<-c()
  obsDICoutBias<-c()
  obsDICoutPrecision<-c()
  obsDOCout<-c()
  obsDOCoutBias<-c()
  obsDOCoutPrecision<-c()
  mdfobsDOCout<-c()
  mdfobsDOCoutBias<-c()
  mdfobsDOCoutPrecision<-c()
  mdfobsDICout<-c()
  mdfobsDICoutBias<-c()
  mdfobsDICoutPrecision<-c()
  mdfDOCcv<-c()
  mdfDICcv<-c()
  for(m in 1:n){
    fracLabile0<-0.01 # fraction of initial DOC pool that is labile
    fracLabile<-0.08 # fraction loaded t-DOC that is labile; estimate from Berggren et al. 2010 Eco Letts & ISME
    
    # setting up initial parameter values for all ensemble members 
    rVec<-rVec_slow
    rVec_fast<-rVec_fast
    fracVec<-fracVec
    fracInVec<-fracInVec
    
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
    
    
    # initial B transition matrix for each ensemble  
    B<-array(NA,dim=c(4,4,nStep,nEn)) # array of transition matrix B[a,b,c,d]; 
    # where [a,b] is transition matrix DIC & DOCr & DOCl, c=timeStep, and d=ensemble member 
    # initial parameters of B at timestep 1 
    for(i in 1:nEn){
      B[,,1,i]<-matrix(c(1-data2$QoutInt[1]/data2$epiVol[1]-data2$kCO2[1]/data2$thermo.depth[1]+
                           (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                         rFunc(rVec[i],data2$wtr[1]),rFunc(rVec_fast[i],data2$wtr[1]),0,
                         0,1-rFunc(rVec[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],0,0,
                         0,0,1-rFunc(rVec_fast[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],0,
                         0,1-rFunc(rVec[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                         1-rFunc(rVec_fast[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],0),
                       nrow=4,byrow=T)
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
      h[1,5,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic 
      h[2,8,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) #doc total (we only have data on total DOC pool)
    }
    
    P <- array(0,dim=c(4,4,nStep))
    
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
    Y[4,1,1,]<-fracInVec # fraction of stream that goes into epi 
    Y[5,1,1,]<-z$X[1,1,1,] # DIC state  
    Y[6,1,1,]<-z$X[2,1,1,] # DOC recalcitrant state  
    Y[7,1,1,]<-z$X[3,1,1,] # DOC labile state  
    Y[8,1,1,]<-z$X[4,1,1,] # DOC total state  
    
    
    #Iterate through time
    for(t in 2:nStep){
      for(i in 1:nEn){
        # Forecasting; need to update parameters, 
        z$pars[1:4,1,t,i]<-Y[1:4,1,t-1,i] # r20
        z$X[1:4,1,t-1,i]<-Y[5:8,1,t-1,i] # updating state variables from Y 
        
        #Predictions
        z$X[,,t,i]<-z$B[,,t-1,i]%*%z$X[,,t-1,i] + z$C[,,t-1,i]%*%z$ut[,,t-1,i] # forecasting state variable predictions 
        z$B[,,t,i]<-matrix(c(1-data2$QoutInt[t]/data2$epiVol[t]-data2$kCO2[t]/data2$thermo.depth[t]+
                               (data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t], # zmix is fraction of pool available for exchange with atm; parameters estimates are the same as previous time step
                             rFunc(z$pars[1,1,t,i],data2$wtr[t]),rFunc(z$pars[2,1,t,i],data2$wtr[t]),0,
                             0,1-rFunc(z$pars[1,1,t,i],data2$wtr[t])-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],0,0,
                             0,0,1-rFunc(z$pars[2,1,t,i],data2$wtr[t])-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],0,
                             0,1-rFunc(z$pars[1,1,t,i],data2$wtr[t])-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],
                             1-rFunc(z$pars[2,1,t,i],data2$wtr[t])-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],0),
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
    DOCout<-apply(Y[8,1,,],MARGIN = 1,FUN=mean)
    DOCout<-DOCout/data2$epiVol*12
    TrueDOC<-true[8,1,,1]/data2$epiVol*12 
    DOCrmse<-sqrt(mean((DOCout-TrueDOC)^2))
    DOCbias<-mean(DOCout-TrueDOC,na.rm=T)^2
    DOCprecision<-summary(lm(DOCout~TrueDOC))$r.squared
    docOut<-rbind(docOut,DOCrmse)
    docOutBias<-rbind(docOutBias,DOCbias)
    docOutPrecision<-rbind(docOutPrecision,DOCprecision)
    obsDOC<-z$y[2,1,,1]/data2$epiVol*12
    obsDOCrmse<-sqrt(mean((obsDOC-TrueDOC)^2,na.rm = T))
    obsDOCBias<-mean(obsDOC-TrueDOC,na.rm = T)^2
    obsDOCPrecision<-summary(lm(obsDOC~TrueDOC))$r.squared
    obsDOCout<-rbind(obsDOCout,obsDOCrmse)
    obsDOCoutBias<-rbind(obsDOCoutBias,obsDOCBias)
    obsDOCoutPrecision<-rbind(obsDOCoutPrecision,obsDOCPrecision)
    mdfobsDOCrmse<-sqrt(mean((DOCout-obsDOC)^2,na.rm = T))
    mdfobsDOCbias<-mean(DOCout-obsDOC,na.rm = T)^2
    mdfobsDOCprecision<-summary(lm(DOCout~obsDOC))$r.squared
    mdfobsDOCout<-rbind(mdfobsDOCout,mdfobsDOCrmse)
    mdfobsDOCoutBias<-rbind(mdfobsDOCoutBias,mdfobsDOCbias)
    mdfobsDOCoutPrecision<-rbind(mdfobsDOCoutPrecision,mdfobsDOCprecision)
    DOCcv<-mean(apply(Y[8,1,,],MARGIN = 1,FUN=sd)/apply(Y[8,1,,],MARGIN = 1,FUN=mean))
    mdfDOCcv<-rbind(mdfDOCcv,DOCcv)
    
    DICout<-apply(Y[5,1,,],MARGIN = 1,FUN=mean)
    DICout<-DICout/data2$epiVol*12
    TrueDIC<-true[5,1,,1]/data2$epiVol*12 
    DICrmse<-sqrt(mean((DICout-TrueDIC)^2))
    DICbias<-mean(DICout-TrueDIC,na.rm=T)^2
    DICprecision<-summary(lm(DICout~TrueDIC))$r.squared
    dicOut<-rbind(dicOut,DICrmse)
    dicOutBias<-rbind(dicOutBias,DICbias)
    dicOutPrecision<-rbind(dicOutPrecision,DICprecision)
    obsDIC<-z$y[1,1,,1]/data2$epiVol*12
    obsDICrmse<-sqrt(mean((obsDIC-TrueDIC)^2,na.rm = T))
    obsDICBias<-mean(obsDIC-TrueDIC,na.rm = T)^2
    obsDICPrecision<-summary(lm(obsDIC~TrueDIC))$r.squared
    obsDICout<-rbind(obsDICout,obsDICrmse)
    obsDICoutBias<-rbind(obsDICoutBias,obsDICBias)
    obsDICoutPrecision<-rbind(obsDICoutPrecision,obsDICPrecision)
    mdfobsDICrmse<-sqrt(mean((DICout-obsDIC)^2,na.rm = T))
    mdfobsDICbias<-mean(DICout-obsDIC,na.rm = T)^2
    mdfobsDICprecision<-summary(lm(DICout~obsDIC))$r.squared
    mdfobsDICout<-rbind(mdfobsDICout,mdfobsDICrmse)
    mdfobsDICoutBias<-rbind(mdfobsDICoutBias,mdfobsDICbias)
    mdfobsDICoutPrecision<-rbind(mdfobsDICoutPrecision,mdfobsDICprecision)
    DICcv<-mean(apply(Y[5,1,,],MARGIN = 1,FUN=sd)/apply(Y[5,1,,],MARGIN = 1,FUN=mean))
    mdfDICcv<-rbind(mdfDICcv,DICcv)
    
  }
  DOCrmse<-mean(docOut,na.rm = T)
  DOCbias<-mean(docOutBias,na.rm=T)
  DOCprecision<-mean(docOutPrecision,na.rm=T)
  obsDOCrmse<-mean(obsDOCout,na.rm=T)
  obsDOCbias<-mean(obsDOCoutBias,na.rm = T)
  obsDOCprecision<-mean(obsDOCoutPrecision,na.rm = T)
  mdfobsDOCrmse<-mean(mdfobsDOCout,na.rm=T)
  mdfobsDOCbias<-mean(mdfobsDOCoutBias,na.rm = T)
  mdfobsDOCprecision<-mean(mdfobsDOCoutPrecision,na.rm = T)
  mdfDOCcv<-mean(mdfDOCcv,na.rm = T)
  
  DICrmse<-mean(dicOut,na.rm = T)
  DICbias<-mean(dicOutBias,na.rm=T)
  DICprecision<-mean(dicOutPrecision,na.rm=T)
  obsDICrmse<-mean(obsDICout,na.rm=T)
  obsDICbias<-mean(obsDICoutBias,na.rm = T)
  obsDICprecision<-mean(obsDICoutPrecision,na.rm = T)
  mdfobsDICrmse<-mean(mdfobsDICout,na.rm=T)
  mdfobsDICbias<-mean(mdfobsDICoutBias,na.rm = T)
  mdfobsDICprecision<-mean(mdfobsDICoutPrecision,na.rm = T)
  mdfDICcv<-mean(mdfDICcv,na.rm = T)
  
  return(list("DOCrmse"=DOCrmse,"DOCbias"=DOCbias,'DOCprecision'=DOCprecision,
              "DICrmse"=DICrmse,'DICbias'=DICbias,'DICprecision'=DICprecision,
              "obsDOCrmse"=obsDOCrmse,'obsDOCbias'=obsDOCbias,'obsDOCprecision'=obsDOCprecision,
              "obsDICrmse"=obsDICrmse,'obsDICbias'=obsDICbias,'obsDICprecision'=obsDICprecision,
              'mdfobsDOCrmse'=mdfobsDOCrmse,'mdfobsDOCbias'=mdfobsDOCbias,'mdfobsDOCprecision'=mdfobsDOCprecision,
              'mdfobsDICrmse'=mdfobsDICrmse,'mdfobsDICbias'=mdfobsDICbias,'mdfobsDICprecision'=mdfobsDICprecision,
              'mdfDOCcv'=mdfDOCcv,'mdfDICcv'=mdfDICcv))
}

twoDOCpools_noMM_noTemp<-function(obs,freq,reps,rVec_slow,rVec_fast,fracVec,fracInVec,n){
  docOut<-c()
  docOutBias<-c()
  docOutPrecision<-c()
  dicOut<-c()
  dicOutBias<-c()
  dicOutPrecision<-c()
  obsDICout<-c()
  obsDICoutBias<-c()
  obsDICoutPrecision<-c()
  obsDOCout<-c()
  obsDOCoutBias<-c()
  obsDOCoutPrecision<-c()
  mdfobsDOCout<-c()
  mdfobsDOCoutBias<-c()
  mdfobsDOCoutPrecision<-c()
  mdfobsDICout<-c()
  mdfobsDICoutBias<-c()
  mdfobsDICoutPrecision<-c()
  mdfDOCcv<-c()
  mdfDICcv<-c()
  for(m in 1:n){
    fracLabile0<-0.01 # fraction of initial DOC pool that is labile
    fracLabile<-0.08 # fraction loaded t-DOC that is labile; estimate from Berggren et al. 2010 Eco Letts & ISME
    
    # setting up initial parameter values for all ensemble members 
    rVec<-rVec_slow
    rVec_fast<-rVec_fast
    fracVec<-fracVec
    fracInVec<-fracInVec
    
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
    
    
    # initial B transition matrix for each ensemble  
    B<-array(NA,dim=c(4,4,nStep,nEn)) # array of transition matrix B[a,b,c,d]; 
    # where [a,b] is transition matrix DIC & DOCr & DOCl, c=timeStep, and d=ensemble member 
    # initial parameters of B at timestep 1 
    for(i in 1:nEn){
      B[,,1,i]<-matrix(c(1-data2$QoutInt[1]/data2$epiVol[1]-data2$kCO2[1]/data2$thermo.depth[1]+
                           (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                         rVec[i],rVec_fast[i],0,
                         0,1-rVec[i]-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],0,0,
                         0,0,1-rVec_fast[i]-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],0,
                         0,1-rVec[i]-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                         1-rVec_fast[i]-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],0),
                       nrow=4,byrow=T)
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
      h[1,5,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic 
      h[2,8,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) #doc total (we only have data on total DOC pool)
    }
    
    P <- array(0,dim=c(4,4,nStep))
    
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
    Y[4,1,1,]<-fracInVec # fraction of stream that goes into epi 
    Y[5,1,1,]<-z$X[1,1,1,] # DIC state  
    Y[6,1,1,]<-z$X[2,1,1,] # DOC recalcitrant state  
    Y[7,1,1,]<-z$X[3,1,1,] # DOC labile state  
    Y[8,1,1,]<-z$X[4,1,1,] # DOC total state  
    
    
    #Iterate through time
    for(t in 2:nStep){
      for(i in 1:nEn){
        # Forecasting; need to update parameters, 
        z$pars[1:4,1,t,i]<-Y[1:4,1,t-1,i] # r20
        z$X[1:4,1,t-1,i]<-Y[5:8,1,t-1,i] # updating state variables from Y 
        
        #Predictions
        z$X[,,t,i]<-z$B[,,t-1,i]%*%z$X[,,t-1,i] + z$C[,,t-1,i]%*%z$ut[,,t-1,i] # forecasting state variable predictions 
        z$B[,,t,i]<-matrix(c(1-data2$QoutInt[t]/data2$epiVol[t]-data2$kCO2[t]/data2$thermo.depth[t]+
                               (data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t], # zmix is fraction of pool available for exchange with atm; parameters estimates are the same as previous time step
                             z$pars[1,1,t,i],z$pars[2,1,t,i],0,
                             0,1-z$pars[1,1,t,i]-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],0,0,
                             0,0,1-z$pars[2,1,t,i]-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],0,
                             0,1-z$pars[1,1,t,i]-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],
                             1-z$pars[2,1,t,i]-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],0),
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
    DOCout<-apply(Y[8,1,,],MARGIN = 1,FUN=mean)
    DOCout<-DOCout/data2$epiVol*12
    TrueDOC<-true[8,1,,1]/data2$epiVol*12 
    DOCrmse<-sqrt(mean((DOCout-TrueDOC)^2))
    DOCbias<-mean(DOCout-TrueDOC,na.rm=T)^2
    DOCprecision<-summary(lm(DOCout~TrueDOC))$r.squared
    docOut<-rbind(docOut,DOCrmse)
    docOutBias<-rbind(docOutBias,DOCbias)
    docOutPrecision<-rbind(docOutPrecision,DOCprecision)
    obsDOC<-z$y[2,1,,1]/data2$epiVol*12
    obsDOCrmse<-sqrt(mean((obsDOC-TrueDOC)^2,na.rm = T))
    obsDOCBias<-mean(obsDOC-TrueDOC,na.rm = T)^2
    obsDOCPrecision<-summary(lm(obsDOC~TrueDOC))$r.squared
    obsDOCout<-rbind(obsDOCout,obsDOCrmse)
    obsDOCoutBias<-rbind(obsDOCoutBias,obsDOCBias)
    obsDOCoutPrecision<-rbind(obsDOCoutPrecision,obsDOCPrecision)
    mdfobsDOCrmse<-sqrt(mean((DOCout-obsDOC)^2,na.rm = T))
    mdfobsDOCbias<-mean(DOCout-obsDOC,na.rm = T)^2
    mdfobsDOCprecision<-summary(lm(DOCout~obsDOC))$r.squared
    mdfobsDOCout<-rbind(mdfobsDOCout,mdfobsDOCrmse)
    mdfobsDOCoutBias<-rbind(mdfobsDOCoutBias,mdfobsDOCbias)
    mdfobsDOCoutPrecision<-rbind(mdfobsDOCoutPrecision,mdfobsDOCprecision)
    DOCcv<-mean(apply(Y[8,1,,],MARGIN = 1,FUN=sd)/apply(Y[8,1,,],MARGIN = 1,FUN=mean))
    mdfDOCcv<-rbind(mdfDOCcv,DOCcv)
    
    DICout<-apply(Y[5,1,,],MARGIN = 1,FUN=mean)
    DICout<-DICout/data2$epiVol*12
    TrueDIC<-true[5,1,,1]/data2$epiVol*12 
    DICrmse<-sqrt(mean((DICout-TrueDIC)^2))
    DICbias<-mean(DICout-TrueDIC,na.rm=T)^2
    DICprecision<-summary(lm(DICout~TrueDIC))$r.squared
    dicOut<-rbind(dicOut,DICrmse)
    dicOutBias<-rbind(dicOutBias,DICbias)
    dicOutPrecision<-rbind(dicOutPrecision,DICprecision)
    obsDIC<-z$y[1,1,,1]/data2$epiVol*12
    obsDICrmse<-sqrt(mean((obsDIC-TrueDIC)^2,na.rm = T))
    obsDICBias<-mean(obsDIC-TrueDIC,na.rm = T)^2
    obsDICPrecision<-summary(lm(obsDIC~TrueDIC))$r.squared
    obsDICout<-rbind(obsDICout,obsDICrmse)
    obsDICoutBias<-rbind(obsDICoutBias,obsDICBias)
    obsDICoutPrecision<-rbind(obsDICoutPrecision,obsDICPrecision)
    mdfobsDICrmse<-sqrt(mean((DICout-obsDIC)^2,na.rm = T))
    mdfobsDICbias<-mean(DICout-obsDIC,na.rm = T)^2
    mdfobsDICprecision<-summary(lm(DICout~obsDIC))$r.squared
    mdfobsDICout<-rbind(mdfobsDICout,mdfobsDICrmse)
    mdfobsDICoutBias<-rbind(mdfobsDICoutBias,mdfobsDICbias)
    mdfobsDICoutPrecision<-rbind(mdfobsDICoutPrecision,mdfobsDICprecision)
    DICcv<-mean(apply(Y[5,1,,],MARGIN = 1,FUN=sd)/apply(Y[5,1,,],MARGIN = 1,FUN=mean))
    mdfDICcv<-rbind(mdfDICcv,DICcv)
  }
  DOCrmse<-mean(docOut,na.rm = T)
  DOCbias<-mean(docOutBias,na.rm=T)
  DOCprecision<-mean(docOutPrecision,na.rm=T)
  obsDOCrmse<-mean(obsDOCout,na.rm=T)
  obsDOCbias<-mean(obsDOCoutBias,na.rm = T)
  obsDOCprecision<-mean(obsDOCoutPrecision,na.rm = T)
  mdfobsDOCrmse<-mean(mdfobsDOCout,na.rm=T)
  mdfobsDOCbias<-mean(mdfobsDOCoutBias,na.rm = T)
  mdfobsDOCprecision<-mean(mdfobsDOCoutPrecision,na.rm = T)
  mdfDOCcv<-mean(mdfDOCcv,na.rm = T)
  
  DICrmse<-mean(dicOut,na.rm = T)
  DICbias<-mean(dicOutBias,na.rm=T)
  DICprecision<-mean(dicOutPrecision,na.rm=T)
  obsDICrmse<-mean(obsDICout,na.rm=T)
  obsDICbias<-mean(obsDICoutBias,na.rm = T)
  obsDICprecision<-mean(obsDICoutPrecision,na.rm = T)
  mdfobsDICrmse<-mean(mdfobsDICout,na.rm=T)
  mdfobsDICbias<-mean(mdfobsDICoutBias,na.rm = T)
  mdfobsDICprecision<-mean(mdfobsDICoutPrecision,na.rm = T)
  mdfDICcv<-mean(mdfDICcv,na.rm = T)
  
  return(list("DOCrmse"=DOCrmse,"DOCbias"=DOCbias,'DOCprecision'=DOCprecision,
              "DICrmse"=DICrmse,'DICbias'=DICbias,'DICprecision'=DICprecision,
              "obsDOCrmse"=obsDOCrmse,'obsDOCbias'=obsDOCbias,'obsDOCprecision'=obsDOCprecision,
              "obsDICrmse"=obsDICrmse,'obsDICbias'=obsDICbias,'obsDICprecision'=obsDICprecision,
              'mdfobsDOCrmse'=mdfobsDOCrmse,'mdfobsDOCbias'=mdfobsDOCbias,'mdfobsDOCprecision'=mdfobsDOCprecision,
              'mdfobsDICrmse'=mdfobsDICrmse,'mdfobsDICbias'=mdfobsDICbias,'mdfobsDICprecision'=mdfobsDICprecision,
              'mdfDOCcv'=mdfDOCcv,'mdfDICcv'=mdfDICcv))
}

twoDOCpools_MM_noTemp<-function(obs,freq,reps,rVec_slow,rVec_fast,fracVec,fracInVec,n){
  docOut<-c()
  docOutBias<-c()
  docOutPrecision<-c()
  dicOut<-c()
  dicOutBias<-c()
  dicOutPrecision<-c()
  obsDICout<-c()
  obsDICoutBias<-c()
  obsDICoutPrecision<-c()
  obsDOCout<-c()
  obsDOCoutBias<-c()
  obsDOCoutPrecision<-c()
  mdfobsDOCout<-c()
  mdfobsDOCoutBias<-c()
  mdfobsDOCoutPrecision<-c()
  mdfobsDICout<-c()
  mdfobsDICoutBias<-c()
  mdfobsDICoutPrecision<-c()
  mdfDOCcv<-c()
  mdfDICcv<-c()
  for(m in 1:n){
    fracLabile0<-0.01 # fraction of initial DOC pool that is labile
    fracLabile<-0.08 # fraction loaded t-DOC that is labile; estimate from Berggren et al. 2010 Eco Letts & ISME
    
    # setting up initial parameter values for all ensemble members 
    rVec<-rVec_slow
    rVec_fast<-rVec_fast
    fracVec<-fracVec
    fracInVec<-fracInVec
    
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
    
    
    # initial B transition matrix for each ensemble  
    B<-array(NA,dim=c(4,4,nStep,nEn)) # array of transition matrix B[a,b,c,d]; 
    # where [a,b] is transition matrix DIC & DOCr & DOCl, c=timeStep, and d=ensemble member 
    # initial parameters of B at timestep 1 
    for(i in 1:nEn){
      B[,,1,i]<-matrix(c(1-data2$QoutInt[1]/data2$epiVol[1]-data2$kCO2[1]/data2$thermo.depth[1]+
                           (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                         mm(rVec[i],docVec[1],halfSat,data2$epiVol[1]),mm(rVec_fast[i],docVec[1],halfSat,data2$epiVol[1]),0,
                         0,1-mm(rVec[i],docVec[1],halfSat,data2$epiVol[1])-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],0,0,
                         0,0,1-mm(rVec_fast[i],docVec[1],halfSat,data2$epiVol[1])-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],0,
                         0,1-mm(rVec[i],docVec[1],halfSat,data2$epiVol[1])-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                         1-mm(rVec_fast[i],docVec[1],halfSat,data2$epiVol[1])-data2$QoutInt[1]/data2$epiVol[1]+(data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],0),
                       nrow=4,byrow=T)
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
      h[1,5,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic 
      h[2,8,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) #doc total (we only have data on total DOC pool)
    }
    
    P <- array(0,dim=c(4,4,nStep))
    
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
    Y[4,1,1,]<-fracInVec # fraction of stream that goes into epi 
    Y[5,1,1,]<-z$X[1,1,1,] # DIC state  
    Y[6,1,1,]<-z$X[2,1,1,] # DOC recalcitrant state  
    Y[7,1,1,]<-z$X[3,1,1,] # DOC labile state  
    Y[8,1,1,]<-z$X[4,1,1,] # DOC total state  
    
    
    #Iterate through time
    for(t in 2:nStep){
      for(i in 1:nEn){
        # Forecasting; need to update parameters, 
        z$pars[1:4,1,t,i]<-Y[1:4,1,t-1,i] # r20
        z$X[1:4,1,t-1,i]<-Y[5:8,1,t-1,i] # updating state variables from Y 
        
        #Predictions
        z$X[,,t,i]<-z$B[,,t-1,i]%*%z$X[,,t-1,i] + z$C[,,t-1,i]%*%z$ut[,,t-1,i] # forecasting state variable predictions 
        z$B[,,t,i]<-matrix(c(1-data2$QoutInt[t]/data2$epiVol[t]-data2$kCO2[t]/data2$thermo.depth[t]+
                               (data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t], # zmix is fraction of pool available for exchange with atm; parameters estimates are the same as previous time step
                             mm(z$pars[1,1,t,i],z$X[2,1,t,i],halfSat,data2$epiVol[t]),mm(z$pars[2,1,t,i],z$X[2,1,t,i],halfSat,data2$epiVol[t]),0,
                             0,1-mm(z$pars[1,1,t,i],z$X[2,1,t,i],halfSat,data2$epiVol[t])-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],0,0,
                             0,0,1-mm(z$pars[2,1,t,i],z$X[2,1,t,i],halfSat,data2$epiVol[t])-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],0,
                             0,1-mm(z$pars[1,1,t,i],z$X[2,1,t,i],halfSat,data2$epiVol[t])-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],
                             1-mm(z$pars[2,1,t,i],z$X[2,1,t,i],halfSat,data2$epiVol[t])-data2$QoutInt[t]/data2$epiVol[t]+(data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[4,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],0),
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
    DOCout<-apply(Y[8,1,,],MARGIN = 1,FUN=mean)
    DOCout<-DOCout/data2$epiVol*12
    TrueDOC<-true[8,1,,1]/data2$epiVol*12 
    DOCrmse<-sqrt(mean((DOCout-TrueDOC)^2))
    DOCbias<-mean(DOCout-TrueDOC,na.rm=T)^2
    DOCprecision<-summary(lm(DOCout~TrueDOC))$r.squared
    docOut<-rbind(docOut,DOCrmse)
    docOutBias<-rbind(docOutBias,DOCbias)
    docOutPrecision<-rbind(docOutPrecision,DOCprecision)
    obsDOC<-z$y[2,1,,1]/data2$epiVol*12
    obsDOCrmse<-sqrt(mean((obsDOC-TrueDOC)^2,na.rm = T))
    obsDOCBias<-mean(obsDOC-TrueDOC,na.rm = T)^2
    obsDOCPrecision<-summary(lm(obsDOC~TrueDOC))$r.squared
    obsDOCout<-rbind(obsDOCout,obsDOCrmse)
    obsDOCoutBias<-rbind(obsDOCoutBias,obsDOCBias)
    obsDOCoutPrecision<-rbind(obsDOCoutPrecision,obsDOCPrecision)
    mdfobsDOCrmse<-sqrt(mean((DOCout-obsDOC)^2,na.rm = T))
    mdfobsDOCbias<-mean(DOCout-obsDOC,na.rm = T)^2
    mdfobsDOCprecision<-summary(lm(DOCout~obsDOC))$r.squared
    mdfobsDOCout<-rbind(mdfobsDOCout,mdfobsDOCrmse)
    mdfobsDOCoutBias<-rbind(mdfobsDOCoutBias,mdfobsDOCbias)
    mdfobsDOCoutPrecision<-rbind(mdfobsDOCoutPrecision,mdfobsDOCprecision)
    DOCcv<-mean(apply(Y[8,1,,],MARGIN = 1,FUN=sd)/apply(Y[8,1,,],MARGIN = 1,FUN=mean))
    mdfDOCcv<-rbind(mdfDOCcv,DOCcv)
    
    DICout<-apply(Y[5,1,,],MARGIN = 1,FUN=mean)
    DICout<-DICout/data2$epiVol*12
    TrueDIC<-true[5,1,,1]/data2$epiVol*12 
    DICrmse<-sqrt(mean((DICout-TrueDIC)^2))
    DICbias<-mean(DICout-TrueDIC,na.rm=T)^2
    DICprecision<-summary(lm(DICout~TrueDIC))$r.squared
    dicOut<-rbind(dicOut,DICrmse)
    dicOutBias<-rbind(dicOutBias,DICbias)
    dicOutPrecision<-rbind(dicOutPrecision,DICprecision)
    obsDIC<-z$y[1,1,,1]/data2$epiVol*12
    obsDICrmse<-sqrt(mean((obsDIC-TrueDIC)^2,na.rm = T))
    obsDICBias<-mean(obsDIC-TrueDIC,na.rm = T)^2
    obsDICPrecision<-summary(lm(obsDIC~TrueDIC))$r.squared
    obsDICout<-rbind(obsDICout,obsDICrmse)
    obsDICoutBias<-rbind(obsDICoutBias,obsDICBias)
    obsDICoutPrecision<-rbind(obsDICoutPrecision,obsDICPrecision)
    mdfobsDICrmse<-sqrt(mean((DICout-obsDIC)^2,na.rm = T))
    mdfobsDICbias<-mean(DICout-obsDIC,na.rm = T)^2
    mdfobsDICprecision<-summary(lm(DICout~obsDIC))$r.squared
    mdfobsDICout<-rbind(mdfobsDICout,mdfobsDICrmse)
    mdfobsDICoutBias<-rbind(mdfobsDICoutBias,mdfobsDICbias)
    mdfobsDICoutPrecision<-rbind(mdfobsDICoutPrecision,mdfobsDICprecision)
    DICcv<-mean(apply(Y[5,1,,],MARGIN = 1,FUN=sd)/apply(Y[5,1,,],MARGIN = 1,FUN=mean))
    mdfDICcv<-rbind(mdfDICcv,DICcv)
  }
  DOCrmse<-mean(docOut,na.rm = T)
  DOCbias<-mean(docOutBias,na.rm=T)
  DOCprecision<-mean(docOutPrecision,na.rm=T)
  obsDOCrmse<-mean(obsDOCout,na.rm=T)
  obsDOCbias<-mean(obsDOCoutBias,na.rm = T)
  obsDOCprecision<-mean(obsDOCoutPrecision,na.rm = T)
  mdfobsDOCrmse<-mean(mdfobsDOCout,na.rm=T)
  mdfobsDOCbias<-mean(mdfobsDOCoutBias,na.rm = T)
  mdfobsDOCprecision<-mean(mdfobsDOCoutPrecision,na.rm = T)
  mdfDOCcv<-mean(mdfDOCcv,na.rm = T)
  
  DICrmse<-mean(dicOut,na.rm = T)
  DICbias<-mean(dicOutBias,na.rm=T)
  DICprecision<-mean(dicOutPrecision,na.rm=T)
  obsDICrmse<-mean(obsDICout,na.rm=T)
  obsDICbias<-mean(obsDICoutBias,na.rm = T)
  obsDICprecision<-mean(obsDICoutPrecision,na.rm = T)
  mdfobsDICrmse<-mean(mdfobsDICout,na.rm=T)
  mdfobsDICbias<-mean(mdfobsDICoutBias,na.rm = T)
  mdfobsDICprecision<-mean(mdfobsDICoutPrecision,na.rm = T)
  mdfDICcv<-mean(mdfDICcv,na.rm = T)
  
  return(list("DOCrmse"=DOCrmse,"DOCbias"=DOCbias,'DOCprecision'=DOCprecision,
              "DICrmse"=DICrmse,'DICbias'=DICbias,'DICprecision'=DICprecision,
              "obsDOCrmse"=obsDOCrmse,'obsDOCbias'=obsDOCbias,'obsDOCprecision'=obsDOCprecision,
              "obsDICrmse"=obsDICrmse,'obsDICbias'=obsDICbias,'obsDICprecision'=obsDICprecision,
              'mdfobsDOCrmse'=mdfobsDOCrmse,'mdfobsDOCbias'=mdfobsDOCbias,'mdfobsDOCprecision'=mdfobsDOCprecision,
              'mdfobsDICrmse'=mdfobsDICrmse,'mdfobsDICbias'=mdfobsDICbias,'mdfobsDICprecision'=mdfobsDICprecision,
              'mdfDOCcv'=mdfDOCcv,'mdfDICcv'=mdfDICcv))
}

twoDOCpools_MM_Temp<-function(obs,freq,reps,rVec_slow,rVec_fast,fracVec,fracInVec,n){
  docOut<-c()
  docOutBias<-c()
  docOutPrecision<-c()
  dicOut<-c()
  dicOutBias<-c()
  dicOutPrecision<-c()
  obsDICout<-c()
  obsDICoutBias<-c()
  obsDICoutPrecision<-c()
  obsDOCout<-c()
  obsDOCoutBias<-c()
  obsDOCoutPrecision<-c()
  mdfobsDOCout<-c()
  mdfobsDOCoutBias<-c()
  mdfobsDOCoutPrecision<-c()
  mdfobsDICout<-c()
  mdfobsDICoutBias<-c()
  mdfobsDICoutPrecision<-c()
  mdfDOCcv<-c()
  mdfDICcv<-c()
  for(m in 1:n){
    fracLabile0<-0.01 # fraction of initial DOC pool that is labile
    fracLabile<-0.08 # fraction loaded t-DOC that is labile; estimate from Berggren et al. 2010 Eco Letts & ISME
    
    # setting up initial parameter values for all ensemble members 
    rVec<-rVec_slow
    rVec_fast<-rVec_fast
    fracVec<-fracVec
    fracInVec<-fracInVec
    
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
      h[1,5,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic 
      h[2,8,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) #doc total (we only have data on total DOC pool)
    }
    
    P <- array(0,dim=c(4,4,nStep))
    
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
    Y[4,1,1,]<-fracInVec # fraction of stream that goes into epi 
    Y[5,1,1,]<-z$X[1,1,1,] # DIC state  
    Y[6,1,1,]<-z$X[2,1,1,] # DOC recalcitrant state  
    Y[7,1,1,]<-z$X[3,1,1,] # DOC labile state  
    Y[8,1,1,]<-z$X[4,1,1,] # DOC total state  
    
    
    #Iterate through time
    for(t in 2:nStep){
      for(i in 1:nEn){
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
    DOCout<-apply(Y[8,1,,],MARGIN = 1,FUN=mean)
    DOCout<-DOCout/data2$epiVol*12
    TrueDOC<-true[8,1,,1]/data2$epiVol*12 
    DOCrmse<-sqrt(mean((DOCout-TrueDOC)^2))
    DOCbias<-mean(DOCout-TrueDOC,na.rm=T)^2
    DOCprecision<-summary(lm(DOCout~TrueDOC))$r.squared
    docOut<-rbind(docOut,DOCrmse)
    docOutBias<-rbind(docOutBias,DOCbias)
    docOutPrecision<-rbind(docOutPrecision,DOCprecision)
    obsDOC<-z$y[2,1,,1]/data2$epiVol*12
    obsDOCrmse<-sqrt(mean((obsDOC-TrueDOC)^2,na.rm = T))
    obsDOCBias<-mean(obsDOC-TrueDOC,na.rm = T)^2
    obsDOCPrecision<-summary(lm(obsDOC~TrueDOC))$r.squared
    obsDOCout<-rbind(obsDOCout,obsDOCrmse)
    obsDOCoutBias<-rbind(obsDOCoutBias,obsDOCBias)
    obsDOCoutPrecision<-rbind(obsDOCoutPrecision,obsDOCPrecision)
    mdfobsDOCrmse<-sqrt(mean((DOCout-obsDOC)^2,na.rm = T))
    mdfobsDOCbias<-mean(DOCout-obsDOC,na.rm = T)^2
    mdfobsDOCprecision<-summary(lm(DOCout~obsDOC))$r.squared
    mdfobsDOCout<-rbind(mdfobsDOCout,mdfobsDOCrmse)
    mdfobsDOCoutBias<-rbind(mdfobsDOCoutBias,mdfobsDOCbias)
    mdfobsDOCoutPrecision<-rbind(mdfobsDOCoutPrecision,mdfobsDOCprecision)
    DOCcv<-mean(apply(Y[8,1,,],MARGIN = 1,FUN=sd)/apply(Y[8,1,,],MARGIN = 1,FUN=mean))
    mdfDOCcv<-rbind(mdfDOCcv,DOCcv)
    
    DICout<-apply(Y[5,1,,],MARGIN = 1,FUN=mean)
    DICout<-DICout/data2$epiVol*12
    TrueDIC<-true[5,1,,1]/data2$epiVol*12 
    DICrmse<-sqrt(mean((DICout-TrueDIC)^2))
    DICbias<-mean(DICout-TrueDIC,na.rm=T)^2
    DICprecision<-summary(lm(DICout~TrueDIC))$r.squared
    dicOut<-rbind(dicOut,DICrmse)
    dicOutBias<-rbind(dicOutBias,DICbias)
    dicOutPrecision<-rbind(dicOutPrecision,DICprecision)
    obsDIC<-z$y[1,1,,1]/data2$epiVol*12
    obsDICrmse<-sqrt(mean((obsDIC-TrueDIC)^2,na.rm = T))
    obsDICBias<-mean(obsDIC-TrueDIC,na.rm = T)^2
    obsDICPrecision<-summary(lm(obsDIC~TrueDIC))$r.squared
    obsDICout<-rbind(obsDICout,obsDICrmse)
    obsDICoutBias<-rbind(obsDICoutBias,obsDICBias)
    obsDICoutPrecision<-rbind(obsDICoutPrecision,obsDICPrecision)
    mdfobsDICrmse<-sqrt(mean((DICout-obsDIC)^2,na.rm = T))
    mdfobsDICbias<-mean(DICout-obsDIC,na.rm = T)^2
    mdfobsDICprecision<-summary(lm(DICout~obsDIC))$r.squared
    mdfobsDICout<-rbind(mdfobsDICout,mdfobsDICrmse)
    mdfobsDICoutBias<-rbind(mdfobsDICoutBias,mdfobsDICbias)
    mdfobsDICoutPrecision<-rbind(mdfobsDICoutPrecision,mdfobsDICprecision)
    DICcv<-mean(apply(Y[5,1,,],MARGIN = 1,FUN=sd)/apply(Y[5,1,,],MARGIN = 1,FUN=mean))
    mdfDICcv<-rbind(mdfDICcv,DICcv)
  }
  DOCrmse<-mean(docOut,na.rm = T)
  DOCbias<-mean(docOutBias,na.rm=T)
  DOCprecision<-mean(docOutPrecision,na.rm=T)
  obsDOCrmse<-mean(obsDOCout,na.rm=T)
  obsDOCbias<-mean(obsDOCoutBias,na.rm = T)
  obsDOCprecision<-mean(obsDOCoutPrecision,na.rm = T)
  mdfobsDOCrmse<-mean(mdfobsDOCout,na.rm=T)
  mdfobsDOCbias<-mean(mdfobsDOCoutBias,na.rm = T)
  mdfobsDOCprecision<-mean(mdfobsDOCoutPrecision,na.rm = T)
  mdfDOCcv<-mean(mdfDOCcv,na.rm = T)
  
  DICrmse<-mean(dicOut,na.rm = T)
  DICbias<-mean(dicOutBias,na.rm=T)
  DICprecision<-mean(dicOutPrecision,na.rm=T)
  obsDICrmse<-mean(obsDICout,na.rm=T)
  obsDICbias<-mean(obsDICoutBias,na.rm = T)
  obsDICprecision<-mean(obsDICoutPrecision,na.rm = T)
  mdfobsDICrmse<-mean(mdfobsDICout,na.rm=T)
  mdfobsDICbias<-mean(mdfobsDICoutBias,na.rm = T)
  mdfobsDICprecision<-mean(mdfobsDICoutPrecision,na.rm = T)
  mdfDICcv<-mean(mdfDICcv,na.rm = T)
  
  return(list("DOCrmse"=DOCrmse,"DOCbias"=DOCbias,'DOCprecision'=DOCprecision,
              "DICrmse"=DICrmse,'DICbias'=DICbias,'DICprecision'=DICprecision,
              "obsDOCrmse"=obsDOCrmse,'obsDOCbias'=obsDOCbias,'obsDOCprecision'=obsDOCprecision,
              "obsDICrmse"=obsDICrmse,'obsDICbias'=obsDICbias,'obsDICprecision'=obsDICprecision,
              'mdfobsDOCrmse'=mdfobsDOCrmse,'mdfobsDOCbias'=mdfobsDOCbias,'mdfobsDOCprecision'=mdfobsDOCprecision,
              'mdfobsDICrmse'=mdfobsDICrmse,'mdfobsDICbias'=mdfobsDICbias,'mdfobsDICprecision'=mdfobsDICprecision,
              'mdfDOCcv'=mdfDOCcv,'mdfDICcv'=mdfDICcv))
}

oneDOCpool_noMM_Temp<-function(obs,freq,reps,rVec,fracInVec,n){
  docOut<-c()
  docOutBias<-c()
  docOutPrecision<-c()
  dicOut<-c()
  dicOutBias<-c()
  dicOutPrecision<-c()
  obsDICout<-c()
  obsDICoutBias<-c()
  obsDICoutPrecision<-c()
  obsDOCout<-c()
  obsDOCoutBias<-c()
  obsDOCoutPrecision<-c()
  mdfobsDOCout<-c()
  mdfobsDOCoutBias<-c()
  mdfobsDOCoutPrecision<-c()
  mdfobsDICout<-c()
  mdfobsDICoutBias<-c()
  mdfobsDICoutPrecision<-c()
  mdfDOCcv<-c()
  mdfDICcv<-c()
  for(m in 1:n){
  rVec<-rVec
  fracInVec<-fracInVec

  # initial B transition matrix for each ensemble  
  B<-array(t(c(rep(NA,nEn),rep(NA,nEn),rep(NA,nEn),rep(NA,nEn))),dim=c(2,2,nStep,nEn)) # array of transition matrix B[a,b,c,d]; 
  # where [a,b] is transition matrix DIC & DOC, c=timeStep, and d=ensemble member 
  
  # initial parameters of B at timestep 1 
  for(i in 1:nEn){
    B[,,1,i]<-matrix(c(1-data2$QoutInt[1]/data2$epiVol[1]-data2$kCO2[1]/data2$thermo.depth[1]+
                         (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                       rFunc(rVec[i],data2$wtr[1]),0,1-rFunc(rVec[i],data2$wtr[1])-data2$QoutInt[1]/data2$epiVol[1]+
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
  h<-array(0,dim=c(2,4,nStep))
  for(i in 1:nStep){
    h[1,3,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic 
    h[2,4,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) #doc
  }
  
  
  P <- array(0,dim=c(2,2,nStep))
  
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
  
  pars<-array(rep(NA,nEn),dim=c(2,1,nStep,nEn)) # parameters: r20
  pars[1,1,1,]<-rVec
  pars[2,1,1,]<-fracInVec
  
  # set up a list for all matrices 
  z=list(B=B,y=y,X=X,C=C,ut=ut,pars=pars)
  
  # i is ensemble member; t is timestep 
  i=1
  t=2
  
  # set up Y vector for which we concatonate parameters, states, and observed data 
  Y<-array(NA,c(4,1,nStep,nEn))
  Y[1,1,1,]<-rVec # r20 parameter 
  Y[2,1,1,]<-fracInVec # fraction that goes into epi 
  Y[3,1,1,]<-z$X[1,1,1,] # DIC state  
  Y[4,1,1,]<-z$X[2,1,1,] # DOC state 
  
  #Iterate through time
  for(t in 2:nStep){
    for(i in 1:nEn){
      # Forecasting; need to update parameters, 
      z$pars[1:2,1,t,i]<-Y[1:2,1,t-1,i] # r20
      z$X[1:2,1,t-1,i]<-Y[3:4,1,t-1,i] # updating state variables from Y 
      
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
      Y[1:2,1,t,i]<-Y[1:2,1,t-1,i] #r20,same as previous timestep
      Y[3:4,1,t,i]<-z$X[1:2,1,t,i] #forecasted states 
      
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
  
  DOCout<-apply(Y[4,1,,],MARGIN = 1,FUN=mean)
  DOCout<-DOCout/data2$epiVol*12
  TrueDOC<-true[8,1,,1]/data2$epiVol*12 
  DOCrmse<-sqrt(mean((DOCout-TrueDOC)^2))
  DOCbias<-mean(DOCout-TrueDOC,na.rm=T)^2
  DOCprecision<-summary(lm(DOCout~TrueDOC))$r.squared
  docOut<-rbind(docOut,DOCrmse)
  docOutBias<-rbind(docOutBias,DOCbias)
  docOutPrecision<-rbind(docOutPrecision,DOCprecision)
  obsDOC<-z$y[2,1,,1]/data2$epiVol*12
  obsDOCrmse<-sqrt(mean((obsDOC-TrueDOC)^2,na.rm = T))
  obsDOCBias<-mean(obsDOC-TrueDOC,na.rm = T)^2
  obsDOCPrecision<-summary(lm(obsDOC~TrueDOC))$r.squared
  obsDOCout<-rbind(obsDOCout,obsDOCrmse)
  obsDOCoutBias<-rbind(obsDOCoutBias,obsDOCBias)
  obsDOCoutPrecision<-rbind(obsDOCoutPrecision,obsDOCPrecision)
  mdfobsDOCrmse<-sqrt(mean((DOCout-obsDOC)^2,na.rm = T))
  mdfobsDOCbias<-mean(DOCout-obsDOC,na.rm = T)^2
  mdfobsDOCprecision<-summary(lm(DOCout~obsDOC))$r.squared
  mdfobsDOCout<-rbind(mdfobsDOCout,mdfobsDOCrmse)
  mdfobsDOCoutBias<-rbind(mdfobsDOCoutBias,mdfobsDOCbias)
  mdfobsDOCoutPrecision<-rbind(mdfobsDOCoutPrecision,mdfobsDOCprecision)
  DOCcv<-mean(apply(Y[4,1,,],MARGIN = 1,FUN=sd)/apply(Y[4,1,,],MARGIN = 1,FUN=mean))
  mdfDOCcv<-rbind(mdfDOCcv,DOCcv)
  
  DICout<-apply(Y[3,1,,],MARGIN = 1,FUN=mean)
  DICout<-DICout/data2$epiVol*12
  TrueDIC<-true[5,1,,1]/data2$epiVol*12 
  DICrmse<-sqrt(mean((DICout-TrueDIC)^2))
  DICbias<-mean(DICout-TrueDIC,na.rm=T)^2
  DICprecision<-summary(lm(DICout~TrueDIC))$r.squared
  dicOut<-rbind(dicOut,DICrmse)
  dicOutBias<-rbind(dicOutBias,DICbias)
  dicOutPrecision<-rbind(dicOutPrecision,DICprecision)
  obsDIC<-z$y[1,1,,1]/data2$epiVol*12
  obsDICrmse<-sqrt(mean((obsDIC-TrueDIC)^2,na.rm = T))
  obsDICBias<-mean(obsDIC-TrueDIC,na.rm = T)^2
  obsDICPrecision<-summary(lm(obsDIC~TrueDIC))$r.squared
  obsDICout<-rbind(obsDICout,obsDICrmse)
  obsDICoutBias<-rbind(obsDICoutBias,obsDICBias)
  obsDICoutPrecision<-rbind(obsDICoutPrecision,obsDICPrecision)
  mdfobsDICrmse<-sqrt(mean((DICout-obsDIC)^2,na.rm = T))
  mdfobsDICbias<-mean(DICout-obsDIC,na.rm = T)^2
  mdfobsDICprecision<-summary(lm(DICout~obsDIC))$r.squared
  mdfobsDICout<-rbind(mdfobsDICout,mdfobsDICrmse)
  mdfobsDICoutBias<-rbind(mdfobsDICoutBias,mdfobsDICbias)
  mdfobsDICoutPrecision<-rbind(mdfobsDICoutPrecision,mdfobsDICprecision)
  DICcv<-mean(apply(Y[3,1,,],MARGIN = 1,FUN=sd)/apply(Y[3,1,,],MARGIN = 1,FUN=mean))
  mdfDICcv<-rbind(mdfDICcv,DICcv)
  
  }
  DOCrmse<-mean(docOut,na.rm = T)
  DOCbias<-mean(docOutBias,na.rm=T)
  DOCprecision<-mean(docOutPrecision,na.rm=T)
  obsDOCrmse<-mean(obsDOCout,na.rm=T)
  obsDOCbias<-mean(obsDOCoutBias,na.rm = T)
  obsDOCprecision<-mean(obsDOCoutPrecision,na.rm = T)
  mdfobsDOCrmse<-mean(mdfobsDOCout,na.rm=T)
  mdfobsDOCbias<-mean(mdfobsDOCoutBias,na.rm = T)
  mdfobsDOCprecision<-mean(mdfobsDOCoutPrecision,na.rm = T)
  mdfDOCcv<-mean(mdfDOCcv,na.rm = T)
  
  DICrmse<-mean(dicOut,na.rm = T)
  DICbias<-mean(dicOutBias,na.rm=T)
  DICprecision<-mean(dicOutPrecision,na.rm=T)
  obsDICrmse<-mean(obsDICout,na.rm=T)
  obsDICbias<-mean(obsDICoutBias,na.rm = T)
  obsDICprecision<-mean(obsDICoutPrecision,na.rm = T)
  mdfobsDICrmse<-mean(mdfobsDICout,na.rm=T)
  mdfobsDICbias<-mean(mdfobsDICoutBias,na.rm = T)
  mdfobsDICprecision<-mean(mdfobsDICoutPrecision,na.rm = T)
  mdfDICcv<-mean(mdfDICcv,na.rm = T)
  
  return(list("DOCrmse"=DOCrmse,"DOCbias"=DOCbias,'DOCprecision'=DOCprecision,
              "DICrmse"=DICrmse,'DICbias'=DICbias,'DICprecision'=DICprecision,
              "obsDOCrmse"=obsDOCrmse,'obsDOCbias'=obsDOCbias,'obsDOCprecision'=obsDOCprecision,
              "obsDICrmse"=obsDICrmse,'obsDICbias'=obsDICbias,'obsDICprecision'=obsDICprecision,
              'mdfobsDOCrmse'=mdfobsDOCrmse,'mdfobsDOCbias'=mdfobsDOCbias,'mdfobsDOCprecision'=mdfobsDOCprecision,
              'mdfobsDICrmse'=mdfobsDICrmse,'mdfobsDICbias'=mdfobsDICbias,'mdfobsDICprecision'=mdfobsDICprecision,
              'mdfDOCcv'=mdfDOCcv,'mdfDICcv'=mdfDICcv))
}

oneDOCpool_noMM_noTemp<-function(obs,freq,reps,rVec,fracInVec,n){
  docOut<-c()
  docOutBias<-c()
  docOutPrecision<-c()
  dicOut<-c()
  dicOutBias<-c()
  dicOutPrecision<-c()
  obsDICout<-c()
  obsDICoutBias<-c()
  obsDICoutPrecision<-c()
  obsDOCout<-c()
  obsDOCoutBias<-c()
  obsDOCoutPrecision<-c()
  mdfobsDOCout<-c()
  mdfobsDOCoutBias<-c()
  mdfobsDOCoutPrecision<-c()
  mdfobsDICout<-c()
  mdfobsDICoutBias<-c()
  mdfobsDICoutPrecision<-c()
  mdfDOCcv<-c()
  mdfDICcv<-c()
  for(m in 1:n){
    rVec<-rVec
    fracInVec<-fracInVec
    
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
    h<-array(0,dim=c(2,4,nStep))
    for(i in 1:nStep){
      h[1,3,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic 
      h[2,4,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) #doc
    }
    
    
    P <- array(0,dim=c(2,2,nStep))
    
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
    
    pars<-array(rep(NA,nEn),dim=c(2,1,nStep,nEn)) # parameters: r20
    pars[1,1,1,]<-rVec
    pars[2,1,1,]<-fracInVec
    
    # set up a list for all matrices 
    z=list(B=B,y=y,X=X,C=C,ut=ut,pars=pars)
    
    # i is ensemble member; t is timestep 
    i=1
    t=2
    
    # set up Y vector for which we concatonate parameters, states, and observed data 
    Y<-array(NA,c(4,1,nStep,nEn))
    Y[1,1,1,]<-rVec # r20 parameter 
    Y[2,1,1,]<-fracInVec # fraction that goes into epi 
    Y[3,1,1,]<-z$X[1,1,1,] # DIC state  
    Y[4,1,1,]<-z$X[2,1,1,] # DOC state 
    
    #Iterate through time
    for(t in 2:nStep){
      for(i in 1:nEn){
        # Forecasting; need to update parameters, 
        z$pars[1:2,1,t,i]<-Y[1:2,1,t-1,i] # r20
        z$X[1:2,1,t-1,i]<-Y[3:4,1,t-1,i] # updating state variables from Y 
        
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
        Y[1:2,1,t,i]<-Y[1:2,1,t-1,i] #r20,same as previous timestep
        Y[3:4,1,t,i]<-z$X[1:2,1,t,i] #forecasted states 
        
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
    
    DOCout<-apply(Y[4,1,,],MARGIN = 1,FUN=mean)
    DOCout<-DOCout/data2$epiVol*12
    TrueDOC<-true[8,1,,1]/data2$epiVol*12 
    DOCrmse<-sqrt(mean((DOCout-TrueDOC)^2))
    DOCbias<-mean(DOCout-TrueDOC,na.rm=T)^2
    DOCprecision<-summary(lm(DOCout~TrueDOC))$r.squared
    docOut<-rbind(docOut,DOCrmse)
    docOutBias<-rbind(docOutBias,DOCbias)
    docOutPrecision<-rbind(docOutPrecision,DOCprecision)
    obsDOC<-z$y[2,1,,1]/data2$epiVol*12
    obsDOCrmse<-sqrt(mean((obsDOC-TrueDOC)^2,na.rm = T))
    obsDOCBias<-mean(obsDOC-TrueDOC,na.rm = T)^2
    obsDOCPrecision<-summary(lm(obsDOC~TrueDOC))$r.squared
    obsDOCout<-rbind(obsDOCout,obsDOCrmse)
    obsDOCoutBias<-rbind(obsDOCoutBias,obsDOCBias)
    obsDOCoutPrecision<-rbind(obsDOCoutPrecision,obsDOCPrecision)
    mdfobsDOCrmse<-sqrt(mean((DOCout-obsDOC)^2,na.rm = T))
    mdfobsDOCbias<-mean(DOCout-obsDOC,na.rm = T)^2
    mdfobsDOCprecision<-summary(lm(DOCout~obsDOC))$r.squared
    mdfobsDOCout<-rbind(mdfobsDOCout,mdfobsDOCrmse)
    mdfobsDOCoutBias<-rbind(mdfobsDOCoutBias,mdfobsDOCbias)
    mdfobsDOCoutPrecision<-rbind(mdfobsDOCoutPrecision,mdfobsDOCprecision)
    DOCcv<-mean(apply(Y[4,1,,],MARGIN = 1,FUN=sd)/apply(Y[4,1,,],MARGIN = 1,FUN=mean))
    mdfDOCcv<-rbind(mdfDOCcv,DOCcv)
    
    DICout<-apply(Y[3,1,,],MARGIN = 1,FUN=mean)
    DICout<-DICout/data2$epiVol*12
    TrueDIC<-true[5,1,,1]/data2$epiVol*12 
    DICrmse<-sqrt(mean((DICout-TrueDIC)^2))
    DICbias<-mean(DICout-TrueDIC,na.rm=T)^2
    DICprecision<-summary(lm(DICout~TrueDIC))$r.squared
    dicOut<-rbind(dicOut,DICrmse)
    dicOutBias<-rbind(dicOutBias,DICbias)
    dicOutPrecision<-rbind(dicOutPrecision,DICprecision)
    obsDIC<-z$y[1,1,,1]/data2$epiVol*12
    obsDICrmse<-sqrt(mean((obsDIC-TrueDIC)^2,na.rm = T))
    obsDICBias<-mean(obsDIC-TrueDIC,na.rm = T)^2
    obsDICPrecision<-summary(lm(obsDIC~TrueDIC))$r.squared
    obsDICout<-rbind(obsDICout,obsDICrmse)
    obsDICoutBias<-rbind(obsDICoutBias,obsDICBias)
    obsDICoutPrecision<-rbind(obsDICoutPrecision,obsDICPrecision)
    mdfobsDICrmse<-sqrt(mean((DICout-obsDIC)^2,na.rm = T))
    mdfobsDICbias<-mean(DICout-obsDIC,na.rm = T)^2
    mdfobsDICprecision<-summary(lm(DICout~obsDIC))$r.squared
    mdfobsDICout<-rbind(mdfobsDICout,mdfobsDICrmse)
    mdfobsDICoutBias<-rbind(mdfobsDICoutBias,mdfobsDICbias)
    mdfobsDICoutPrecision<-rbind(mdfobsDICoutPrecision,mdfobsDICprecision)
    DICcv<-mean(apply(Y[3,1,,],MARGIN = 1,FUN=sd)/apply(Y[3,1,,],MARGIN = 1,FUN=mean))
    mdfDICcv<-rbind(mdfDICcv,DICcv)
  }
  DOCrmse<-mean(docOut,na.rm = T)
  DOCbias<-mean(docOutBias,na.rm=T)
  DOCprecision<-mean(docOutPrecision,na.rm=T)
  obsDOCrmse<-mean(obsDOCout,na.rm=T)
  obsDOCbias<-mean(obsDOCoutBias,na.rm = T)
  obsDOCprecision<-mean(obsDOCoutPrecision,na.rm = T)
  mdfobsDOCrmse<-mean(mdfobsDOCout,na.rm=T)
  mdfobsDOCbias<-mean(mdfobsDOCoutBias,na.rm = T)
  mdfobsDOCprecision<-mean(mdfobsDOCoutPrecision,na.rm = T)
  mdfDOCcv<-mean(mdfDOCcv,na.rm = T)
  
  DICrmse<-mean(dicOut,na.rm = T)
  DICbias<-mean(dicOutBias,na.rm=T)
  DICprecision<-mean(dicOutPrecision,na.rm=T)
  obsDICrmse<-mean(obsDICout,na.rm=T)
  obsDICbias<-mean(obsDICoutBias,na.rm = T)
  obsDICprecision<-mean(obsDICoutPrecision,na.rm = T)
  mdfobsDICrmse<-mean(mdfobsDICout,na.rm=T)
  mdfobsDICbias<-mean(mdfobsDICoutBias,na.rm = T)
  mdfobsDICprecision<-mean(mdfobsDICoutPrecision,na.rm = T)
  mdfDICcv<-mean(mdfDICcv,na.rm = T)
  
  return(list("DOCrmse"=DOCrmse,"DOCbias"=DOCbias,'DOCprecision'=DOCprecision,
              "DICrmse"=DICrmse,'DICbias'=DICbias,'DICprecision'=DICprecision,
              "obsDOCrmse"=obsDOCrmse,'obsDOCbias'=obsDOCbias,'obsDOCprecision'=obsDOCprecision,
              "obsDICrmse"=obsDICrmse,'obsDICbias'=obsDICbias,'obsDICprecision'=obsDICprecision,
              'mdfobsDOCrmse'=mdfobsDOCrmse,'mdfobsDOCbias'=mdfobsDOCbias,'mdfobsDOCprecision'=mdfobsDOCprecision,
              'mdfobsDICrmse'=mdfobsDICrmse,'mdfobsDICbias'=mdfobsDICbias,'mdfobsDICprecision'=mdfobsDICprecision,
              'mdfDOCcv'=mdfDOCcv,'mdfDICcv'=mdfDICcv))
}

oneDOCpool_MM_noTemp<-function(obs,freq,reps,rVec,fracInVec,n){
  docOut<-c()
  docOutBias<-c()
  docOutPrecision<-c()
  dicOut<-c()
  dicOutBias<-c()
  dicOutPrecision<-c()
  obsDICout<-c()
  obsDICoutBias<-c()
  obsDICoutPrecision<-c()
  obsDOCout<-c()
  obsDOCoutBias<-c()
  obsDOCoutPrecision<-c()
  mdfobsDOCout<-c()
  mdfobsDOCoutBias<-c()
  mdfobsDOCoutPrecision<-c()
  mdfobsDICout<-c()
  mdfobsDICoutBias<-c()
  mdfobsDICoutPrecision<-c()
  mdfDOCcv<-c()
  mdfDICcv<-c()
  for(m in 1:n){
    rVec<-rVec
    fracInVec<-fracInVec
    
    # initial B transition matrix for each ensemble  
    B<-array(t(c(rep(NA,nEn),rep(NA,nEn),rep(NA,nEn),rep(NA,nEn))),dim=c(2,2,nStep,nEn)) # array of transition matrix B[a,b,c,d]; 
    # where [a,b] is transition matrix DIC & DOC, c=timeStep, and d=ensemble member 
    
    # initial parameters of B at timestep 1 
    for(i in 1:nEn){
      B[,,1,i]<-matrix(c(1-data2$QoutInt[1]/data2$epiVol[1]-data2$kCO2[1]/data2$thermo.depth[1]+
                           (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                         mm(rVec[i],docVec[1],halfSat,data2$epiVol[1]),0,1-mm(rVec[i],docVec[1],halfSat,data2$epiVol[1])-data2$QoutInt[1]/data2$epiVol[1]+
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
    h<-array(0,dim=c(2,4,nStep))
    for(i in 1:nStep){
      h[1,3,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic 
      h[2,4,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) #doc
    }
    
    
    P <- array(0,dim=c(2,2,nStep))
    
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
    
    pars<-array(rep(NA,nEn),dim=c(2,1,nStep,nEn)) # parameters: r20
    pars[1,1,1,]<-rVec
    pars[2,1,1,]<-fracInVec
    
    # set up a list for all matrices 
    z=list(B=B,y=y,X=X,C=C,ut=ut,pars=pars)
    
    # i is ensemble member; t is timestep 
    i=1
    t=2
    
    # set up Y vector for which we concatonate parameters, states, and observed data 
    Y<-array(NA,c(4,1,nStep,nEn))
    Y[1,1,1,]<-rVec # r20 parameter 
    Y[2,1,1,]<-fracInVec # fraction that goes into epi 
    Y[3,1,1,]<-z$X[1,1,1,] # DIC state  
    Y[4,1,1,]<-z$X[2,1,1,] # DOC state 
    
    #Iterate through time
    for(t in 2:nStep){
      for(i in 1:nEn){
        # Forecasting; need to update parameters, 
        z$pars[1:2,1,t,i]<-Y[1:2,1,t-1,i] # r20
        z$X[1:2,1,t-1,i]<-Y[3:4,1,t-1,i] # updating state variables from Y 
        
        #Predictions
        z$X[,,t,i]<-z$B[,,t-1,i]%*%z$X[,,t-1,i] + z$C[,,t-1,i]%*%z$ut[,,t-1,i] # forecasting state variable predictions 
        z$B[,,t,i]<-matrix(c(1-data2$QoutInt[t]/data2$epiVol[t]-data2$kCO2[t]/data2$thermo.depth[t]+
                               (data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[2,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],
                             mm(z$pars[1,1,t,i],z$X[2,1,t,i],halfSat,data2$epiVol[t]),0,1-mm(z$pars[1,1,t,i],z$X[2,1,t,i],halfSat,data2$epiVol[t])-data2$QoutInt[t]/data2$epiVol[t]+
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
        Y[1:2,1,t,i]<-Y[1:2,1,t-1,i] #r20,same as previous timestep
        Y[3:4,1,t,i]<-z$X[1:2,1,t,i] #forecasted states 
        
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
    
    DOCout<-apply(Y[4,1,,],MARGIN = 1,FUN=mean)
    DOCout<-DOCout/data2$epiVol*12
    TrueDOC<-true[8,1,,1]/data2$epiVol*12 
    DOCrmse<-sqrt(mean((DOCout-TrueDOC)^2))
    DOCbias<-mean(DOCout-TrueDOC,na.rm=T)^2
    DOCprecision<-summary(lm(DOCout~TrueDOC))$r.squared
    docOut<-rbind(docOut,DOCrmse)
    docOutBias<-rbind(docOutBias,DOCbias)
    docOutPrecision<-rbind(docOutPrecision,DOCprecision)
    obsDOC<-z$y[2,1,,1]/data2$epiVol*12
    obsDOCrmse<-sqrt(mean((obsDOC-TrueDOC)^2,na.rm = T))
    obsDOCBias<-mean(obsDOC-TrueDOC,na.rm = T)^2
    obsDOCPrecision<-summary(lm(obsDOC~TrueDOC))$r.squared
    obsDOCout<-rbind(obsDOCout,obsDOCrmse)
    obsDOCoutBias<-rbind(obsDOCoutBias,obsDOCBias)
    obsDOCoutPrecision<-rbind(obsDOCoutPrecision,obsDOCPrecision)
    mdfobsDOCrmse<-sqrt(mean((DOCout-obsDOC)^2,na.rm = T))
    mdfobsDOCbias<-mean(DOCout-obsDOC,na.rm = T)^2
    mdfobsDOCprecision<-summary(lm(DOCout~obsDOC))$r.squared
    mdfobsDOCout<-rbind(mdfobsDOCout,mdfobsDOCrmse)
    mdfobsDOCoutBias<-rbind(mdfobsDOCoutBias,mdfobsDOCbias)
    mdfobsDOCoutPrecision<-rbind(mdfobsDOCoutPrecision,mdfobsDOCprecision)
    DOCcv<-mean(apply(Y[4,1,,],MARGIN = 1,FUN=sd)/apply(Y[4,1,,],MARGIN = 1,FUN=mean))
    mdfDOCcv<-rbind(mdfDOCcv,DOCcv)
    
    DICout<-apply(Y[3,1,,],MARGIN = 1,FUN=mean)
    DICout<-DICout/data2$epiVol*12
    TrueDIC<-true[5,1,,1]/data2$epiVol*12 
    DICrmse<-sqrt(mean((DICout-TrueDIC)^2))
    DICbias<-mean(DICout-TrueDIC,na.rm=T)^2
    DICprecision<-summary(lm(DICout~TrueDIC))$r.squared
    dicOut<-rbind(dicOut,DICrmse)
    dicOutBias<-rbind(dicOutBias,DICbias)
    dicOutPrecision<-rbind(dicOutPrecision,DICprecision)
    obsDIC<-z$y[1,1,,1]/data2$epiVol*12
    obsDICrmse<-sqrt(mean((obsDIC-TrueDIC)^2,na.rm = T))
    obsDICBias<-mean(obsDIC-TrueDIC,na.rm = T)^2
    obsDICPrecision<-summary(lm(obsDIC~TrueDIC))$r.squared
    obsDICout<-rbind(obsDICout,obsDICrmse)
    obsDICoutBias<-rbind(obsDICoutBias,obsDICBias)
    obsDICoutPrecision<-rbind(obsDICoutPrecision,obsDICPrecision)
    mdfobsDICrmse<-sqrt(mean((DICout-obsDIC)^2,na.rm = T))
    mdfobsDICbias<-mean(DICout-obsDIC,na.rm = T)^2
    mdfobsDICprecision<-summary(lm(DICout~obsDIC))$r.squared
    mdfobsDICout<-rbind(mdfobsDICout,mdfobsDICrmse)
    mdfobsDICoutBias<-rbind(mdfobsDICoutBias,mdfobsDICbias)
    mdfobsDICoutPrecision<-rbind(mdfobsDICoutPrecision,mdfobsDICprecision)
    DICcv<-mean(apply(Y[3,1,,],MARGIN = 1,FUN=sd)/apply(Y[3,1,,],MARGIN = 1,FUN=mean))
    mdfDICcv<-rbind(mdfDICcv,DICcv)
  }
  DOCrmse<-mean(docOut,na.rm = T)
  DOCbias<-mean(docOutBias,na.rm=T)
  DOCprecision<-mean(docOutPrecision,na.rm=T)
  obsDOCrmse<-mean(obsDOCout,na.rm=T)
  obsDOCbias<-mean(obsDOCoutBias,na.rm = T)
  obsDOCprecision<-mean(obsDOCoutPrecision,na.rm = T)
  mdfobsDOCrmse<-mean(mdfobsDOCout,na.rm=T)
  mdfobsDOCbias<-mean(mdfobsDOCoutBias,na.rm = T)
  mdfobsDOCprecision<-mean(mdfobsDOCoutPrecision,na.rm = T)
  mdfDOCcv<-mean(mdfDOCcv,na.rm = T)
  
  DICrmse<-mean(dicOut,na.rm = T)
  DICbias<-mean(dicOutBias,na.rm=T)
  DICprecision<-mean(dicOutPrecision,na.rm=T)
  obsDICrmse<-mean(obsDICout,na.rm=T)
  obsDICbias<-mean(obsDICoutBias,na.rm = T)
  obsDICprecision<-mean(obsDICoutPrecision,na.rm = T)
  mdfobsDICrmse<-mean(mdfobsDICout,na.rm=T)
  mdfobsDICbias<-mean(mdfobsDICoutBias,na.rm = T)
  mdfobsDICprecision<-mean(mdfobsDICoutPrecision,na.rm = T)
  mdfDICcv<-mean(mdfDICcv,na.rm = T)
  
  return(list("DOCrmse"=DOCrmse,"DOCbias"=DOCbias,'DOCprecision'=DOCprecision,
              "DICrmse"=DICrmse,'DICbias'=DICbias,'DICprecision'=DICprecision,
              "obsDOCrmse"=obsDOCrmse,'obsDOCbias'=obsDOCbias,'obsDOCprecision'=obsDOCprecision,
              "obsDICrmse"=obsDICrmse,'obsDICbias'=obsDICbias,'obsDICprecision'=obsDICprecision,
              'mdfobsDOCrmse'=mdfobsDOCrmse,'mdfobsDOCbias'=mdfobsDOCbias,'mdfobsDOCprecision'=mdfobsDOCprecision,
              'mdfobsDICrmse'=mdfobsDICrmse,'mdfobsDICbias'=mdfobsDICbias,'mdfobsDICprecision'=mdfobsDICprecision,
              'mdfDOCcv'=mdfDOCcv,'mdfDICcv'=mdfDICcv))
}


oneDOCpool_MM_Temp<-function(obs,freq,reps,rVec,fracInVec,n){
  docOut<-c()
  docOutBias<-c()
  docOutPrecision<-c()
  dicOut<-c()
  dicOutBias<-c()
  dicOutPrecision<-c()
  obsDICout<-c()
  obsDICoutBias<-c()
  obsDICoutPrecision<-c()
  obsDOCout<-c()
  obsDOCoutBias<-c()
  obsDOCoutPrecision<-c()
  mdfobsDOCout<-c()
  mdfobsDOCoutBias<-c()
  mdfobsDOCoutPrecision<-c()
  mdfobsDICout<-c()
  mdfobsDICoutBias<-c()
  mdfobsDICoutPrecision<-c()
  mdfDOCcv<-c()
  mdfDICcv<-c()
  for(m in 1:n){
    rVec<-rVec
    fracInVec<-fracInVec
    
    # initial B transition matrix for each ensemble  
    B<-array(t(c(rep(NA,nEn),rep(NA,nEn),rep(NA,nEn),rep(NA,nEn))),dim=c(2,2,nStep,nEn)) # array of transition matrix B[a,b,c,d]; 
    # where [a,b] is transition matrix DIC & DOC, c=timeStep, and d=ensemble member 
    
    # initial parameters of B at timestep 1 
    for(i in 1:nEn){
      B[,,1,i]<-matrix(c(1-data2$QoutInt[1]/data2$epiVol[1]-data2$kCO2[1]/data2$thermo.depth[1]+
                           (data2$entrainVol[1]-data2$streamWaterdisch[1]*(1-splitFunc(data2$epiDens[1],data2$streamDens[1],fracInVec[i])))*data2$entrainHypo[1]/data2$epiVol[1],
                         mm(rFunc(rVec[i],data2$wtr[1]),docVec[1],halfSat,data2$epiVol[1]),0,1-mm(rFunc(rVec[i],data2$wtr[1]),docVec[1],halfSat,data2$epiVol[1])-data2$QoutInt[1]/data2$epiVol[1]+
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
    h<-array(0,dim=c(2,4,nStep))
    for(i in 1:nStep){
      h[1,3,i]<-ifelse(!is.na(y[1,1,i,1]),1,0) #dic 
      h[2,4,i]<-ifelse(!is.na(y[2,1,i,1]),1,0) #doc
    }
    
    
    P <- array(0,dim=c(2,2,nStep))
    
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
    
    pars<-array(rep(NA,nEn),dim=c(2,1,nStep,nEn)) # parameters: r20
    pars[1,1,1,]<-rVec
    pars[2,1,1,]<-fracInVec
    
    # set up a list for all matrices 
    z=list(B=B,y=y,X=X,C=C,ut=ut,pars=pars)
    
    # i is ensemble member; t is timestep 
    i=1
    t=2
    
    # set up Y vector for which we concatonate parameters, states, and observed data 
    Y<-array(NA,c(4,1,nStep,nEn))
    Y[1,1,1,]<-rVec # r20 parameter 
    Y[2,1,1,]<-fracInVec # fraction that goes into epi 
    Y[3,1,1,]<-z$X[1,1,1,] # DIC state  
    Y[4,1,1,]<-z$X[2,1,1,] # DOC state 
    
    #Iterate through time
    for(t in 2:nStep){
      for(i in 1:nEn){
        # Forecasting; need to update parameters, 
        z$pars[1:2,1,t,i]<-Y[1:2,1,t-1,i] # r20
        z$X[1:2,1,t-1,i]<-Y[3:4,1,t-1,i] # updating state variables from Y 
        
        #Predictions
        z$X[,,t,i]<-z$B[,,t-1,i]%*%z$X[,,t-1,i] + z$C[,,t-1,i]%*%z$ut[,,t-1,i] # forecasting state variable predictions 
        z$B[,,t,i]<-matrix(c(1-data2$QoutInt[t]/data2$epiVol[t]-data2$kCO2[t]/data2$thermo.depth[t]+
                               (data2$entrainVol[t]-data2$streamWaterdisch[t]*(1-splitFunc(data2$epiDens[t],data2$streamDens[t],z$pars[2,1,t,i])))*data2$entrainHypo[t]/data2$epiVol[t],
                             mm(rFunc(z$pars[1,1,t,i],data2$wtr[t]),z$X[2,1,t,i],halfSat,data2$epiVol[t]),0,1-mm(rFunc(z$pars[1,1,t,i],data2$wtr[t]),z$X[2,1,t,i],halfSat,data2$epiVol[t])-data2$QoutInt[t]/data2$epiVol[t]+
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
        Y[1:2,1,t,i]<-Y[1:2,1,t-1,i] #r20,same as previous timestep
        Y[3:4,1,t,i]<-z$X[1:2,1,t,i] #forecasted states 
        
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
    
    DOCout<-apply(Y[4,1,,],MARGIN = 1,FUN=mean)
    DOCout<-DOCout/data2$epiVol*12
    TrueDOC<-true[8,1,,1]/data2$epiVol*12 
    DOCrmse<-sqrt(mean((DOCout-TrueDOC)^2))
    DOCbias<-mean(DOCout-TrueDOC,na.rm=T)^2
    DOCprecision<-summary(lm(DOCout~TrueDOC))$r.squared
    docOut<-rbind(docOut,DOCrmse)
    docOutBias<-rbind(docOutBias,DOCbias)
    docOutPrecision<-rbind(docOutPrecision,DOCprecision)
    obsDOC<-z$y[2,1,,1]/data2$epiVol*12
    obsDOCrmse<-sqrt(mean((obsDOC-TrueDOC)^2,na.rm = T))
    obsDOCBias<-mean(obsDOC-TrueDOC,na.rm = T)^2
    obsDOCPrecision<-summary(lm(obsDOC~TrueDOC))$r.squared
    obsDOCout<-rbind(obsDOCout,obsDOCrmse)
    obsDOCoutBias<-rbind(obsDOCoutBias,obsDOCBias)
    obsDOCoutPrecision<-rbind(obsDOCoutPrecision,obsDOCPrecision)
    mdfobsDOCrmse<-sqrt(mean((DOCout-obsDOC)^2,na.rm = T))
    mdfobsDOCbias<-mean(DOCout-obsDOC,na.rm = T)^2
    mdfobsDOCprecision<-summary(lm(DOCout~obsDOC))$r.squared
    mdfobsDOCout<-rbind(mdfobsDOCout,mdfobsDOCrmse)
    mdfobsDOCoutBias<-rbind(mdfobsDOCoutBias,mdfobsDOCbias)
    mdfobsDOCoutPrecision<-rbind(mdfobsDOCoutPrecision,mdfobsDOCprecision)
    DOCcv<-mean(apply(Y[4,1,,],MARGIN = 1,FUN=sd)/apply(Y[4,1,,],MARGIN = 1,FUN=mean))
    mdfDOCcv<-rbind(mdfDOCcv,DOCcv)
    
    DICout<-apply(Y[3,1,,],MARGIN = 1,FUN=mean)
    DICout<-DICout/data2$epiVol*12
    TrueDIC<-true[5,1,,1]/data2$epiVol*12 
    DICrmse<-sqrt(mean((DICout-TrueDIC)^2))
    DICbias<-mean(DICout-TrueDIC,na.rm=T)^2
    DICprecision<-summary(lm(DICout~TrueDIC))$r.squared
    dicOut<-rbind(dicOut,DICrmse)
    dicOutBias<-rbind(dicOutBias,DICbias)
    dicOutPrecision<-rbind(dicOutPrecision,DICprecision)
    obsDIC<-z$y[1,1,,1]/data2$epiVol*12
    obsDICrmse<-sqrt(mean((obsDIC-TrueDIC)^2,na.rm = T))
    obsDICBias<-mean(obsDIC-TrueDIC,na.rm = T)^2
    obsDICPrecision<-summary(lm(obsDIC~TrueDIC))$r.squared
    obsDICout<-rbind(obsDICout,obsDICrmse)
    obsDICoutBias<-rbind(obsDICoutBias,obsDICBias)
    obsDICoutPrecision<-rbind(obsDICoutPrecision,obsDICPrecision)
    mdfobsDICrmse<-sqrt(mean((DICout-obsDIC)^2,na.rm = T))
    mdfobsDICbias<-mean(DICout-obsDIC,na.rm = T)^2
    mdfobsDICprecision<-summary(lm(DICout~obsDIC))$r.squared
    mdfobsDICout<-rbind(mdfobsDICout,mdfobsDICrmse)
    mdfobsDICoutBias<-rbind(mdfobsDICoutBias,mdfobsDICbias)
    mdfobsDICoutPrecision<-rbind(mdfobsDICoutPrecision,mdfobsDICprecision)
    DICcv<-mean(apply(Y[3,1,,],MARGIN = 1,FUN=sd)/apply(Y[3,1,,],MARGIN = 1,FUN=mean))
    mdfDICcv<-rbind(mdfDICcv,DICcv)
  }
  DOCrmse<-mean(docOut,na.rm = T)
  DOCbias<-mean(docOutBias,na.rm=T)
  DOCprecision<-mean(docOutPrecision,na.rm=T)
  obsDOCrmse<-mean(obsDOCout,na.rm=T)
  obsDOCbias<-mean(obsDOCoutBias,na.rm = T)
  obsDOCprecision<-mean(obsDOCoutPrecision,na.rm = T)
  mdfobsDOCrmse<-mean(mdfobsDOCout,na.rm=T)
  mdfobsDOCbias<-mean(mdfobsDOCoutBias,na.rm = T)
  mdfobsDOCprecision<-mean(mdfobsDOCoutPrecision,na.rm = T)
  mdfDOCcv<-mean(mdfDOCcv,na.rm = T)
  
  DICrmse<-mean(dicOut,na.rm = T)
  DICbias<-mean(dicOutBias,na.rm=T)
  DICprecision<-mean(dicOutPrecision,na.rm=T)
  obsDICrmse<-mean(obsDICout,na.rm=T)
  obsDICbias<-mean(obsDICoutBias,na.rm = T)
  obsDICprecision<-mean(obsDICoutPrecision,na.rm = T)
  mdfobsDICrmse<-mean(mdfobsDICout,na.rm=T)
  mdfobsDICbias<-mean(mdfobsDICoutBias,na.rm = T)
  mdfobsDICprecision<-mean(mdfobsDICoutPrecision,na.rm = T)
  mdfDICcv<-mean(mdfDICcv,na.rm = T)
  
  return(list("DOCrmse"=DOCrmse,"DOCbias"=DOCbias,'DOCprecision'=DOCprecision,
              "DICrmse"=DICrmse,'DICbias'=DICbias,'DICprecision'=DICprecision,
              "obsDOCrmse"=obsDOCrmse,'obsDOCbias'=obsDOCbias,'obsDOCprecision'=obsDOCprecision,
              "obsDICrmse"=obsDICrmse,'obsDICbias'=obsDICbias,'obsDICprecision'=obsDICprecision,
              'mdfobsDOCrmse'=mdfobsDOCrmse,'mdfobsDOCbias'=mdfobsDOCbias,'mdfobsDOCprecision'=mdfobsDOCprecision,
              'mdfobsDICrmse'=mdfobsDICrmse,'mdfobsDICbias'=mdfobsDICbias,'mdfobsDICprecision'=mdfobsDICprecision,
              'mdfDOCcv'=mdfDOCcv,'mdfDICcv'=mdfDICcv))
}

startingValues<-function(){
  # draws from priors to create ensemble 
  parGuess <- c(0.007) #r20; units: day-1
  min<-c(0.001)
  max<-c(0.02)
  
  rPDF<-abs(rnorm(n=nEn,mean = parGuess[1],sd = (max[1]-min[1])/5)) # max-min / 5 is rule of thumb sd; forcing positive for negative draws 
  # hist(rPDF)
  
  # setting up initial parameter values for all ensemble members 
  rVec<-matrix(rPDF) # each row is an ensemble member 
  
  parGuess <- c(0.004,0.3,0.3,0.1) #r20; units: day-1; fraction labile of loaded DOC; fraction inlet that goes into epi; turnover rate parameters set constant frac labile estimated 
  min<-c(0.001,0.06,0.005,0.3)
  max<-c(0.001,0.06,0.5,0.5)

  rPDF_slow<-abs(rnorm(n=nEn,mean = parGuess[1],sd = (max[1]-min[1])/5)) # max-min / 5 is rule of thumb sd; forcing positive for negative draws 
  rPDF_fast<-abs(rnorm(n=nEn,mean = parGuess[2],sd = (max[2]-min[2])/5))
  fracPDF<-abs(rnorm(n=nEn,mean=parGuess[3],sd=(max[3]-min[3])/5))
  fracInPDF<-abs(rnorm(n=nEn,mean=parGuess[4],sd=(max[4]-min[4])/5))
  # fracInPDF<-rbeta(n = nEn,shape1 = 4,shape2 = 1)
  
  # setting up initial parameter values for all ensemble members 
  rVec_slow<-matrix(rPDF_slow) # each row is an ensemble member 
  rVec_fast<-matrix(rPDF_fast)
  fracVec<-matrix(fracPDF)
  fracInVec<-matrix(fracInPDF)
  
  return(list('rVec'=rVec,'rVec_slow'=rVec_slow,'rVec_fast'=rVec_fast,'fracVec'=fracVec,'fracInVec'=fracInVec))
}



