
EnKF_2pools <- function(Y, z, i, t){
  #Iterate through time
  for(t in 2:nStep){
    for(i in 1:nEn){
      # Forecasting; need to update parameters, 
      z$pars[1:5,1,t,i]<-Y[1:5,1,t-1,i] # r20
      z$X[1:4,1,t-1,i]<-Y[6:9,1,t-1,i] # updating state variables from Y 
      
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
      Y[1:5,1,t,i]<-Y[1:5,1,t-1,i] #r20, r20_fast, fraction labile loaded DOC parameters same as previous timestep
      Y[6:9,1,t,i]<-z$X[1:4,1,t,i] #forecasted states 
      
    } # end forecast for each ensemble for timestep t
    
    #begin data assimilation if there are any observations and obs is an assimilation obs 
    if(any(!is.na(z$y[,,t,i]))==TRUE & (assim_number %% 2) == 0){ # update vector as long as there is one observation of state 
      assim_number <- assim_number+1
      assim_obs[t] <- 1
      # logging values to eliminate negatives 
      
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
      
      
      covar_inflat <- mean(Y[5,1,t,])
      
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
      
      # state convergence 
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
  return(list(Y = Y, assim_obs = assim_obs))
}
