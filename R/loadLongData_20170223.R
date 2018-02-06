# readinging in and formating data for MDF long 
# JAZ ; 2016-03-09 

library(LakeMetabolizer)
library(rLakeAnalyzer)
source('C:/Users/Jake/Desktop/R functions/DOY.r')

lake<-'EL'
site<-c('Inlet1')
years<-c(2013,2014) # changed to only inlude 2014 for now; inclduing 2013 on 2016-07-11
years<-c(2014,2015)
dateRange<-c('2012-05-23','2015-10-25')

# evap function for CR 
evapCalc<-function(airT,jDay,lat){
  #saturated Vapor Density 
  svd<-5.018+(.32321*airT)+(0.0081847*airT^2)+(0.00031243*airT^3) 
  
  #Daylength 
  degToRad<-2*pi/360
  radToDeg<-180/pi
  
  #day angle gamma (radians) 
  dayAngle<-2*pi*(jDay-1)/365 
  
  #declination of the sun 'delta' (radians)
  dec<-0.006918-0.399912*cos(dayAngle)+0.070257*sin(dayAngle)-
    0.006758*cos(2*dayAngle)+0.000907*sin(2*dayAngle)-0.002697*
    cos(3*dayAngle)+0.00148*sin(3*dayAngle)
  
  # sunrise hour angle 'omega' (degrees) 
  latRad<-lat*degToRad
  sunriseHourAngle<-acos(-tan(latRad)*tan(dec))*radToDeg 
  
  #sunrise and sunset times (decimal hours, relative to solar time) 
  sunrise<-12-sunriseHourAngle/15
  sunset<-12+sunriseHourAngle/15
  dayLength<-sunset-sunrise #day length in hours
  
  evap = 0.55*((dayLength/12)^2)*(svd/100)*25.4 #calculates evaporation for each jDay (units are mm/day)
  return(evap)
}

source('/Users/Jake/Desktop/R functions/dbTable.R')

bathy<-dbTable('lake_bathymetry',lakeID=lake)
lakeArea<-bathy$area_m2[as.numeric(bathy$depth_m)==0]
if(lake=='MO'){
  lakeVol<-as.numeric(bathy$volumeToBottom_m3[as.numeric(bathy$depth_m)==0])
}else{
  lakeVol<-as.numeric(bathy$volumeToBottom3D_m3[as.numeric(bathy$depth_m)==0])
}

staff<-dbTable('staff_gauges',lakeID=lake)
staff<-staff[staff$siteName=='WholeLake',]
staff$dateTimeSample<-strftime(strptime(staff$dateTimeSample,'%Y-%m-%d %H:%M:%S'),'%Y-%m-%d')
colnames(staff)[colnames(staff)=='dateTimeSample']<-'datetime'
staff<-staff[,c(6,12)]
staff2014<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Data/staffGaugeLog.csv',stringsAsFactor=F)
staff2014$datetime<-strftime(strptime(paste(staff2014$Date,staff2014$Time),'%m/%d/%Y %H:%S'),'%Y-%m-%d')
staff2014<-staff2014[staff2014$Lake.ID==lake&staff2014$Site=='Whole Lake',]
staff2014$Water.Height<-staff2014$Water.Height*.3048 # converting to meters 
staff2014<-staff2014[,c('datetime','Water.Height')]
colnames(staff2014)<-c('datetime','waterHeight_m')
staff<-rbind(staff,staff2014)
staff2015<-read.csv('/Users/Jake/Documents/Jake/UNDERC 2015/Data/StaffGauges2015_1.csv',stringsAsFactors = F)
staff2015$dateTimeSample<-strftime(strptime(staff2015$dateTimeSample,'%Y/%m/%d %H:%M'),'%Y-%m-%d')
staff2015<-staff2015[staff2015$lakeID==lake&staff2015$siteName=='WholeLake',]
staff2015<-staff2015[,c('dateTimeSample','waterHeight_m')]
colnames(staff2015)<-c('datetime','waterHeight_m')
staff<-rbind(staff,staff2015)
staff<-staff[sort.list(staff$datetime),]
doc<-dbTable('doc',lakeID=lake,depthClass='PML')
doc<-doc[doc$flag==0,] #grabbing PML doc data for CR - change in storage term (back calc inflows)
doc<-aggregate(doc$DOC,by=list(doc$dateSample),FUN=mean)
doc<-data.frame(datetime=strftime(strptime(doc$Group.1,'%Y-%m-%d %H:%M:%S'),'%Y-%m-%d'),doc=doc$x)
doc2014<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Data/DOC2014.csv',stringsAsFactor=F)
doc2014<-doc2014[doc2014$Lake.ID==lake&doc2014$Depth.Class=='PML'&doc2014$flag==0,]
doc2014$Date.Sample<-strftime(strptime(doc2014$Date.Sample,'%m/%d/%Y'),'%Y-%m-%d')
doc2014<-aggregate(doc2014$DOC_mgL,by=list(doc2014$Date.Sample),FUN=mean)
colnames(doc2014)<-c('datetime','doc')
doc<-rbind(doc,doc2014)
doc2015<-read.csv('/Users/Jake/Documents/Jake/UNDERC 2015/Data/Limnology/DOC2015.csv',stringsAsFactors = F)
doc2015<-doc2015[doc2015$Lake.ID==lake&doc2015$Depth.Class=='PML'&doc2015$flag==0,]
doc2015$Date.Sample<-strftime(strptime(doc2015$Date.Sample,'%m/%d/%Y'),'%Y-%m-%d') 
doc2015<-aggregate(doc2015$DOC_mgL,by=list(doc2015$Date.Sample),FUN=mean)
colnames(doc2015)<-c('datetime','doc')
doc<-rbind(doc,doc2015)
airT<-read.csv('C:/Users/Jake/Documents/Jake/UNDERC 2013/B3/Output/WL/WL2013/WL2013.csv',header=T)
datetime<-strftime(strptime(paste(airT[,1],airT[,2]),'%d/%m/%Y %H:%M'),'%Y-%m-%d')
airT<-data.frame(datetime=datetime,airT=airT[,grep('airTemp',colnames(airT))])
airT<-na.omit(aggregate(airT$airT,by=list(airT$datetime),FUN=mean,na.rm=T))
colnames(airT)<-c('datetime','airT') #airTemp in degrees C
airT$datetime<-as.Date(airT$datetime)
airT2014<-read.table('/Users/Jake/Documents/Jake/UNDERC 2014/LakeMetabolism/InputFiles/WL/airTemp.txt',stringsAsFactor=F,
                     sep='\t',header=T)
airT2014<-na.omit(aggregate(airT2014$airT,by=list(as.Date(airT2014$datetime)),FUN=mean,na.rm=T))
colnames(airT2014)<-c('datetime','airT')
airT<-rbind(airT,airT2014)
airT2015<-read.table('/Users/Jake/Documents/Jake/UNDERC 2015/LakeMetabolism/InputFiles/Pelagic/WL/airTemp.txt',stringsAsFactors = F,sep='\t',header=T)
airT2015<-na.omit(aggregate(airT2015$airT,by=list(as.Date(airT2015$datetime)),FUN=mean,na.rm=T))
colnames(airT2015)<-c('datetime','airT')
airT<-rbind(airT,airT2015)
airT$datetime<-strftime(airT$datetime,'%Y-%m-%d')
evap<-data.frame(datetime=airT$datetime,evap=evapCalc(airT$airT,DOY(airT$datetime),46.15)*lakeArea/1000) # evap in m^3 day-1 
# loading TP data for TP loading 2015-11-02
source('/Users/Jake/Desktop/R functions/dbTable.R')
tp<-dbTable('nutrients',lakeID = lake) #function that grabs table from database 
tp2014<-read.csv('/Users/Jake/Documents/Jake/UNDERC 2014/Data/Limnology/TP//TP2014.csv',stringsAsFactors=F) #read in 2014 data 
tp<-tp[,c(1,3,4,6:10,16,17)]
tp$dateTimeSample<-as.POSIXct(tp$dateTimeSample)
tp2014$Date.Sample<-as.POSIXct(paste(tp2014$Date.Sample,tp2014$Time.Sample),format='%m/%d/%Y %H:%M')
tp2014<-tp2014[,c(2:5,7:9,12:14)]
colnames(tp2014)<-colnames(tp)
tp2014<-tp2014[tp2014$lakeID==lake,]
tp<-rbind(tp,tp2014)
rm(tp2014)
tp$siteName<-ifelse(tp$siteName=='Inlet','Inlet1',tp$siteName)
tpPML<-tp[tp$depthClass=='PML'&tp$TP>0&tp$TP<60,c('dateTimeSample','TP')]

### stream doc contribution 
dir<-'C:/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Data/StreamDischarge/'
data<-data.frame()
years<-2014
for(i in 1:length(years)){
  for(j in 1:length(site)){
    cur<-read.csv(paste(dir,years[i],'/',lake,'/',site[j],'/',lake,site[j],'Discharge.csv',sep=''))
    colnames(cur)<-tolower(colnames(cur))
    cur$datetime<-as.POSIXct(cur$datetime)
    cur<-cur[,c('datetime','docdischarge','tpdischarge','dischargetotal','inlet1.temp')]
    data<-rbind(data,cur)
  }
}
years<-c(2014,2015)
cur2015<-read.csv('/Users/Jake/Documents/Jake/UNDERC 2015/StreamDischarge/EL/Inlet1/ELInlet1Discharge.csv',stringsAsFactors = F)
colnames(cur2015)<-tolower(colnames(cur2015))
cur2015$datetime<-as.POSIXct(cur2015$datetime)
cur2015<-cur2015[,c('datetime','docdischarge','tpdischarge','dischargetotal','inlet1.temp')]
data<-rbind(data,cur2015)

data<-na.omit(data)
dataDOC<-aggregate(data$docdischarge,by=list(data$datetime),FUN=sum) #aggregate all inlets doc discharge
dataTP<-aggregate(data$tpdischarge,by=list(data$datetime),FUN=sum) #aggregate all inlets tp discharge
dataWater<-aggregate(data$dischargetotal,by=list(data$datetime),FUN=sum) #aggregate all inlets water discharge m3
dataTemp<-aggregate(data$inlet1.temp,by=list(data$datetime),FUN=mean)
colnames(dataDOC)<-c('datetime','streamDOCdisch')
colnames(dataTP)<-c('datetime','streamTPdisch')
colnames(dataWater)<-c('datetime','streamWaterdisch')
colnames(dataTemp)<-c('datetime','streamTemp')
dataDOC<-dataDOC[sort.list(dataDOC$datetime),] #tells how much DOC was discharged over 10 minutes 
dataTP<-dataTP[sort.list(dataTP$datetime),] #tells how much DOC was discharged over 10 minutes 
dataWater<-dataWater[sort.list(dataWater$datetime),]
dataTemp<-dataTemp[sort.list(dataTemp$datetime),]
#aggregating loads to daily means in mg C/m2/day 
dataDOC$datetime<-strftime(strptime(dataDOC$datetime,'%Y-%m-%d %H:%M:%S'),'%Y-%m-%d')
dataDOC<-aggregate(dataDOC$streamDOCdisch,by=list(dataDOC$datetime),FUN=sum)
colnames(dataDOC)<-c('datetime','streamDOCdisch')
dataTP$datetime<-strftime(strptime(dataTP$datetime,'%Y-%m-%d %H:%M:%S'),'%Y-%m-%d')
dataTP<-aggregate(dataTP$streamTPdisch,by=list(dataTP$datetime),FUN=sum)
colnames(dataTP)<-c('datetime','streamTPdisch')
dataWater$datetime<-strftime(strptime(dataWater$datetime,'%Y-%m-%d %H:%M:%S'),'%Y-%m-%d')
dataWater<-aggregate(dataWater$streamWaterdisch,by=list(dataWater$datetime),FUN=sum)
colnames(dataWater)<-c('datetime','streamWaterdisch')
dataTemp$datetime<-strftime(strptime(dataTemp$datetime,'%Y-%m-%d %H:%M:%S'),'%Y-%m-%d')
dataTemp<-aggregate(dataTemp$streamTemp,by=list(dataTemp$datetime),FUN=mean)
colnames(dataTemp)<-c('datetime','streamTemp')
data<-merge(dataWater,dataDOC,by='datetime')
data<-merge(data,dataTP,by='datetime')
data<-merge(data,dataTemp,by='datetime')
data$streamDens<-water.density(data$streamTemp)

#read in precipitation data 
precip<-read.csv('C:/Users/Jake/Documents/Jake/UNDERC 2013/B3/Output/UNDERC/UNDERC2013/UNDERC2013.csv')
datetime<-strftime(strptime(paste(precip[,1],precip[,2]),'%d/%m/%Y %H:%M'),'%Y-%m-%d')
precip<-data.frame(datetime=datetime,precip=precip[,grep('rain',colnames(precip))])
precip<-na.omit(aggregate(precip$precip,by=list(precip$datetime),FUN=sum,na.rm=T))
colnames(precip)<-c('datetime','precip') #precip from UNDERC met station in mm/day 
## including 2014 precip
april2014<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Data/2014Precip/April2014_KingsLandOLakesAirport.csv',
                    stringsAsFactor=F)
may2014<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Data/2014Precip/May2014_KingsLandOLakesAirport.csv',
                  stringsAsFactor=F)
june2014<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Data/2014Precip/June2014_KingsLandOLakesAirport.csv',
                   stringsAsFactor=F)
july2014<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Data/2014Precip/July2014_KingsLandOLakesAirport.csv',
                   stringsAsFactor=F)
kingsPrecip<-rbind(april2014,may2014,june2014,july2014)
kingsPrecip<-kingsPrecip[,c('CDT','PrecipitationIn')]
colnames(kingsPrecip)<-c('datetime','precip')
kingsPrecip$precip<-kingsPrecip$precip*25.4 #inches to mm conversion 
kingsPrecip$datetime<-strftime(strptime(kingsPrecip$datetime,'%m/%d/%Y'))
# 2014 precip data from WL tipping bucket
wlPrecip<-read.csv('/Users/Jake/Documents/Jake/UNDERC 2014/WL_Precip.csv',stringsAsFactor=F)
colnames(wlPrecip)<-c('datetime','precip')
wlPrecip$datetime<-strftime(strptime(wlPrecip$datetime,'%Y-%m-%d %H:%M:%S'),'%Y-%m-%d')
wlPrecip<-aggregate(wlPrecip$precip,by=list(wlPrecip$datetime),FUN=sum) # units already in mm 
colnames(wlPrecip)<-c('datetime','precip')
precip2014<-rbind(wlPrecip,kingsPrecip)
precip2014<-precip2014[!duplicated(precip2014$datetime),] #keeping non-duplicated dates, wl gauge when datapoints overlap 
precip2014<-precip2014[sort.list(precip2014$datetime),]
precip2015<-read.table('/Users/Jake/Documents/Jake/UNDERC 2015/WL_Precip.txt',sep='\t',stringsAsFactors = F,header=T)
precip2015$datetime<-strftime(strptime(precip2015$datetime,'%Y-%m-%d %H:%M:%S'),'%Y-%m-%d')
precip2015<-aggregate(precip2015$precip,by=list(precip2015$datetime),FUN=sum)
colnames(precip2015)<-c('datetime','precip')
precip<-rbind(precip,precip2014)
precip<-rbind(precip,precip2015)
precip$datetime<-as.character(precip$datetime)
staff$datetime<-as.character(staff$datetime)
doc$datetime<-as.character(doc$datetime)
evap$datetime<-as.character(evap$datetime)
data<-merge(data,precip,by='datetime',all.x=T)
data<-merge(data,staff,by='datetime',all.x=T)
data<-merge(data,doc,by='datetime',all.x=T)
data<-merge(data,evap,by='datetime',all.x=T)
colnames(tpPML)<-c('datetime','tp')
tpPML$datetime<-strftime(tpPML$datetime,'%Y-%m-%d')
data<-merge(data,tpPML,by='datetime',all.x=T)
# 
# data$waterHeight_m[1]<-data$waterHeight_m[!is.na(data$waterHeight_m)][1] # making first staff gauge first non-na 
# data$waterHeight_m[length(data$datetime)]<-
#   data$waterHeight_m[!is.na(data$waterHeight_m)][length(data$waterHeight_m[!is.na(data$waterHeight_m)])] # making last staff gauge last non-na 
# data$waterHeight_m<-approx(x=as.Date(data$datetime),y=data$waterHeight_m,xout=as.Date(data$datetime))$y
#change in storage term for DOC 

data$precip<-data$precip*lakeArea/1000 #volume of precip in m^3/day
data$precipDOC<-data$precip*1000*3.195 #total amount of DOC from precip in mg/day - DOC is average of precip in our database
data$precipTP<-data$precip*1000*1.3 #total amount of TP from precip in ug/day - TP is average of precip from Grimshaw and Dolske 2001 

data$gwDisch<-rep(20,length(data$datetime)) #m3/day
data$gwDischTP<-data$gwDisch*data$tp*1000
data$gwDischDOC<-data$gwDisch*data$doc*1000 #total amount of DOC leaving in GW mg/day


if(lake=='MO'){
  data$totalAlloC<-(data$streamDOCdisch+data$precipDOC+data$gwDischDOC)/lakeVol #total allochthonous carbon flux in mg C m-3 
  data$totalAlloP<-(data$streamTPdisch+data$precipTP+data$gwDischTP)/lakeVol # total allochthonous P flux in ug P m-3 
}else if(lake=='EL'){
  data$totalAlloC<-(data$streamDOCdisch+data$precipDOC)/lakeVol #total allochthonous carbon flux in mg C m-3 
  data$totalAlloP<-(data$streamTPdisch+data$precipTP)/lakeVol # total allochthonous P flux in ug P m-3 
}else if(lake=='CR'){
  tempData<-data
  data2<-data.frame()
  for(i in 1:length(years)){
    curData<-tempData[as.POSIXct(tempData$datetime)>as.POSIXct(paste(years[i],'01-01',sep='-'),format='%Y-%m-%d')&
                        as.POSIXct(tempData$datetime)<as.POSIXct(paste(years[i],'12-30',sep='-'),format='%Y-%m-%d'),]
    data<-curData
  # balancing water budget on a daily scale so that In = Out. Precip + Overland = GWout + Outlet + Evap
  overland<-data$gwDisch+data$streamWaterdisch+data$evap-data$precip # some of overland flow is negative 
  overlandNeg<--1*overland[which(overland<0)] #all days for which overland is negative 
  overland<-ifelse(overland<0,overland*-1,overland) # forcing negative overland values to 0
  # take away overland water flux from days on which no rain occured 
  
  #   noPrecip<-length(which(data$precip==0)) # days for with precip 
  overlandNeg<-sum(overlandNeg*2) #need to subtract this much from the overland flow (preferably on days for which there was no precip)
  #   evenSubtraction<-overlandNeg/noPrecip # amount of water subtracted from each day if there excess was evenly spread over all days 
  #   overland<-ifelse(data$precip==0,overland-evenSubtraction,overland) #will create some negative values
  #   newNeg<-sum(-1*overland[which(overland<0)]) # amount of water still needed to subtract from overland 
  #   overland<-ifelse(overland<0,0,overland) #forcing overland to 0 if negative 
  evenSubtraction<-overlandNeg/length(which(data$precip==0)) # amount of water subtracted from each day if there excess was evenly spread over all days with no precip 
  subtracted<-sum(overland[which(overland<evenSubtraction&data$precip==0)]) # amount of water forcing overland less than evenSubtraction to 0
  overland<-ifelse(overland<evenSubtraction&data$precip==0,0,overland) # subtracting water  from overland less than evenSubtraction 
  needToSubtract<-overlandNeg-subtracted #still need to subtract this much water 
  subtracted<-length(overland[which(overland>evenSubtraction&data$precip==0)])*evenSubtraction
  overland<-ifelse(overland>evenSubtraction&data$precip==0,overland-evenSubtraction,overland) 
  needToSubtract<-needToSubtract-subtracted
  evenSubtraction<-needToSubtract/length(which(data$precip==0&overland>0)) # amount of water subtracted from each day if there excess was evenly spread over all days with no precip 
  subtracted<-sum(overland[which(overland<evenSubtraction&data$precip==0)]) # amount of water forcing overland less than evenSubtraction to 0
  overland<-ifelse(overland<evenSubtraction&data$precip==0,0,overland) # subtracting water  from overland less than evenSubtraction 
  needToSubtract<-needToSubtract-subtracted #still need to subtract this much water 
  evenSubtraction<-needToSubtract/length(which(data$precip==0&overland>evenSubtraction))
  subtracted<-length(overland[which(overland>evenSubtraction&data$precip==0)])*evenSubtraction
  overland<-ifelse(overland>evenSubtraction&data$precip==0,overland-evenSubtraction,overland) 
  needToSubtract<-needToSubtract-subtracted
  evenSubtraction<-needToSubtract/length(which(data$precip==0&overland>0)) # amount of water subtracted from each day if there excess was evenly spread over all days with no precip 
  subtracted<-sum(overland[which(overland<evenSubtraction&data$precip==0)]) # amount of water forcing overland less than evenSubtraction to 0
  overland<-ifelse(overland<evenSubtraction&data$precip==0,0,overland) # subtracting water  from overland less than evenSubtraction 
  needToSubtract<-needToSubtract-subtracted #still need to subtract this much water 
  evenSubtraction<-needToSubtract/length(which(data$precip==0&overland>evenSubtraction))
  subtracted<-length(overland[which(overland>evenSubtraction&data$precip==0)])*evenSubtraction
  overland<-ifelse(overland>evenSubtraction&data$precip==0,overland-evenSubtraction,overland)  # this 
  needToSubtract<-needToSubtract-subtracted #still need to subtract this much water - should be 0 now 
  
  data$overland<-overland #m3 day-1
  # sebestyen lab water chem for overland DOC 
  ssLog<-read.csv('/Users/Jake/Documents/Jake/UNDERC 2015/Sebestyen Lab Water Chem Log.csv',
                  stringsAsFactor=F)
  ssData1<-read.csv('/Users/Jake/Documents/Jake/UNDERC 2014/Data/Limnology/2015.09.15_101.2014Data_toZwartJ_sebestyen.csv',
                    stringsAsFactor=F)#2014 data 
  ssData2<-read.csv('/Users/Jake/Documents/Jake/UNDERC 2014/Data/Limnology/2015.09.16_101.2013Data_toZwartJ_sebestyen.csv',
                    stringsAsFactor=F)# 2013 data 
  ssData3<-read.csv('/Users/Jake/Documents/Jake/UNDERC 2014/Data/Limnology/2012.10.02_101000-101068_chem.data_to.Zwart.J_sebestyen.csv',
                    stringsAsFactor=F)# 2011-2012 data 
  colnames(ssData3)<-ssData3[1,]
  ssData3<-ssData3[c(4:87),c(1:44)]
  colnames(ssData2)<-ssData2[1,]
  ssData2<-ssData2[c(4:66),c(1:44)]
  colnames(ssData1)<-ssData1[1,]
  ssData1<-ssData1[c(4:70),c(1:41)]
  ssData1$idnbr<-ssData1$idnbr.2
  ssData2$idnbr<-ssData2$idnbr.1
  vars<-c('idnbr','pH','ANC','K25','Cl','NO3N','DIP','SO4','UV vis','NH3N','NO3+NO2N','TP','TN','P as PO4',
          'Al','Ca','Fe','K','Mg','Mn','Na','Si','Sr','DOC ')
  ssData1<-ssData1[,colnames(ssData1)%in%vars]
  ssData2<-ssData2[,colnames(ssData2)%in%vars]
  ssData3<-ssData3[,colnames(ssData3)%in%vars]
  ssData<-rbind(ssData1,ssData2,ssData3)
  colnames(ssLog)<-c('idnbr','projID','lakeID','site','dateSample','timeSample','depthClass','depthTop',
                     'depthBot','comments')
  colnames(ssData)[24]<-'DOC'
  ssData<-merge(ssLog,ssData,by='idnbr',all.x=T)
  ssData$datetime<-as.POSIXct(ssData$dateSample,format='%m/%d/%Y')
  colnames(ssData)[grep('P as PO4',colnames(ssData))]<-'PO4'
  overlandDOC<-ssData[ssData$site%in%c('S1','S4','S5','S3','P1','P19','P10','P18','P16','P13','P15','P11','P9','P7','P5','P2','P14','P4','P3'),]
  overlandDOCmean<-mean(as.numeric(overlandDOC$DOC),na.rm=T) # mean DOC concnetration from shallow sub-surface well collectors 
  overlandTPmean<-mean(as.numeric(overlandDOC$TP),na.rm=T) # mg L-1 PO4 
  
  data$overlandDOCdisch<-data$overland*overlandDOCmean*1000 # mg C day-1
  data$overlandTPdisch<-data$overland*overlandTPmean*1000*1000 # ug P day-1 
  
  data$totalAlloC<-(data$overlandDOCdisch+data$precipDOC)/lakeVol #total allochthonous carbon flux in mg C m-3 
  data$totalAlloP<-(data$overlandTPdisch+data$precipTP)/lakeVol
  # if precip adds more DOC to lake than leaves on a given day in outlet+gw, then totalAlloC = precip doc 
  #   data$totalAlloC<-ifelse(data$totalAlloC<(data$precipDOC/lakeArea),(data$precipDOC/lakeArea),data$totalAlloC)
  data$totalAlloC<-ifelse(data$totalAlloC<0,0,data$totalAlloC)
  data$totalAlloP<-ifelse(data$totalAlloP<0,0,data$totalAlloP)
  
  data2<-rbind(data2,data)
  }
  data<-data2
  data$overland<-ifelse(is.na(data$overland),0,data$overland)
  data$totalAlloC<-ifelse(is.na(data$totalAlloC),0,data$totalAlloC)
  data$totalAlloP<-ifelse(is.na(data$totalAlloP),0,data$totalAlloP)
  data$overlandDOCdisch<-ifelse(is.na(data$overlandDOCdisch),0,data$overlandDOCdisch)
  data$overlandTPdisch<-ifelse(is.na(data$overlandTPdisch),0,data$overlandTPdisch)
}
# should CR DOC load be outlet + GW out? 

# outlet discharge for EL based on salt slugs ~ stage relationship 
if(lake=='EL'){
  data$outletDischarge<-5e-13*exp(31.127*data$waterHeight_m) # discharge in m3 sec-1 - see excel spreadsheet in Pulse / lag folder (****change to 1/2 of excel spreadsheet value cuz water wasn't balanced)
  # updated on 2016-12-08
  data$outletDischarge<-3.1889*data$waterHeight_m^21.396
  data$outletDischarge<-data$outletDischarge*60*60*24 #discharge in m3 day-1
}else if(lake=='MO'){ 
  # unable to get a sucessful salt slug because outlet is so diffuse; precip is about equal to evap, gw is neglible, so streamIN = streamOut
  
  Stage0<-data$waterHeight_m[1] #initial stage height in m (all water added to the lake will be in reference to this height)
  V0<-lakeVol # initial lake Volume - assuming lake stays a pretty constant volume 
  
  waterIn<-data$streamWaterdisch+data$precip+data$gwDisch # all water that goes into MO
  data$delta_storage<-ifelse(abs(data$delta_storage>2000),0,data$delta_storage)
  waterOut<-data$evap-data$delta_storage #water (except outlet) that leaves MO or change in storage 
  outlet<-waterIn-waterOut #difference between waterIn and evap+storage should be outlet discharge 
  data$outletDischarge<-outlet
}


metab<-data.frame()
years<-c(2014)
for(i in 1:length(years)){
  files<-list.files(file.path('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Data/Bootstrapped/',
                              years[i],lake))
  if(length(files[grep('results',files)])>0){
    files<-files[-grep('results',files)]
  }
  for(j in 1:length(files)){
    cur<-read.table(file.path('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Data/Bootstrapped/',
                              years[i],lake,files[j]),sep='\t',header=T,stringsAsFactor=F)
    metab<-rbind(metab,cur)
  }
}
years<-c(2014,2015)
aggGPP<-aggregate(metab$GPP,by = list(metab$solarDay),FUN=mean)
aggR<-aggregate(metab$rho,by=list(metab$solarDay),FUN=mean)
aggIota<-aggregate(metab$iota,by = list(metab$solarDay),FUN=mean)
sdIota<-aggregate(metab$iota,by = list(metab$solarDay),FUN=sd)
sdGPP<-aggregate(metab$GPP,by=list(metab$solarDay),FUN=sd)
sdR<-aggregate(metab$rho,by=list(metab$solarDay),FUN=sd)
metabAgg<-merge(aggGPP,aggR,by='Group.1')
metabAgg<-merge(metabAgg,aggIota,by='Group.1')
metabAgg<-merge(metabAgg,sdGPP,by='Group.1')
metabAgg<-merge(metabAgg,sdR,by='Group.1')
metabAgg<-merge(metabAgg,sdIota,by='Group.1')
colnames(metabAgg)<-c('datetime','GPP','R','Iota','sdGPP','sdR','sdIota')
metabAgg$NEP<-metabAgg$GPP-metabAgg$R
metab2015<-read.table('/Users/Jake/Documents/Jake/UNDERC 2015/LakeMetabolism/Results/Pelagic/EL/EL_2015 optimOut.txt',
                      sep='',stringsAsFactors = F,header=T)
gpp2015<-read.table('/Users/Jake/Documents/Jake/UNDERC 2015/LakeMetabolism/Results/Pelagic/EL/EL_2015 GPPFitOut.txt',
                    sep='',stringsAsFactors = F, header = T)
metab2015<-merge(metab2015,gpp2015,by='solarDay',all.x=T)
metab2015<-metab2015[,c('solarDay','GPP','rhoEst','iotaEst','GPPSd','rhoSd','iotaSd')]
metab2015$NEP<-metab2015$GPP-metab2015$rhoEst
colnames(metab2015)<-colnames(metabAgg)
metabAgg<-rbind(metabAgg,metab2015)
data<-merge(data,metabAgg,by='datetime',all.x=T)
data<-data[data$datetime>=dateRange[1]&data$datetime<=dateRange[2],]
data<-data[sort.list(data$datetime),]
source('C:/Users/Jake/Desktop/R functions/MA.R')
source('C:/Users/Jake/Desktop/R functions/ma.weight.R')

# adding in overland flow estimates based on CR estimates 
if(lake=='CR'){
  crData<-data
  save(crData,file='/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/R Data/CRdata.Rda')
}
if(lake%in%c('EL','MO')){
  load('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/R Data/CRdata.Rda')
  crData<-crData[,c('datetime','overland','overlandDOCdisch','overlandTPdisch')]
  # trying to estimate overland flow for EL and MO assuming it's about the same perimeter rate as in CR 
  perim<-data.frame(lakeID=c('MO','EL','CR'),perim=c(1050.957,800.90,3455.0459),vol=c(142454,123734,1302506),
                    wa=c(1122138,72003,319951))
  crData$overland<-crData$overland/perim$perim[perim$lakeID=='CR']*perim$perim[perim$lakeID==lake] # m3 / m shoreline / day 
  crData$overlandDOCdisch<-crData$overlandDOCdisch/perim$perim[perim$lakeID=='CR']*perim$perim[perim$lakeID==lake] # m3 / m shoreline / day 
  crData$overlandTPdisch<-crData$overlandTPdisch/perim$perim[perim$lakeID=='CR']*perim$perim[perim$lakeID==lake] # m3 / m shoreline / day 
  colnames(crData)[1]<-c('datetime')
  
  data<-merge(data,crData,by='datetime',all.x=T)
  # we don't want to double count the inlet streams so overland flow is overland flow minus inlet stream discharge 
  data$overland<-data$overland-data$streamWaterdisch
  data$overland<-ifelse(data$overland<0,0,data$overland)
  data$overlandDOCdisch<-data$overland*24.2*1000
  data$overlandTPdisch<-data$overland*1000*1000*0.073333333
}

if(lake=='MO'){
  data$totalAlloC<-(data$streamDOCdisch+data$precipDOC+data$gwDischDOC+data$overlandDOCdisch)/lakeVol #total allochthonous carbon flux in mg C m-3 
  data$totalAlloP<-(data$streamTPdisch+data$precipTP+data$gwDischTP+data$overlandTPdisch)/lakeVol # total allochthonous P flux in ug P m-3 
}else if(lake=='EL'){
  data$totalAlloC<-(data$streamDOCdisch+data$precipDOC+data$overlandDOCdisch)/lakeVol #total allochthonous carbon flux in mg C m-3 
  data$totalAlloP<-(data$streamTPdisch+data$precipTP+data$overlandTPdisch)/lakeVol # total allochthonous P flux in ug P m-3 
}


par13<-read.table('/Users/Jake/Documents/Jake/MyPapers/FINISHED/Long Lake Metabolism/Data/Pelagic/2013/EL/PAR.txt',stringsAsFactor=F,
                  header=T,sep='\t')
par14<-read.table('/Users/Jake/Documents/Jake/MyPapers/FINISHED/Long Lake Metabolism/Data/Pelagic/2014/EL/PAR.txt',stringsAsFactor=F,
                  header=T,sep='\t')
par15<-read.table('/Users/Jake/Documents/Jake/UNDERC 2015/LakeMetabolism/InputFiles/Pelagic/EL/PAR.txt',
                  stringsAsFactors = F,header=T,sep='\t')
parData<-rbind(par13,par14,par15)
rm(par13,par14,par15)
parData$datetime<-as.character(as.Date(parData$datetime))
parData<-aggregate(parData$PAR,by=list(parData$datetime),FUN=sum)
colnames(parData)<-c('datetime','PAR')
parData$PAR<-parData$PAR*10*60/1000 # converting par from umol photons m-2 s-1 (discrete measure) to mmol photons m-2 day-1 
ws13<-read.table('/Users/Jake/Documents/Jake/MyPapers/FINISHED/Long Lake Metabolism/Data/Pelagic/2013/EL/WS.txt',stringsAsFactor=F,
                  header=T,sep='\t')
ws14<-read.table('/Users/Jake/Documents/Jake/MyPapers/FINISHED/Long Lake Metabolism/Data/Pelagic/2014/EL/WS.txt',stringsAsFactor=F,
                  header=T,sep='\t')
ws15<-read.table('/Users/Jake/Documents/Jake/UNDERC 2015/LakeMetabolism/InputFiles/Pelagic/EL/WS.txt',stringsAsFactors = F,
                 sep='\t',header=T)
wsData<-rbind(ws13,ws14,ws15)
rm(ws13,ws14,ws15)
wsData$datetime<-as.character(as.Date(wsData$datetime))
wsData<-aggregate(wsData$WS,by=list(wsData$datetime),FUN=mean) # mean windspeed day-1 for k estimates 
colnames(wsData)<-c('datetime','wnd')

thermo.depth<-data.frame()
zmix<-data.frame()
epiTemp<-data.frame()
for(i in 1:length(years)){
  cur<-read.table(file.path('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Data/',years[i],lake,'TEMP_PROF.txt'),
                  stringsAsFactor=F,header=T,sep='\t')
  if(length(cur[1,grep('Outlet',colnames(cur))])>0){
    cur<-cur[,-grep('Outlet',colnames(cur))]
  }
  cols<-colnames(cur)
  cols<-gsub('temp','wtr_',cols)
  colnames(cur)<-cols
  curZall<-ts.thermo.depth(cur,na.rm=T,seasonal=F)
  curMixall<-ts.meta.depths(cur,na.rm = T)
  curZ<-aggregate(curZall$thermo.depth,by = list(as.Date(curZall$datetime)),FUN=max,na.rm=T)
  curMix<-aggregate(curMixall$top,by = list(as.Date(curMixall$datetime)),FUN=mean,na.rm=T)
  colnames(curZ)<-c('datetime','thermo.depth')
  colnames(curMix)<-c('datetime','z.mix')
  curZ$thermo.depth<-ifelse(curZ$thermo.depth=='-Inf',max(get.offsets(cur)),curZ$thermo.depth) # if INF, set thermo to deepest temp chain depth 
  curMix$z.mix<-ifelse(curMix$z.mix=='-Inf',max(get.offsets(cur)),curMix$z.mix)
  thermo.depth<-rbind(thermo.depth,curZ)
  zmix<-rbind(zmix,curMix)
  curBathy<-bathy[,c('area_m2','depth_m')]
  colnames(curBathy)<-c('areas','depths')
  curEpiTemp<-ts.layer.temperature(wtr=cur,top = 0,bottom = curZall$thermo.depth,bathy = curBathy,na.rm = T)
  curEpiTemp<-aggregate(curEpiTemp$layer.temp,by=list(as.Date(curEpiTemp$datetime)),FUN=max,na.rm=T)
  colnames(curEpiTemp)<-c('datetime','epiTemp')
  epiTemp<-rbind(epiTemp,curEpiTemp)
}

thermo.depth<-thermo.depth[!duplicated(thermo.depth$datetime),]
thermo.depth$datetime<-as.character(thermo.depth$datetime)
zmix<-zmix[!duplicated(zmix$datetime),]
zmix$datetime<-as.character(zmix$datetime)
epiTemp<-epiTemp[!duplicated(epiTemp$datetime),]
epiTemp$datetime<-as.character(epiTemp$datetime)

data<-merge(data,thermo.depth,by='datetime',all.x=T)
data<-merge(data,zmix,by='datetime',all.x=T)
data<-merge(data,epiTemp,by='datetime',all.x=T)
source('/Users/Jake/Desktop/R functions/layer.volume.R')

data$epiVol<-rep(NA,length(data$datetime))
for(j in 1:length(data$datetime)){
  data$epiVol[j]<-layer.volume(bthA = as.numeric(bathy$area_m2),bthD=as.numeric(bathy$depth_m),top=0,bot=data$thermo.depth[j]) # volume of water (m3) in mixed layer
}

data$V0<-bathy$volumeToBottom_m3[1] #m3 lake Volume original from database 
data$A0<-bathy$area_m2[1] #m2 lake area original from database 
data<-merge(data,parData,by='datetime',all.x=T)
data<-merge(data,wsData,by='datetime',all.x=T)

gc2013<-read.table('/Users/Jake/Documents/Jake/UNDERC 2013/GC data/Results/2013_UNDERC_database_2015-01-20.txt',
                   stringsAsFactor=F,sep='\t',header=T)
gc2014<-read.table('/Users/Jake/Documents/Jake/UNDERC 2014/Data/GC/GC2014.txt',stringsAsFactor=F,
                   sep='\t',header=T)
gc15<-read.table('/Users/Jake/Documents/Jake/UNDERC 2015/Data/UNDERC 2015 GC/GC2015.txt',stringsAsFactors = F,
                 sep='\t',header=T)
gc15$dateSample<-strftime(strptime(gc15$dateSample,'%m/%d/%Y'),'%Y-%m-%d')
gc15$dateTimeSample<-strftime(strptime(gc15$dateTimeSample,'%m/%d/%Y %H:%M'),'%Y-%m-%d %H:%M:%S')
gc15$subsampleDateTime<-strftime(strptime(gc15$subsampleDateTime,'%m/%d/%Y %H:%M'),'%Y-%m-%d %H:%M:%S')
gc<-rbind(gc2013,gc2014)
gc<-gc[,-grep('location',colnames(gc))]
gc<-rbind(gc,gc15)
gc$siteName<-ifelse(gc$siteName=='Inlet','Inlet1',gc$siteName)
gc$depthClass<-ifelse(gc$depthClass=='surface','point',gc$depthClass)
rm(gc2013,gc2014,gc15)
gc$depthTop<-as.numeric(gc$depthTop)
gc$depthBottom<-as.numeric(gc$depthBottom)
# assuming that DIC = pCO2 in EL for now (ignoring any carbonate equilibrium modeling)
elInletCO2<-gc[paste(gc$lakeID,gc$siteName,gc$subsampleClass)=='EL Inlet1 pCO2',]
elDeepholeCO2<-gc[paste(gc$lakeID,gc$siteName,gc$depthClass,gc$subsampleClass)==
                    'EL deepHole point pCO2'&gc$depthTop<0.9,]
elHypoDIC<-gc[paste(gc$lakeID,gc$siteName,gc$depthClass,gc$depthTop,gc$subsampleClass)==
                'EL deepHole point 8 DIC',]

meanInletCO2<-mean(elInletCO2$CO2original_umolL[elInletCO2$CO2original_umolL>0],na.rm=T)
meanInletCO2<-meanInletCO2/1000/1000*1000 # converting from umol L-1 to mol C m-3
elInletCO2<-aggregate(elInletCO2$CO2original_umolL~elInletCO2$dateSample,FUN=mean)
colnames(elInletCO2)<-c('datetime','Inlet_dic')
elInletCO2$Inlet_dic<-elInletCO2$Inlet_dic/1000/1000*1000 #converting from umol L-1 to mol C m-3 
elDeepholeCO2<-aggregate(elDeepholeCO2$CO2original_umolL~elDeepholeCO2$dateSample,FUN=mean)
colnames(elDeepholeCO2)<-c('datetime','dic') 
elDeepholeCO2$dic<-elDeepholeCO2$dic/1000/1000*1000 #converting from umol L-1 to mol C m-3
elHypoDIC<-aggregate(elHypoDIC$CO2original_umolL~elHypoDIC$dateSample,FUN=mean)
colnames(elHypoDIC)<-c('datetime','hypo_dic') 
elHypoDIC$hypo_dic<-elHypoDIC$hypo_dic/1000/1000*1000 #converting from umol L-1 to mol C m-3
data<-merge(data,elDeepholeCO2,by='datetime',all.x=T)
data<-merge(data,elInletCO2,by='datetime',all.x=T)
data<-merge(data,elHypoDIC,by='datetime',all.x=T)

# linearly interpolating inlet dic to get loads 
data$Inlet_dic[1]<-data$Inlet_dic[min(which(!is.na(data$Inlet_dic)))]
data$Inlet_dic[length(data$datetime)]<-data$Inlet_dic[max(which(!is.na(data$Inlet_dic)))]
data$Inlet_dic<-approx(as.Date(data$datetime),data$Inlet_dic,xout = as.Date(data$datetime))$y

# find lake DIC data ************************ 
data$docIn<-data$totalAlloC/12/1000*data$epiVol #mg C m-3 day-1 to mol C day-1 into epi
data$docIn<-data$totalAlloC/12/1000*data$V0 #mol C day-1 
data$tpIn<-data$totalAlloP/31/1000/1000*data$epiVol #ug P m-3 day-1 to mol P day-1 
data$tpIn<-data$totalAlloP/31/1000/1000*data$V0 #ug P m-3 day-1 to mol P day-1 
data$Qout<-(data$outletDischarge+data$gwDisch)/data$V0*data$epiVol # water discharge out m3 day-1 from the epi (partitioned equally from Epi and hypo )
data$Qout<-(data$outletDischarge+data$gwDisch) # water discharge out m3 day-1 
data$dicIn<-(data$streamWaterdisch+data$overland)*data$Inlet_dic/data$V0*data$epiVol # mol C day-1 into epi
data$dicIn<-(data$streamWaterdisch+data$overland)*data$Inlet_dic # mol C day-1 
data$DICeq<-rep(0.0136,length(data$datetime)) #DIC equilibrium from Chris;#Concentration of CO2 in water (mol m-3) at equilibrium with atmosphere, assume atmosphere is 400 ppm = 400 uatm, use Henry's law and gas constant 29.41 L atm mol-1
# put in precip DIC ???? 

# k estimates based on cole & caraco 1998; m day-1 
data$wnd_10<-wind.scale(data,wnd.z = 2)$wnd_10 # scaling wind speed to 10m 
data$k600<-k.cole(data)$k600 # k600 based on cole & caraco 1998; m day-1
data$k600<-k.vachon(data,lake.area = lakeArea)$k600 # k600 from Vachon and Prairie 2013 ; m day-1
data$wtr<-data$epiTemp # creating a wtr column for kCO2 
data$kCO2<-k600.2.kGAS(data,gas = 'CO2')$k.gas # m day-1

# pressure sensor water height for more resolved outlet dischage estimates 
highFreqStaff<-read.table('/Users/Jake/Documents/Jake/UNDERC 2014/PressureSensorWaterHeight/EL/StaffGauge/EL_StaffGauge_waterHeight.txt',
                          stringsAsFactors = F,sep='\t',header=T)
highFreqStaff15<-read.table('/Users/Jake/Documents/Jake/UNDERC 2015/PressureSensorWaterHeight/EL/WholeLake/EL_WholeLake_waterHeight.txt',
                            stringsAsFactors = F,sep='\t',header=T)
highFreqStaff$datetime<-as.POSIXct(highFreqStaff$datetime)
temp<-data[,c('datetime','waterHeight_m')]
temp$datetime<-as.POSIXct(temp$datetime,format='%Y-%m-%d')
highFreqStaff<-merge(highFreqStaff,temp,all.x=T)
temp2<-highFreqStaff[which(!is.na(highFreqStaff$waterHeight_m)),]
mean(temp2$waterHeight-temp2$waterHeight_m) # correction factor for actual water height in m 
highFreqStaff$waterHeightCorr<-highFreqStaff$waterHeight-mean(temp2$waterHeight-temp2$waterHeight_m)
highFreqStaff<-aggregate(highFreqStaff$waterHeightCorr,by=list(as.Date(highFreqStaff$datetime)),FUN=mean,na.rm=T)
colnames(highFreqStaff)<-c('datetime','highFreqWaterHeight')
highFreqStaff$datetime<-as.character(highFreqStaff$datetime)

highFreqStaff15$datetime<-as.POSIXct(highFreqStaff15$datetime)
temp<-data[,c('datetime','waterHeight_m')]
temp$datetime<-as.POSIXct(temp$datetime,format='%Y-%m-%d')
highFreqStaff15<-merge(highFreqStaff15,temp,all.x=T)
temp2<-highFreqStaff15[which(!is.na(highFreqStaff15$waterHeight_m)),]
mean(temp2$waterHeight-temp2$waterHeight_m,na.rm=T) # correction factor for actual water height in m 
highFreqStaff15$waterHeightCorr<-highFreqStaff15$waterHeight-mean(temp2$waterHeight-temp2$waterHeight_m,na.rm = T)
highFreqStaff15<-aggregate(highFreqStaff15$waterHeightCorr,by=list(as.Date(highFreqStaff15$datetime)),FUN=mean,na.rm=T)
colnames(highFreqStaff15)<-c('datetime','highFreqWaterHeight')
highFreqStaff15$datetime<-as.character(highFreqStaff15$datetime)
highFreqStaff<-rbind(highFreqStaff,highFreqStaff15)

data<-merge(data,highFreqStaff,by='datetime',all.x=T)
# adding in highfrequency water height if available 
data$waterHeight_m<-ifelse(is.na(data$highFreqWaterHeight),data$waterHeight_m,data$highFreqWaterHeight)

# linearly interpolating water height for missing data (especially when there is no high frequency data)
# data$waterHeight_m[1]<-data$waterHeight_m[min(which(!is.na(data$waterHeight_m)))]
# data$waterHeight_m<-approx(as.Date(data$datetime),data$waterHeight_m,xout = as.Date(data$datetime))$y

data$outletDischarge<-5e-13*exp(31.127*data$waterHeight_m) # discharge in m3 sec-1 - see excel spreadsheet in Pulse / lag folder 
data$outletDischarge<-3.1889*data$waterHeight_m^21.396 # updated on 2016-12-08
data$outletDischarge<-data$outletDischarge*60*60*24 #discharge in m3 day-1

data$Qout<-(data$outletDischarge+data$gwDisch) # water discharge out m3 day-1 
data$epiDens<-water.density(data$epiTemp) # water density of epilimnion 

# entrainment of hypo water into epi or loss of epi water into hypo 
data$epiDiff<-NA
# using negative to indicate epi water going into hypo; take into account evap / water loads / losses 
data$epiDiff[2:length(data$epiDiff)]<-diff(data$epiVol) 
data$overland<-ifelse(is.na(data$overland),0,data$overland) # overland water flux m3 day-1 
# # quick fix for missing outlet discharge data
data$QoutInt<-data$Qout # Qout is outlet + gw discharge out 
data$QoutInt[1]<-data$Qout[min(which(!is.na(data$Qout)))]
data$QoutInt[length(data$datetime)]<-data$Qout[max(which(!is.na(data$Qout)))]
data$timeStep<-seq(1:length(data$datetime))
data$QoutInt<-approx(data$timeStep,data$QoutInt,data$timeStep)$y
# water loads and losses 
data$waterLoad<-data$streamWaterdisch+data$precip+data$overland # m3 day-1
# water load - water loss / area = delta stage; we have uncertainty in outlet discharge 
data$outletDischarge2<-NA 
data$stageDiff<-NA
data$stageDiff[2:length(data$stageDiff)]<-diff(data$waterHeight_m)
for(i in 2:length(data$outletDischarge)){
  data$outletDischarge2[i-1]=data$waterLoad[i-1]-data$gwDisch[i-1]-data$evap[i-1]-data$stageDiff[i]*data$A0[i-1]
}
data$outletDischarge2<-ifelse(data$outletDischarge2<0,0,data$outletDischarge2)

#crest spillway - have to figure out height of outlet weir 
stageOut=0.68 # water height in m when outlet does not discharge 
C=(2/3)^1.5*9.806^0.5  
L=1 # length of spillway - see optimization from NHLD model (this is roughly Crampton's optim L - CR was L=2)
H=data$waterHeight_m-stageOut
H=ifelse(H<0,0,H)
Q=C*L*H^1.5  #m3 s-1
Q=Q*24*60*60 #m3 day-1 


plot(Q)
points(data$outletDischarge,col='red')
plot(Q~data$outletDischarge)
abline(0,1)

plot(Q~data$outletDischarge2)
abline(0,1)

plot(data$outletDischarge~data$outletDischarge2)
abline(0,1)

data$QoutInt<-Q # Qout is outlet + gw discharge out 
data$QoutInt[1]<-data$Qout[min(which(!is.na(data$Qout)))]
data$QoutInt[length(data$datetime)]<-data$Qout[max(which(!is.na(data$Qout)))]
data$timeStep<-seq(1:length(data$datetime))
data$QoutInt<-approx(data$timeStep,data$QoutInt,data$timeStep)$y

data$waterLoss<-data$QoutInt+data$evap # m3 day-1
# data$entrainVol<-data$waterLoad-data$waterLoss+data$epiDiff # Entrain volume balances the water budget for the Epilimnion, negative means loss to hypo, positive means gain from hypo
data$entrainVol<-NA
for(i in 2:length(data$entrainVol)){
  data$entrainVol[i-1]<-data$epiDiff[i]-data$waterLoad[i-1]+data$waterLoss[i-1]
}

# hypolimnion DOC and DIC entrainment 
doc2014<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Data/DOC2014.csv',stringsAsFactor=F)
doc2015<-read.csv('/Users/Jake/Documents/Jake/UNDERC 2015/Data/Limnology/DOC2015.csv',stringsAsFactors = F)
doc2014<-rbind(doc2014,doc2015)
doc2014<-doc2014[doc2014$Lake.ID==lake&doc2014$Depth.Class=='Hypo'&doc2014$flag==0,]
doc2014$Date.Sample<-strftime(strptime(doc2014$Date.Sample,'%m/%d/%Y'),'%Y-%m-%d')
doc2014<-aggregate(doc2014$DOC_mgL,by=list(doc2014$Date.Sample),FUN=mean)
colnames(doc2014)<-c('datetime','hypo_doc')
data<-merge(data,doc2014,all.x=T)
data$hypo_docInt<-data$hypo_doc #interpolating hypo DOC and DIC for entrainment load estimates 
data$hypo_docInt[1]<-data$hypo_docInt[min(which(!is.na(data$hypo_docInt)))]
data$hypo_docInt[length(data$hypo_docInt)]<-data$hypo_docInt[max(which(!is.na(data$hypo_docInt)))]
data$hypo_docInt<-approx(data$timeStep,data$hypo_docInt,data$timeStep)$y
data$hypo_dicInt<-data$hypo_dic #interpolating hypo DOC and DIC for entrainment load estimates 
data$hypo_dicInt[1]<-data$hypo_dicInt[min(which(!is.na(data$hypo_dicInt)))]
data$hypo_dicInt[length(data$hypo_dicInt)]<-data$hypo_dicInt[max(which(!is.na(data$hypo_dicInt)))]
data$hypo_dicInt<-approx(data$timeStep,data$hypo_dicInt,data$timeStep)$y

keep<-c('data')
toRm<-ls()
toRm<-toRm[toRm!=keep]
rm(list=toRm)




