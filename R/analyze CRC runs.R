# JAZ 2016-12-11 

dir<-'/Users/jzwart/LakeCarbonEnKF/Data/CRC/'
files<-list.files(dir)

models<-data.frame()
for(i in 1:length(files)){
  load(file.path(dir,files[i]))
  curOut$paramsRemoved<-0
  curOut$paramsRemoved<-curOut$paramsRemoved+length(grep('one',curOut$model))+length(grep('noMM',curOut$model))+
    length(grep('noTemp',curOut$model))
  models<-rbind(models,curOut)
}

windows()
plot(models$DOCrmse~models$paramsRemoved)

plot(models$obsDOCrmse~models$freq) # should be no relationship 
plot(models$obsDOCrmse~models$obs) # should be positive relationship 
plot(models$obsDOCrmse~models$reps) # should be a negative relationship 

plot(models$DICrmse~models$paramsRemoved)
points(models$obsDICrmse~models$paramsRemoved,col='red')

# when is DIC rmse and/or DOC rmse better than observed rmse? 

models$mdfDOCbetter<-as.numeric(models$DOCrmse<models$obsDOCrmse) # 1 means yes mdf is closer to truth 
models$mdfDICbetter<-as.numeric(models$DICrmse<models$obsDICrmse) # 1 means yes mdf is closer to truth 

models$mdfMinusObsDOC<-models$DOCrmse-models$obsDOCrmse
models$mdfMinusObsDIC<-models$DICrmse-models$obsDICrmse

models$obsMinusMdfDOC<-models$obsDOCrmse-models$DOCrmse
models$obsMinusMdfDIC<-models$obsDICrmse-models$DICrmse


sum(models$mdfDOCbetter)
sum(models$mdfDICbetter)
hist(models$mdfMinusObsDOC)
hist(models$mdfMinusObsDIC)

mean(models$mdfMinusObsDOC)
mean(models$mdfMinusObsDIC)

summary(lm(models$mdfMinusObsDOC~models$obs*models$freq*models$reps))

plot(models$mdfMinusObsDOC~models$paramsRemoved)
windows()
plot(models$mdfMinusObsDOC[models$obs==1&models$reps==3&models$freq==8]~
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==8],pch=16,cex=2)
text(models$mdfMinusObsDOC[models$obs==1&models$reps==3&models$freq==8]~
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==8],
     labels=models$model[models$obs==1&models$reps==3&models$freq==8],pos=1)

plot(models$mdfMinusObsDOC[models$obs==1&models$reps==3&models$freq==15]~
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==15],pch=16,cex=2)
text(models$mdfMinusObsDOC[models$obs==1&models$reps==3&models$freq==15]~
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==15],
     labels=models$model[models$obs==1&models$reps==3&models$freq==15],pos=1)

plot(models$mdfMinusObsDOC[models$obs==1&models$reps==3&models$freq==1]~
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==1],pch=16,cex=2)
text(models$mdfMinusObsDOC[models$obs==1&models$reps==3&models$freq==1]~
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==1],
     labels=models$model[models$obs==1&models$reps==3&models$freq==1],pos=1)

plot(models$mdfMinusObsDOC[models$obs==.5&models$reps==9&models$freq==1]~
       models$paramsRemoved[models$obs==.5&models$reps==9&models$freq==1],pch=16,cex=2)
text(models$mdfMinusObsDOC[models$obs==.5&models$reps==9&models$freq==1]~
       models$paramsRemoved[models$obs==.5&models$reps==9&models$freq==1],
     labels=models$model[models$obs==.5&models$reps==9&models$freq==1],pos=1)

plot(models$mdfMinusObsDOC[models$obs==.5&models$reps==9&models$freq==57]~
       models$paramsRemoved[models$obs==.5&models$reps==9&models$freq==57],pch=16,cex=2)
text(models$mdfMinusObsDOC[models$obs==.5&models$reps==9&models$freq==57]~
       models$paramsRemoved[models$obs==.5&models$reps==9&models$freq==57],
     labels=models$model[models$obs==.5&models$reps==9&models$freq==57],pos=1)

# impact of frequency of sampling on MDF rmse 
plot(models$mdfMinusObsDOC[models$obs==1&models$reps==3]~
       models$freq[models$obs==1&models$reps==3],pch=16,cex=2) # more frequent sampling makes models closer to truth over obs 

plot(models$DOCrmse[models$obs==1&models$reps==3]~
       models$freq[models$obs==1&models$reps==3],pch=16,cex=2) # more frequent sampling makes models closer to truth  

plot(models$obsDOCrmse[models$obs==1&models$reps==3]~
       models$freq[models$obs==1&models$reps==3],pch=16,cex=2) # more frequent sampling makes obs further from truth 
plot(models$obsDOCrmse[models$obs==1&models$reps==9]~
       models$freq[models$obs==1&models$reps==9],pch=16,cex=2) # more frequent sampling makes obs further from truth across all replicates 


plot(models$mdfMinusObsDIC[models$obs==1&models$reps==3]~
       models$freq[models$obs==1&models$reps==3],pch=16,cex=2) # more frequent sampling makes models closer to truth over data 

plot(models$DICrmse[models$obs==1&models$reps==3]~
       models$freq[models$obs==1&models$reps==3],pch=16,cex=2) # more frequent sampling makes models closer to truth  

plot(models$obsDICrmse[models$obs==1&models$reps==3]~
       models$freq[models$obs==1&models$reps==3],pch=16,cex=2) # more frequent sampling makes obs further from truth
plot(models$obsDICrmse[models$obs==1&models$reps==5]~
       models$freq[models$obs==1&models$reps==5],pch=16,cex=2) # more frequent sampling makes obs further from truth across all replicates 

plot(models$DOCrmse[models$obs==1&models$reps==3]~
       models$freq[models$obs==1&models$reps==3],pch=16,cex=2,ylim=c(0.1,0.9)) # more frequent sampling makes models closer to truth  
points(models$obsDOCrmse[models$obs==1&models$reps==3]~
       models$freq[models$obs==1&models$reps==3],pch=16,cex=2,col='red')
plot(models$DOCrmse[models$obs==1&models$reps==3&models$model=='twoDOCpools_MM_Temp']~
       models$freq[models$obs==1&models$reps==3&models$model=='twoDOCpools_MM_Temp'],pch=16,cex=2,ylim=c(0.1,0.9)) # more frequent sampling makes models closer to truth  
points(models$obsDOCrmse[models$obs==1&models$reps==3&models$model=='twoDOCpools_MM_Temp']~
         models$freq[models$obs==1&models$reps==3&models$model=='twoDOCpools_MM_Temp'],pch=16,cex=2,col='red')
plot(models$mdfMinusObsDOC[models$obs==1&models$reps==3&models$model=='twoDOCpools_MM_Temp']~
       models$freq[models$obs==1&models$reps==3&models$model=='twoDOCpools_MM_Temp'],pch=16,cex=2)

plot(models$DICrmse[models$obs==1&models$reps==3]~
       models$freq[models$obs==1&models$reps==3],pch=16,cex=2,ylim=c(0.01,0.04)) # more frequent sampling makes models closer to truth  
points(models$obsDICrmse[models$obs==1&models$reps==3]~
         models$freq[models$obs==1&models$reps==3],pch=16,cex=2,col='red')

# impact of observation error on MDF rmse 
plot(models$mdfMinusObsDOC[models$freq==8&models$reps==3]~
       models$obs[models$freq==8&models$reps==3],pch=16,cex=2) # increased observation error makes models closer to truth over data 
plot(models$DOCrmse[models$freq==8&models$reps==3]~
       models$obs[models$freq==8&models$reps==3],pch=16,cex=2) # increased observation error makes models further from truth  
plot(models$obsDOCrmse[models$freq==8&models$reps==3]~
       models$obs[models$freq==8&models$reps==3],pch=16,cex=2) # increased observation error makes obs further from truth 

plot(models$mdfMinusObsDIC[models$freq==8&models$reps==3]~
       models$obs[models$freq==8&models$reps==3],pch=16,cex=2) 

# impact of replicates on MDF rmse 
plot(models$mdfMinusObsDOC[models$freq==8&models$obs==1]~
       models$reps[models$freq==8&models$obs==1],pch=16,cex=2) # increased replicates makes models less important when trying to capture truth 
plot(models$DOCrmse[models$freq==8&models$obs==1]~
       models$reps[models$freq==8&models$obs==1],pch=16,cex=2) # increased reps decreases model rmse 
plot(models$obsDOCrmse[models$freq==8&models$obs==1]~
       models$reps[models$freq==8&models$obs==1],pch=16,cex=2) # increased reps decreases obs rmse 


# plotting on 3 dimensions 
library(plot3D)
library(scatterplot3d)
curMod<-models[models$model=='twoDOCpools_MM_Temp',]
x=curMod$obs
y=curMod$reps
z=curMod$mdfMinusObsDOC
scatterplot3d(x,y,z,pch = 16)

summary(lm(models$mdfMinusObsDOC~models$obs*models$freq*models$reps))

plot(models$mdfMinusObsDIC[models$obs==1&models$reps==3&models$freq==8]~
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==8],pch=16,cex=2)
text(models$mdfMinusObsDIC[models$obs==1&models$reps==3&models$freq==8]~
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==8],
     labels=models$model[models$obs==1&models$reps==3&models$freq==8],pos=1)
abline(0,0,lwd=2,lty=2)

plot(models$mdfMinusObsDIC[models$obs==1&models$reps==1&models$freq==8]~
       models$paramsRemoved[models$obs==1&models$reps==1&models$freq==8],pch=16,cex=2)
text(models$mdfMinusObsDIC[models$obs==1&models$reps==1&models$freq==8]~
       models$paramsRemoved[models$obs==1&models$reps==1&models$freq==8],
     labels=models$model[models$obs==1&models$reps==1&models$freq==8],pos=1)
abline(0,0,lwd=2,lty=2)

plot(models$mdfMinusObsDIC[models$obs==1&models$reps==3&models$freq==57]~
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==57],pch=16,cex=2)
text(models$mdfMinusObsDIC[models$obs==1&models$reps==3&models$freq==57]~
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==57],
     labels=models$model[models$obs==1&models$reps==3&models$freq==57],pos=1)
abline(0,0,lwd=2,lty=2)

plot(models$mdfMinusObsDIC[models$obs==1&models$reps==9&models$freq==1]~ ## you can take 9 samples every day and you'll still be closer to truth with a true process; and not that far off with a slight wrong model 
       models$paramsRemoved[models$obs==1&models$reps==9&models$freq==1],pch=16,cex=2)
text(models$mdfMinusObsDIC[models$obs==1&models$reps==9&models$freq==1]~
       models$paramsRemoved[models$obs==1&models$reps==9&models$freq==1],
     labels=models$model[models$obs==1&models$reps==9&models$freq==1],pos=1)
abline(0,0,lwd=2,lty=2)

plot(models$mdfMinusObsDIC[models$obs==.5&models$reps==9&models$freq==1]~ # even 
       models$paramsRemoved[models$obs==.5&models$reps==9&models$freq==1],pch=16,cex=2)
text(models$mdfMinusObsDIC[models$obs==.5&models$reps==9&models$freq==1]~
       models$paramsRemoved[models$obs==.5&models$reps==9&models$freq==1],
     labels=models$model[models$obs==.5&models$reps==9&models$freq==1],pos=1)
abline(0,0,lwd=2,lty=2)

plot(models$mdfMinusObsDIC[models$obs==1&models$reps==3&models$freq==15]~ ## this is probably what most people sample at
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==15],pch=16,cex=2)
text(models$mdfMinusObsDIC[models$obs==1&models$reps==3&models$freq==15]~
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==15],
     labels=models$model[models$obs==1&models$reps==3&models$freq==15],pos=1)
abline(0,0,lwd=2,lty=2)

plot(models$mdfMinusObsDIC[models$obs==1&models$reps==3&models$freq==29]~
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==29],pch=16,cex=2)
text(models$mdfMinusObsDIC[models$obs==1&models$reps==3&models$freq==29]~
       models$paramsRemoved[models$obs==1&models$reps==3&models$freq==29],
     labels=models$model[models$obs==1&models$reps==3&models$freq==29],pos=1)
abline(0,0,lwd=2,lty=2)

windows()
colors<-heat.colors(n=10)
plot(models$DOCrmse[models$obs==1]~models$freq[models$obs==1],cex=0)
for(i in 1:length(unique(models$reps))){
  lines(models$DOCrmse[models$obs==1&models$reps==unique(models$reps)[i]]~models$freq[models$obs==1&models$reps==unique(models$reps)[i]],ylab='DOC RMSE',cex.lab=1.5,pch=16,
        type='o',col=colors[i])
}
legend('topleft',col=colors[1:9],fill=colors[1:9],legend = unique(models$reps))

windows()
colors<-heat.colors(n=10)
plot(models$DICrmse[models$obs==1]~models$freq[models$obs==1],cex=0)
for(i in 1:length(unique(models$reps))){
  lines(models$DICrmse[models$obs==1&models$reps==unique(models$reps)[i]]~models$freq[models$obs==1&models$reps==unique(models$reps)[i]],ylab='DIC RMSE',cex.lab=1.5,pch=16,
        type='o',col=colors[i])
}
legend('topleft',col=colors[1:9],fill=colors[1:9],legend = unique(models$reps))


# what explains obsDOCrmse? 
summary(lm(models$obsDOCrmse~models$obs*models$freq*models$reps))

summary(lm(models$obsDOCrmse[models$obs==1]~models$freq[models$obs==1]*models$reps[models$obs==1]))
summary(lm(models$obsDICrmse[models$obs==1]~models$freq[models$obs==1]*models$reps[models$obs==1]))

summary(lm(models$obsMinusMdfDOC[models$obs==1]~models$freq[models$obs==1]*models$reps[models$obs==1]*models$model[models$obs==1]))
summary(lm(models$obsMinusMdfDIC[models$obs==1]~models$freq[models$obs==1]*models$reps[models$obs==1]*models$model[models$obs==1]))

summary(lm(models$obsMinusMdfDOC[models$obs==1]~models$freq[models$obs==1]+models$reps[models$obs==1]+factor(models$model[models$obs==1])))
summary(lm(models$obsMinusMdfDIC[models$obs==1]~models$freq[models$obs==1]+models$reps[models$obs==1]+factor(models$model[models$obs==1])))

summary(lm(models$obsMinusMdfDOC~models$freq+models$reps+factor(models$model)))
summary(lm(models$obsMinusMdfDIC~models$freq+models$reps+factor(models$model)))

# changing reference model to true data process 'twoDOCpools_MM_Temp' 
table(models$model)
models$model<-relevel(models$model,ref = 'twoDOCpools_MM_Temp')
table(models$model)

summary(lm(models$obsMinusMdfDOC[models$obs==1]~models$freq[models$obs==1]+models$reps[models$obs==1]+factor(models$model[models$obs==1])+
             models$freq[models$obs==1]*models$reps[models$obs==1]))
summary(lm(models$obsMinusMdfDIC[models$obs==1]~models$freq[models$obs==1]+models$reps[models$obs==1]+factor(models$model[models$obs==1])+
             models$freq[models$obs==1]*models$reps[models$obs==1]))

summary(lm(models$obsMinusMdfDOC~models$freq+models$reps+factor(models$model)))
summary(lm(models$obsMinusMdfDIC~models$freq+models$reps+factor(models$model)))

summary(lm(models$obsMinusMdfDOC[models$obs==1]~models$freq[models$obs==1]*models$reps[models$obs==1]*factor(models$model[models$obs==1])))
summary(lm(models$obsMinusMdfDIC[models$obs==1]~models$freq[models$obs==1]*models$reps[models$obs==1]*factor(models$model[models$obs==1])))

summary(lm(models$obsMinusMdfDOC[models$obs==1]~factor(models$model[models$obs==1])))
summary(lm(models$obsMinusMdfDIC[models$obs==1]~factor(models$model[models$obs==1])))


models$mdfDICbetter[models$reps<=3&models$obs==1]


# some other things to include in the CRC runs - Bias, Precision, CV. also compare model to obs in addition to obs & model to true? 
plot(models$DOCbias~models$obs)
plot(models$DOCprecision~models$obs)



# Fig 4. DOC and pCO2 RMSE for all model runs 
library(vioplot)

windows()
plot(models$obsDOCrmse~models$model)
plot(models$DOCrmse~models$model)
plot(models$mdfMinusObsDOC~models$model)

plot(models$mdfMinusObsDIC[models$obs==1]~models$model[models$obs==1])
abline(0,0,lty=2)

library(ggplot2) # ggplot is kickass!! first time using it on 2017-01-19 

#assigning models a number in MS table so it's easier to organize 
unique(models$model)
models$modelNumber<-NA
models$modelNumber<-ifelse(models$model=='oneDOCpool_MM_noTemp',5,models$modelNumber)
models$modelNumber<-ifelse(models$model=='oneDOCpool_MM_Temp',6,models$modelNumber)
models$modelNumber<-ifelse(models$model=='oneDOCpool_noMM_noTemp',3,models$modelNumber)
models$modelNumber<-ifelse(models$model=='oneDOCpool_noMM_Temp',4,models$modelNumber)
models$modelNumber<-ifelse(models$model=='twoDOCpools_MM_Temp',10,models$modelNumber)
models$modelNumber<-ifelse(models$model=='twoDOCpools_MM_noTemp',9,models$modelNumber)
models$modelNumber<-ifelse(models$model=='twoDOCpools_noMM_noTemp',7,models$modelNumber)
models$modelNumber<-ifelse(models$model=='twoDOCpools_noMM_Temp',8,models$modelNumber)

models<-models[sort.list(models$modelNumber),]

models<-models[models$reps<10&models$freq<60,]

png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig4_DOC_RMSE.png', 
    res=300, width=7, height=7, units = 'in')
# cex=2
# bp.at<-c(2,4,6,8,10,12,14,16)
# bp.cols<-c('white')
# par(mar=c(5,6,4,2))
# ylim=c(-1.1,0.1)
# boxplot(models$mdfMinusObsDOC[models$obs==1]~models$model[models$obs==1],outline=T,ylim=ylim,xlab='Model',
#         lwd=6, pars=list(boxfill=bp.cols, medcol='black', staplecol='black'),at=bp.at,xaxt='n',
#         boxwex=1,cex.axis=cex,ylab=expression(DOC~RMSE~(Modeled~-~Observed)),cex.lab=cex)
# abline(0,0,lty=2,lwd=4,col='black')
# box('plot',lwd=2)
# axis(1, at=bp.at, labels=c(3,4,5,6,7,8,9,10),cex.axis=cex)
p <- ggplot(models[models$obs==1,], aes(x=as.character(modelNumber), y=obsMinusMdfDOC)) + ylim(low=-3,high=2) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Model), y = expression(Obs~DOC~RMSE~-~DA~DOC~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(unique(models$modelNumber)))+
  scale_fill_brewer(palette="Greys")+theme(legend.position="none",axis.text=element_text(size=16),
                                           axis.title=element_text(size=16,face="bold"),
                                           axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
dev.off()

####################practicing 

# Basic violin plot

# 
# p+stat_summary(fun.y=median, geom="point", size=2, color="red")+geom_boxplot(width=0.1)
# p+stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="red")
# p+ geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)


png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig4_DIC_RMSE.png', 
    res=300, width=7, height=7, units = 'in')
# cex=2
# bp.at<-c(2,4,6,8,10,12,14,16)
# bp.cols<-c('white')
# par(mar=c(5,6,4,2))
# boxplot(models$obsMinusMdfDIC[models$obs==1]~models$model[models$obs==1],outline=T,xlab='Model',
#         lwd=6, pars=list(boxfill=bp.cols, medcol='black', staplecol='black'),at=bp.at,xaxt='n',
#         boxwex=1,cex.axis=cex,ylab=expression(CO[2]~RMSE~(Modeled~-~Observed)),cex.lab=cex)
# abline(0,0,lty=2,lwd=4,col='black')
# box('plot',lwd=2)
# axis(1, at=bp.at, labels=c(3,4,5,6,7,8,9,10),cex.axis=cex)
p <- ggplot(models[models$obs==1,], aes(x=as.character(modelNumber), y=obsMinusMdfDIC)) + ylim(low=-.3,high=.2) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Model), y = expression(Obs~CO[2]~RMSE~-~DA~CO[2]~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(unique(models$modelNumber)))+
  scale_fill_brewer(palette="Greys")+theme(legend.position="none",axis.text=element_text(size=16),
                                           axis.title=element_text(size=16,face="bold"),
                                           axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
dev.off()

# plot(models$mdfobsDICrmse~models$model)
# plot(models$mdfobsDICrmse~models$obsMinusMdfDIC)
# plot(models$mdfobsDOCrmse~models$obsMinusMdfDOC)
# plot(models$obsMinusMdfDOC[models$obs==1]~models$freq[models$obs==1])
# plot(models$obsMinusMdfDOC[models$obs==1]~models$reps[models$obs==1])
# plot(models$obsMinusMdfDIC[models$obs==1]~models$reps[models$obs==1])

png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig4_DIC_repsInfluence.png', 
    res=300, width=7, height=7, units = 'in')
# cex=2
# par(mar=c(5,6,4,2))
# plot(models$obsMinusMdfDIC[models$obs==1]~models$reps[models$obs==1],xlab='Replicate Samples (n)',
#         lwd=6,cex.axis=cex,ylab=expression(CO[2]~RMSE~(Modeled~-~Observed)),cex.lab=cex)
# abline(0,0,lty=2,lwd=4,col='black')
# box('plot',lwd=2)
p <- ggplot(models[models$obs==1,], aes(x=as.character(reps), y=obsMinusMdfDIC)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Replicate~Samples~(n)), y = expression(Obs~CO[2]~RMSE~-~DA~CO[2]~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(c(1,2,3,4,5,6)))+
    theme(legend.position="none",axis.text=element_text(size=16),
                                           axis.title=element_text(size=16,face="bold"),
                                           axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig4_DOC_repsInfluence.png', 
    res=300, width=7, height=7, units = 'in')
# cex=2
# par(mar=c(5,6,4,2))
# ylim=c(-1.1,0.1)
# plot(models$obsMinusMdfDOC[models$obs==1]~models$reps[models$obs==1],xlab='Replicate Samples (n)',
#      lwd=6,cex.axis=cex,ylab=expression(DOC~RMSE~(Modeled~-~Observed)),cex.lab=cex,ylim=ylim)
# abline(0,0,lty=2,lwd=4,col='black')
# box('plot',lwd=2)
p <- ggplot(models[models$obs==1,], aes(x=as.character(reps), y=obsMinusMdfDOC)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Replicate~Samples~(n)), 
                                                            y = expression(Obs~DOC~RMSE~-~DA~DOC~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(c(1,2,3,4,5,6)))+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig4_DIC_freqInfluence.png', 
    res=300, width=7, height=7, units = 'in')
# cex=2
# par(mar=c(5,6,4,2))
# plot(models$obsMinusMdfDIC[models$obs==1]~models$freq[models$obs==1],xlab='Sampling Frequency (days)',
#      lwd=6,cex.axis=cex,ylab=expression(CO[2]~RMSE~(Modeled~-~Observed)),cex.lab=cex)
# abline(0,0,lty=2,lwd=4,col='black')
# box('plot',lwd=2)
p <- ggplot(models[models$obs==1,], aes(x=as.character(freq), y=obsMinusMdfDIC)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Sampling~Interval~(days)), 
                                                            y = expression(Obs~CO[2]~RMSE~-~DA~CO[2]~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(c(1,7,14,21,28,35)))+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Model Data Fusion/Figures/Fig4_DOC_freqInfluence.png', 
    res=300, width=7, height=7, units = 'in')
# cex=2
# par(mar=c(5,6,4,2))
# ylim=c(-1.1,0.1)
# plot(models$obsMinusMdfDOC[models$obs==1]~models$freq[models$obs==1],xlab='Sampling Frequency (days)',
#      lwd=6,cex.axis=cex,ylab=expression(DOC~RMSE~(Modeled~-~Observed)),cex.lab=cex,ylim=ylim)
# abline(0,0,lty=2,lwd=4,col='black')
# box('plot',lwd=2)
p <- ggplot(models[models$obs==1,], aes(x=as.character(freq), y=obsMinusMdfDOC)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Sampling~Interval~(days)),
                                                            y = expression(Obs~DOC~RMSE~-~DA~DOC~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(c(1,7,14,21,28,35)))+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
dev.off()

# mdf compared to obs rmse given our sampling frequency and replicates 
windows()
p <- ggplot(models[models$obs==1&models$freq==7&models$reps==2,], aes(x=as.character(modelNumber), y=mdfobsDICrmse)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Model), y = expression(Mod~CO[2]~Obs~CO[2]~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(unique(models$modelNumber)))+
  scale_fill_brewer(palette="Greys")+theme(legend.position="none",axis.text=element_text(size=16),
                                           axis.title=element_text(size=16,face="bold"),
                                           axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
windows()
p <- ggplot(models[models$obs==1&models$freq==7&models$reps==2,], aes(x=as.character(modelNumber), y=mdfobsDOCrmse)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Model), y = expression(Mod~DOC~Obs~DOC~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(unique(models$modelNumber)))+
  scale_fill_brewer(palette="Greys")+theme(legend.position="none",axis.text=element_text(size=16),
                                           axis.title=element_text(size=16,face="bold"),
                                           axis.line = element_line(colour = "black",size = 1, linetype = "solid"))

# sampling frequency effect on MDF RMSE compared to truth 
windows()
p <- ggplot(models[models$obs==1,], aes(x=as.character(freq), y=DOCrmse)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Sampling~Frequency~(days)), y = expression(MDF~DOC~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(c(1,7,14,21,28,35,63)))+
  scale_fill_brewer(palette="Greys")+theme(legend.position="none",axis.text=element_text(size=16),
                                           axis.title=element_text(size=16,face="bold"),
                                           axis.line = element_line(colour = "black",size = 1, linetype = "solid"))

windows()
p <- ggplot(models[models$obs==1,], aes(x=as.character(freq), y=DICrmse)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Sampling~Frequency~(days)), y = expression(MDF~CO[2]~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(c(1,8,15,22,29,36,43,50,57)))+
  scale_fill_brewer(palette="Greys")+theme(legend.position="none",axis.text=element_text(size=16),
                                           axis.title=element_text(size=16,face="bold"),
                                           axis.line = element_line(colour = "black",size = 1, linetype = "solid"))

# sample replicates effect on MDF RMSE compared to truth 
windows()
p <- ggplot(models[models$obs==1,], aes(x=as.character(reps), y=DOCrmse)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Replicate~Samples~(n)), y = expression(MDF~DOC~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(colour = "black",size = 1, linetype = "solid"))

windows()
p <- ggplot(models[models$obs==1,], aes(x=as.character(reps), y=DICrmse)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Replicate~Samples~(n)), y = expression(MDF~CO[2]~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(colour = "black",size = 1, linetype = "solid"))


# sample replicates effect on MDF CV  
windows()
p <- ggplot(models[models$obs==1,], aes(x=as.character(reps), y=mdfDOCcv)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Replicate~Samples~(n)), y = expression(MDF~DOC~CV))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(colour = "black",size = 1, linetype = "solid"))

windows()
p <- ggplot(models[models$obs==1,], aes(x=as.character(reps), y=mdfDICcv)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Replicate~Samples~(n)), y = expression(MDF~CO[2]~CV))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(colour = "black",size = 1, linetype = "solid"))

# sampling frequency effect on MDF CV  
windows()
p <- ggplot(models[models$obs==1,], aes(x=as.character(freq), y=mdfDOCcv)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Sampling~Frequency~(days)), y = expression(MDF~DOC~CV))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(c(1,8,15,22,29,36,43,50,57)))+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(colour = "black",size = 1, linetype = "solid"))

windows()
p <- ggplot(models[models$obs==1,], aes(x=as.character(freq), y=mdfDICcv)) + 
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Sampling~Frequency~(days)), y = expression(MDF~CO[2]~CV))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(c(1,8,15,22,29,36,43,50,57)))+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(colour = "black",size = 1, linetype = "solid"))




