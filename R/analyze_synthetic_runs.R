# JAZ 2016-12-11

dir<-'Results/Results_20170319/'
files<-list.files(dir)

models<-data.frame()
for(i in 1:length(files)){
  load(file.path(dir,files[i]))
  curOut$paramsRemoved<-0
  curOut$paramsRemoved<-curOut$paramsRemoved+length(grep('one',curOut$model))+length(grep('noMM',curOut$model))+
    length(grep('noTemp',curOut$model))
  models<-rbind(models,curOut)
}

models$mdfDOCbetter<-as.numeric(models$DOCrmse<models$obsDOCrmse) # 1 means yes mdf is closer to truth
models$mdfDICbetter<-as.numeric(models$DICrmse<models$obsDICrmse) # 1 means yes mdf is closer to truth

models$mdfMinusObsDOC<-models$DOCrmse-models$obsDOCrmse
models$mdfMinusObsDIC<-models$DICrmse-models$obsDICrmse

models$obsMinusMdfDOC<-models$obsDOCrmse-models$DOCrmse
models$obsMinusMdfDIC<-models$obsDICrmse-models$DICrmse



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

library(ggplot2) # ggplot is kickass!! first time using it on 2017-01-19

#assigning models a number in MS table so it's easier to organize
unique(models$model)
models$modelNumber<-NA
models$modelNumber<-ifelse(models$model=='oneDOCpool_MM_noTemp',4,models$modelNumber)
models$modelNumber<-ifelse(models$model=='oneDOCpool_MM_Temp',5,models$modelNumber)
models$modelNumber<-ifelse(models$model=='oneDOCpool_noMM_noTemp',2,models$modelNumber)
models$modelNumber<-ifelse(models$model=='oneDOCpool_noMM_Temp',3,models$modelNumber)
models$modelNumber<-ifelse(models$model=='twoDOCpools_MM_Temp',9,models$modelNumber)
models$modelNumber<-ifelse(models$model=='twoDOCpools_MM_noTemp',8,models$modelNumber)
models$modelNumber<-ifelse(models$model=='twoDOCpools_noMM_noTemp',6,models$modelNumber)
models$modelNumber<-ifelse(models$model=='twoDOCpools_noMM_Temp',7,models$modelNumber)

models<-models[sort.list(models$modelNumber),]

models<-models[models$reps<10&models$freq<60,]

ylim_doc = c(-.2,1.3)
ylim_dic = c(-0.05, 0.075)

png('/Users/jzwart/LakeCarbonEnKF/Figures/Fig4_DOC_RMSE.png',
    res=300, width=7, height=7, units = 'in')

p <- ggplot(models[models$obs==1,], aes(x=as.character(modelNumber), y=obsMinusMdfDOC)) + ylim(low=ylim_doc[1],high=ylim_doc[2]) +
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Model), y = expression(Obs~DOC~RMSE~-~DA~DOC~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(unique(models$modelNumber)))+
  scale_fill_brewer(palette="Greys")+theme(legend.position="none",axis.text=element_text(size=16),
                                           axis.title=element_text(size=16,face="bold"),
                                           axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
dev.off()


png('/Users/jzwart/LakeCarbonEnKF/Figures/Fig4_DIC_RMSE.png',
    res=300, width=7, height=7, units = 'in')

p <- ggplot(models[models$obs==1,], aes(x=as.character(modelNumber), y=obsMinusMdfDIC)) + ylim(low=ylim_dic[1],high=ylim_dic[2]) +
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Model), y = expression(Obs~CO[2]~RMSE~-~DA~CO[2]~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(unique(models$modelNumber)))+
  scale_fill_brewer(palette="Greys")+theme(legend.position="none",axis.text=element_text(size=16),
                                           axis.title=element_text(size=16,face="bold"),
                                           axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
dev.off()

png('/Users/jzwart/LakeCarbonEnKF/Figures/Fig4_DIC_repsInfluence.png',
    res=300, width=7, height=7, units = 'in')

p <- ggplot(models[models$obs==1,], aes(x=as.character(reps), y=obsMinusMdfDIC)) + ylim(low=ylim_dic[1],high=ylim_dic[2]) +
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Replicate~Samples~(n)), y = expression(Obs~CO[2]~RMSE~-~DA~CO[2]~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(c(1,2,3,4,5,6)))+
    theme(legend.position="none",axis.text=element_text(size=16),
                                           axis.title=element_text(size=16,face="bold"),
                                           axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
dev.off()

png('/Users/jzwart/LakeCarbonEnKF/Figures/Fig4_DOC_repsInfluence.png',
    res=300, width=7, height=7, units = 'in')

p <- ggplot(models[models$obs==1,], aes(x=as.character(reps), y=obsMinusMdfDOC)) + ylim(low=ylim_doc[1],high=ylim_doc[2]) +
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Replicate~Samples~(n)),
                                                            y = expression(Obs~DOC~RMSE~-~DA~DOC~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(c(1,2,3,4,5,6)))+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
dev.off()

png('/Users/jzwart/LakeCarbonEnKF/Figures/Fig4_DIC_freqInfluence.png',
    res=300, width=7, height=7, units = 'in')

p <- ggplot(models[models$obs==1,], aes(x=as.character(freq), y=obsMinusMdfDIC)) + ylim(low=ylim_dic[1],high=ylim_dic[2]) +
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Sampling~Interval~(days)),
                                                            y = expression(Obs~CO[2]~RMSE~-~DA~CO[2]~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(c(1,7,14,21,28,35)))+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
dev.off()

png('/Users/jzwart/LakeCarbonEnKF/Figures/Fig4_DOC_freqInfluence.png',
    res=300, width=7, height=7, units = 'in')

p <- ggplot(models[models$obs==1,], aes(x=as.character(freq), y=obsMinusMdfDOC)) + ylim(low=ylim_doc[1],high=ylim_doc[2]) +
  geom_violin(trim = F,fill='grey')+ theme_classic() + labs(x=expression(Sampling~Interval~(days)),
                                                            y = expression(Obs~DOC~RMSE~-~DA~DOC~RMSE))
p + geom_jitter(shape=16, position=position_jitter(0.2),col='black')+geom_hline(yintercept=0,lty=2,lwd=1.5)+
  scale_x_discrete(limits=as.character(c(1,7,14,21,28,35)))+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
dev.off()



