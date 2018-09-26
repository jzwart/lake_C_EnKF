
out1 # CO2 / DOC model
out # DIC / CO2 / DOC model

Y1 = out1$Y
Y = out$Y

co2_doc_mod = apply(Y1[6,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12

dic_co2_doc_mod = apply(Y[7,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12


cor(co2_doc_mod, dic_co2_doc_mod)
sqrt(mean((co2_doc_mod - dic_co2_doc_mod)^2))

ylim = c(0.2,0.7)
xlim = ylim
plot(co2_doc_mod~ dic_co2_doc_mod, xlim = xlim, ylim =ylim, ylab= 'Old Model', xlab='New Model')
abline(0,1, lty = 2, lwd = 2)
# abline(lm(co2_doc_mod~dic_co2_doc_mod))


AIC(logLik(lm(co2_doc_mod~dic_co2_doc_mod)))


t_col <- function(color, percent = 50, name = NULL) {
  #	  color = color name
  #	percent = % transparency
  #	   name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  ## Save the color
  invisible(t.col)

}

windows()
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
CO2out<-apply(Y[7,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
ylim=range(c(Y[7,1,,]/data2$epiVol*12,z$y[2,1,,1]/data2$epiVol*12+co2PoolSD/data2$epiVol*12,z$y[2,1,,1]/data2$epiVol*12-co2PoolSD/data2$epiVol*12),na.rm=T)
plot(CO2out~as.POSIXct(data2$datetime),type='l',ylim=ylim,lwd=lwd,cex.axis=cex.axis, col= t_col('red',80),
     ylab=expression(CO[2]~(mg~C~L^-1)),xlab='',cex.lab=cex.lab)
for(i in 1:nEn){
  lines(Y[7,1,,i]/data2$epiVol*12~as.POSIXct(data2$datetime),
        col=t_col('red', 80),ylab='',lwd=lwd)
}
lines(CO2out~as.POSIXct(data2$datetime),ylab='',lwd=3,col=t_col('red',80))
par(new=T)
CO2out<-apply(Y1[6,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
plot(CO2out~as.POSIXct(data2$datetime),type='l',ylim=ylim,lwd=lwd,cex.axis=cex.axis,col=t_col('blue',80),
     ylab=expression(CO[2]~(mg~C~L^-1)),xlab='',cex.lab=cex.lab)
for(i in 1:nEn){
  lines(Y1[6,1,,i]/data2$epiVol*12~as.POSIXct(data2$datetime),
        col=t_col('light blue',80),ylab='',lwd=lwd)
}
lines(CO2out~as.POSIXct(data2$datetime),ylab='',lwd=3,col=t_col('blue',80))

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

