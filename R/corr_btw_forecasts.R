


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

out1 # CO2 / DOC model
out # DIC / CO2 / DOC model

Y1 = out1$Y
Y = out$Y

# co2
co2_doc_mod = apply(Y1[6,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12

dic_co2_doc_mod = apply(Y[7,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12


cor(co2_doc_mod, dic_co2_doc_mod)
sqrt(mean((co2_doc_mod - dic_co2_doc_mod)^2))


ylim = c(0.2,0.7)
xlim = ylim
plot(co2_doc_mod~ dic_co2_doc_mod, xlim = xlim, ylim =ylim, pch=16, ylab= 'CO2 Only Model', xlab='Carbonate Model')
abline(0,1, lty = 2, lwd = 2)

windows()
old = Y1[6,1,,1]/data2$epiVol*12
new = Y[7,1,,1]/data2$epiVol*12
ylim = c(0.2,0.7)
xlim = ylim
plot(new~old, pch =16, col = t_col('black', 80), ylim = ylim, xlim= xlim)
for(i in 2:nEn){
  old = Y1[6,1,,i]/data2$epiVol*12
  new = Y[7,1,,i]/data2$epiVol*12

  points(new~old, pch =16, col = t_col('black', 80), ylim = ylim, xlim= xlim)
}
abline(0,1, lwd=2, lty=2, col ='red')

# AIC(logLik(lm(co2_doc_mod~dic_co2_doc_mod)))

# doc
co2_doc_mod = apply(Y1[9,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12

dic_co2_doc_mod = apply(Y[10,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12


cor(co2_doc_mod, dic_co2_doc_mod)
sqrt(mean((co2_doc_mod - dic_co2_doc_mod)^2))

ylim = c(9,11)
xlim = ylim
plot(co2_doc_mod~ dic_co2_doc_mod, xlim = xlim, ylim =ylim, pch=16, ylab= 'CO2 Only Model', xlab='Carbonate Model')
abline(0,1, lty = 2, lwd = 2)


#co2
windows()
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels
leg = 0.05 # distance from corner for panel label

cex=4
cex.lab=2
cex.axis=2
lwd=3
ylim=c(0.5,3)
xlim=range(as.POSIXct(data2$datetime))
par(mar = c(5,6,1,1))
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


# doC
windows()
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels
leg = 0.05 # distance from corner for panel label

cex=4
cex.lab=2
cex.axis=2
lwd=3
ylim=c(4,12)
xlim=range(as.POSIXct(data2$datetime))
par(mar = c(5,6,1,1))
DOCout<-apply(Y[10,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
ylim=range(c(Y[10,1,,]/data2$epiVol*12,z$y[3,1,,1]/data2$epiVol*12+docPoolSD/data2$epiVol*12,z$y[3,1,,1]/data2$epiVol*12-docPoolSD/data2$epiVol*12),na.rm=T)
plot(DOCout~as.POSIXct(data2$datetime),type='l',ylim=ylim,lwd=lwd,cex.axis=cex.axis, col= t_col('red',80),
     ylab=expression(DOC~(mg~C~L^-1)),xlab='',cex.lab=cex.lab)
for(i in 1:nEn){
  lines(Y[10,1,,i]/data2$epiVol*12~as.POSIXct(data2$datetime),
        col=t_col('red', 80),ylab='',lwd=lwd)
}
lines(DOCout~as.POSIXct(data2$datetime),ylab='',lwd=3,col=t_col('red',80))
par(new=T)
DOCout<-apply(Y1[9,1,,]/data2$epiVol*12,MARGIN = 1,FUN=mean)
plot(DOCout~as.POSIXct(data2$datetime),type='l',ylim=ylim,lwd=lwd,cex.axis=cex.axis,col=t_col('blue',80),
     ylab=expression(DOC~(mg~C~L^-1)),xlab='',cex.lab=cex.lab)
for(i in 1:nEn){
  lines(Y1[9,1,,i]/data2$epiVol*12~as.POSIXct(data2$datetime),
        col=t_col('light blue',80),ylab='',lwd=lwd)
}
lines(DOCout~as.POSIXct(data2$datetime),ylab='',lwd=3,col=t_col('blue',80))

par(new=T)
plot(z$y[3,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12~as.POSIXct(data2$datetime[assim_obs==0]),cex=cex,
     ylim=ylim,col='black',pch=21,ylab = '',xlab='',cex.lab=cex.lab,cex.axis=cex.axis,xlim=xlim)
arrows(as.POSIXct(data2$datetime[assim_obs==0]),z$y[3,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12-docPoolSD[assim_obs==0]/data2$epiVol[assim_obs==0]*12,
       as.POSIXct(data2$datetime[assim_obs==0]),z$y[3,1,assim_obs==0,1]/data2$epiVol[assim_obs==0]*12+docPoolSD[assim_obs==0]/data2$epiVol[assim_obs==0]*12,
       code=3,length=0.1,angle=90,col='black',lwd=3)
par(new=T)
plot(z$y[3,1,assim_obs==1,1]/data2$epiVol[assim_obs==1]*12~as.POSIXct(data2$datetime[assim_obs==1]),cex=cex,
     ylim=ylim,col='black',pch=16,ylab='',xlab='',cex.lab=cex.lab,cex.axis=cex.axis,xlim=xlim)
arrows(as.POSIXct(data2$datetime[assim_obs==1]),z$y[3,1,assim_obs==1,1]/data2$epiVol[assim_obs==1]*12-docPoolSD[assim_obs==1]/data2$epiVol[assim_obs==1]*12,
       as.POSIXct(data2$datetime[assim_obs==1]),z$y[3,1,assim_obs==1,1]/data2$epiVol[assim_obs==1]*12+docPoolSD[assim_obs==1]/data2$epiVol[assim_obs==1]*12,
       code=3,length=0.1,angle=90,col='black',lwd=3)













# diff in pH
data2$pH_sensor_corr
data2$ph_diff = c(NA,diff(data2$pH_sensor_corr))
windows()
plot(data2$ph_diff, pch =16, ylab= 'Day-to-day change in pH', cex = 2, cex.lab =1.4, xlab= '')
abline(0,0, lty=2, lwd= 2)


# how much CO2 flux occurs with day-to-day change in pH?
# pH_vals = seq(5,8.5, by = .1)
# pH_range = 0.01

data2$CO2_diff <- NA

for(i in 2:nrow(data2)){
  sal = 0
  temp = data2$epiTemp[i]
  DIC_pool = mean(data2$dic, na.rm = T) # mol C; don't have data every day so taking average
  Vepi = mean(data2$epiVol,na.rm = T) # m^3
  today_pH = data2$pH_sensor_corr[i]
  yesterday_pH = data2$pH_sensor_corr[i-1]

  today_pH = AquaEnv::aquaenv(S = sal, # salinity
                            t = temp, # water temp
                            d = 0, # depth
                            lat = 46, # latitude
                            SumCO2 = DIC_pool/Vepi/rLakeAnalyzer::water.density(temp,sal),
                            pH = today_pH)

  today_CO2=today_pH$CO2[1]*water.density(temp,sal)*Vepi
  today_fracCO2=today_CO2/DIC_pool


  yesterday_pH = AquaEnv::aquaenv(S = sal, # salinity
                             t = temp, # water temp
                             d = 0, # depth
                             lat = 46, # latitude
                             SumCO2 = DIC_pool/Vepi/rLakeAnalyzer::water.density(temp,sal),
                             pH = yesterday_pH)

  yesterday_CO2=yesterday_pH$CO2[1]*water.density(temp,sal)*Vepi
  yesterday_fracCO2=yesterday_CO2/DIC_pool

  data2$CO2_diff[i] <- today_CO2-yesterday_CO2
}

plot(data2$CO2_diff/Vepi*12 ~ as.POSIXct(data2$datetime), pch=16, cex = 2, xlab = '',
     ylab = 'day-to-day CO2 change due to pH (mg C / L / day)')
abline(0,0, lty =2 , lwd= 2)

