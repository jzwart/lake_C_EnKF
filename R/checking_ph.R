
# looking at ph in EL
library(dplyr)

minDate = '2014-05-01 00:00'
maxDate = '2014-07-10 00:00'

ph = dbTable(table = 'limno_profiles', lakeID = 'EL', minDate = minDate, maxDate = maxDate, dateFormat = '%Y-%m-%d %H:%M') %>%
  as_tibble() %>%
  dplyr::filter(depthBottom < 2, pH<7)


plot(ph$pH ~ ph$dateTimeSample)

sensordbTableList()

ph_sensor = sensordbTable(table = 'ysi_corr', lakeID = 'EL', minDate = minDate, maxDate = maxDate, dateFormat = '%Y-%m-%d %H:%M') %>%
  as_tibble() %>%
  dplyr::filter(!is.na(pH))

ph_sensor

plot(ph_sensor$pH ~ ph_sensor$dateTime, type = 'l')

windows()
plot(ph$pH ~ ph$dateTimeSample, pch = 16, ylim = c(4.5,6))
lines(ph_sensor$pH~ph_sensor$dateTime, col ='grey40')

# showing the time window of when we sample DIC / CO2 (10-11:40 am)
minDate = '2014-07-01 00:00'
maxDate = '2014-07-06 00:00'
ph_sensor = sensordbTable(table = 'ysi_corr', lakeID = 'EL', minDate = minDate, maxDate = maxDate, dateFormat = '%Y-%m-%d %H:%M') %>%
  as_tibble() %>%
  dplyr::filter(!is.na(pH))

ave_daily <- ph_sensor %>%
  mutate(date = as.Date(dateTime)) %>%
  group_by(date) %>%
  summarise(ave_pH = mean(pH, na.rm=T))

windows()
plot(ph_sensor$pH ~ ph_sensor$dateTime, type ='l')
abline(v = as.POSIXct(c('2014-07-01 10:00', '2014-07-01 11:40',
                        '2014-07-02 10:00', '2014-07-02 11:40',
                        '2014-07-03 10:00', '2014-07-03 11:40',
                        '2014-07-04 10:00', '2014-07-04 11:40',
                        '2014-07-05 10:00', '2014-07-05 11:40'), tz = 'GMT'), col='red', lwd= 2 )
segments(x0 = as.POSIXct(paste(ave_daily$date[1],'00:00'),tz='GMT'), y0 = ave_daily$ave_pH[1],
         x1 = as.POSIXct(paste(ave_daily$date[1],'23:50'),tz='GMT'), y1 = ave_daily$ave_pH[1], lty = 2, lwd=2)
segments(x0 = as.POSIXct(paste(ave_daily$date[2],'00:00'),tz='GMT'), y0 = ave_daily$ave_pH[2],
         x1 = as.POSIXct(paste(ave_daily$date[2],'23:50'),tz='GMT'), y1 = ave_daily$ave_pH[2], lty = 2, lwd=2)
segments(x0 = as.POSIXct(paste(ave_daily$date[3],'00:00'),tz='GMT'), y0 = ave_daily$ave_pH[3],
         x1 = as.POSIXct(paste(ave_daily$date[3],'23:50'),tz='GMT'), y1 = ave_daily$ave_pH[3], lty = 2, lwd=2)
segments(x0 = as.POSIXct(paste(ave_daily$date[4],'00:00'),tz='GMT'), y0 = ave_daily$ave_pH[4],
         x1 = as.POSIXct(paste(ave_daily$date[4],'23:50'),tz='GMT'), y1 = ave_daily$ave_pH[4], lty = 2, lwd=2)
segments(x0 = as.POSIXct(paste(ave_daily$date[5],'00:00'),tz='GMT'), y0 = ave_daily$ave_pH[5],
         x1 = as.POSIXct(paste(ave_daily$date[5],'23:50'),tz='GMT'), y1 = ave_daily$ave_pH[5], lty = 2, lwd=2)

# day to day changes in pH
minDate = '2014-05-01 00:00'
maxDate = '2014-07-10 00:00'
ph_sensor = sensordbTable(table = 'ysi_corr', lakeID = 'EL', minDate = minDate, maxDate = maxDate, dateFormat = '%Y-%m-%d %H:%M') %>%
  as_tibble() %>%
  dplyr::filter(!is.na(pH))

daily_pH <- ph_sensor %>%
  mutate(date = as.Date(dateTime)) %>%
  group_by(date) %>%
  summarise(ave_pH = mean(pH, na.rm=T)) %>%
  mutate(diff_pH = c(NA,diff(ave_pH)))

plot(daily_pH$diff_pH, pch=16, cex =2, ylab= 'day-to-day change in ave pH')


# how much does fraction of DIC pool as CO2 change within range during day?

# carbonate speciation
# notes: alkalinity is ANC of filtered water; ANC is ANC of unfiltered water - some people assume ANC = Alkalinity
# ANC is in microEq/L; convert to mol / kg of solution by dividing by 1e6 multiplying by 1000 (mol/m3) and dividing by water density (kg/m3) to get to mol/kg
# calculating speciations, fraction of DIC pool that is CO2, and pH by using water temperature and total alkalinity
sal = 0
temp = 20
DIC_pool = 1800 # mol C
Vepi = 51000 # m^3
pH= 6

aquaOut = AquaEnv::aquaenv(S = sal, # salinity
                           t = temp, # water temp
                           d = 0, # depth
                           lat = 46, # latitude
                           SumCO2 = DIC_pool/Vepi/rLakeAnalyzer::water.density(temp,sal),
                           pH = pH)

CO2=aquaOut$CO2[1]*water.density(temp,sal)*Vepi
HCO3=aquaOut$HCO3[1]*water.density(temp,sal)*Vepi
CO3=aquaOut$CO3[1]*water.density(temp,sal)*Vepi
fracCO2=CO2/DIC_pool

# if our range during day is 0.2 pH, how much CO2 flux occurs across the various pH's?
pH_vals = seq(5,8.5, by = .1)
pH_range = 0.01

out<- data.frame()
for(i in 1:length(pH_vals)){
  sal = 0
  temp = 20
  DIC_pool = 1800 # mol C
  Vepi = 51000 # m^3
  low_pH= pH_vals[i] - pH_range/2
  high_pH = pH_vals[i] + pH_range/2

  low_ph = AquaEnv::aquaenv(S = sal, # salinity
                             t = temp, # water temp
                             d = 0, # depth
                             lat = 46, # latitude
                             SumCO2 = DIC_pool/Vepi/rLakeAnalyzer::water.density(temp,sal),
                             pH = low_pH)

  low_CO2=low_ph$CO2[1]*water.density(temp,sal)*Vepi
  low_fracCO2=low_CO2/DIC_pool


  high_ph = AquaEnv::aquaenv(S = sal, # salinity
                            t = temp, # water temp
                            d = 0, # depth
                            lat = 46, # latitude
                            SumCO2 = DIC_pool/Vepi/rLakeAnalyzer::water.density(temp,sal),
                            pH = high_pH)

  high_CO2=high_ph$CO2[1]*water.density(temp,sal)*Vepi
  high_fracCO2=high_CO2/DIC_pool

  cur <- data.frame(pH = pH_vals[i], pH_range = pH_range, low_CO2 = low_CO2, high_CO2 = high_CO2, CO2_range = low_CO2 - high_CO2,
                    low_fracCO2 = low_fracCO2, high_fracCO2 = high_fracCO2, fracCO2_range = low_fracCO2 - high_fracCO2)

  out <- rbind(out, cur)
}

out

plot(out$CO2_range/Vepi*12 ~ out$pH, ylab = 'CO2 range (mg C / L)')



# how does this sub-daily flux compare to other daily fluxes?
# loaded CO2 from streams: (had to run first part of code in 'run_lake_C_EnKF.R')
mean(data2$dicIn)


# range in CO2 efflux throughout day
co2_efflux_range = mean(data2$kCO2/data2$epiVol*12)*out$CO2_range
plot(co2_efflux_range~ out$pH, ylab = 'CO2 Efflux range (mg C / L / day)')
