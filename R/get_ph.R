# fetch pH data

# looking at phi n EL
library(dplyr)
source('/Users/jzwart/Documents/Jake/Database/R Code/dbTable.R')
source('/Users/jzwart/Documents/Jake/Database/R Code/sensordbTable.R')

minDate = '2014-05-01 00:00'
maxDate = '2014-07-22 00:00' # dates after this seem pretty odd for manual YSI

ph = dbTable(table = 'limno_profiles', lakeID = 'EL', minDate = minDate, maxDate = maxDate, dateFormat = '%Y-%m-%d %H:%M') %>%
  as_tibble() %>%
  dplyr::filter(depthBottom <= 2, depthBottom >0, !is.na(pH)) %>%
  group_by(dateSample) %>%
  summarize(pH = mean(pH)) %>%
  mutate(dateTimeSample = as.POSIXct(paste(dateSample, '11:00'),tz = 'GMT'))

plot(ph$pH ~ ph$dateTimeSample)

maxDate = '2014-10-10 00:00'
ph_sensor = sensordbTable(table = 'ysi_corr', lakeID = 'EL', minDate = minDate, maxDate = maxDate, dateFormat = '%Y-%m-%d %H:%M') %>%
  as_tibble() %>%
  dplyr::filter(!is.na(pH))

offset = left_join(ph_sensor,ph, by = c('dateTime' = 'dateTimeSample'), suffix = c('sensor', 'manual')) %>%
  dplyr::filter(!is.na(pHsensor)) %>%
  mutate(offset = pHmanual - pHsensor)

offset
ave_offset = mean(offset$offset, na.rm = T)

plot(offset$offset ~ offset$dateTime)

offset = mutate(offset, pH_sensor_corr = pHsensor + ave_offset)

plot(offset$pH_sensor_corr ~ offset$dateTime, type = 'l', ylab = 'pH', xlab ='')
points(offset$pHmanual ~ offset$dateTime, col = 'red', pch=16,cex =2)

# pH to be used in DA
out <- offset %>%
  select(dateTime, pH_sensor_corr) %>%
  mutate(date = as.Date(dateTime)) %>%
  group_by(date) %>%
  summarise(pH_sensor_corr = mean(pH_sensor_corr))

plot(out$pH_sensor_corr, type ='l')

saveRDS(out, 'Data/EL_ph.rds')
