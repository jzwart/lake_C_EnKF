# fetch pH data

# looking at phi n EL
library(dplyr)
source('/Users/jzwart/Documents/Jake/Database/R Code/dbTable.R')
source('/Users/jzwart/Documents/Jake/Database/R Code/sensordbTable.R')

minDate = '2014-05-01 00:00'
maxDate = '2014-10-10 00:00'

ph = dbTable(table = 'limno_profiles', lakeID = 'EL', minDate = minDate, maxDate = maxDate, dateFormat = '%Y-%m-%d %H:%M') %>%
  as_tibble() %>%
  dplyr::filter(depthBottom < .5, pH<6.3)


plot(ph$pH ~ ph$dateTimeSample)

ph_sensor = sensordbTable(table = 'ysi_corr', lakeID = 'EL', minDate = minDate, maxDate = maxDate, dateFormat = '%Y-%m-%d %H:%M') %>%
  as_tibble() %>%
  dplyr::filter(!is.na(pH))

offset = left_join(ph_sensor,ph, by = c('dateTime' = 'dateTimeSample'), suffix = c('sensor', 'manual')) %>%
  dplyr::filter(!is.na(pHmanual))

offset
mean(offset$pHmanual - offset$pHsensor)

plot(offset$pHmanual - offset$pHsensor)


