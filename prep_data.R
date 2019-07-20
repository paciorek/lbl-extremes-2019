library(dplyr)
library(tidyr)

## Dec-Feb (winter) precip for stations in California, 1949-2017
load('CA_precip.Rda')  

dat <- CA_data[!is.na(CA_data$precipitation), ]
tbl <- table(dat$stationID)
## only stations with at least 6000 daily values (out of a max of 6120, ignoring leap years)
sub <- tbl[tbl > 6000]

## northern California
avail <- dat[dat$stationID %in% names(sub) & dat$latitude > 36, ]

## summarize by station and year
summ <- avail %>% group_by(stationID, latitude, longitude, stationName, season.year) %>%
    summarize(total = sum(precipitation), mx = max(precipitation),
              cnt = n())

totals <- summ %>% select(stationID, season.year, longitude, latitude, stationName, total) %>% spread(season.year, total)
maxs <- summ %>% select(stationID, season.year, longitude, latitude, stationName, mx) %>% spread(season.year, mx)
cnts <- summ %>% select(stationID, season.year, longitude, latitude, stationName, cnt) %>% spread(season.year, cnt)

save(totals, maxs, cnts, file = 'precip.Rda')

