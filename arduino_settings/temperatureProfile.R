#Making arduino Settings.ini file

library(stringr)
library(tidyverse)
library(lubridate)

startTemp <- 28.0 #Starting temperature
holdTemp <- c(28, 31, 31, 28) #Target temperature make list for multiple for low tank of 30 set to 31 and then it will switch back later 
rampStart <- '12:00' #Time HH:MM
rampEnd <- '12:00'
rampTime <- 0 #in hours
holdTime <- 0 #in hours
downStart <- '18:00' #Time HH:MM
downTime <-  0 #in hours
downEnd  <- '18:00'

rampRate <- (holdTemp-startTemp)/(rampTime*30) # 2 minute increase

rampProfile <- NA
for(i in 1:length(holdTemp)){
  print(holdTemp[i])
  ramp <- data.frame(i=round(seq(startTemp, holdTemp[i], rampRate[i]), 3))
  names(ramp)[1] <- paste0('Tank', i)
  rampProfile <- cbind(rampProfile, ramp)
}
rampProfile


rampProfile <- rampProfile %>%
  dplyr::select(-rampProfile)


rampTimes <- expand_grid(hour=seq(as.numeric(str_sub(rampStart, 1,2)), 
                                  as.numeric(str_sub(rampStart, 1,2))+rampTime, 1), 
                         minute=seq(0, 58, 2)) %>%
            arrange(hour) %>%
  rbind(., c(hour= as.numeric(str_sub(rampStart, 1,2))+rampTime, minute=as.numeric(str_sub(rampStart, 4, 5))))%>%
  mutate(minute=ifelse(minute<10, paste0('0', minute), as.character(minute)), Time=paste(sprintf("%02d", hour), minute, sep=':')) %>%
  dplyr::filter(Time > rampStart, Time <= rampEnd)
  #dplyr::select(Time) 

ramp <- bind_cols(rampTimes, rampProfile)

hold <- expand_grid(hour=seq(as.numeric(str_sub(rampStart, 1,2))+rampTime, 
                                  as.numeric(str_sub(rampStart, 1,2))+rampTime+holdTime, 1), minute=seq(0, 50, 5)) %>%
  arrange(hour) %>%
  mutate(minute=ifelse(minute<10, '00', as.character(minute)), Time=paste(sprintf("%02d", hour), minute, sep=':')) 


hold <- hold %>%
  mutate('Tank1'=holdTemp[1], 'Tank2'=holdTemp[2], Tank3=holdTemp[3], Tank4=holdTemp[4]) %>%
  filter(Time >= rampEnd, Time <= downStart)
                         

base <- expand_grid(hour=seq(0, 23), minute=seq(0, 50, 10)) %>%
  arrange(hour) %>%
  mutate(minute=ifelse(minute<10, paste0('0', minute), as.character(minute)), Time=paste(sprintf("%02d", hour), minute, sep=':'))
  #dplyr::select(Time) 
  
base <- base %>%
  mutate('Tank1'=startTemp, 'Tank2'=startTemp, Tank3=startTemp, Tank4=startTemp) %>%
  filter(Time < rampStart | Time > downEnd)

rampDownRate <- (holdTemp-startTemp)/(downTime*30) #2 minute increase

downProfile <- NA
for(i in 1:length(holdTemp)){
  print(holdTemp[i])
  rampD <- data.frame(i=round(seq(holdTemp[i], startTemp, -rampDownRate[i]), 3))
  names(rampD)[1] <- paste0('Tank', i)
  downProfile <- cbind(downProfile, rampD)
}

downTimes <- expand_grid(hour=seq(as.numeric(str_sub(downStart, 1,2)), 
                                  as.numeric(str_sub(downStart, 1,2))+downTime, 1), minute=seq(0, 58, 2)) %>%
  arrange(hour) %>%
  rbind(., c(hour= as.numeric(str_sub(downStart, 1,2))+downTime, minute=as.numeric(str_sub(downStart, 4, 5))))%>%
  mutate(minute=ifelse(minute<10, paste0('0', minute), as.character(minute)), Time=paste(sprintf("%02d", hour), minute, sep=':')) %>%
  filter(Time > downStart, Time <= downEnd)
#dplyr::select(Time) 

rampDown <- bind_cols(downTimes, downProfile) %>%
  dplyr::select(-downProfile)

profile <- rbind(base, ramp, hold, rampDown) %>%
  arrange(desc(Tank1)) %>%
  filter(!duplicated(Time)) %>%
  arrange(as.numeric(hour), as.numeric(minute)) %>%
  dplyr::select(Time, Tank1, Tank2, Tank3, Tank4) 
  #dplyr::select(Time, Tank1, Tank2, Tank3)
#profile$Tank1 <- startTemp
#profile$Tank2 <- startTemp #Uncomment if more than two tanks are staying at control

#remove leading zeros
profile <- profile %>%
  mutate(Time = sub("^0", "", Time))

write_delim(profile, '20221201_E5pulchra/20221201_E5pulchra_28.5_35_34_37_profiles.txt', delim = '\t', col_names = FALSE)
#Set to Setttings.ini on sd card
