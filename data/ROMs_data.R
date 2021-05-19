#ROMs Data from Mike Jacox
library(tidyverse)

# Averaging ROMS data from July-June
dat<-read.csv("roms_data_south.csv")
south.dat <- dat[dat$month %in% c(7:12,1:6),]
south.dat$yr <- ifelse(south.dat$month %in% 7:12, south.dat$year+1, south.dat$year)
avgROMs.south<-south.dat %>%
  group_by(yr) %>%
  summarise(south.sst=mean(sst,na.rm=T),
            south.ssh=mean(ssh,na.rm=T),
            south.ild=mean(ild,na.rm=T),
            south.BV=mean(BV,na.rm=T), #stratification
            south.cuti=mean(CUTI,na.rm=T), #upwelling
            south.beuti=mean(BEUTI,na.rm=T)) #nitrate flux


dat<-read.csv("roms_data_central.csv")
central.dat <- dat[dat$month %in% c(7:12,1:6),]
central.dat$yr <- ifelse(central.dat$month %in% 7:12, central.dat$year+1, central.dat$year)
avgROMs.central<-central.dat %>%
  group_by(yr) %>%
  summarise(central.sst=mean(sst,na.rm=T),
            central.ssh=mean(ssh,na.rm=T),
            central.ild=mean(ild,na.rm=T),
            central.BV=mean(BV,na.rm=T), #stratification
            central.cuti=mean(CUTI,na.rm=T), #upwelling
            central.beuti=mean(BEUTI,na.rm=T)) #nitrate flux

roms.july_june.mean=cbind(avgROMs.central,avgROMs.south)
roms.july_june.mean<-roms.july_june.mean[-c(8)]
write.csv(roms.july_june.mean,"roms.july_june.mean.csv")


# Averaging ROMS data across spring season (Mar-May)
avgROMs<-dat %>%
  group_by(year) %>%
  summarise(spr.sst=mean(sst[month%in%c(3:5)],na.rm=T),
            spr.ssh=mean(ssh[month%in%c(3:5)],na.rm=T),
            spr.ild=mean(ild[month%in%c(3:5)],na.rm=T),
            spr.BV=mean(BV[month%in%c(3:5)],na.rm=T), #stratification
            spr.cuti=mean(CUTI[month%in%c(3:5)],na.rm=T), #upwelling
            spr.beuti=mean(BEUTI[month%in%c(3:5)],na.rm=T)) #nitrate flux

colnames(avgROMs)<-c("year","spr.sst","spr.ssh","spr.ild","spr.BV","spr.upwelling","spr.nitrate_flux",
                         "win.sst","win.ssh","win.ild","win.BV","win.upwelling","win.nitrate_flux")

write.csv(avgROMs,"ROMs_south_spr.csv")

#For winter, modift this code (non-som repo) to use nov and dec data from previous year

win.pdo <- pdo[pdo$month %in% c("NOV", "DEC", "JAN", "FEB", "MAR"),]
win.pdo$win.yr <- ifelse(win.pdo$month %in% c("NOV", "DEC"), win.pdo$YEAR+1, win.pdo$YEAR)
win.pdo <- tapply(win.pdo$value, win.pdo$win.yr, mean)

win.npgo <- npgo[npgo$month %in% c(11,12,1:3),]
win.npgo$win.yr <- ifelse(win.npgo$month %in% 11:12, win.npgo$Year+1, win.npgo$Year)

win.npgo <- tapply(win.npgo$value, win.npgo$win.yr, mean)


colnames(avgROMs)<-c("year","spr.sst","spr.ssh","spr.ild","spr.BV","spr.upwelling","spr.nitrate_flux",
                     "win.sst","win.ssh","win.ild","win.BV","win.upwelling","win.nitrate_flux")

