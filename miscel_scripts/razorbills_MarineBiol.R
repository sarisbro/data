# root name
basedir <- "/Users/stephane/Documents/data/akiko/131004/"

# load processed data
load("/Users/stephane/Documents/data/akiko/131004/whoisbaka_all.RData")
load("/Users/stephane/Documents/data/akiko/131004/tdr1.2/tdr_new1.2.RData")

# colony locations: uncomment as appropriate...
col.loc1 <- c(-5.304686, 51.732764) # 2011, 2013
col.loc2 <- c(-5.300254, 51.745299) # 2012

# loading up the libs
library(diveMove)
library(fields)
library(geosphere)
library(ggplot2)
library(gridExtra)
library(ggmap)
library(maptools)
library(genefilter)
library(RColorBrewer)
library(maps)
library(mapdata)
library(robust)
library(rgl)
library(reshape)
library(zoo)
library(gdata)
library(hopach)
library(cluster)
library(rgdal)
library(aspace)
library(geneplotter)
library(pvclust)
library(graphics)
library(ncdf)
library(lme4)
library(marmap)
library(brainwaver)
library(lubridate)
library(nlme)

mytdr_d.thresh <- 1 # threshold depth (m) below which an underwater phase should be considered a dive.
mytdr_t.thresh <- 2 # threshold time (s) above which an underwater phase should be considered a dive.

####################################################################
#                                  WARNING
#
# IF YOU DON'T WANT TO WAIT FOR 72 HOURS, LOAD THE EXISTING RData 
# FILE WITH THE LINE BELOW...
#
####################################################################

####################################################################
#load(paste(basedir,"tdr", mytdr_d.thresh,".", mytdr_t.thresh,"/","tdr", mytdr_d.thresh,".", mytdr_t.thresh, ".RData", sep=""))
####################################################################

online <- 1#1 				# if not online (option 0), not maps
get.weather.data <- 1		# if not, the RData files are read for weather2 and weather3
							# WARNING: only valid for the 1/150 TDR offset & threshold!!!
interpolate.weather <- 0	# to interpolate or not weather data: WARNING, VERY SLOW...
normalize_speeds <- 1		# otherwise, recomputes GPS speeds (innacurate)

dir.create(paste(basedir,"tdr", mytdr_d.thresh,".", mytdr_t.thresh, sep=""))

##############################################################################
##############################################################################
# read and process TDR data
# don't forget to cleanup TDR the file!
# myfilelist[i] -> tdrfile
proc_tdr_data <- function(tdrfile, t.thresh=3, d.thresh=350){
	setwd(basedir)
	srcfile <- paste(basedir,"razor_press/",tdrfile,".txt",sep="")
	# read file
	rawdat <- read.table(srcfile, skip=0, sep=",",header=T)
	subsamp <- 1 #5
	concurrentCols <- 4:6
	speed <- F
	dateCol = 1; timeCol = 1; depthCol = 2
	names(rawdat) <-tolower(names(rawdat))
	rawdat.ncol <- seq(ncol(rawdat))
	dtpasted <- paste(rawdat[, dateCol], rawdat[, timeCol])
	# the TDRs in 2012 were an hour off...
	if((length(grep("2012",rawdat[1,1]))>0)|(length(grep("2013",rawdat[1,1]))>0)){
		datetime <- as.POSIXct(strptime(dtpasted, format = "%d/%m/%Y %H:%M:%S"), tz = "GMT") - 3600
	}else{
		datetime <- as.POSIXct(strptime(dtpasted, format = "%d/%m/%Y %H:%M:%S"), tz = "GMT")
	}
	
	# read and interpolate temperatures (cubic spline)
	srcfile_t <- paste(basedir,"razor_temp/",tdrfile,"_tmp.csv",sep="")
	rawtemp <- read.table(srcfile_t, skip=0, sep=",",header=T)
	if((length(grep("2012",rawdat[1,1]))>0)|(length(grep("2013",rawdat[1,1]))>0)){
		datetime.temp <- as.POSIXct(strptime(rawtemp[,1], format = "%d/%m/%Y %H:%M:%S"), tz = "GMT") - 3600
	}else{
		datetime.temp <- as.POSIXct(strptime(rawtemp[,1], format = "%d/%m/%Y %H:%M:%S"), tz = "GMT")
	}
	x <- datetime.temp
	y <- rawtemp[,2]
	#plot(x, y, pch=".", xlim=c(datetime.temp[1],datetime.temp[1000]))
	z <- spline(x, y, n=15+difftime(datetime.temp[length(datetime.temp)], datetime.temp[1], units="secs"))
	#z <- spline(x, y, n=15+difftime(max(datetime[length(datetime)],datetime.temp[length(datetime.temp)]), datetime.temp[1], units="secs"))
	#lines(z, col="red")
	interpol.temp <- z$y
	if(length(rawdat[,1])>length(interpol.temp)){
		interpol.temp <- c(interpol.temp,rep(NA,length(rawdat[,1])-length(interpol.temp)))
	}else{
		diff_t <- 1+length(interpol.temp)-length(rawdat[,1])
		interpol.temp <- interpol.temp[diff_t:length(interpol.temp)]
	}
	rawdat <- cbind(rawdat, interpol.temp)
	rawdat.ncol <- seq(ncol(rawdat))
	write.csv(rawdat,paste(basedir,"razor_temp/",tdrfile,"_alltemps.csv",sep=""))
	rm(x,y,z)
	
	# sync start times with GPS times
	gpsfile <- paste(basedir,"razor_gps/",tdrfile,".csv",sep="")
	gpsdata <- read.csv(gpsfile)
	datetime.gps <- as.POSIXct(strptime(gpsdata$datetime, format = "%Y-%m-%dT%H:%M:%SZ"), tz = "GMT")
	lag <- difftime(datetime[1],datetime.gps[1])
	#if(abs(lag)<60){
		datetime <- datetime - lag
	#}
	#else{
	#	if(lag>0){
	#	}
	#}
	
	
	# back to TDR data
	origint <- diveMove:::.getInterval(datetime) 
		if(!identical(all.equal(origint,subsamp), TRUE)){
		steptim <- as.numeric((subsamp)/origint)
		stepind <- seq(from = 1, to = length(datetime), by = round(steptim))
		datetime <- datetime[stepind]
		rawdat <- rawdat[stepind, ]
	}
	goodcc <- concurrentCols[is.finite(concurrentCols)]
	okconcurCols <-- goodcc %in% rawdat.ncol
	allbadcc <- all(!okconcurCols)
	somebadcc <- any(!okconcurCols) && !allbadcc
	if(somebadcc){
		warning(paste("Colums", concurrentCols[!okconcurCols], "given as concurrentCols could not be found\n"))
	}
	if(allbadcc && !is.null(concurrentCols)){
		warning("None of the columns given as concurrentCols exist\n")
		temp.df <- data.frame(rawdat[,3])
		colnames(temp.df) <- "temperature";
		tdr <- new("TDR", file = srcfile, time = datetime, depth = rawdat[, depthCol], concurrentData = temp.df, dtime = diveMove:::.getInterval(datetime))
	}
	# calibrate
	dcalib <- calibrateDepth(tdr, dive.thr=3, wet.thr = 3610, zoc.method="offset", offset=3, descent.crit.q=0.01, ascent.crit.q=0, knot.factor=20)
	
	pdf(paste(basedir,"tdr", mytdr_d.thresh,".", mytdr_t.thresh,"/",tdrfile,"_plotTDR.pdf",sep=""), width=10.5, height=6)                          
	plotTDR(dcalib, concurVars=c("temperature"),surface=T,interact=F)
	plotTDR(dcalib, diveNo=1:10, what="phases",interact=F)
	plot(tdr@time,tdr@depth,type="l")
	dev.off()
	
	# remove dives of duration < t.thresh
	tdrSumm1 <- diveStats(dcalib)
	tdrSumm1 <- subset(tdrSumm1, tdrSumm1$divetim > t.thresh)

	# IPQ computation: follows Elliott et al., 2008
	#res <- lm(log(tdrSumm1$postdive.dur[tdrSumm1$postdive.dur<t.thresh & tdrSumm1$postdive.dur>0])~tdrSumm1$divetim[tdrSumm1$postdive.dur<t.thresh & tdrSumm1$postdive.dur>0])
	res <- lmRob(log(tdrSumm1$postdive.dur[tdrSumm1$postdive.dur>t.thresh])~tdrSumm1$divetim[tdrSumm1$postdive.dur>t.thresh])
	b <- exp(as.numeric(res$coefficients[1]))
	c <- as.numeric(res$coefficients[2])
	#b <- 3.18; c <- 0.0191
	pdf(paste(basedir,"tdr", mytdr_d.thresh,".", mytdr_t.thresh,"/",tdrfile,"_logPostDivDur.vs.TimeDiv.pdf",sep=""), width=10.5, height=6)
	plot(log(tdrSumm1$postdive.dur[tdrSumm1$postdive.dur>t.thresh])~tdrSumm1$divetim[tdrSumm1$postdive.dur>t.thresh],pch=16,xlab="Dive time (s)",ylab="log PostDivDur")
	text(max(tdrSumm1$divetim[tdrSumm1$postdive.dur>t.thresh],na.rm=T)/2,max(log(tdrSumm1$postdive.dur[tdrSumm1$postdive.dur>t.thresh]),na.rm=T)-2,paste("b = ",format(b, scientific = F, digits=4),"; c = ",format(c, scientific = F, digits=4),sep=""))
	dev.off()
	
	# dive duration (s)
	u <- tdrSumm1$divetim
	# max depth
	d <- tdrSumm1$maxdep
	
	# the most *biutiful* formula
	transit_time <- tdrSumm1$desctim + tdrSumm1$asctim
	#res_tau <- lmRob(transit_time ~ -1 + tdrSumm1$maxdep)
	#tau <- res_tau$coefficients[[1]] * d # so that now, each track has its own coefficient
	#tau <- 1.321 * d
	tau <- transit_time
	ipq <- ((1 + b * c * exp(c * u)) * (u - tau))/(b * exp(c * u) + u)
	#print(paste("IPQ = ",ipq,"; tau = ",tau,"; depth = ",d,"; divDur = ",u,sep=""))
	
	# add *begdesc* temps to tdrSumm1
	# that is, temp at begasc
	print("Adding temperatures to tdrSumm1: this should be quick now!")
	ndives <- length(tdrSumm1[,1])
	temperature <- rep(NA,ndives)
	for(k in 1:ndives){
		temperature[k] <- rawdat[which(difftime(tdrSumm1$begdesc[k], datetime, units="secs") == 0),3]
	}
	# binding tdrSumm1 and tmp.temp
	tdrSumm1 <- cbind(tdrSumm1, temperature)
	
	return(list(tdrdata=tdr, tdrcalib=dcalib, tdrsum=tdrSumm1, ipq=ipq, tau=tau))
}
##############################################################################
##############################################################################
# read and process GPS data
# myfilelist[1] -> myfile
proc_gps_data <- function(myfile,tdrSumm1,curipq){
	setwd(basedir)
	gpsfile <- paste(basedir,"razor_gps/",myfile,".csv",sep="")
	# read file
	gpsdata <- read.csv(gpsfile)
	datetime.gps<- as.POSIXct(strptime(gpsdata$datetime, format = "%Y-%m-%dT%H:%M:%SZ"), tz = "GMT")	
	
	ndives <- length(tdrSumm1$begasc)
	t_inf <- numeric(ndives) # lower bound on time
	t_sup <- numeric(ndives) # upper bound
	
	for(i in 1:ndives){
		t_inf[i] <- sum(datetime.gps < tdrSumm1$begasc[i])
		t_sup[i] <- t_inf[i] + 1
	}
	
	lat_inf <- numeric(ndives); lon_inf <- numeric(ndives)
	lat_sup <- numeric(ndives); lon_sup <- numeric(ndives)
	
	lat_inf <- gpsdata$Latitude[t_inf]; lon_inf <- gpsdata$Longitude[t_inf]
	lat_sup <- gpsdata$Latitude[t_sup]; lon_sup <- gpsdata$Longitude[t_sup]
	gpsdata$latitude <- gpsdata$Latitude; gpsdata$longitude <- gpsdata$Longitude;
	gpsdata$distance_km <- gpsdata$Distance
	#if((length(grep("2013", myfile))>0) | (length(grep("2012", myfile))>0)){
	#	gpsdata$speed <- gpsdata$Speed/1000; 
	#}else{
	#	gpsdata$speed <- gpsdata$Speed/3600 # in km/h now
	#}
	if(length(grep("2013", datetime.gps[1]))>0){
		gpsdata$speed <- gpsdata$Speed/3600; 
	}
	if(length(grep("2011", datetime.gps[1]))>0){
		gpsdata$speed <- gpsdata$Speed/3600; 
	}
	if(length(grep("2012", datetime.gps[1]))>0){
		gpsdata$speed <- gpsdata$Speed/1000; 
	}
	if(length(grep("m93992", myfile))>0){
		gpsdata$speed <- gpsdata$speed*3.6; 
	}
	
	
	# interpolate assuming constant speed and linear splines along lon and lat independtly
	tot_time <- as.numeric(difftime(datetime.gps[t_sup],datetime.gps[t_inf], units="secs"))
	tinf.tmaxdepth <- as.numeric(difftime(tdrSumm1$begasc,datetime.gps[t_inf], units="secs"))
	proptimeduration <- tinf.tmaxdepth/tot_time
	# for lat
	interpol.lat <- lat_inf + proptimeduration * (lat_sup - lat_inf)
	interpol.lon <- lon_inf + proptimeduration * (lon_sup - lon_inf)
	
	# just to check, and guess what!?
	#plot(gpsdata$longitude,gpsdata$latitude,type="l", col="blue")
	#points(interpol.lon,interpol.lat, pch=3, col="red")
	
	# dist from colony to interpolated max depth dive positions -- in meters
	mydist2col <- numeric(ndives)
	sampled.year <- as.numeric(substring(datetime.gps[1],1,4))
	for(k in 1:ndives){
		if((sampled.year == 2011)|(sampled.year == 2013)){
			mydist2col[k] <- distMeeus(c(interpol.lon[k],interpol.lat[k]),col.loc1)
			#mydist2col[i] <- distVincentySphere(c(interpol.lon[i],interpol.lat[i]),col.loc1) #check if changing distance measure changes anything in results
		}else{
			mydist2col[k] <- distMeeus(c(interpol.lon[k],interpol.lat[k]),col.loc2)
			#mydist2col[i] <- distVincentySphere(c(interpol.lon[i],interpol.lat[i]),col.loc2) #check if changing distance measure changes anything in results
		}
	}
	

	if(normalize_speeds){
		return(list(lat=gpsdata$latitude, lon=gpsdata$longitude, ilat=interpol.lat, ilon=interpol.lon, dist2col=mydist2col, year=sampled.year, speed= gpsdata$speed, datetime= datetime.gps, traveldistance = gpsdata$distance_km, ipq=curipq))
	}else{
		# recompute speeds and travelled distances
		cur_speed <- c()
		dist_trav <- c()
		pairwiseD <- c()
		cur_speed[1] <- dist_trav[1] <- 0
		for(l in 2:length(gpsdata$latitude)){
			pairwiseD <- distMeeus(c(gpsdata$longitude[l],gpsdata$latitude[l]),c(gpsdata$longitude[l-1],gpsdata$latitude[l-1])) / 1000
			cur_speed[l] <- pairwiseD / difftime(datetime.gps[l],datetime.gps[l-1],units = "hours")[[1]]
			dist_trav[l] <- dist_trav[l-1] + pairwiseD
		}
		#plot(cur_speed,type="l")
		#plot(dist_trav,type="l")
		return(list(lat=gpsdata$latitude, lon=gpsdata$longitude, ilat=interpol.lat, ilon=interpol.lon, dist2col=mydist2col, year=sampled.year, speed= cur_speed, datetime= datetime.gps, traveldistance = dist_trav, ipq=curipq))
	}
}

#
#plot(as.numeric(coord_all[,8]),type="l")
#tapply(as.numeric(coord_all[,6]),coord_all[,5],mean)
#tapply(as.numeric(coord_all[,6]),coord_all[,1],mean)

##############################################################################
##############################################################################
# to convert unix time to POSIX time
unix2POSIXct <- function (time)   as.POSIXlt(structure(time, class = c("POSIXt", "POSIXct")) , "GMT")
##############################################################################
##############################################################################


mytdrfilelist <- list.files(paste(basedir,"razor_press/",sep=""), pattern=".txt")
myfilelist <- sub(".txt","",mytdrfilelist)
res_all <- matrix(nrow = 1, ncol = 45, byrow=T)
colnames(res_all) <- c("name","ID","year","ipq","tau","dist2col","ilon","ilat","begdesc","enddesc","begasc","desctim","botttim","asctim","divetim","descdist","bottdist","ascdist","bottdep.mean","bottdep.median","bottdep.sd","maxdep","postdive.dur","descD.min","descD.1stqu","descD.median","descD.mean","descD.3rdqu","descD.max","descD.sd","bottD.min","bottD.1stqu","bottD.median","bottD.mean","bottD.3rdqu","bottD.max","bottD.sd","ascD.min","ascD.1stqu","ascD.median","ascD.mean","ascD.3rdqu","ascD.max","ascD.sd","temperature")
coord_all <- matrix(nrow = 1, ncol = 9, byrow=T) 
colnames(coord_all) <- c("name","ID","lon","lat","year","speed","datetime","traveldistance","ipq")

for(i in 1:length(mytdrfilelist)){
	print(paste("Processing ",myfilelist[i]," (",i,"/",length(mytdrfilelist),")","!",sep=""))
	tmp.tdr <- proc_tdr_data(myfilelist[i], mytdr_t.thresh, mytdr_d.thresh) # this is here that you should change arguments such as offset and depth
	tdrsum <- tmp.tdr$tdrsum
	curipq <- tmp.tdr$ipq
	curtau <- tmp.tdr$tau
	tmp.gps <- proc_gps_data(myfilelist[i], tdrsum, curipq)
	
	setwd(paste(basedir,"tdr", mytdr_d.thresh,".", mytdr_t.thresh, sep=""))
	
	# plotting individual dives -- default sunsire/sunset
	pdf(paste(myfilelist[i],".diveprofile.default.pdf",sep=""), width=10.5, height=11)
	par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,1))
	plotTDR(tmp.tdr$tdrdata, interact=F)
	plotTDR(tmp.tdr$tdrcalib, interact=F)
	par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
	dev.off()
	
	# plotting individual dives -- custom sunsire/sunset
	pdf(paste(myfilelist[i],".diveprofile.pdf",sep=""), width=10.5, height=11)
	par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,1))
	myloc <- matrix(col.loc1, nrow=1)
	mylocSP <- SpatialPoints(myloc, proj4string=CRS("+proj=longlat +datum=WGS84"))
	seq_days <- seq(from=min(tmp.tdr$tdrsum$begdesc), to=max(tmp.tdr$tdrsum$begdesc), by="days")
	sunup    <- sunriset(mylocSP, as.POSIXct(seq_days), direction="sunrise", POSIXct.out=T)
	sundown  <- sunriset(mylocSP, seq_days, direction="sunset", POSIXct.out=T)
	plotTDR(tmp.tdr$tdrdata, interact=F, sunrise.time= format(mean(sundown$time), "%H:%M:%S"), sunset.time= format(mean(sunup$time), "%H:%M:%S"))
	plotTDR(tmp.tdr$tdrcalib, interact=F, sunrise.time= format(mean(sundown$time), "%H:%M:%S"), sunset.time= format(mean(sunup$time), "%H:%M:%S"))
	par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
	dev.off()
	
	if(online){
		# plotting individual tracks
		pdf(paste(myfilelist[i],".trax.pdf",sep=""))
		plot(tmp.gps$lon,tmp.gps$lat,type="l", col="blue", xlab="Longtitude", ylab="Latitude", ylim=c(51.2, 52.4), xlim=c(-5.9,-5.1), main=myfilelist[i])
		points(tmp.gps$ilon,tmp.gps$ilat, pch=3, col="red")
		dev.off()
		SHmap <- qmap(c(lon=mean(tmp.gps$lon), lat=mean(tmp.gps$lat)), maptype = 'satellite', zoom=9)
		p <- SHmap + 
			geom_path(data=data.frame(tmp.gps$lon,tmp.gps$lat), aes(tmp.gps.lon, tmp.gps.lat), lineend="round", col="red") + 
			geom_point(data=data.frame(tmp.gps$ilon,tmp.gps$ilat), aes(tmp.gps.ilon, tmp.gps.ilat), col="black", pch=3, cex=.9)
		p
		ggsave(paste(myfilelist[i],".google.trax.pdf",sep=""), plot=p, width=11, height=11)
	}
	
	# output for res_all
	ndives <- length(tmp.tdr$tdrsum[,1])
	name <- matrix(rep(myfilelist[i],ndives),nrow=ndives)
	year <- matrix(rep(tmp.gps$year,ndives),nrow=ndives)
	id <- matrix(rep(i,ndives),nrow=ndives)
	res_tmp <- cbind(name,id, year, tmp.tdr$ipq, tmp.tdr$tau, tmp.gps$dist2col, tmp.gps$ilon,tmp.gps$ilat,tmp.tdr$tdrsum)
	colnames(res_tmp) <- c("name","ID","year","ipq","tau","dist2col","ilon","ilat","begdesc","enddesc","begasc","desctim","botttim","asctim","divetim","descdist","bottdist","ascdist","bottdep.mean","bottdep.median","bottdep.sd","maxdep","postdive.dur","descD.min","descD.1stqu","descD.median","descD.mean","descD.3rdqu","descD.max","descD.sd","bottD.min","bottD.1stqu","bottD.median","bottD.mean","bottD.3rdqu","bottD.max","bottD.sd","ascD.min","ascD.1stqu","ascD.median","ascD.mean","ascD.3rdqu","ascD.max","ascD.sd","temperature")
	res_all <- rbind(res_all,res_tmp)
	
	# output for coord_all
	nrec <- length(tmp.gps$lon)
	name <- matrix(rep(myfilelist[i], nrec),nrow= nrec)
	id <- matrix(rep(i, nrec),nrow= nrec)
	year <- matrix(rep(tmp.gps$year,nrec),nrow=nrec)
	coord_tmp <- cbind(name,id,tmp.gps$lon,tmp.gps$lat,year,tmp.gps$speed,tmp.gps$datetime,tmp.gps$traveldistance,tmp.gps$ipq)
	colnames(coord_tmp) <- c("name","ID","lon","lat","year","speed","datetime","traveldistance","ipq")
	coord_all <- rbind(coord_all, coord_tmp)
	
}

# eliminate first row (that contains NA)
res_all <- res_all[2:length(res_all[,1]),]
coord_all <- coord_all[2:length(coord_all[,1]),]
#format(head(unix2POSIXct(coord_all[,7])), "%H:%M:%S")
#format(head(unix2POSIXct(coord_all[,7])), "%Y/%m/%d")

# eliminates IPQ 0 dives [those with no bottom time]
res_all_allIPQs <- res_all
coord_all_allIPQs <- coord_all
res_all <- subset(res_all, res_all[,4]>0)
coord_all <- subset(coord_all, coord_all[,9]>0)

# takes log IPQ [normalization]
res_all[,4] <- log(res_all[,4])
coord_all[,9] <- log(as.numeric(coord_all[,9]))


#plot(res_all$divetim,res_all$maxdep,pch=16)
pdf("Divetim_Maxdep.pdf",width=6,height=4)
smoothScatter(res_all$divetim,res_all$maxdep,xlab="Dive duration (s)",ylab="Maximum depth (m)")
#x <- seq(0,80,.1)
#y <- 1.321 * x
#points(x,y,type="l")
dev.off()



#############
# read mass #
#############

mass <- read.csv(paste(basedir,"razor_mass.csv",sep=""),header=T)
mass_early <- numeric(length(res_all[,1]))
mass_late <- numeric(length(res_all[,1]))
mass_years <- numeric(length(mass[,1]))

for(i in 1:length(mass[,1])){
	mass_early[which(res_all$name == mass[i,1])] <- mass[i,2]
	mass_late[which(res_all$name == mass[i,1])] <- mass[i,3]
	mass_years[i] <- res_all$year[res_all$name == mass[i,1]][1]
}

mass_tmp <- cbind(mass_early, mass_late)
res_all <- cbind(res_all, mass_tmp)
rm(mass_tmp)

mean(res_all$mass_early[res_all$year==2011] - res_all$mass_late[res_all$year==2011])
mean(res_all$mass_early[res_all$year==2012] - res_all$mass_late[res_all$year==2012])

diffmass2011 <- mass[,3][mass_years == 2011] - mass[,2][mass_years == 2011]
diffmass2012 <- mass[,3][mass_years == 2012] - mass[,2][mass_years == 2012]
mean(diffmass2011)
mean(diffmass2012)
t.test(diffmass2011, diffmass2012)
# Mann-Whitney U Test
wilcox.test(diffmass2011, diffmass2012)
wilcox.test(diffmass2011, diffmass2012, alternative = "greater", correct=F)

####################################################################
save.image(paste("tdr_new", mytdr_d.thresh,".", mytdr_t.thresh, ".RData", sep=""))
####################################################################

# eliminates tau from res_all [was used for debug purposes only]
res_all <- cbind(res_all[,1:4], res_all[,6:length(res_all[1,])])




# creates stage columns
stage <- res_all$name
stage[grep("2011", res_all$year)] <- "chick"
stage[grep("2012", res_all$year)] <- "egg"
stage[grep("chick", res_all$name)] <- "chick"
stage[grep("egg", res_all$name)] <- "egg"
res_all <- cbind(res_all, stage)

stage <- coord_all[,1]
stage[grep("2011", coord_all[,5])] <- "chick"
stage[grep("2012", coord_all[,5])] <- "egg"
stage[grep("chick", coord_all[,1])] <- "chick"
stage[grep("egg", coord_all[,1])] <- "egg"
coord_all <- cbind(coord_all, stage)

# creates bird ID columns
birdIDcol <- res_all$name
birdIDcol <- sub("2013_chick_", "", birdIDcol)
birdIDcol <- sub("2013_egg_", "", birdIDcol)
res_all <- cbind(res_all, birdIDcol)

birdIDcol <- coord_all[,1]
birdIDcol <- sub("2013_chick_", "", birdIDcol)
birdIDcol <- sub("2013_egg_", "", birdIDcol)
coord_all <- cbind(coord_all, birdIDcol)

# creates the stage ID column (only in coord_all of course, since all res_all entries are about dives...)
speedlimit <- 5 # km/h
#plot(density(as.numeric(coord_all[,6])),xlim=c(0,5))
hist(as.numeric(coord_all[,6])[as.numeric(coord_all[,6])>5],100)
hist(as.numeric(coord_all[,6])[as.numeric(coord_all[,6])>5],100,ylim=c(0,50))
stateID <- rep("rest",length(coord_all[,1]))
stateID[coord_all[,6]>speedlimit] <- "fly"
for(i in 1:length(coord_all[,7])){
	if( sum((coord_all[i,7]>res_all$begdesc[res_all$name==coord_all[i,1]]) & (coord_all[i,7]<res_all$begdesc[res_all$name==coord_all[i,1]]+res_all$divetim[res_all$name==coord_all[i,1]]))>0 ){
		stateID[i] <- "dive"
	}
}
stateID
coord_all <- cbind(coord_all, stateID)


ks.test(cumsum(as.numeric(coord_all[,6])[as.numeric(coord_all[,6])>5]), cumsum(as.numeric(coord_all[,6])[as.numeric(coord_all[,6])>10]))

#stateIDn <- rep(0,length(stateID))
#stateIDn[stateID=="fly"] <- 1
#stateIDn[stateID=="dive"] <- -1
#plot(stateIDn,type="l",col=factor(birdIDcol))


tmp_res_all <- res_all[complete.cases(res_all),]
tmp_coord_all <- coord_all[complete.cases(coord_all),]
pdf("activity_maps.pdf",width=10, height=15)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(3,2))
# rest
bathy <- getNOAA.bathy(lon1 = min(tmp_res_all$ilon,na.rm=T), lon2 = max(tmp_res_all$ilon,na.rm=T), lat1 = min(tmp_res_all$ilat,na.rm=T), lat2 = max(tmp_res_all$ilat,na.rm=T),resolution = 1)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="rest"], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="rest"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: rest")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="rest"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="rest"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="rest"], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="rest"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: rest")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="rest"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="rest"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
# fly
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="fly"], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="fly"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: fly")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="fly"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="fly"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="fly"], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="fly"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: fly")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="fly"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="fly"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
# dive
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: dive")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: dive")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
# all
#smoothScatter(tmp_res_all$ilon[tmp_res_all$stage=="egg"], tmp_res_all$ilat[tmp_res_all$stage=="egg"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: TDR dives")
#z <- kde2d(tmp_res_all$ilon[tmp_res_all$stage=="egg"], tmp_res_all$ilat[tmp_res_all$stage=="egg"], n=50)
#contour(z, drawlabels=T, nlevels=4, add=TRUE)
#lines(tmp_res_all$ilon[tmp_res_all$stage=="egg"], tmp_res_all$ilat[tmp_res_all$stage=="egg"],type="l",col="black",lwd=.5)
#map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
#points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
#smoothScatter(tmp_res_all$ilon[tmp_res_all$stage=="chick"], tmp_res_all$ilat[tmp_res_all$stage=="chick"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: TDR dives")
#z <- kde2d(tmp_res_all$ilon[tmp_res_all$stage=="chick"], tmp_res_all$ilat[tmp_res_all$stage=="chick"], n=50)
#contour(z, drawlabels=T, nlevels=4, add=TRUE)
#lines(tmp_res_all$ilon[tmp_res_all$stage=="chick"], tmp_res_all$ilat[tmp_res_all$stage=="chick"],type="l",col="black",lwd=.5)
#map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
#points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()
rm(z)


#########################
# environmental factors #
#########################

# NB. all are taken at begdesc (ie diving time)
if(get.weather.data){
	if(interpolate.weather){
		############################
		# weather2: from bouy data #
		############################
		weather2 <- matrix(nrow =length(res_all[,1]), ncol = 6, byrow=T) 
		colnames(weather2) <- c("sea_temp","windspeed","winddirection","WindDir_SD1","Gust_3s_Max","AirTC_Avg")
	
		# seaTemps
		seaTemps <- read.table(paste(basedir,"razor_temp/skomer_buoy_seaTemp.csv",sep=""),sep=",",head=T)
		seaTemps[,1] <- as.POSIXct(strptime(seaTemps[,1], format = "%m/%d/%Y %H:%M"), tz = "GMT")
		# First: interpolate unsampled data
		datetime_t <- seaTemps[,1]
		x <- datetime_t
		y <- seaTemps[,2]
		z <- spline(x, y, n=difftime(datetime_t[length(datetime_t)], datetime_t[1], units="secs"))
		timeseries_t <- z$x
		interpol.temp <- z$y
		
		# wind data
		wind <- read.table(paste(basedir,"razor_temp/skomer_buoy_wind.csv",sep=""),sep=",",head=T)
		wind[,1] <- as.POSIXct(strptime(wind[,1], format = "%m/%d/%Y %H:%M"), tz = "GMT")
		# First: interpolate unsampled data
		datetime_w <- wind[,1]
		x <- datetime_w
		y <- wind[,2]
		z <- spline(x, y, n=difftime(datetime_w[length(datetime_w)], datetime_w[1], units="secs"))
		interpol.w2 <- z$y
		y <- wind[,3]
		z <- spline(x, y, n=difftime(datetime_w[length(datetime_w)], datetime_w[1], units="secs"))
		interpol.w3 <- z$y
		y <- wind[,4]
		z <- spline(x, y, n=difftime(datetime_w[length(datetime_w)], datetime_w[1], units="secs"))
		interpol.w4 <- z$y
		y <- wind[,5]
		z <- spline(x, y, n=difftime(datetime_w[length(datetime_w)], datetime_w[1], units="secs"))
		interpol.w5 <- z$y
		y <- wind[,6]
		z <- spline(x, y, n=difftime(datetime_w[length(datetime_w)], datetime_w[1], units="secs"))
		interpol.w6 <- z$y
		timeseries_w <- z$x
	
		for(i in 1:length(res_all[,1])){
			# seaTemps data
			j <- which(abs(timeseries_t - res_all[i,8])<1)[1]
			weather2[i,1] <- interpol.temp[j]
			print(paste("Weather for ",i,"/",length(res_all[,1])," (",unix2POSIXct(res_all[i,8]),"): ",weather2[i,1]," deg C",sep=""))
			# wind data
			j <- which(abs(timeseries_w - res_all[i,8])<1)[1]
			weather2[i,2] <- interpol.w2[j]
			weather2[i,3] <- interpol.w3[j]
			weather2[i,4] <- interpol.w4[j]
			weather2[i,5] <- interpol.w5[j]
			weather2[i,6] <- interpol.w6[j]
			print(paste("Weather for ",i,"/",length(res_all[,1]),": ",weather2[i,6]," deg C",sep=""))
		}
		
		################################
		# weather3: oceanographic data #
		################################
		seaTemps3 <- read.table(paste(basedir,"razor_temp/oceanograph_all.csv",sep=""),sep=",",head=T)
		seaTemps3[,10] <- as.POSIXct(strptime(seaTemps3[,10], format = "%m/%d/%Y %H:%M"), tz = "GMT")
		weather3 <- matrix(nrow =length(res_all[,1]), ncol = length(colnames(seaTemps3)[5:16]), byrow=T) 
		colnames(weather3) <- colnames(seaTemps3)[5:16]
		# First: interpolate unsampled data
		datetime_st <- seaTemps3[,10]
		x <- datetime_st
		y <- seaTemps3[,5]
		z <- spline(x, y, n=difftime(datetime_st[length(datetime_st)], datetime_st[1], units="secs"))
		interpol.st5 <- z$y
		y <- seaTemps3[,6]
		z <- spline(x, y, n=difftime(datetime_st[length(datetime_st)], datetime_st[1], units="secs"))
		interpol.st6 <- z$y
		y <- seaTemps3[,7]
		z <- spline(x, y, n=difftime(datetime_st[length(datetime_st)], datetime_st[1], units="secs"))
		interpol.st7 <- z$y
		y <- seaTemps3[,8]
		z <- spline(x, y, n=difftime(datetime_st[length(datetime_st)], datetime_st[1], units="secs"))
		interpol.st8 <- z$y
		y <- seaTemps3[,9]
		z <- spline(x, y, n=difftime(datetime_st[length(datetime_st)], datetime_st[1], units="secs"))
		interpol.st9 <- floor(z$y)
		y <- seaTemps3[,10]
		z <- spline(x, y, n=difftime(datetime_st[length(datetime_st)], datetime_st[1], units="secs"))
		interpol.st10 <- z$y
		y <- seaTemps3[,11]
		z <- spline(x, y, n=difftime(datetime_st[length(datetime_st)], datetime_st[1], units="secs"))
		interpol.st11 <- z$y
		y <- seaTemps3[,12]
		z <- spline(x, y, n=difftime(datetime_st[length(datetime_st)], datetime_st[1], units="secs"))
		interpol.st12 <- z$y
		y <- seaTemps3[,13]
		z <- spline(x, y, n=difftime(datetime_st[length(datetime_st)], datetime_st[1], units="secs"))
		interpol.st13 <- z$y
		y <- seaTemps3[,14]
		z <- spline(x, y, n=difftime(datetime_st[length(datetime_st)], datetime_st[1], units="secs"))
		interpol.st14 <- z$y
		y <- seaTemps3[,15]
		z <- spline(x, y, n=difftime(datetime_st[length(datetime_st)], datetime_st[1], units="secs"))
		interpol.st15 <- z$y
		y <- seaTemps3[,16]
		z <- spline(x, y, n=difftime(datetime_st[length(datetime_st)], datetime_st[1], units="secs"))
		interpol.st16 <- z$y
		timeseries_st <- z$x
	
		for(i in 1:length(res_all[,1])){
			# seaTemps data
			j <- which(abs(timeseries_st - res_all[i,8])<1)[1]
			weather3[i,1] <- interpol.st5[j]
			weather3[i,2] <- interpol.st6[j]
			weather3[i,3] <- interpol.st7[j]
			weather3[i,4] <- interpol.st8[j]
			weather3[i,5] <- interpol.st9[j]
			weather3[i,6] <- interpol.st10[j]
			weather3[i,7] <- interpol.st11[j]
			weather3[i,8] <- interpol.st12[j]
			weather3[i,9] <- interpol.st13[j]
			weather3[i,10] <- interpol.st14[j]
			weather3[i,11] <- interpol.st15[j]
			weather3[i,12] <- interpol.st16[j]
			print(paste("Weather for ",i,"/",length(res_all[,1])," (",unix2POSIXct(res_all[i,8]),"): ",weather3[i,12]," deg C",sep=""))
		}
	
		# final saves
		save(weather2,file=paste(basedir,"razor_temp/weather2.i.RData",sep=""))
		save(weather3,file=paste(basedir,"razor_temp/weather3.i.RData",sep=""))	
		
		rm(x,y,z, seaTemps, datetime_st, datetime_t, timeseries_w, timeseries_t, interpol.temp, wind, datetime_w, interpol.w2, interpol.w3, interpol.w4, interpol.w5, interpol.w6, interpol.st10, interpol.st11, interpol.st12, interpol.st13, interpol.st14, interpol.st15, interpol.st16, interpol.st5, interpol.st6, interpol.st7, interpol.st8, interpol.st9)
		
	}else{ # no interpolations for weather
	
		############################
		# weather2: from bouy data #
		############################
		weather2 <- matrix(nrow =length(res_all[,1]), ncol = 6, byrow=T) 
		colnames(weather2) <- c("sea_temp","windspeed","winddirection","WindDir_SD1","Gust_3s_Max","AirTC_Avg")
	
		# seaTemps
		seaTemps <- read.table(paste(basedir,"razor_temp/skomer_buoy_seaTemp.csv",sep=""),sep=",",head=T)
		seaTemps[,1] <- as.POSIXct(strptime(seaTemps[,1], format = "%m/%d/%Y %H:%M"), tz = "GMT")
		
		# wind data
		wind <- read.table(paste(basedir,"razor_temp/skomer_buoy_wind.csv",sep=""),sep=",",head=T)
		wind[,1] <- as.POSIXct(strptime(wind[,1], format = "%m/%d/%Y %H:%M"), tz = "GMT")
	
		for(i in 1:length(res_all[,1])){
			# seaTemps data
			j <- which(abs(difftime(seaTemps[,1], unix2POSIXct(res_all[i,8]), units="hours"))<.75)[[1]]
			weather2[i,1] <- seaTemps[j,2]
			if(!(i %% 500)){
				print(paste("Weather for ",i,"/",length(res_all[,1])," (",unix2POSIXct(res_all[i,8]),"): ",weather2[i,1]," deg C",sep=""))
			}
			# wind data
			j <- which(abs(difftime(wind[,1], unix2POSIXct(res_all[i,8]), units="hours"))<.75)[[1]]
			weather2[i,2] <- wind[j,2]
			weather2[i,3] <- wind[j,3]
			weather2[i,4] <- wind[j,4]
			weather2[i,5] <- wind[j,5]
			weather2[i,6] <- wind[j,6]
			if(!(i %% 500)){
				print(paste("Weather for ",i,"/",length(res_all[,1]),": ",weather2[i,6]," deg C",sep=""))
			}
		}
		
		################################
		# weather3: oceanographic data #
		################################
		seaTemps3 <- read.table(paste(basedir,"razor_temp/oceanograph_all.csv",sep=""),sep=",",head=T)
		seaTemps3[,10] <- as.POSIXct(strptime(seaTemps3[,10], format = "%m/%d/%Y %H:%M"), tz = "GMT")
		weather3 <- matrix(nrow =length(res_all[,1]), ncol = length(colnames(seaTemps3)[5:16]), byrow=T) 
		colnames(weather3) <- colnames(seaTemps3)[5:16]
		for(i in 1:length(res_all[,1])){
			#print(i)
			j <- which(abs(difftime(seaTemps3[,10], unix2POSIXct(res_all[i,8]), units="hours"))<.75)[[1]]
			weather3[i,1] <- seaTemps3[j,5]
			weather3[i,2] <- seaTemps3[j,6]
			weather3[i,3] <- seaTemps3[j,7]
			weather3[i,4] <- seaTemps3[j,8]
			weather3[i,5] <- as.numeric(seaTemps3[j,9])
			weather3[i,6] <- seaTemps3[j,10]
			weather3[i,7] <- seaTemps3[j,11]
			weather3[i,8] <- seaTemps3[j,12]
			weather3[i,9] <- seaTemps3[j,13]
			weather3[i,10] <- seaTemps3[j,14]
			weather3[i,11] <- seaTemps3[j,15]
			weather3[i,12] <- seaTemps3[j,16]
			if(!(i %% 500)){
				print(paste("Weather for ",i,"/",length(res_all[,1])," (",unix2POSIXct(res_all[i,8]),"): ",weather3[i,12]," deg C",sep=""))
			}
		}
		
		# final saves
		save(weather2,file=paste(basedir,"razor_temp/weather2.ni.RData",sep=""))
		save(weather3,file=paste(basedir,"razor_temp/weather3.ni.RData",sep=""))	
		
		rm(x,y,z, seaTemps, wind, seaTemps3)
	}
	
	write.csv(weather2,paste(basedir,"razor_temp/weather2.csv",sep=""))
	save(weather2,file=paste(basedir,"razor_temp/weather2.RData",sep=""))
	
	write.csv(weather3,paste(basedir,"razor_temp/weather3.csv",sep=""))
	save(weather3,file=paste(basedir,"razor_temp/weather3.RData",sep=""))
}else{
	load(paste(basedir,"weather2.i.RData",sep=""))
	load(paste(basedir,"weather3.i.RData",sep=""))
}


res_all_B4_weather2 <- res_all
res_all <- cbind(res_all,weather2)
res_all_B4_weather3 <- res_all
res_all <- cbind(res_all,weather3)
	

####################################################################
save.image(paste("tdr_new", mytdr_d.thresh,".", mytdr_t.thresh, ".RData", sep=""))
####################################################################





####################################################################
# filter dives: omit those "too" close to colony
####################################################################
threshold <- 1000 # unit is meters; 1000 seems to help "iron out" false positives

# compute dist to colony: for res_all
ncoords <- length(res_all[,1])
dist2col <- numeric(ncoords)
for(i in 1:ncoords){
	if((as.numeric(res_all[i,3]) == 2011)|(as.numeric(res_all[i,3]) == 2013)){
		dist2col[i] <- distMeeus(c(as.numeric(res_all[i,6]),as.numeric(res_all[i,7])),col.loc1)
	}else{
		dist2col[i] <- distMeeus(c(as.numeric(res_all[i,6]),as.numeric(res_all[i,7])),col.loc2)
	}
}
res_all_ori <- res_all	# saves a copy, as original object gets overwritten
res_all <- subset(res_all, dist2col>threshold)

# compute dist to colony: for coord_all
ncoords_coordall <- length(coord_all[,1])
dist2col_coordall <- numeric(ncoords_coordall)
for(i in 1:ncoords_coordall){
	if((as.numeric(coord_all[i,5]) == 2011)|(as.numeric(coord_all[i,5]) == 2013)){
		dist2col_coordall[i] <- distMeeus(c(as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4])),col.loc1)
	}else{
		dist2col_coordall[i] <- distMeeus(c(as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4])),col.loc2)
	}
}
coord_all_ori <- coord_all	# saves a copy, as original object gets overwritten
coord_all <- subset(coord_all, dist2col_coordall>threshold)




tmp_res_all <- res_all[complete.cases(res_all),]
tmp_coord_all <- coord_all[complete.cases(coord_all),]
pdf("activity_maps_colony_trimmed.pdf",width=10, height=15)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(3,2))
# rest
bathy <- getNOAA.bathy(lon1 = min(tmp_res_all$ilon,na.rm=T), lon2 = max(tmp_res_all$ilon,na.rm=T), lat1 = min(tmp_res_all$ilat,na.rm=T), lat2 = max(tmp_res_all$ilat,na.rm=T),resolution = 1)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="rest"], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="rest"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: rest")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="rest"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="rest"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="rest"], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="rest"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: rest")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="rest"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="rest"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
# fly
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="fly"], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="fly"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: fly")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="fly"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="fly"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="fly"], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="fly"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: fly")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="fly"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="fly"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
# dive
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: dive")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: dive")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
# all
#smoothScatter(tmp_res_all$ilon[tmp_res_all$stage=="egg"], tmp_res_all$ilat[tmp_res_all$stage=="egg"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: TDR dives")
#z <- kde2d(tmp_res_all$ilon[tmp_res_all$stage=="egg"], tmp_res_all$ilat[tmp_res_all$stage=="egg"], n=50)
#contour(z, drawlabels=T, nlevels=4, add=TRUE)
#lines(tmp_res_all$ilon[tmp_res_all$stage=="egg"], tmp_res_all$ilat[tmp_res_all$stage=="egg"],type="l",col="black",lwd=.5)
#map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
#points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
#smoothScatter(tmp_res_all$ilon[tmp_res_all$stage=="chick"], tmp_res_all$ilat[tmp_res_all$stage=="chick"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: TDR dives")
#z <- kde2d(tmp_res_all$ilon[tmp_res_all$stage=="chick"], tmp_res_all$ilat[tmp_res_all$stage=="chick"], n=50)
#contour(z, drawlabels=T, nlevels=4, add=TRUE)
#lines(tmp_res_all$ilon[tmp_res_all$stage=="chick"], tmp_res_all$ilat[tmp_res_all$stage=="chick"],type="l",col="black",lwd=.5)
#map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
#points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()
rm(z)

tmp_res_all <- res_all[complete.cases(res_all),]
tmp_coord_all <- coord_all[complete.cases(coord_all),]
pdf("activity_maps_colony_trimmed2.pdf",width=10, height=10)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2))
# rest
bathy <- getNOAA.bathy(lon1 = min(tmp_res_all$ilon,na.rm=T), lon2 = max(tmp_res_all$ilon,na.rm=T), lat1 = min(tmp_res_all$ilat,na.rm=T), lat2 = max(tmp_res_all$ilat,na.rm=T),resolution = 1)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: fly/rest")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: fly/rest")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
# dive
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: dive")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: dive")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
# all
#smoothScatter(tmp_res_all$ilon[tmp_res_all$stage=="egg"], tmp_res_all$ilat[tmp_res_all$stage=="egg"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: TDR dives")
#z <- kde2d(tmp_res_all$ilon[tmp_res_all$stage=="egg"], tmp_res_all$ilat[tmp_res_all$stage=="egg"], n=50)
#contour(z, drawlabels=T, nlevels=4, add=TRUE)
#lines(tmp_res_all$ilon[tmp_res_all$stage=="egg"], tmp_res_all$ilat[tmp_res_all$stage=="egg"],type="l",col="black",lwd=.5)
#map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
#points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
#smoothScatter(tmp_res_all$ilon[tmp_res_all$stage=="chick"], tmp_res_all$ilat[tmp_res_all$stage=="chick"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: TDR dives")
#z <- kde2d(tmp_res_all$ilon[tmp_res_all$stage=="chick"], tmp_res_all$ilat[tmp_res_all$stage=="chick"], n=50)
#contour(z, drawlabels=T, nlevels=4, add=TRUE)
#lines(tmp_res_all$ilon[tmp_res_all$stage=="chick"], tmp_res_all$ilat[tmp_res_all$stage=="chick"],type="l",col="black",lwd=.5)
#map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
#points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()
rm(z)


######################################################
# relationship between bathymetry and dive intensity #
######################################################
pdf("activity_maps_fitted_colony_trimmed.pdf",width=6, height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,1))
mean_bathy <- mean(bathy)
# egg -- dive
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), n=50)
dive_density <- c()
depth <- c()
counter <- 1
for(i in 1:length(z$x)){
	for(j in 1:length(z$y)){
		dive_density[counter] <- z$z[i,j]
		depth[counter] <- bathy[floor(dim(bathy)[1]*i/length(z$x)),floor(dim(bathy)[2]*j/length(z$y))]
		counter <- counter + 1
	}
}
counter
dd <- data.frame(dive_density,depth)
dd$depth[dd$depth>0] <- NA
x <- dd$dive_density
y <- dd$depth
plot(x,y,xlab="Dive density",ylab="Bathymetry (m)",pch=20,main="Incubating bird: dive")
dd1.lm <- lm(y ~ x)
dd2.lm <- lm(y ~ x +I(x^2))
summary(dd1.lm)
summary(dd2.lm)
anova(dd1.lm,dd2.lm)
xv <- seq(0,7,.01)
yv <- predict(dd2.lm,list(x=xv))
points(xv, yv,pch=".",col="red")
abline(h= mean_bathy, lty=2, lwd=2, col="orange")
# chick -- dive
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), n=50)
dive_density <- c()
depth <- c()
counter <- 1
for(i in 1:length(z$x)){
	for(j in 1:length(z$y)){
		dive_density[counter] <- z$z[i,j]
		depth[counter] <- bathy[floor(dim(bathy)[1]*i/length(z$x)),floor(dim(bathy)[2]*j/length(z$y))]
		counter <- counter + 1
	}
}
counter
dd <- data.frame(dive_density,depth)
dd$depth[dd$depth>0] <- NA
x <- dd$dive_density
y <- dd$depth
plot(x,y,xlab="Dive density",ylab="Bathymetry (m)",pch=20,main="Chick rearing bird: dive")
dd1.lm <- lm(y ~ x)
dd2.lm <- lm(y ~ x +I(x^2))
dd3.lm <- nls(y ~ 1/(1+a*exp(-b*x)),start = list(a = 1, b = 1),algorithm = "plinear")
summary(dd1.lm)
summary(dd2.lm)
summary(dd3.lm)
anova(dd1.lm,dd2.lm)
xv <- seq(0,160,.01)
yv <- (predict(dd3.lm,list(x=xv)))
points(xv, yv,pch=".",col="red")
abline(h= mean_bathy, lty=2, lwd=2, col="orange")

par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()






###################################
# Plots all tracks on same figure #
###################################

if(online){
	all_tracks <- data.frame(coord_all[,1],as.numeric(coord_all[,2]),as.numeric(coord_all[,3]),as.numeric(coord_all[,4]),as.numeric(coord_all[,5]))
	colnames(all_tracks) <- c("name","ID","lon","lat","stage")
	SHmap <- qmap(c(lon=mean(as.numeric(coord_all[,3]), na.rm=T), lat=mean(as.numeric(coord_all[,4]), na.rm=T)), maptype = 'satellite', zoom=9)
	p <- SHmap + 
		geom_path(data= all_tracks, aes(lon, lat, col=all_tracks$ID), lineend="round") #+ 
		#opts(legend.position="left", legend.text = theme_text(colour="black", size = 6)) +
		#geom_point(data=data.frame(res_all$ilon,res_all$ilat), aes(res_all$ilon,res_all$ilat), col="black", pch=3, cex=.9)
	#p
	ggsave("all.google.trax.pdf", plot=p, width=11, height=11)

	p2 <- SHmap + 
		geom_path(data= all_tracks, aes(lon, lat, col=all_tracks$stage), lineend="round") #+ 
		#opts(legend.position="left", legend.text = theme_text(colour="black", size = 6)) +
		#geom_point(data=data.frame(res_all$ilon,res_all$ilat), aes(res_all$ilon,res_all$ilat), col="black", pch=3, cex=.9)
	#p
	ggsave("all.google.trax_year.pdf", plot=p2, width=11, height=11)
}

# the interactive plot3d version: not easy to read tho -- or to save!
plot3d(as.numeric(coord_all[,3]), as.numeric(coord_all[,4]), as.numeric(coord_all[,2]), type="l", col=rainbow(1000), xlab="Longitude", ylab="Latitude", zlab="Bird ID")
rgl.postscript("persp3dd.pdf","pdf")

coord_all[which(coord_all[,4]>53),] # line 6711
coord_all[which(coord_all[,4]>53),c(3,4,6,8)] <- NA

# the interactive plot3d version: not easy to read tho -- or to save!
plot3d(as.numeric(coord_all[,3]), as.numeric(coord_all[,4]), as.numeric(coord_all[,2]), type="l", col=rainbow(1000), xlab="Longitude", ylab="Latitude", zlab="Bird ID")
rgl.postscript("persp3dd.pdf","pdf")

## ~~~~ do not run until line 1420  ~~~~ ##

#####################################
# temperature / tide / IQP heatmaps #
#####################################

lon.min <- min(res_all$ilon,na.rm=T)
lon.max <- max(res_all$ilon,na.rm=T)
lat.min <- min(res_all$ilat,na.rm=T)
lat.max <- max(res_all$ilat,na.rm=T)

grid.size <- 100
lon.grid <- (lon.max-lon.min)/grid.size
lat.grid <- (lat.max-lat.min)/grid.size

temp.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
lon.breaks <- seq(lon.min, lon.max, lon.grid)
lat.breaks <- seq(lat.min, lat.max, lat.grid)

for(i in 1:grid.size){
	for(j in 1:grid.size){
		temp.grid[i,j] <- mean(res_all$temperature[res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]], na.rm=T)
	}
}
temp.grid[temp.grid=="NaN"] <- NA

# these lines are adapted from: http://www.r-bloggers.com/controlling-heatmap-colors-with-ggplot2/
quantile_range <- quantile(temp.grid, probs = seq(0, 1, 0.2),na.rm=T)
color_palette <- colorRampPalette(c("#3794bf", "#FFFFFF", "#df8640"))(length(quantile_range) - 1)
label_text <- rollapply(round(quantile_range, 2), width = 2, by = 1, FUN = function(i) paste(i, collapse = " : "))
mod_mat <- matrix(findInterval(temp.grid, quantile_range, all.inside = TRUE), nrow = nrow(temp.grid))

p <- ggplot(melt(mod_mat), aes(x = X1, y = X2, fill = factor(value))) +
geom_tile(color = "grey70") + theme_bw() +
scale_fill_manual(values = color_palette, name = "", labels = label_text) +
scale_x_continuous(name=paste("Longitude (rescaled from ",format(lon.min,digits=3)," to ", format(lon.max,digits=3),")",sep="")) + 
scale_y_continuous(name=paste("Latitude (rescaled from ",format(lat.min,digits=3)," to ", format(lat.max,digits=3),")",sep=""))
ggsave(paste("temp.map.",grid.size,"x",grid.size,".pdf",sep=""), plot= p, width=7, height=4.5)


temp.threshold <- 20 # deg C
### for 2011 ###
temp.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
for(i in 1:grid.size){
	for(j in 1:grid.size){
		temp.grid[i,j] <- mean(res_all$temperature[res_all$year==2011 & res_all$temperature<temp.threshold & res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]], na.rm=T)
	}
}
temp.grid[temp.grid=="NaN"] <- NA
# this is a workaround for the ggplot2 bug... (will create boggus pixels in top left corner)
temp.grid[1:length(quantile_range),grid.size] <- quantile_range
mod_mat <- matrix(findInterval(temp.grid, quantile_range, all.inside = TRUE), nrow = nrow(temp.grid))
t_2011 <- ggplot(melt(mod_mat), aes(x = X1, y = X2, fill = factor(value))) +
geom_tile(color = "grey70") + theme_bw() +
scale_fill_manual(values = color_palette, name = "", labels = label_text) +
scale_x_continuous(name=paste("Longitude (rescaled from ",format(lon.min,digits=3)," to ", format(lon.max,digits=3),")",sep="")) + 
scale_y_continuous(name=paste("Latitude (rescaled from ",format(lat.min,digits=3)," to ", format(lat.max,digits=3),")",sep="")) +
ggtitle("Temperature in 2011") +
theme(panel.background = element_rect(fill = "grey70", colour = NA), panel.grid.major = element_line(colour = "grey70"), panel.grid.minor = element_line(colour = "grey70"))
ggsave(paste("temp.map.2011.",grid.size,"x",grid.size,".pdf",sep=""), plot= t_2011, width=7, height=4.5)
### for 2012 ###
temp.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
for(i in 1:grid.size){
	for(j in 1:grid.size){
		temp.grid[i,j] <- mean(res_all$temperature[res_all$year==2012 & res_all$temperature<temp.threshold & res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]], na.rm=T)
	}
}
temp.grid[temp.grid=="NaN"] <- NA
# this is a workaround for the ggplot2 bug... (will create boggus pixels in top left corner)
temp.grid[1:length(quantile_range),grid.size] <- quantile_range
mod_mat <- matrix(findInterval(temp.grid, quantile_range, all.inside = TRUE), nrow = nrow(temp.grid))
t_2012 <- ggplot(melt(mod_mat), aes(x = X1, y = X2, fill = factor(value))) +
geom_tile(color = "grey70") + theme_bw() +
scale_fill_manual(values = color_palette, name = "", labels = label_text) +
scale_x_continuous(name=paste("Longitude (rescaled from ",format(lon.min,digits=3)," to ", format(lon.max,digits=3),")",sep="")) + 
scale_y_continuous(name=paste("Latitude (rescaled from ",format(lat.min,digits=3)," to ", format(lat.max,digits=3),")",sep="")) +
ggtitle("Temperature in 2012") +
theme(panel.background = element_rect(fill = "grey70", colour = NA), panel.grid.major = element_line(colour = "grey70"), panel.grid.minor = element_line(colour = "grey70"))
ggsave(paste("temp.map.2012.",grid.size,"x",grid.size,".pdf",sep=""), plot= t_2012, width=7, height=4.5)
### for 2013 ###
temp.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
for(i in 1:grid.size){
	for(j in 1:grid.size){
		temp.grid[i,j] <- mean(res_all$temperature[res_all$year==2013 & res_all$temperature<temp.threshold & res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]], na.rm=T)
	}
}
temp.grid[temp.grid=="NaN"] <- NA
# this is a workaround for the ggplot2 bug... (will create boggus pixels in top left corner)
temp.grid[1:length(quantile_range),grid.size] <- quantile_range
mod_mat <- matrix(findInterval(temp.grid, quantile_range, all.inside = TRUE), nrow = nrow(temp.grid))
t_2013 <- ggplot(melt(mod_mat), aes(x = X1, y = X2, fill = factor(value))) +
geom_tile(color = "grey70") + theme_bw() +
scale_fill_manual(values = color_palette, name = "", labels = label_text) +
scale_x_continuous(name=paste("Longitude (rescaled from ",format(lon.min,digits=3)," to ", format(lon.max,digits=3),")",sep="")) + 
scale_y_continuous(name=paste("Latitude (rescaled from ",format(lat.min,digits=3)," to ", format(lat.max,digits=3),")",sep="")) +
ggtitle("Temperature in 2013") +
theme(panel.background = element_rect(fill = "grey70", colour = NA), panel.grid.major = element_line(colour = "grey70"), panel.grid.minor = element_line(colour = "grey70"))
ggsave(paste("temp.map.2013.",grid.size,"x",grid.size,".pdf",sep=""), plot= t_2013, width=7, height=4.5)

### for incubating birds ###
temp.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
for(i in 1:grid.size){
	for(j in 1:grid.size){
		temp.grid[i,j] <- mean(res_all$temperature[res_all$stage=="egg" & res_all$temperature<temp.threshold & res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]], na.rm=T)
	}
}
temp.grid[temp.grid=="NaN"] <- NA
# this is a workaround for the ggplot2 bug... (will create boggus pixels in top left corner)
temp.grid[1:length(quantile_range),grid.size] <- quantile_range
mod_mat <- matrix(findInterval(temp.grid, quantile_range, all.inside = TRUE), nrow = nrow(temp.grid))
t_egg <- ggplot(melt(mod_mat), aes(x = X1, y = X2, fill = factor(value))) +
geom_tile(color = "grey70") + theme_bw() +
scale_fill_manual(values = color_palette, name = "", labels = label_text) +
scale_x_continuous(name=paste("Longitude (rescaled from ",format(lon.min,digits=3)," to ", format(lon.max,digits=3),")",sep="")) + 
scale_y_continuous(name=paste("Latitude (rescaled from ",format(lat.min,digits=3)," to ", format(lat.max,digits=3),")",sep="")) +
ggtitle("Temperature during incubation") +
theme(panel.background = element_rect(fill = "grey70", colour = NA), panel.grid.major = element_line(colour = "grey70"), panel.grid.minor = element_line(colour = "grey70"))
ggsave(paste("temp.map.egg.",grid.size,"x",grid.size,".pdf",sep=""), plot= t_egg, width=7, height=4.5)
### for chich-rearing birds ###
temp.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
for(i in 1:grid.size){
	for(j in 1:grid.size){
		temp.grid[i,j] <- mean(res_all$temperature[res_all$stage=="chick" & res_all$temperature<temp.threshold & res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]], na.rm=T)
	}
}
temp.grid[temp.grid=="NaN"] <- NA
# this is a workaround for the ggplot2 bug... (will create boggus pixels in top left corner)
temp.grid[1:length(quantile_range),grid.size] <- quantile_range
mod_mat <- matrix(findInterval(temp.grid, quantile_range, all.inside = TRUE), nrow = nrow(temp.grid))
t_chick <- ggplot(melt(mod_mat), aes(x = X1, y = X2, fill = factor(value))) +
geom_tile(color = "grey70") + theme_bw() +
scale_fill_manual(values = color_palette, name = "", labels = label_text) +
scale_x_continuous(name=paste("Longitude (rescaled from ",format(lon.min,digits=3)," to ", format(lon.max,digits=3),")",sep="")) + 
scale_y_continuous(name=paste("Latitude (rescaled from ",format(lat.min,digits=3)," to ", format(lat.max,digits=3),")",sep="")) +
ggtitle("Temperature during chich rearing") +
theme(panel.background = element_rect(fill = "grey70", colour = NA), panel.grid.major = element_line(colour = "grey70"), panel.grid.minor = element_line(colour = "grey70"))
ggsave(paste("temp.map.chick.",grid.size,"x",grid.size,".pdf",sep=""), plot= t_chick, width=7, height=4.5)




# same grid treatment for IPQ
ipq.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
for(i in 1:grid.size){
	for(j in 1:grid.size){
		ipq.grid[i,j] <- mean(exp(res_all$ipq[res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]]), na.rm=T)
	}
}
ipq.grid[ipq.grid=="NaN"] <- NA
quantile_range <- quantile(ipq.grid, probs = seq(0, 1, 0.2),na.rm=T)
color_palette <- colorRampPalette(c("#3794bf", "#FFFFFF", "#df8640"))(length(quantile_range) - 1)
label_text <- rollapply(round(quantile_range, 2), width = 2, by = 1, FUN = function(i) paste(i, collapse = " : "))
### 2011 ###
ipq.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
for(i in 1:grid.size){
	for(j in 1:grid.size){
		ipq.grid[i,j] <- mean(exp(res_all$ipq[res_all$year==2011 & res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]]), na.rm=T)
	}
}
ipq.grid[ipq.grid=="NaN"] <- NA
# this is a workaround for the ggplot2 bug... (will create boggus pixels in top left corner)
ipq.grid[1:length(quantile_range),grid.size] <- quantile_range
mod_mat <- matrix(findInterval(ipq.grid, quantile_range, all.inside = TRUE), nrow = nrow(ipq.grid))
ipq_2011 <- ggplot(melt(mod_mat), aes(x = X1, y = X2, fill = factor(value))) +
geom_tile(color = "grey70") + theme_bw() +
scale_fill_manual(values = color_palette, name = "", labels = label_text) +
scale_x_continuous(name=paste("Longitude (rescaled from ",format(lon.min,digits=3)," to ", format(lon.max,digits=3),")",sep="")) + 
scale_y_continuous(name=paste("Latitude (rescaled from ",format(lat.min,digits=3)," to ", format(lat.max,digits=3),")",sep="")) +
ggtitle("IPQ in 2011 (not on log scale)") +
theme(panel.background = element_rect(fill = "grey70", colour = NA), panel.grid.major = element_line(colour = "grey70"), panel.grid.minor = element_line(colour = "grey70"))
### 2012 ###
ipq.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
for(i in 1:grid.size){
	for(j in 1:grid.size){
		ipq.grid[i,j] <- mean(exp(res_all$ipq[res_all$year==2012 & res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]]), na.rm=T)
	}
}
ipq.grid[ipq.grid=="NaN"] <- NA
# this is a workaround for the ggplot2 bug... (will create boggus pixels in top left corner)
ipq.grid[1:length(quantile_range),grid.size] <- quantile_range
mod_mat <- matrix(findInterval(ipq.grid, quantile_range, all.inside = TRUE), nrow = nrow(ipq.grid))
ipq_2012 <- ggplot(melt(mod_mat), aes(x = X1, y = X2, fill = factor(value))) +
geom_tile(color = "grey70") + theme_bw() +
scale_fill_manual(values = color_palette, name = "", labels = label_text) +
scale_x_continuous(name=paste("Longitude (rescaled from ",format(lon.min,digits=3)," to ", format(lon.max,digits=3),")",sep="")) + 
scale_y_continuous(name=paste("Latitude (rescaled from ",format(lat.min,digits=3)," to ", format(lat.max,digits=3),")",sep="")) +
ggtitle("IPQ in 2012 (not on log scale)") +
theme(panel.background = element_rect(fill = "grey70", colour = NA), panel.grid.major = element_line(colour = "grey70"), panel.grid.minor = element_line(colour = "grey70"))
### 2013 ###
ipq.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
for(i in 1:grid.size){
	for(j in 1:grid.size){
		ipq.grid[i,j] <- mean(exp(res_all$ipq[res_all$year==2013 & res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]]), na.rm=T)
	}
}
ipq.grid[ipq.grid=="NaN"] <- NA
# this is a workaround for the ggplot2 bug... (will create boggus pixels in top left corner)
ipq.grid[1:length(quantile_range),grid.size] <- quantile_range
mod_mat <- matrix(findInterval(ipq.grid, quantile_range, all.inside = TRUE), nrow = nrow(ipq.grid))
ipq_2013 <- ggplot(melt(mod_mat), aes(x = X1, y = X2, fill = factor(value))) +
geom_tile(color = "grey70") + theme_bw() +
scale_fill_manual(values = color_palette, name = "", labels = label_text) +
scale_x_continuous(name=paste("Longitude (rescaled from ",format(lon.min,digits=3)," to ", format(lon.max,digits=3),")",sep="")) + 
scale_y_continuous(name=paste("Latitude (rescaled from ",format(lat.min,digits=3)," to ", format(lat.max,digits=3),")",sep="")) +
ggtitle("IPQ in 2013 (not on log scale)") +
theme(panel.background = element_rect(fill = "grey70", colour = NA), panel.grid.major = element_line(colour = "grey70"), panel.grid.minor = element_line(colour = "grey70"))

# putting temp and IPQ together into same figure
p <- arrangeGrob(t_2011, t_2012, t_2013, ipq_2011, ipq_2012, ipq_2013, ncol=3)
#grid.draw(p)
ggsave("temp_ipq_year_grid.pdf", plot=p, width=20, height=8)


# same grid treatment for IPQ: stage effects
ipq.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
for(i in 1:grid.size){
	for(j in 1:grid.size){
		ipq.grid[i,j] <- mean(exp(res_all$ipq[res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]]), na.rm=T)
	}
}
ipq.grid[ipq.grid=="NaN"] <- NA
quantile_range <- quantile(ipq.grid, probs = seq(0, 1, 0.2),na.rm=T)
color_palette <- colorRampPalette(c("#3794bf", "#FFFFFF", "#df8640"))(length(quantile_range) - 1)
label_text <- rollapply(round(quantile_range, 2), width = 2, by = 1, FUN = function(i) paste(i, collapse = " : "))
### egg ###
ipq.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
for(i in 1:grid.size){
	for(j in 1:grid.size){
		ipq.grid[i,j] <- mean(exp(res_all$ipq[res_all$stage=="egg" & res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]]), na.rm=T)
	}
}
ipq.grid[ipq.grid=="NaN"] <- NA
# this is a workaround for the ggplot2 bug... (will create boggus pixels in top left corner)
ipq.grid[1:length(quantile_range),grid.size] <- quantile_range
mod_mat <- matrix(findInterval(ipq.grid, quantile_range, all.inside = TRUE), nrow = nrow(ipq.grid))
ipq_egg <- ggplot(melt(mod_mat), aes(x = X1, y = X2, fill = factor(value))) +
geom_tile(color = "grey70") + theme_bw() +
scale_fill_manual(values = color_palette, name = "", labels = label_text) +
scale_x_continuous(name=paste("Longitude (rescaled from ",format(lon.min,digits=3)," to ", format(lon.max,digits=3),")",sep="")) + 
scale_y_continuous(name=paste("Latitude (rescaled from ",format(lat.min,digits=3)," to ", format(lat.max,digits=3),")",sep="")) +
ggtitle("IPQ during incubation (not on log scale)") +
theme(panel.background = element_rect(fill = "grey70", colour = NA), panel.grid.major = element_line(colour = "grey70"), panel.grid.minor = element_line(colour = "grey70"))
### chick ###
ipq.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
for(i in 1:grid.size){
	for(j in 1:grid.size){
		ipq.grid[i,j] <- mean(exp(res_all$ipq[res_all$stage=="chick" & res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]]), na.rm=T)
	}
}
ipq.grid[ipq.grid=="NaN"] <- NA
# this is a workaround for the ggplot2 bug... (will create boggus pixels in top left corner)
ipq.grid[1:length(quantile_range),grid.size] <- quantile_range
mod_mat <- matrix(findInterval(ipq.grid, quantile_range, all.inside = TRUE), nrow = nrow(ipq.grid))
ipq_chick <- ggplot(melt(mod_mat), aes(x = X1, y = X2, fill = factor(value))) +
geom_tile(color = "grey70") + theme_bw() +
scale_fill_manual(values = color_palette, name = "", labels = label_text) +
scale_x_continuous(name=paste("Longitude (rescaled from ",format(lon.min,digits=3)," to ", format(lon.max,digits=3),")",sep="")) + 
scale_y_continuous(name=paste("Latitude (rescaled from ",format(lat.min,digits=3)," to ", format(lat.max,digits=3),")",sep="")) +
ggtitle("IPQ during chick rearing (not on log scale)") +
theme(panel.background = element_rect(fill = "grey70", colour = NA), panel.grid.major = element_line(colour = "grey70"), panel.grid.minor = element_line(colour = "grey70"))


# putting temp and IPQ together into same figure
pp <- arrangeGrob(t_egg, t_chick, ipq_egg, ipq_chick, ncol=2)
#grid.draw(p)
ggsave("temp_ipq_stage_grid.pdf", plot=pp, width=16, height=8)



vect.temp <- unmatrix(temp.grid)
vect.ipq <- unmatrix(ipq.grid)
ti.lm <- lmRob(vect.ipq ~ vect.temp)
summary(ti.lm)
#Call: lmRob(formula = vect.ipq ~ vect.temp)
#
#Residuals:
#         Min           1Q       Median           3Q          Max 
#-0.171080449 -0.051604404 -0.003117383  0.062016587  0.374656449 
#
#Coefficients:
#            Value         Std. Error    t value       Pr(>|t|)     
#(Intercept)  0.2483156112  0.0650744648  3.8158686673  0.0001483587
#vect.temp   -0.0060919342  0.0056925275 -1.0701633406  0.2849340575
#
#Residual standard error: 0.0866475 on 665 degrees of freedom
#9333 observations deleted due to missingness 
#Multiple R-Squared: 0.00156161 
#
#Test for Bias:
#            statistic      p-value
#M-estimate   7.631335 2.202301e-02
#LS-estimate 36.515707 1.176833e-08
plot(vect.temp, vect.ipq, xlim=c(10,16), pch=16, cex=.7)
abline(ti.lm, col="red")

vect.temp <- unmatrix(temp.grid)
vect.ipq <- unmatrix(log(ipq.grid))
ti.lm <- lmRob(vect.ipq ~ vect.temp)
summary(ti.lm)
#Call: lmRob(formula = vect.ipq ~ vect.temp)
#
#Residuals:
#       Min         1Q     Median         3Q        Max 
#-2.6651814 -0.3414172 -0.0119659  0.3019333  1.1745307 
#
#Coefficients:
#            Value         Std. Error    t value       Pr(>|t|)     
#(Intercept) -1.4951477232  0.3826150128 -3.9077079392  0.0001027129
#vect.temp   -0.0200389638  0.0334323651 -0.5993881611  0.5491182183
#
#Residual standard error: 0.4879 on 665 degrees of freedom
#9333 observations deleted due to missingness 
#Multiple R-Squared: 0.000321126 
#
#Test for Bias:
#            statistic    p-value
#M-estimate   6.051431 0.04852308
#LS-estimate 80.916012 0.00000000


#######################################################################
# plots by year of individual response to temperature in terms of IPQ #
#######################################################################
indiv.ipq <- function(birdID){
	# temps
	temp.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
	for(i in 1:grid.size){
		for(j in 1:grid.size){
			temp.grid[i,j] <- mean(res_all$temperature[res_all$ID==birdID & res_all$temperature < temp.threshold & res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]], na.rm=T)
		}
	}
	temp.grid[temp.grid=="NaN"] <- NA
	# IPQs
	ipq.grid <- matrix(rep(NA, grid.size^2),nrow= grid.size)
	for(i in 1:grid.size){
		for(j in 1:grid.size){
			ipq.grid[i,j] <- mean(res_all$ipq[res_all$ID==birdID & res_all$ipq>-1 & res_all$ipq < 1 & res_all$ilon>=lon.breaks[i] & res_all$ilon<lon.breaks[i+1] & res_all$ilat>=lat.breaks[j] & res_all$ilat<lat.breaks[j+1]], na.rm=T)
		}
	}
	ipq.grid[ipq.grid=="NaN"] <- NA
	# convert to vectors
	vect.temp <- unmatrix(temp.grid)
	vect.ipq <- unmatrix(ipq.grid)
	#ti.lm <- lmRob(vect.ipq ~ vect.temp)
	#summary(ti.lm)
	birdid <- rep(birdID, length(vect.temp))
	#return(birdid, vect.temp, vect.ipq, ti.lm)
	return(data.frame(birdid, vect.ipq, vect.temp))
}
# what bird what year
birds.2011 <- as.numeric(levels(as.factor(res_all$ID[res_all$year==2011])))
ipq.temp2011 <- matrix(nrow = 1, ncol = 3, byrow=T) 
colnames(ipq.temp2011) <- c("ID","ipq","temp")
for(i in birds.2011){
	print(paste("Doing bird ",i,sep=""))
	tmp_df <- indiv.ipq(i)
	colnames(tmp_df) <- c("ID","ipq","temp")
	ipq.temp2011  <- rbind(ipq.temp2011, tmp_df)
	# eliminates NA values
	ipq.temp2011  <- ipq.temp2011[which(!is.na(ipq.temp2011$ipq)),]
}
# 2011
ipqTemp_2011 <- ggplot(ipq.temp2011, aes(x = temp, y = ipq, col=factor(ID), fill = factor(ID))) + 
geom_point() + geom_smooth(method="rlm") +
theme_bw() +
scale_x_continuous(name="Temperature") + 
scale_y_continuous(name="IPQ") + 
ggtitle("log IPQ vs. temperature in 2011")

birds.2012 <- as.numeric(levels(as.factor(res_all$ID[res_all$year==2012])))
ipq.temp2012 <- matrix(nrow = 1, ncol = 3, byrow=T) 
colnames(ipq.temp2012) <- c("ID","ipq","temp")
for(i in birds.2012){
	print(paste("Doing bird ",i,sep=""))
	tmp_df <- indiv.ipq(i)
	colnames(tmp_df) <- c("ID","ipq","temp")
	ipq.temp2012  <- rbind(ipq.temp2012, tmp_df)
	# eliminates NA values
	ipq.temp2012  <- ipq.temp2012[which(!is.na(ipq.temp2012$ipq)),]
}
# 2012
ipqTemp_2012 <- ggplot(ipq.temp2012, aes(x = temp, y = ipq, col=factor(ID), fill = factor(ID))) + 
geom_point() + geom_smooth(method="rlm") +
theme_bw() +
scale_x_continuous(name="Temperature") + 
scale_y_continuous(name="IPQ") + 
ggtitle("log IPQ vs. temperature in 2012")

birds.2013 <- as.numeric(levels(as.factor(res_all$ID[res_all$year==2013])))
ipq.temp2013 <- matrix(nrow = 1, ncol = 3, byrow=T) 
colnames(ipq.temp2013) <- c("ID","ipq","temp")
for(i in birds.2013){
	print(paste("Doing bird ",i,sep=""))
	tmp_df <- indiv.ipq(i)
	colnames(tmp_df) <- c("ID","ipq","temp")
	ipq.temp2013  <- rbind(ipq.temp2013, tmp_df)
	# eliminates NA values
	ipq.temp2013  <- ipq.temp2013[which(!is.na(ipq.temp2013$ipq)),]
}
# 2013
ipqTemp_2013 <- ggplot(ipq.temp2013, aes(x = temp, y = ipq, col=factor(ID), fill = factor(ID))) + 
geom_point() + geom_smooth(method="rlm") +
theme_bw() +
scale_x_continuous(name="Temperature") + 
scale_y_continuous(name="IPQ") + 
ggtitle("log IPQ vs. temperature in 2013")

# the ANOVAs
aov2011 <- aov(ipq.temp2011$ipq ~ ipq.temp2011$ID + ipq.temp2011$temp)
summary.lm(aov2011)
aov2012 <- aov(ipq.temp2012$ipq ~ ipq.temp2012$ID + ipq.temp2012$temp)
summary.lm(aov2012)
aov2013 <- aov(ipq.temp2013$ipq ~ ipq.temp2013$ID + ipq.temp2013$temp)
summary.lm(aov2013)
summary(lmRob(ipq.temp2011$ipq ~ ipq.temp2011$temp))
summary(lmRob(ipq.temp2012$ipq ~ ipq.temp2012$temp))
summary(lmRob(ipq.temp2013$ipq ~ ipq.temp2013$temp))
summary(lmRob(c(ipq.temp2011$ipq,ipq.temp2012$ipq,ipq.temp2013$ipq) ~ c(ipq.temp2011$temp,ipq.temp2012$temp,ipq.temp2013$temp)))
# putting temp and IPQ together into same figure
ppp <- arrangeGrob(ipqTemp_2011, ipqTemp_2012, ipqTemp_2013, ncol=3)
ggsave("temp_regress_year.pdf", plot=ppp, width=20, height=4)


birds.egg <- as.numeric(levels(as.factor(res_all$ID[res_all$stage=="egg"])))
ipq.tempegg <- matrix(nrow = 1, ncol = 3, byrow=T) 
colnames(ipq.tempegg) <- c("ID","ipq","temp")
for(i in birds.egg){
	print(paste("Doing bird ",i,sep=""))
	tmp_df <- indiv.ipq(i)
	colnames(tmp_df) <- c("ID","ipq","temp")
	ipq.tempegg  <- rbind(ipq.tempegg, tmp_df)
	# eliminates NA values
	ipq.tempegg  <- ipq.tempegg[which(!is.na(ipq.tempegg$ipq)),]
}
# egg
ipqTemp_egg <- ggplot(ipq.tempegg, aes(x = temp, y = ipq, col=factor(ID), fill = factor(ID))) + 
geom_point() + geom_smooth(method="rlm") +
theme_bw() +
scale_x_continuous(name="Temperature") + 
scale_y_continuous(name="IPQ") + 
ggtitle("log IPQ vs. temperature during incubation")

birds.chick <- as.numeric(levels(as.factor(res_all$ID[res_all$stage=="chick"])))
ipq.tempchick <- matrix(nrow = 1, ncol = 3, byrow=T) 
colnames(ipq.tempchick) <- c("ID","ipq","temp")
for(i in birds.chick){
	print(paste("Doing bird ",i,sep=""))
	tmp_df <- indiv.ipq(i)
	colnames(tmp_df) <- c("ID","ipq","temp")
	ipq.tempchick  <- rbind(ipq.tempchick, tmp_df)
	# eliminates NA values
	ipq.tempchick  <- ipq.tempchick[which(!is.na(ipq.tempchick$ipq)),]
}
# chick
ipqTemp_chick <- ggplot(ipq.tempchick, aes(x = temp, y = ipq, col=factor(ID), fill = factor(ID))) + 
geom_point() + geom_smooth(method="rlm") +
theme_bw() +
scale_x_continuous(name="Temperature") + 
scale_y_continuous(name="IPQ") + 
ggtitle("log IPQ vs. temperature during chick rearing")

# putting temp and IPQ together into same figure
p5 <- arrangeGrob(ipqTemp_egg, ipqTemp_chick, ncol=2)
ggsave("temp_regress_stage.pdf", plot=p5, width=16, height=4)
aovegg <- aov(ipq.tempegg$ipq ~ ipq.tempegg$ID + ipq.tempegg$temp)
summary.lm(aovegg)
aovchick <- aov(ipq.tempchick$ipq ~ ipq.tempchick$ID + ipq.tempchick$temp)
summary.lm(aovchick)


		######################################
		# START AGAIN RUNNING CODE FROM HERE #
		#            13/10/01                #
		######################################


####################################################
# Passage IDs  -- the last one needs to be removed #
####################################################
#threshold <- 1000 # unit is meters; 1000 seems to help "iron out" false positives
# this is now defined above...

threshold_plus <- threshold + 2000 #1200 #700

# basic aov
aov.dist <- aov(res_all$dist2col ~ as.factor(res_all$ID) + res_all$ilon + res_all$ilat + res_all$bottdist + res_all$bottdep.mean)
summary(aov.dist)


# compute dist to colony
ncoords <- length(coord_all[,1])
dist2col <- numeric(ncoords)
for(i in 1:ncoords){
	if((as.numeric(coord_all[i,5]) == 2011)|(as.numeric(coord_all[i,5]) == 2013)){
		dist2col[i] <- distMeeus(c(as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4])),col.loc1)
	}else{
		dist2col[i] <- distMeeus(c(as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4])),col.loc2)
	}
	#if(dist2col[i]>100000){
	#	warning(i)
	#}
}
#plot(dist2col, type="l")
#plot(dist2col, type="l", xlim=c(0,1000))
#plot(dist2col, pch=".", xlim=c(0,200))

#which(res_all$dist2col>100000) # 10313 10314 10315 10316
#res_all$dist2col[which(res_all$dist2col>100000)] <- NA

# for res_all
ncoords_resall <- length(res_all[,1])
dist2col_resall <- numeric(ncoords_resall)
for(i in 1:ncoords_resall){
	if((as.numeric(res_all[i,3]) == 2011)|(as.numeric(res_all[i,3]) == 2013)){
		dist2col_resall[i] <- distMeeus(c(as.numeric(res_all[i,6]),as.numeric(res_all[i,7])),col.loc1)
	}else{
		dist2col_resall[i] <- distMeeus(c(as.numeric(res_all[i,6]),as.numeric(res_all[i,7])),col.loc2)
	}
	#if(dist2col[i]>100000){
	#	warning(i)
	#}
}

# transform in km
res_all$dist2col <- dist2col_resall/1000
#res_all$dist2col <- res_all$dist2col*1000
threshold_plus_km <- threshold_plus/1000

which(is.na(res_all$dist2col))
which(is.na(dist2col))
dist2col[which(is.na(dist2col))-1]
dist2col[which(is.na(dist2col))+1]
dist2col[which(is.na(dist2col))] <- mean(c(dist2col[which(is.na(dist2col))-1],dist2col[which(is.na(dist2col))+1]))

# prelim visualization
pdf("passageIDs1.pdf", width=6, height=4)
plot(res_all$dist2col, type="l", col="gray", ylab="Distance to colony (km)")
#points(res_all$dist2col, pch=".", col=as.numeric(res_all$year))
points(res_all$dist2col, pch=".", col=as.numeric(res_all$ID))
abline(h=threshold_plus_km, lty=2)

hist(res_all$dist2col[res_all$dist2col> threshold_plus_km],50)
plot(density(res_all$dist2col[res_all$dist2col> threshold_plus_km],na.rm=T,bw = 2.5),main="",xlab="Distance to colony (km)")
abline(v=quantile(dist2col/1000, probs=.95,na.rm=T),lty=2)
dev.off()


bird_counter <- 1
counter   <- 0 # this assumes that the logger got turned on outside of colony AND outside of threshold_plus
cur.in.col <- T
passageID <- rep(NA,ncoords) 					# records IDs of passage for each bird
#passageDirection <- rep(NA,ncoords) 			# records out- (1) and in- (-1) bound trips
passageCounter <- rep(0,length(mytdrfilelist))	# records # of trips each bird took - last "trip" *included*
curPassageID <- 0
for(i in 1:ncoords){
	# check bird ID
	passageCounter[bird_counter] <- counter
	if(!(as.numeric(coord_all[i,2]) == bird_counter)){
		passageCounter[bird_counter] <- counter
		bird_counter <- as.numeric(coord_all[i,2])
		counter <- 0
		cur.in.col <- T
	}
#	if(passageID[i] < curPassageID){
#		passageCounter <- curPassageID
#	}else{
#		curPassageID <- passageID[i]
#	}
	# get passage IDs
	if(dist2col[i]>threshold_plus & cur.in.col){ # we're getting out of colony
		passageID[i] <- counter
		cur.in.col <- F
	}else{
		if(dist2col[i]>threshold_plus & !cur.in.col){
			passageID[i] <- counter
		}else{
			if(dist2col[i]<threshold_plus & !cur.in.col){
				cur.in.col <- T
				counter <- counter + 1
				
			}else{
				if(dist2col[i]>threshold_plus & cur.in.col){
				passageID[i] <- counter
				}
			}
		}
	}
}

passageID
#passageCounter

passageCounter  <- tapply(passageID,coord_all[,2],max,na.rm=T)
passageCounter
#rownames(passageCounter)
#ddd <- data.frame(tapply(passageID,coord_all[,2],max,na.rm=T))
#rownames(ddd) <- rownames(tapply(passageID,coord_all[,2],max,na.rm=T))


##############
# plot speed #
##############

pdf("speed.pdf", width=11, height=6.5)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,1))
plot(coord_all[,6], type="l", col="gray", ylab="speed (km/h)")
points(coord_all[,6], pch=".", col=as.numeric(coord_all[,2]))
abline(h=5,col="red",lwd=2,lty=2)
# hist of speed for nonzero passages -- includes the last passage
# which should probably be removed but too complicated right now, 
# so this is overestimating slightly low speed frequencies
#hist(as.numeric(coord_all[,6][passageID >0 & !is.na(passageID)]), nclass=50, xlab=paste("speed (km/h); dist2col threshold_plus = ", threshold_plus, "m", sep=""), main="")
#abline(v=5,col="red",lwd=2,lty=2)
hist(log(as.numeric(coord_all[,6][passageID >0 & !is.na(passageID)]),10), nclass=50, xlab="log10 speed (km/h)", main="")
abline(v=log(5,10),col="red",lwd=2,lty=2)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()

#speedlimit <- 20 # was defined on line 467
pdf("speed.5.pdf", width=11, height=6.5)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,1))
plot(coord_all[,6], type="l", col="gray", ylab="speed (km/h)")
points(coord_all[,6], pch=".", col=as.numeric(coord_all[,2]))
abline(h= speedlimit, lty=2)
# hist of speed for nonzero passages -- includes the last passage
# which should probably be removed but too complicated right now, 
# so this is overestimating slightly low speed frequencies
hist(as.numeric(coord_all[,6][passageID >0 & !is.na(passageID)])[as.numeric(coord_all[,6][passageID >0 & !is.na(passageID)])>speedlimit], nclass=50, xlab=paste("speed (km/h) > ", speedlimit,"; dist2col threshold_plus = ", threshold_plus, "m", sep=""), main="")
#plot(density(as.numeric(coord_all[,6][passageID >0 & !is.na(passageID)])[as.numeric(coord_all[,6][passageID >0 & !is.na(passageID)])>speedlimit]))
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()


mean_speeds <- tapply(as.numeric(coord_all[,6]), coord_all[,2], mean, na.rm=T)
sd_speeds <- tapply(as.numeric(coord_all[,6]), coord_all[,2], sd, na.rm=T)

#speedlimit <- .0 #20 # km/h
coord_all.df <- data.frame(cbind(coord_all[,1][as.numeric(coord_all[,6])>speedlimit],coord_all[,5][as.numeric(coord_all[,6])>speedlimit],coord_all[,6][as.numeric(coord_all[,6])>speedlimit],coord_all[,10][as.numeric(coord_all[,6])>speedlimit]),stringsAsFactors=F)
coord_all.df <- coord_all.df[complete.cases(coord_all.df),]
colnames(coord_all.df) <- c("bird","year","speed","stage")
p1 <- qplot(factor(bird),as.numeric(speed), fill=factor(year), data=coord_all.df, geom="boxplot", position="dodge")+theme_bw()+scale_x_discrete(name="Bird IDs")+scale_y_continuous(name="Speed (km/h)") + theme(axis.text.x=element_text(angle=-90))
p2 <- qplot(factor(year),as.numeric(speed), fill=factor(year), data=coord_all.df, geom="boxplot", position="dodge")+theme_bw()+scale_x_discrete(name="Year")+scale_y_continuous(name="Speed (km/h)")
p3 <- qplot(factor(stage),as.numeric(speed), fill=factor(stage), data=coord_all.df, geom="boxplot", position="dodge")+theme_bw()+scale_x_discrete(name="Stage")+scale_y_continuous(name="Speed (km/h)")

mp <- grid.arrange(p1, p2, p3, nrow=3)
pdf("speed.by.year.by.stage.pdf",width=5,height=12)
print(grid.arrange(p1, p2, p3, ncol=1))
dev.off()



################
# in/out times #
################
# this depends on threshold_plus as defined above (thru passageID)
# NB the last "trip" is never counted, as it should!
# all durations are in seconds

in.out.times <- matrix(nrow =sum(passageCounter), ncol = 25, byrow=T) 
colnames(in.out.times) <- c("name","birdID","passageID","year","goingouttime","comingbacktime","maxDist2col","TraveledDist_out","TraveledDist_back","totTravDist","meanIPQ","maxIPQ","timeTotal","timeDiving","timeFlying","timeFloating","i_out","i_backin","meanSpeed","nb_stops","nb_flights","nb_dives","TD_SPD","stage","birdName")
ipq_time <- data.frame(c1 <- numeric(0),c2 <- numeric(0),c3 <- numeric(0),c4 <- numeric(0),c5 <- numeric(0),c6 <- numeric(0),c7 <- numeric(0),c8 <- numeric(0),c9 <- numeric(0),c10 <- numeric(0),c11 <- numeric(0),c12 <- numeric(0))
colnames(ipq_time) <- c("BirdID","year","tripID","time","duration","ipq","lon","lat","dist2col","inbound","theta","speed")
counter <- 1
outindex <- 0
#ncoords <- length(coord_all[,1]) # already defined above!
fly_speed <-  speedlimit #10 # km/h
curBird <- "1"
curYear <- curLon <- curLat <- curSpeed <- curIPQ <- curStage <- NA
curTime <- goingouttime <- comingbacktime <- "" #res_all$begdesc[1]
nb_dives <- 0
start_coming_back <- 0
ncolor <- 7
mypalette<-rev(brewer.pal(ncolor,"BuPu"))
max_speed <- max(as.numeric(coord_all[,6]),na.rm=T)
for(i in 2:ncoords){
	t.time <- coord_all[i,7][[1]]
	# going out
	if((is.na(passageID[i-1])) & !(is.na(passageID[i])) & (passageID[i]>0) & (passageID[i] <= passageCounter[rownames(passageCounter) == curBird][[1]])){
		in.out.times[counter,1] <- coord_all[i,1]
		in.out.times[counter,2] <- curBird <- coord_all[i,2] # coord_all[i,2]
		in.out.times[counter,24] <- curStage<- coord_all[i,10] 
		in.out.times[counter,25] <- coord_all[i,11] 
		in.out.times[counter,3] <- passageID[i]
		in.out.times[counter,4] <- curYear <- as.numeric(coord_all[i,5])
		in.out.times[counter,5] <- goingouttime <- t.time
		in.out.times[counter,17] <- i
		in.out.times[counter,8] <- coord_all[i,8][[1]]
		curLon <- c(curLon,as.numeric(coord_all[i,3]))
		curLat <- c(curLat,as.numeric(coord_all[i,4]))
		curSpeed <- c(curSpeed,as.numeric(coord_all[i,6]))
		curTime <- c(curTime, t.time)
		curIPQ <- c(curIPQ,as.numeric(coord_all[i,9][[1]]))
		outindex <- i
		#tmp_ipq_time <- data.frame(curBird, curYear, counter, t.time, 0, as.numeric(coord_all[i,9][[1]]),as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4]),dist2col[i],1)
		tmp_ipq_time <- data.frame(curBird, curYear, passageID[i], t.time, 0, as.numeric(coord_all[i,9][[1]]),as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4]),dist2col[i],1,atan_d(as.numeric(coord_all[i,4])/as.numeric(coord_all[i,3])),as.numeric(coord_all[i,6]))
		colnames(tmp_ipq_time) <- colnames(ipq_time)
		ipq_time <- rbind(ipq_time, tmp_ipq_time)
		nb_dives <- nb_dives + 1
		#print(paste(counter,i, t.time))
	}else{
		# during the trip
		if(!(is.na(passageID[i-1])) & !(is.na(passageID[i])) & (t.time > goingouttime) & (passageID[i]>0) & (passageID[i] <= passageCounter[rownames(passageCounter) == curBird][[1]])){
			curLon <- c(curLon,as.numeric(coord_all[i,3]))
			curLat <- c(curLat,as.numeric(coord_all[i,4]))
			curSpeed <- c(curSpeed,as.numeric(coord_all[i,6]))
			curTime <- c(curTime, t.time)
			curIPQ <- c(curIPQ,as.numeric(coord_all[i,9][[1]]))
			#print(paste("[",counter,"]",i, t.time))
			if(dist2col[i]<dist2col[i-1] & start_coming_back==0){
				tmp_ipq_time <- data.frame(curBird, curYear, passageID[i], t.time, (as.numeric(t.time)-as.numeric(goingouttime)), as.numeric(coord_all[i,9][[1]]),as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4]),dist2col[i],1,atan_d(as.numeric(coord_all[i,4])/as.numeric(coord_all[i,3])),as.numeric(coord_all[i,6]))
				 start_coming_back <- 1
			}
			if(start_coming_back==0){
				tmp_ipq_time <- data.frame(curBird, curYear, passageID[i], t.time, (as.numeric(t.time)-as.numeric(goingouttime)), as.numeric(coord_all[i,9][[1]]),as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4]),dist2col[i],1,atan_d(as.numeric(coord_all[i,4])/as.numeric(coord_all[i,3])),as.numeric(coord_all[i,6]))
			}else{
				tmp_ipq_time <- data.frame(curBird, curYear, passageID[i], t.time, (as.numeric(t.time)-as.numeric(goingouttime)), as.numeric(coord_all[i,9][[1]]),as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4]),dist2col[i],0,atan_d(as.numeric(coord_all[i,4])/as.numeric(coord_all[i,3])),as.numeric(coord_all[i,6]))
			}
			colnames(tmp_ipq_time) <- colnames(ipq_time)
			ipq_time <- rbind(ipq_time, tmp_ipq_time)
			nb_dives <- nb_dives + 1
		}else{
			# coming back
			if((!(is.na(passageID[i-1]))) & (is.na(passageID[i])) & (passageID[i-1]>0) & (t.time > goingouttime) & (passageID[i-1] <= passageCounter[rownames(passageCounter) == curBird][[1]])){
					in.out.times[counter,6] <- comingbacktime <- coord_all[i-1,7][[1]]
					in.out.times[counter,18] <- (i-1)
					in.out.times[counter,7] <- max(dist2col[outindex:i])
					in.out.times[counter,9] <- coord_all[i-1,8][[1]]
					in.out.times[counter,10] <- as.numeric(coord_all[i-1,8]) - as.numeric(coord_all[outindex,8])
					in.out.times[counter,11] <- mean(curIPQ,na.rm=T) 
					in.out.times[counter,12] <- max(curIPQ,na.rm=T)
					in.out.times[counter,13] <- time_tot <- as.numeric(comingbacktime) - as.numeric(goingouttime)
					in.out.times[counter,14] <- time_div <- sum(res_all$divetim[(res_all$ID == curBird) & (res_all$begdesc >= goingouttime) & ((res_all$begdesc + res_all$divetim)< comingbacktime )])
					in.out.times[counter,19] <- mean(curSpeed,na.rm=T) 
					in.out.times[counter,23] <- time_div + sum(res_all$postdive.dur[(res_all$postdive.dur<100) & (res_all$ID == curBird) & (res_all$begdesc >= goingouttime) & ((res_all$begdesc + res_all$postdive.dur)< comingbacktime )])
					
					# eliminate runs < 10 steps in passageID
					if(as.numeric(in.out.times[counter,18]) - as.numeric(in.out.times[counter,17]) < 10){
						in.out.times[counter,1:16] <- NA
					}else{
						x <- as.numeric(curTime)
						y <- curSpeed
						z <- spline(x, y, method = "natural", n=length(seq(as.numeric(goingouttime),as.numeric(comingbacktime),by=1)))
						time_fly <- (table(z$y>fly_speed)[2] - sum(rle(z$y>fly_speed)$values))
						# number of stops, flights and dives, resp.
						in.out.times[counter,20] <- sum(rle(z$y<fly_speed)$values) 
						in.out.times[counter,21] <- sum(rle(z$y>fly_speed)$values)
						in.out.times[counter,22] <- nb_dives
						
						in.out.times[counter,15] <- time_fly
						in.out.times[counter,16] <- time_tot - (time_div + time_fly)
						
						pdf(paste("spdPrfl_B",curBird,"_T",counter,curStage,".pdf",sep=""), width=12, height=4)
						par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,2))
						plot(coord_all[,3][coord_all[,2]==curBird],coord_all[,4][coord_all[,2]==curBird],type="l", col="gray", xlab="Longtitude", ylab="Latitude", ylim=c(51.2, 52.4), xlim=c(-5.9,-5.1))
						symb <- rep(3,length(curSpeed))
						symb[curSpeed> fly_speed] <- 4
						symb[curSpeed> (2*fly_speed)] <- 5
						points(curLon, curLat,pch=symb,cex=.5)
						s <- seq(length(curLon)-1)
						arrows(curLon[s], curLat[s], curLon[s+1], curLat[s+1],col="gray20", length = 0.05, angle = 15)
						points(curLon, curLat,pch=16,cex=.3,col=mypalette[floor(ncolor*curSpeed/max_speed)+1])
						#curTime[which(diff(curTime,1)<60)] <- NA # to space out markings by n seconds
						curTime[curSpeed<fly_speed] <- NA
						text(curLon+.01, curLat+.01,format(unix2POSIXct(curTime), "%H:%M (%y/%m/%d)"),cex=.3)
						legend("topleft",col=mypalette, cex=.7, pch=3, c("low speed","2","3","4","5","6","high speed"))
						legend("topright",pch=c(3,4,5),cex=.7,c(paste("<= ",fly_speed," km/h",sep=""),paste("> ",fly_speed," km/h",sep=""),paste("> ",2*fly_speed," km/h",sep="")))
						plot(unix2POSIXct(as.numeric(x)),y,ylim=c(0,max_speed),type="l",xlab="Time of day",ylab="Speed (km/h)",main=paste("Bird ",coord_all[i-1,1]," (ID ",curBird,") trip ",counter,sep=""));abline(v=z$x[z$y>fly_speed],col="green")
						lines(unix2POSIXct(as.numeric(x)),y);abline(h= fly_speed,col="blue");lines(z,col="red")
						par(new=T)
						plot(unix2POSIXct(as.numeric(x)), curIPQ, axes=F, type="l", col="orange",xlab="",ylab="")
						axis(3,at=x[curSpeed<fly_speed],labels=F, tck=.02)
						par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
						dev.off()
					}
					
					counter <- counter + 1
					# debug below...
					if(counter>8){
					#	break
					}
	
					all_ipq <- NA
					curLon <- curLat <- curSpeed <- curIPQ <- NA
					curTime <- "" #res_all$begdesc[1]
					nb_dives <- start_coming_back <- 0
			}
		}
	}
}

# clean up passages time that should not be there
rm(x,y,z)
#colnames(ipq_time) <- c("BirdID","year","tripID","time","duration","ipq","lon","lat","dist2col","inbound","theta","speed")
in.out.times <- subset(in.out.times, !(is.na(in.out.times[,1])))
in.out.times <- subset(in.out.times, !(is.na(in.out.times[,16])))
in.out.times <- subset(in.out.times, in.out.times[,16]>0)
in.out.times

plot(ipq_time$dist2col,col=(1+ipq_time$inbound),pch=".")

# % time diving -- includes PDD
as.numeric(in.out.times[,23])/as.numeric(in.out.times[,13])
100*mean(as.numeric(in.out.times[,23])/as.numeric(in.out.times[,13])) # 12.96378%

# % time diving
as.numeric(in.out.times[,14])/as.numeric(in.out.times[,13])
100*mean(as.numeric(in.out.times[,14])/as.numeric(in.out.times[,13])) # 9.095171%
100*mean(as.numeric(in.out.times[,14][in.out.times[,4]=="2011"])/as.numeric(in.out.times[,13][in.out.times[,4]=="2011"])) # 9.094818%
100*mean(as.numeric(in.out.times[,14][in.out.times[,4]=="2012"])/as.numeric(in.out.times[,13][in.out.times[,4]=="2012"])) # 9.096041%
# % time flying
100*mean(as.numeric(in.out.times[,15])/as.numeric(in.out.times[,13])) # 15.67565%
# % time floating
100*mean(as.numeric(in.out.times[,16])/as.numeric(in.out.times[,13])) # 75.22918%
# trip duration (h)
mean(as.numeric(in.out.times[,13]))/3600
sd(as.numeric(in.out.times[,13]))/3600
# dive duration (h)
mean(as.numeric(in.out.times[,14]))/3600
sd(as.numeric(in.out.times[,14]))/3600
# flight duration (h)
mean(as.numeric(in.out.times[,15]))/3600
sd(as.numeric(in.out.times[,15]))/3600
# floating duration (h)
mean(as.numeric(in.out.times[,16]))/3600
sd(as.numeric(in.out.times[,16]))/3600


##########################
# outbound trip analysis #
##########################
# NB. ipq_time$inbound actually describes the outbound trip
# excludes time when speed < speedlimit, so that only flights are considered; issue: some trips become too short and need to be eliminated from matrix
ntrips <- length(in.out.times[,1])
outbound_sumlatlon <- matrix(nrow= ntrips,ncol= ntrips)
rownames(outbound_sumlatlon) <- colnames(outbound_sumlatlon) <- paste("B",in.out.times[,2],"_T",in.out.times[,3],"_Y",in.out.times[,4],"_",in.out.times[,24],sep="")
for(i in 1: ntrips){
	for(j in 1: ntrips){
		#zeros <- numeric(0)
		path1a <- ipq_time$lat[ipq_time$BirdID==in.out.times[i,2] & ipq_time$tripID==in.out.times[i,3] & ipq_time$inbound==1 & ipq_time$speed>speedlimit]
		path1b <- ipq_time$lon[ipq_time$BirdID==in.out.times[i,2] & ipq_time$tripID==in.out.times[i,3] & ipq_time$inbound==1 & ipq_time$speed>speedlimit]
		path2a <- ipq_time$lat[ipq_time$BirdID==in.out.times[j,2] & ipq_time$tripID==in.out.times[j,3] & ipq_time$inbound==1 & ipq_time$speed>speedlimit]
		path2b <- ipq_time$lon[ipq_time$BirdID==in.out.times[j,2] & ipq_time$tripID==in.out.times[j,3] & ipq_time$inbound==1 & ipq_time$speed>speedlimit]
		mycolmeans <- colMeans(matchpt(path1a+path1b,path2a+path2b))[[2]]
		if(!is.na(mycolmeans) & is.finite(mycolmeans)){
			outbound_sumlatlon[i,j] <- mycolmeans
		}else{
			outbound_sumlatlon[i,j] <- 0
		}
	}
}
sym_outbound_sumlatlon <- matrix(0,nrow= ntrips,ncol= ntrips)
rownames(sym_outbound_sumlatlon) <- colnames(sym_outbound_sumlatlon) <- paste("B",in.out.times[,2],"_T",in.out.times[,3],"_Y",in.out.times[,4],"_",in.out.times[,24],sep="")
for(i in 1: ntrips){
	for(j in 1:i){
		sym_outbound_sumlatlon[i,j] <- sym_outbound_sumlatlon[j,i] <- (outbound_sumlatlon[i,j]+ outbound_sumlatlon[j,i])/2
	}
}
zeros <- numeric(0)
for(i in 1:ntrips){
	if(all(sym_outbound_sumlatlon[,i]==0)){
		zeros <- c(zeros,i)
	}
}
sym_outbound_sumlatlon <- sym_outbound_sumlatlon[-zeros,]
sym_outbound_sumlatlon <- sym_outbound_sumlatlon[,-zeros]

pdf("heatmap_outbound_sumlatlon.pdf",width=17,height=20)
heatmap(((sym_outbound_sumlatlon)),col=greenred.colors(100),symm=T,distfun = dist, hclustfun = hclust)
dev.off()
sumlatlon.pv <- pvclust(sym_outbound_sumlatlon, nboot=1000, method.hclust="complete", method.dist="euclidean")
pdf("heatmap_outbound_sumlatlon_pvclust_col.pdf",width=13,height=8)
plot(sumlatlon.pv)
pvrect(sumlatlon.pv)
dev.off()
sumlatlon.pv <- pvclust(t(sym_outbound_sumlatlon), nboot=1000)
pdf("heatmap_outbound_sumlatlon_pvclust_row.pdf",width=13,height=8)
plot(sumlatlon.pv)
pvrect(sumlatlon.pv)
dev.off()


# for incubating birds
sym_outbound_sumlatlon_egg <- sym_outbound_sumlatlon
for(i in 1: dim(sym_outbound_sumlatlon_egg)[1]){
	for(j in 1: dim(sym_outbound_sumlatlon_egg)[2]){
		if( (length(grep("chick",colnames(sym_outbound_sumlatlon_egg)[i]))>0)|(length(grep("chick",rownames(sym_outbound_sumlatlon_egg)[j]))>0) ){
			sym_outbound_sumlatlon_egg[i,j] <- NA
		}
	}
}
mean(sym_outbound_sumlatlon_egg,na.rm=T)
# for chick rearing birds
sym_outbound_sumlatlon_chick <- sym_outbound_sumlatlon
for(i in 1: dim(sym_outbound_sumlatlon_chick)[1]){
	for(j in 1: dim(sym_outbound_sumlatlon_chick)[2]){
		if( (length(grep("egg",colnames(sym_outbound_sumlatlon_chick)[i]))>0)|(length(grep("egg",rownames(sym_outbound_sumlatlon_chick)[j]))>0) ){
			sym_outbound_sumlatlon_chick[i,j] <- NA
		}
	}
}
mean(sym_outbound_sumlatlon_chick,na.rm=T)
pdf("matchpt_outbound.pdf",width=6,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(density(as.vector(sym_outbound_sumlatlon_chick),na.rm=T),col="red",main="",xlab="MATCHPT distance for outbound trips")
lines(density(as.vector(sym_outbound_sumlatlon_egg),na.rm=T),col="blue")
legend(.25,20,c("incubation","chick rearing"),col=c("blue","red"),lty=1)
text(.3,10,paste("P = ",t.test(sym_outbound_sumlatlon_chick, sym_outbound_sumlatlon_egg)$p.value,sep=""))
abline(v=mean(sym_outbound_sumlatlon_egg,na.rm=T),col="blue",lty=2)
abline(v=mean(sym_outbound_sumlatlon_chick,na.rm=T),col="red",lty=2)
dev.off()
t.test(sym_outbound_sumlatlon_chick, sym_outbound_sumlatlon_egg)
ks.test(sym_outbound_sumlatlon_chick, sym_outbound_sumlatlon_egg)


##########################
# inbound trip analysis #
##########################
# NB. ipq_time$inbound actually describes the inbound trip
# excludes time when speed < speedlimit, so that only flights are considered; issue: some trips become too short and need to be eliminated from matrix
ntrips <- length(in.out.times[,1])
inbound_sumlatlon <- matrix(nrow= ntrips,ncol= ntrips)
rownames(inbound_sumlatlon) <- colnames(inbound_sumlatlon) <- paste("B",in.out.times[,2],"_T",in.out.times[,3],"_Y",in.out.times[,4],"_",in.out.times[,24],sep="")
for(i in 1: ntrips){
	for(j in 1: ntrips){
		#zeros <- numeric(0)
		path1a <- ipq_time$lat[ipq_time$BirdID==in.out.times[i,2] & ipq_time$tripID==in.out.times[i,3] & ipq_time$inbound==0 & ipq_time$speed>speedlimit]
		path1b <- ipq_time$lon[ipq_time$BirdID==in.out.times[i,2] & ipq_time$tripID==in.out.times[i,3] & ipq_time$inbound==0 & ipq_time$speed>speedlimit]
		path2a <- ipq_time$lat[ipq_time$BirdID==in.out.times[j,2] & ipq_time$tripID==in.out.times[j,3] & ipq_time$inbound==0 & ipq_time$speed>speedlimit]
		path2b <- ipq_time$lon[ipq_time$BirdID==in.out.times[j,2] & ipq_time$tripID==in.out.times[j,3] & ipq_time$inbound==0 & ipq_time$speed>speedlimit]
		mycolmeans <- colMeans(matchpt(path1a+path1b,path2a+path2b))[[2]]
		if(!is.na(mycolmeans) & is.finite(mycolmeans)){
			inbound_sumlatlon[i,j] <- mycolmeans
		}else{
			inbound_sumlatlon[i,j] <- 0
		}
	}
}
sym_inbound_sumlatlon <- matrix(0,nrow= ntrips,ncol= ntrips)
rownames(sym_inbound_sumlatlon) <- colnames(sym_inbound_sumlatlon) <- paste("B",in.out.times[,2],"_T",in.out.times[,3],"_Y",in.out.times[,4],"_",in.out.times[,24],sep="")
for(i in 1: ntrips){
	for(j in 1:i){
		sym_inbound_sumlatlon[i,j] <- sym_inbound_sumlatlon[j,i] <- (inbound_sumlatlon[i,j]+ inbound_sumlatlon[j,i])/2
	}
}
#zeros <- numeric(0)
#for(i in 1:ntrips){
#	if(all(sym_inbound_sumlatlon[,i]==0)){
#		zeros <- c(zeros,i)
#	}
#}
#sym_inbound_sumlatlon <- sym_inbound_sumlatlon[-zeros,]
#sym_inbound_sumlatlon <- sym_inbound_sumlatlon[,-zeros]

pdf("heatmap_inbound_sumlatlon.pdf",width=17,height=20)
heatmap(((sym_inbound_sumlatlon)),col=greenred.colors(100),symm=T,distfun = dist, hclustfun = hclust)
dev.off()
sumlatlon.pv <- pvclust(sym_inbound_sumlatlon, nboot=1000, method.hclust="complete", method.dist="euclidean")
pdf("heatmap_inbound_sumlatlon_pvclust_col.pdf",width=13,height=8)
plot(sumlatlon.pv)
pvrect(sumlatlon.pv)
dev.off()
sumlatlon.pv <- F(t(sym_inbound_sumlatlon), nboot=1000)
pdf("heatmap_inbound_sumlatlon_pvclust_row.pdf",width=13,height=8)
plot(sumlatlon.pv)
pvrect(sumlatlon.pv)
dev.off()


# for incubating birds
sym_inbound_sumlatlon_egg <- sym_inbound_sumlatlon
for(i in 1: dim(sym_inbound_sumlatlon_egg)[1]){
	for(j in 1: dim(sym_inbound_sumlatlon_egg)[2]){
		if( (length(grep("chick",colnames(sym_inbound_sumlatlon_egg)[i]))>0)|(length(grep("chick",rownames(sym_inbound_sumlatlon_egg)[j]))>0) ){
			sym_inbound_sumlatlon_egg[i,j] <- NA
		}
	}
}
mean(sym_inbound_sumlatlon_egg,na.rm=T)
# for chick rearing birds
sym_inbound_sumlatlon_chick <- sym_inbound_sumlatlon
for(i in 1: dim(sym_inbound_sumlatlon_chick)[1]){
	for(j in 1: dim(sym_inbound_sumlatlon_chick)[2]){
		if( (length(grep("egg",colnames(sym_inbound_sumlatlon_chick)[i]))>0)|(length(grep("egg",rownames(sym_inbound_sumlatlon_chick)[j]))>0) ){
			sym_inbound_sumlatlon_chick[i,j] <- NA
		}
	}
}
mean(sym_inbound_sumlatlon_chick,na.rm=T)
pdf("matchpt_inbound.pdf",width=6,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(density(as.vector(sym_inbound_sumlatlon_chick),na.rm=T),col="red",main="",xlab="MATCHPT distance for inbound trips")
lines(density(as.vector(sym_inbound_sumlatlon_egg),na.rm=T),col="blue")
legend(.25,15,c("incubation","chick rearing"),col=c("blue","red"),lty=1)
text(.3,10,paste("P = ",t.test(sym_inbound_sumlatlon_chick, sym_inbound_sumlatlon_egg)$p.value,sep=""))
abline(v=mean(sym_inbound_sumlatlon_egg,na.rm=T),col="blue",lty=2)
abline(v=mean(sym_inbound_sumlatlon_chick,na.rm=T),col="red",lty=2)
dev.off()
t.test(sym_inbound_sumlatlon_chick, sym_inbound_sumlatlon_egg)
ks.test(sym_inbound_sumlatlon_chick, sym_inbound_sumlatlon_egg)


#########################
# inbound trip analysis #
#########################
# based on ipq_time, testing lengths with dist2col
# summarizing matchpt by its mean 

ntrips <- length(in.out.times[,1])
inbound_dist2col <- matrix(nrow= ntrips,ncol= ntrips)
rownames(inbound_dist2col) <- colnames(inbound_dist2col) <- paste("Brd_",in.out.times[,2],"_Trp_",in.out.times[,3],"_Y",in.out.times[,4],"_",in.out.times[,24],sep="")
for(i in 1: ntrips){
	for(j in 1: ntrips){
		path1 <- ipq_time$dist2col[ipq_time$BirdID==in.out.times[i,2] & ipq_time$tripID==in.out.times[i,3] & ipq_time$inbound==1]
		path2 <- ipq_time$dist2col[ipq_time$BirdID==in.out.times[j,2] & ipq_time$tripID==in.out.times[j,3] & ipq_time$inbound==1]
		inbound_dist2col[i,j] <- colMeans(matchpt(path1,path2))[[2]]
	}
}
pdf("heatmap_dist2col.pdf",width=17,height=20)
heatmap(scale(t(inbound_dist2col)),col=greenred.colors(100))
dev.off()



ntrips <- length(in.out.times[,1])
inbound_rmsd <- matrix(nrow= ntrips,ncol= ntrips)
rownames(inbound_rmsd) <- colnames(inbound_rmsd) <- paste("Brd_",in.out.times[,2],"_Trp_",in.out.times[,3],"_Y",in.out.times[,4],"_",in.out.times[,24],sep="")
for(i in 1: ntrips){
	for(j in 1: ntrips){
		path1a <- ipq_time$lat[ipq_time$BirdID==in.out.times[i,2] & ipq_time$tripID==in.out.times[i,3] & ipq_time$inbound==1]
		path1b <- ipq_time$lon[ipq_time$BirdID==in.out.times[i,2] & ipq_time$tripID==in.out.times[i,3] & ipq_time$inbound==1]
		path2a <- ipq_time$lat[ipq_time$BirdID==in.out.times[j,2] & ipq_time$tripID==in.out.times[j,3] & ipq_time$inbound==1]
		path2b <- ipq_time$lon[ipq_time$BirdID==in.out.times[j,2] & ipq_time$tripID==in.out.times[j,3] & ipq_time$inbound==1]
		path1a <- path1a[1:min(length(path1a),length(path2a))]
		path1b <- path1b[1:min(length(path1b),length(path2b))]
		path2a <- path2a[1:min(length(path1a),length(path2a))]
		path2b <- path2b[1:min(length(path1b),length(path2b))]
		rmsd_sq <- (path1b-path2b)^2 + (path1a-path2a)^2
		inbound_rmsd[i,j] <- sqrt(mean(rmsd_sq))
	}
}
pdf("heatmap_rmsd.pdf",width=17,height=20)
heatmap(scale(t(inbound_rmsd)),col=greenred.colors(100),symm=T)
dev.off()

ntrips <- length(in.out.times[,1])
inbound_dist2col <- matrix(nrow= ntrips,ncol= ntrips)
rownames(inbound_dist2col) <- colnames(inbound_dist2col) <- paste("Brd_",in.out.times[,2],"_Trp_",in.out.times[,3],"_Y",in.out.times[,4],"_",in.out.times[,24],sep="")
for(i in 1: ntrips){
	for(j in 1: ntrips){
		path1 <- ipq_time$dist2col[ipq_time$BirdID==in.out.times[i,2] & ipq_time$tripID==in.out.times[i,3]]
		path2 <- ipq_time$dist2col[ipq_time$BirdID==in.out.times[j,2] & ipq_time$tripID==in.out.times[j,3]]
		inbound_dist2col[i,j] <- colMeans(matchpt(path1,path2))[[2]]
	}
}
pdf("heatmap_dist2col_entireTrack.pdf",width=17,height=20)
heatmap(scale(t(inbound_dist2col)),col=greenred.colors(100))
dev.off()

# based on ipq_time, testing initial directions with theta
# summarizing matchpt by its mean

ntrips <- length(in.out.times[,1])
inbound_theta <- matrix(nrow= ntrips,ncol= ntrips)
rownames(inbound_theta) <- colnames(inbound_theta) <- paste("Brd_",in.out.times[,2],"_Trp_",in.out.times[,3],"_Y",in.out.times[,4],"_",in.out.times[,24],sep="")
for(i in 1: ntrips){
	for(j in 1: ntrips){
		path1 <- mean(ipq_time$theta[ipq_time$BirdID==in.out.times[i,2] & ipq_time$tripID==in.out.times[i,3] & ipq_time$inbound==1])
		path2 <- mean(ipq_time$theta[ipq_time$BirdID==in.out.times[j,2] & ipq_time$tripID==in.out.times[j,3] & ipq_time$inbound==1])
		inbound_theta[i,j] <- abs(path1-path2)
	}
}
pdf("heatmap_theta.pdf",width=17,height=20)
heatmap(scale(t(inbound_theta)),col=greenred.colors(100),symm=T)
dev.off()

#ntrips <- length(in.out.times[,1])
#inbound_theta <- matrix(nrow= ntrips,ncol= ntrips)
#rownames(inbound_theta) <- colnames(inbound_theta) <- paste("Brd_",in.out.times[,2],"_Trp_",in.out.times[,3],"_Yr_",in.out.times[,4],sep="")
#for(i in 1: ntrips){
#	for(j in 1: ntrips){
#		path1 <- ipq_time$theta[ipq_time$BirdID==in.out.times[i,2] & ipq_time$tripID==in.out.times[i,3] & ipq_time$inbound==1]
#		path2 <- ipq_time$theta[ipq_time$BirdID==in.out.times[j,2] & ipq_time$tripID==in.out.times[j,3] & ipq_time$inbound==1]
#		inbound_theta[i,j] <- colMeans(matchpt(path1,path2))[[2]]
#	}
#}
#pdf("heatmap_theta.matchpt.pdf",width=17,height=20)
#heatmap(scale(t(inbound_theta)),col=greenred.colors(100))
#dev.off()

ntrips <- length(in.out.times[,1])
inbound_theta <- matrix(nrow= ntrips,ncol= ntrips)
rownames(inbound_theta) <- colnames(inbound_theta) <- paste("Brd_",in.out.times[,2],"_Trp_",in.out.times[,3],"_Y",in.out.times[,4],"_",in.out.times[,24],sep="")
for(i in 1: ntrips){
	for(j in 1: ntrips){
		path1 <- mean(ipq_time$theta[ipq_time$BirdID==in.out.times[i,2] & ipq_time$tripID==in.out.times[i,3]])
		path2 <- mean(ipq_time$theta[ipq_time$BirdID==in.out.times[j,2] & ipq_time$tripID==in.out.times[j,3]])
		inbound_theta[i,j] <- abs(path1-path2)
	}
}
pdf("heatmap_theta_entireTrack.pdf",width=17,height=20)
heatmap(scale(t(inbound_theta)),col=greenred.colors(100),symm=T)
dev.off()


# based on ipq_time
ntrips <- length(in.out.times[,1])
inbound_sumlatlon <- matrix(nrow= ntrips,ncol= ntrips)
rownames(inbound_sumlatlon) <- colnames(inbound_sumlatlon) <- paste("B",in.out.times[,2],"_T",in.out.times[,3],"_Y",in.out.times[,4],"_",in.out.times[,24],sep="")
for(i in 1: ntrips){
	for(j in 1: ntrips){
		path1a <- ipq_time$lat[ipq_time$BirdID==in.out.times[i,2] & ipq_time$tripID==in.out.times[i,3] & ipq_time$inbound==1]
		path1b <- ipq_time$lon[ipq_time$BirdID==in.out.times[i,2] & ipq_time$tripID==in.out.times[i,3] & ipq_time$inbound==1]
		path2a <- ipq_time$lat[ipq_time$BirdID==in.out.times[j,2] & ipq_time$tripID==in.out.times[j,3] & ipq_time$inbound==1]
		path2b <- ipq_time$lon[ipq_time$BirdID==in.out.times[j,2] & ipq_time$tripID==in.out.times[j,3] & ipq_time$inbound==1]
		inbound_sumlatlon[i,j] <- colMeans(matchpt(path1a+path1b,path2a+path2b))[[2]]
	}
}
pdf("heatmap_sumlatlon.pdf",width=17,height=20)
heatmap(scale(t(inbound_sumlatlon)),col=greenred.colors(100))
dev.off()

sumlatlon.pv <- pvclust(inbound_sumlatlon, nboot=1000)
pdf("heatmap_sumlatlon_pvclust_col.pdf",width=10,height=8)
plot(sumlatlon.pv)
pvrect(sumlatlon.pv)
dev.off()

sumlatlon.pv <- pvclust(t(inbound_sumlatlon), nboot=1000)
pdf("heatmap_sumlatlon_pvclust_row.pdf",width=10,height=8)
plot(sumlatlon.pv)
pvrect(sumlatlon.pv)
dev.off()



#library(MASS)
#data(Boston)
## multiscale bootstrap resampling
#boston.pv <- pvclust(Boston, nboot=100)
## CAUTION: nboot=100 may be too small for actual use.
## We suggest nboot=1000 or larger.
## plot/print functions will be useful for diagnostics.
## plot dendrogram with p-values
#plot(boston.pv)
#ask.bak <- par()$ask
#par(ask=TRUE)
## highlight clusters with high au p-values
#pvrect(boston.pv)
## print the result of multiscale bootstrap resampling
#print(boston.pv, digits=3)
## plot diagnostic for curve fitting
#msplot(sumlatlon.pv, edges=c(2,4,6,7))

#binom.test(8,13,alternative = c("greater"))




##########
# SAHFOS #
# 131122 #
##########

copepod <- read.table("/Users/stephane/Documents/data/akiko/131004/copepod_data.txt",header=T)

pdf("sahfos_data.pdf")
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
bathy <- getNOAA.bathy(lon1 = min(res_all$ilon,na.rm=T), lon2 = max(res_all$ilon,na.rm=T), lat1 = min(res_all$ilat,na.rm=T), lat2 = max(res_all$ilat,na.rm=T),resolution = 1)
smoothScatter(copepod$Sample_longitude, copepod$Sample_latitude,xlim=c(min(copepod$Sample_longitude,na.rm=T),max(copepod$Sample_longitude,na.rm=T)),ylim=c(min(copepod$Sample_latitude),max(copepod$Sample_latitude)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="SAHFOS data")
z <- kde2d(as.numeric(copepod$Sample_longitude), as.numeric(copepod$Sample_latitude), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
points(copepod$Sample_longitude, copepod$Sample_latitude, pch=20, cex=.5, col="blue")

roi_x1 <- -5.75; roi_y1 <- 51.65
roi_x2 <- -5.3; roi_y2 <- 51.825
rect(roi_x1, roi_y1, roi_x2, roi_y2, border ="red")
dev.off()

tmp_res_all <- res_all[complete.cases(res_all),]
tmp_coord_all <- coord_all[complete.cases(coord_all),]
pdf("activity_maps_colony_trimmed2_rect.pdf",width=10, height=10)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2))
# rest
bathy <- getNOAA.bathy(lon1 = min(tmp_res_all$ilon,na.rm=T), lon2 = max(tmp_res_all$ilon,na.rm=T), lat1 = min(tmp_res_all$ilat,na.rm=T), lat2 = max(tmp_res_all$ilat,na.rm=T),resolution = 1)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: fly/rest")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
rect(roi_x1, roi_y1, roi_x2, roi_y2, border ="red")
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: fly/rest")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
rect(roi_x1, roi_y1, roi_x2, roi_y2, border ="red")
# dive
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: dive")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
rect(roi_x1, roi_y1, roi_x2, roi_y2, border ="red")
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: dive")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
rect(roi_x1, roi_y1, roi_x2, roi_y2, border ="red")
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()
rm(z)



###########################################
# flight & dive patterns -- Levy/Gaussian #
# 131020  --  131119                      #
###########################################

# small changes of direction
coord_all_clean <- coord_all
step_len <- rep(NA,length(coord_all[,1]))
step_counter <- rep(0,length(coord_all[,1]))
all_alphas <- c()
# first determines steps
for(i in 2:(length(coord_all_clean[,1])-1)){
	# remove data where time difference betwee 2 consecutive GPS logs > 8 min
	if((unix2POSIXct(as.numeric(coord_all_clean[i,7])) - unix2POSIXct(as.numeric(coord_all_clean[i-1,7]))) > 8){
		coord_all_clean[i,8] <- NA
	}
	# if change of bird_ new step
	if(coord_all[i,2] != coord_all[i-1,2]){
		step_counter[i] <- step_counter[i-1] + 1
	}else{
		# if any location is NA, new trip and hence new step
		if(is.na(as.numeric(coord_all[i-1,3])) | is.na(as.numeric(coord_all[i,3])) | is.na(as.numeric(coord_all[i+1,3]))){
			step_counter[i] <- step_counter[i-1] + 1
		}else{
			# if dist traveled == 0: same step
			if(as.numeric(coord_all[i,3])==0){
				step_counter[i] <- step_counter[i-1]
			}else{
				a <- distMeeus(c(as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4])),c(as.numeric(coord_all[i+1,3]),as.numeric(coord_all[i+1,4])))
				b <- distMeeus(c(as.numeric(coord_all[i-1,3]),as.numeric(coord_all[i-1,4])),c(as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4])))
				c <- distMeeus(c(as.numeric(coord_all[i-1,3]),as.numeric(coord_all[i-1,4])),c(as.numeric(coord_all[i+1,3]),as.numeric(coord_all[i+1,4])))
				d <- (c^2 + b^2 - a^2)/(2 * b * c)
				if((a == 0 & b==0 & c==0) | (b==0) | (c==0)){
					# the bird didn't move
					step_counter[i] <- step_counter[i-1]
				}else{
					if(d < -1){
						d <- -1
					}
					if(d > 1){
						d <- 1
					}
					alpha <- acos(d)
					all_alphas <- c(all_alphas,alpha)
					#if(alpha > 3*pi/4){ # this is a new step
					if(alpha > pi/5){ # this is a new step
						step_counter[i] <- step_counter[i-1] + 1
					}else{
						step_counter[i] <- step_counter[i-1]
					}
				}
			}
		}
	}
}
plot(all_alphas)
# now computes step lengths
index_step_start <- 1
step_len <- c()
step_stage <- c()
for(i in 2:(length(step_counter)-1)){
	if(step_counter[i] != step_counter[i-1]){
		step_len <- c(step_len,distMeeus(c(as.numeric(coord_all[index_step_start,3]),as.numeric(coord_all[index_step_start,4])),c(as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4]))))
		step_stage <- c(step_stage,coord_all[i,10])
		index_step_start <- i
	}
}
step_len <- step_len/1000 # in km now
thresh_step_len <- 200 # km
step_len[step_len > thresh_step_len] <- NA
# histograms
#pdf("step_hist_radical.pdf",width=12,height=4)
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,2))
#hist(step_len[step_stage=="egg"],100,main="(a) Incubating",xlab="Step length (km)",col="red")
#hist(step_len[step_stage=="chick"],100,main="(b) Chick rearing",xlab="Step length (km)",col="blue")
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
#dev.off()
# log step plots
nblocks <- 50
step_len_egg <- hist(step_len[step_stage=="egg" & step_len<7], nblocks)
norm_freq_egg <- step_len_egg$counts / max(step_len_egg$counts)
cum.dist_egg <- 1 - cumsum(norm_freq_egg)/nblocks
step_len_chick <- hist(step_len[step_stage=="chick" & step_len<7], nblocks)
norm_freq_chick <- step_len_chick$counts / max(step_len_chick$counts)
cum.dist_chick <- 1 - cumsum(norm_freq_chick)/nblocks
pdf("step_small.pdf",width=18,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,3))
hist(step_len[step_stage=="egg"],100,main="(a) Incubating",xlab="Step length (km)",col="red")
plot(step_len_egg$mids,log(cum.dist_egg,10),col="red")
points(step_len_chick$mids,log(cum.dist_chick,10),col="blue")
x <- step_len_egg$mids
y <-log(cum.dist_egg,10)
a <- min(step_len[step_stage=="egg" & step_len<7],na.rm=T)
egg_ls <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(mu = 3),algorithm = "plinear")
summary(egg_ls)
pred_egg <- predict(egg_ls, seq(-2,1,.005))
lines(x, pred_egg,col="orange")
x <- step_len_chick$mids
y <- log(cum.dist_chick,10)
a <- min(step_len[step_stage=="chick" & step_len<7],na.rm=T)
chick_ls <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(mu = 3),algorithm = "plinear")
summary(chick_ls)
pred_chick <- predict(chick_ls, seq(0,5,.05))
lines(x, pred_chick,col="violet")
legend(3,-.02,col=c("orange","violet"),lty=1,c(paste("mu = ",format(10^(-summary(egg_ls)$coefficients[1,1]),digits=3),sep=""),paste("mu = ",format(10^(-summary(chick_ls)$coefficients[1,1]),digits=3),sep="")) )
hist(step_len[step_stage=="chick"],100,main="(c) Chick rearing",xlab="Step length (km)",col="blue")
dev.off()

# radical changes of direction
coord_all_clean <- coord_all
step_len <- rep(NA,length(coord_all[,1]))
step_counter <- rep(0,length(coord_all[,1]))
# first determines steps
for(i in 2:(length(coord_all_clean[,1])-1)){
	# remove data where time difference betwee 2 consecutive GPS logs > 8 min
	if((unix2POSIXct(as.numeric(coord_all_clean[i,7])) - unix2POSIXct(as.numeric(coord_all_clean[i-1,7]))) > 8){
		coord_all_clean[i,8] <- NA
	}
	# if change of bird_ new step
	if(coord_all[i,2] != coord_all[i-1,2]){
		step_counter[i] <- step_counter[i-1] + 1
	}else{
		# if any location is NA, new trip and hence new step
		if(is.na(as.numeric(coord_all[i-1,3])) | is.na(as.numeric(coord_all[i,3])) | is.na(as.numeric(coord_all[i+1,3]))){
			step_counter[i] <- step_counter[i-1] + 1
		}else{
			# if dist traveled == 0: same step
			if(as.numeric(coord_all[i,3])==0){
				step_counter[i] <- step_counter[i-1]
			}else{
				a <- distMeeus(c(as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4])),c(as.numeric(coord_all[i+1,3]),as.numeric(coord_all[i+1,4])))
				b <- distMeeus(c(as.numeric(coord_all[i-1,3]),as.numeric(coord_all[i-1,4])),c(as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4])))
				c <- distMeeus(c(as.numeric(coord_all[i-1,3]),as.numeric(coord_all[i-1,4])),c(as.numeric(coord_all[i+1,3]),as.numeric(coord_all[i+1,4])))
				d <- (b^2 + b^2 - a^2)/(2 * b * c)
				if(a == 0 & b==0 & c==0){
					# the bird didn't move
					step_counter[i] <- step_counter[i-1]
				}else{
					if(d < -1){
						d <- -1
					}
					if(d > 1){
						d <- 1
					}
					alpha <- acos(d)
					if(alpha > 7*pi/8){ # this is a new step
					#if(alpha > 3*pi/4){ # this is a new step
					#if(alpha > pi/5){ # this is a new step
						step_counter[i] <- step_counter[i-1] + 1
					}else{
						step_counter[i] <- step_counter[i-1]
					}
				}
			}
		}
	}
}
# now computes step lengths
index_step_start <- 1
step_len <- c()
step_stage <- c()
for(i in 2:(length(step_counter)-1)){
	if(step_counter[i] != step_counter[i-1]){
		step_len <- c(step_len,distMeeus(c(as.numeric(coord_all[index_step_start,3]),as.numeric(coord_all[index_step_start,4])),c(as.numeric(coord_all[i,3]),as.numeric(coord_all[i,4]))))
		step_stage <- c(step_stage,coord_all[i,10])
		index_step_start <- i
	}
}
step_len <- step_len/1000 # in km now
thresh_step_len <- 200 # km
step_len[step_len > thresh_step_len] <- NA
# histograms
#pdf("step_hist_radical.pdf",width=12,height=4)
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,2))
#hist(step_len[step_stage=="egg"],100,main="(a) Incubating",xlab="Step length (km)",col="red")
#hist(step_len[step_stage=="chick"],100,main="(b) Chick rearing",xlab="Step length (km)",col="blue")
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
#dev.off()
# log step plots
nblocks <- 50
step_len_egg <- hist(step_len[step_stage=="egg" & step_len<7], nblocks)
norm_freq_egg <- step_len_egg$counts / max(step_len_egg$counts)
cum.dist_egg <- 1 - cumsum(norm_freq_egg)/nblocks
step_len_chick <- hist(step_len[step_stage=="chick" & step_len<7], nblocks)
norm_freq_chick <- step_len_chick$counts / max(step_len_chick$counts)
cum.dist_chick <- 1 - cumsum(norm_freq_chick)/nblocks
pdf("step_radical.pdf",width=18,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,3))
hist(step_len[step_stage=="egg"],100,main="(a) Incubating",xlab="Step length (km)",col="red")
plot(step_len_egg$mids,log(cum.dist_egg,10),col="red")
points(step_len_chick$mids,log(cum.dist_chick,10),col="blue")
x <- step_len_egg$mids
y <-log(cum.dist_egg,10)
a <- min(step_len[step_stage=="egg" & step_len<7],na.rm=T)
egg_ls <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(mu = 3),algorithm = "plinear")
summary(egg_ls)
pred_egg <- predict(egg_ls, seq(-2,1,.005))
lines(x, pred_egg,col="orange")
x <- step_len_chick$mids
y <- log(cum.dist_chick,10)
a <- min(step_len[step_stage=="chick" & step_len<7],na.rm=T)
chick_ls <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(mu = 3),algorithm = "plinear")
summary(chick_ls)
pred_chick <- predict(chick_ls, seq(0,5,.05))
lines(x, pred_chick,col="violet")
legend(3,-.02,col=c("orange","violet"),lty=1,c(paste("mu = ",format(10^(-summary(egg_ls)$coefficients[1,1]),digits=3),sep=""),paste("mu = ",format(10^(-summary(chick_ls)$coefficients[1,1]),digits=3),sep="")) )
hist(step_len[step_stage=="chick"],100,main="(c) Chick rearing",xlab="Step length (km)",col="blue")
dev.off()


################################################################################################

# 141014
library('VGAM')
pdf("LevyPow_flight_raw_distr.pdf",width=6,height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,1))
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg"  & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
fit_egg <- vglm(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg"  & as.numeric(coord_all_clean[,8])>0]) ~ 1, paretoff)
mu_estimate_egg <- exp(fit_egg@coefficients)
step_lens_chick <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Chick-rearing birds",xlim=c(0,10000))
fit_chick <- vglm(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])>0]) ~ 1, paretoff)
mu_estimate_chick <- exp(fit_chick@coefficients)
dev.off()


# flight steps: "series of contiguous in-flight intervals; all single interval (i.e., *5* min) flight steps were ignored." (Humphries et al. PNAS 2012)
#pdf("flight_steps.pdf",width=12,height=8)
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2))
# egg
nblocks <- 200
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])<500 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
#step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
cum.dist_egg <- 1 - cumsum(norm_freq_egg)/nblocks
#plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(cum.dist_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="Incubating birds")
egg.ls <- loess(log(cum.dist_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
#points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")
# linear regression (power law on log-log scale)
egg.lm.ls <- glm(log(cum.dist_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),family = gaussian)
#abline(egg.lm.ls,col="red")
summary(egg.lm.ls)
logLik(egg.lm.ls)
AIC(egg.lm.ls)
# chick
nblocks <- 200
step_lens_chick <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])<500 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
#step_lens_chick <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
norm_freq_chick <- step_lens_chick$counts[step_lens_chick$counts>1] / max(step_lens_chick$counts[step_lens_chick$counts>1])
cum.dist_chick <- 1 - cumsum(norm_freq_chick)/nblocks
#plot(log(step_lens_chick$mids[step_lens_chick$counts>1],10), log(cum.dist_chick,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="Chick rearing birds")
chick.ls <- loess(log(cum.dist_chick,10) ~ log(step_lens_chick$mids[step_lens_chick$counts>1],10),span = 0.2)
pred_chick <- predict(chick.ls, seq(0,5,.005))
#points(seq(0,5,.005), pred_chick, pch=19, cex=.25,col="orange")
# linear regression (power law on log-log scale)
chick.lm.ls <- glm(log(cum.dist_chick,10) ~ log(step_lens_chick$mids[step_lens_chick$counts>1],10),family = gaussian)
#abline(chick.lm.ls,col="red")
summary(chick.lm.ls)
logLik(chick.lm.ls)
AIC(chick.lm.ls)
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
#dev.off()


pdf("LevyPow_flight_500.pdf",width=8,height=5.5)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
######################### egg #########################
x <- log(step_lens_egg$mids[step_lens_egg$counts>1],10)
y <- log(cum.dist_egg,10)
plot(y ~ x,type="l",,col="red",xlim=c(.5,3),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="")
a <- min(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),na.rm=T)
egg_ls <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(mu = 3),algorithm = "plinear")
summary(egg_ls)
pred_egg <- predict(egg_ls, seq(.5,3,.005))
lines(x, pred_egg,col="orange")
mu <- 1/(sum(norm_freq_egg)*13) # 1/(sum(norm_freq_egg)/length(norm_freq_egg))
exp.trace <- 10^(-mu*10^(seq(0,3,.005)))
#lines(log(10^(seq(0,3,.005)),10),log(exp.trace,10),lty=3,lwd=2,col="orange")
lp1.lm <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(a = 1, mu = -10))
summary(lp1.lm)
lp2.lm <- nls(y ~ lambda*exp(-lambda*(x-a)),start = list(a = .01, lambda = -1))
summary(lp2.lm)
#xv1 <- seq(0,3,.005); yv1 <- (predict(lp1.lm,list(x=xv))); lines(xv1, yv1,,lty=3,lwd=2,col="orange")
#xv2 <- seq(0,3,.005); yv2 <- (predict(lp2.lm,list(x=xv))); lines(xv2, yv2,,lty=3,lwd=2,col="orange")
#egg.lm.ls <- glm(log(cum.dist_egg,10)[log(step_lens_egg$mids[step_lens_egg$counts>1],10)<2] ~ 0+log(step_lens_egg$mids[step_lens_egg$counts>1],10)[log(step_lens_egg$mids[step_lens_egg$counts>1],10)<2],family = gaussian)
#abline(egg.lm.ls,col="yellow")
#summary(egg.lm.ls)
nb <- 13 + 10^x[1] #length(norm_freq_egg[norm_freq_egg > 0])
AIC.exp_egg <- -2 * (nb * log(mu) + nb*mu*10^x[1] - mu * sum(norm_freq_egg)) + 2
gamma <- 1 - nb/(nb*10^x[1] - sum(x))
AIC.pow_egg <- -2 * (nb*log(gamma-1)+nb*(gamma-1)*x[1]-gamma*sum(norm_freq_egg)) + 2
######################### chick #########################
x <- log(step_lens_chick$mids[step_lens_chick$counts>1],10)
y <- log(cum.dist_chick,10)
lines(y ~ x,type="l",,col="blue",xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]")
a <- min(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),na.rm=T)
chick_ls <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(mu = 3),algorithm = "plinear")
summary(chick_ls)
pred_chick <- predict(chick_ls, seq(.5,3,.005))
lines(x, pred_chick,col="purple")
mu <- 1/(sum(norm_freq_chick)*50) # 1/(sum(norm_freq_chick)/length(norm_freq_chick))
exp.trace <- 10^(-mu*10^(seq(0,3,.005)))
#lines(log(10^(seq(0,3,.005)),10),log(exp.trace,10),lty=3,lwd=2,col="purple")
lp1.lm <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(a = 1, mu = -1))
summary(lp1.lm)
lp2.lm <- nls(y ~ lambda*exp(-lambda*(x-a)),start = list(a = .01, lambda = -1))
summary(lp2.lm)
xv1 <- seq(0,3,.005); yv1 <- (predict(lp1.lm,list(x=xv))); lines(xv1, yv1,,lty=3,lwd=2,col="purple")
#xv2 <- seq(0,3,.005); yv2 <- (predict(lp2.lm,list(x=xv))); lines(xv2, yv2,,lty=3,lwd=2,col="purple")
#chick.lm.ls <- glm(log(cum.dist_chick,10)[log(step_lens_chick$mids[step_lens_chick$counts>1],10)<2] ~ 0+log(step_lens_chick$mids[step_lens_chick$counts>1],10)[log(step_lens_chick$mids[step_lens_chick$counts>1],10)<2],family = gaussian)
#abline(chick.lm.ls,col="violet")
#summary(chick.lm.ls)
nb <- 50 + 10^x[1] #length(norm_freq_chick[norm_freq_chick > 0])
AIC.exp_chick <- -2 * (nb * log(mu) + nb*mu*10^x[1] - mu * sum(norm_freq_chick)) + 2
gamma <- 1 - nb/(nb*10^x[1] - sum(x))
AIC.pow_chick <- -2 * (nb*log(gamma-1)+nb*(gamma-1)*x[1]-gamma*sum(norm_freq_chick)) + 2
#legend(.5,-.04,c("Incubation: raw",paste("Incubation: exp (AIC=", format(AIC.exp_egg,digits=5),")",sep=""),paste("Incubation: pow (AIC=", format(AIC.pow_egg,digits=5),")",sep=""),
#legend(.5,-.1,c("Incubation: raw",paste("Incubation: exp (AIC=", format(AIC.exp_egg,digits=5),")",sep=""),paste("Incubation: pow (AIC=", format(AIC.pow_egg,digits=5),")",sep=""),
#                "Chick rearing: raw",paste("Chick rearing: exp (AIC=", format(AIC.exp_chick,digits=5),")",sep=""),paste("Chick rearing: pow (AIC=", format(AIC.pow_chick,digits=5),")",sep="")),col=c("red","orange","orange","blue","purple","purple"),lty=c(1,3,3,1,3,3))
legend(.5,-.1,c("Incubation: raw data",paste("Incubation: fit (mu=", format(-summary(egg_ls)$coefficients[1,1],digits=3),")",sep=""),
#legend(.5,-.04,c("Incubation: raw data",paste("Incubation: fit (mu=", format(-summary(egg_ls)$coefficients[1,1],digits=3),")",sep=""),
                "Chick rearing: raw data",paste("Chick rearing: fit (mu=", format(-summary(chick_ls)$coefficients[1,1],digits=3),")",sep="")),col=c("red","orange","blue","purple"),lty=c(1,3,1,3))
dev.off()

format(AIC.exp_egg,digits=5)
format(AIC.pow_egg,digits=5)
format(AIC.exp_chick,digits=5)
format(AIC.pow_chick,digits=5)

AIC.pow_egg - AIC.exp_egg
AIC.pow_chick - AIC.exp_chick

AIC.exp_chick
AIC.pow_chick


# flight steps: "series of contiguous in-flight intervals; all single interval (i.e., *5* min) flight steps were ignored." (Humphries et al. PNAS 2012)
#pdf("flight_steps.pdf",width=12,height=8)
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2))
# egg
nblocks <- 200
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg"  & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
#step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
cum.dist_egg <- 1 - cumsum(norm_freq_egg)/nblocks
#plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(cum.dist_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="Incubating birds")
egg.ls <- loess(log(cum.dist_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
#points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")
# linear regression (power law on log-log scale)
egg.lm.ls <- glm(log(cum.dist_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),family = gaussian)
#abline(egg.lm.ls,col="red")
summary(egg.lm.ls)
logLik(egg.lm.ls)
AIC(egg.lm.ls)
# chick
nblocks <- 200
step_lens_chick <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
#step_lens_chick <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
norm_freq_chick <- step_lens_chick$counts[step_lens_chick$counts>1] / max(step_lens_chick$counts[step_lens_chick$counts>1])
cum.dist_chick <- 1 - cumsum(norm_freq_chick)/nblocks
#plot(log(step_lens_chick$mids[step_lens_chick$counts>1],10), log(cum.dist_chick,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="Chick rearing birds")
chick.ls <- loess(log(cum.dist_chick,10) ~ log(step_lens_chick$mids[step_lens_chick$counts>1],10),span = 0.2)
pred_chick <- predict(chick.ls, seq(0,5,.005))
#points(seq(0,5,.005), pred_chick, pch=19, cex=.25,col="orange")
# linear regression (power law on log-log scale)
chick.lm.ls <- glm(log(cum.dist_chick,10) ~ log(step_lens_chick$mids[step_lens_chick$counts>1],10),family = gaussian)
#abline(chick.lm.ls,col="red")
summary(chick.lm.ls)
logLik(chick.lm.ls)
AIC(chick.lm.ls)
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
#dev.off()


pdf("LevyPow_flight_all.pdf",width=8,height=5.5)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
######################### egg #########################
x <- log(step_lens_egg$mids[step_lens_egg$counts>1],10)
y <- log(cum.dist_egg,10)
plot(y ~ x,type="l",,col="red",xlim=c(1.25,4),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="")
a <- min(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),na.rm=T)
egg_ls <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(mu = 3),algorithm = "plinear")
summary(egg_ls)
pred_egg <- predict(egg_ls, seq(.5,3,.005))
lines(x, pred_egg,col="orange")
mu <- 1/(sum(norm_freq_egg)*13) # 1/(sum(norm_freq_egg)/length(norm_freq_egg))
exp.trace <- 10^(-mu*10^(seq(0,3,.005)))
#lines(log(10^(seq(0,3,.005)),10),log(exp.trace,10),lty=3,lwd=2,col="orange")
lp1.lm <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(a = 1, mu = -10))
summary(lp1.lm)
lp2.lm <- nls(y ~ lambda*exp(-lambda*(x-a)),start = list(a = .01, lambda = -1))
summary(lp2.lm)
#xv1 <- seq(0,3,.005); yv1 <- (predict(lp1.lm,list(x=xv))); lines(xv1, yv1,,lty=3,lwd=2,col="orange")
#xv2 <- seq(0,3,.005); yv2 <- (predict(lp2.lm,list(x=xv))); lines(xv2, yv2,,lty=3,lwd=2,col="orange")
#egg.lm.ls <- glm(log(cum.dist_egg,10)[log(step_lens_egg$mids[step_lens_egg$counts>1],10)<2] ~ 0+log(step_lens_egg$mids[step_lens_egg$counts>1],10)[log(step_lens_egg$mids[step_lens_egg$counts>1],10)<2],family = gaussian)
#abline(egg.lm.ls,col="yellow")
#summary(egg.lm.ls)
nb <- 13 + 10^x[1] #length(norm_freq_egg[norm_freq_egg > 0])
AIC.exp_egg <- -2 * (nb * log(mu) + nb*mu*10^x[1] - mu * sum(norm_freq_egg)) + 2
gamma <- 1 - nb/(nb*10^x[1] - sum(x))
AIC.pow_egg <- -2 * (nb*log(gamma-1)+nb*(gamma-1)*x[1]-gamma*sum(norm_freq_egg)) + 2
######################### chick #########################
x <- log(step_lens_chick$mids[step_lens_chick$counts>1],10)
y <- log(cum.dist_chick,10)
lines(y ~ x,type="l",,col="blue",xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]")
a <- min(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),na.rm=T)
chick_ls <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(mu = 3),algorithm = "plinear")
summary(chick_ls)
pred_chick <- predict(chick_ls, seq(.5,3,.005))
lines(x, pred_chick,col="purple")
mu <- 1/(sum(norm_freq_chick)*50) # 1/(sum(norm_freq_chick)/length(norm_freq_chick))
exp.trace <- 10^(-mu*10^(seq(0,3,.005)))
#lines(log(10^(seq(0,3,.005)),10),log(exp.trace,10),lty=3,lwd=2,col="purple")
lp1.lm <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(a = 1, mu = -1))
summary(lp1.lm)
lp2.lm <- nls(y ~ lambda*exp(-lambda*(x-a)),start = list(a = .01, lambda = -1))
summary(lp2.lm)
xv1 <- seq(0,3,.005); yv1 <- (predict(lp1.lm,list(x=xv))); lines(xv1, yv1,,lty=3,lwd=2,col="purple")
#xv2 <- seq(0,3,.005); yv2 <- (predict(lp2.lm,list(x=xv))); lines(xv2, yv2,,lty=3,lwd=2,col="purple")
#chick.lm.ls <- glm(log(cum.dist_chick,10)[log(step_lens_chick$mids[step_lens_chick$counts>1],10)<2] ~ 0+log(step_lens_chick$mids[step_lens_chick$counts>1],10)[log(step_lens_chick$mids[step_lens_chick$counts>1],10)<2],family = gaussian)
#abline(chick.lm.ls,col="violet")
#summary(chick.lm.ls)
nb <- 50 + 10^x[1] #length(norm_freq_chick[norm_freq_chick > 0])
AIC.exp_chick <- -2 * (nb * log(mu) + nb*mu*10^x[1] - mu * sum(norm_freq_chick)) + 2
gamma <- 1 - nb/(nb*10^x[1] - sum(x))
AIC.pow_chick <- -2 * (nb*log(gamma-1)+nb*(gamma-1)*x[1]-gamma*sum(norm_freq_chick)) + 2
#legend(.5,-.04,c("Incubation: raw",paste("Incubation: exp (AIC=", format(AIC.exp_egg,digits=5),")",sep=""),paste("Incubation: pow (AIC=", format(AIC.pow_egg,digits=5),")",sep=""),
#legend(.5,-.1,c("Incubation: raw",paste("Incubation: exp (AIC=", format(AIC.exp_egg,digits=5),")",sep=""),paste("Incubation: pow (AIC=", format(AIC.pow_egg,digits=5),")",sep=""),
#                "Chick rearing: raw",paste("Chick rearing: exp (AIC=", format(AIC.exp_chick,digits=5),")",sep=""),paste("Chick rearing: pow (AIC=", format(AIC.pow_chick,digits=5),")",sep="")),col=c("red","orange","orange","blue","purple","purple"),lty=c(1,3,3,1,3,3))
legend(1.225,-.015,c("Incubation: raw data",paste("Incubation: fit (mu=", format(-summary(egg_ls)$coefficients[1,1],digits=3),")",sep=""),
#legend(.5,-.04,c("Incubation: raw data",paste("Incubation: fit (mu=", format(-summary(egg_ls)$coefficients[1,1],digits=3),")",sep=""),
                "Chick rearing: raw data",paste("Chick rearing: fit (mu=", format(-summary(chick_ls)$coefficients[1,1],digits=3),")",sep="")),col=c("red","orange","blue","purple"),lty=c(1,3,1,3))
dev.off()

format(AIC.exp_egg,digits=5)
format(AIC.pow_egg,digits=5)
format(AIC.exp_chick,digits=5)
format(AIC.pow_chick,digits=5)

AIC.pow_egg - AIC.exp_egg
AIC.pow_chick - AIC.exp_chick

AIC.exp_chick
AIC.pow_chick


# flight steps: "series of contiguous in-flight intervals; all single interval (i.e., *5* min) flight steps were ignored." (Humphries et al. PNAS 2012)
#pdf("flight_steps.pdf",width=12,height=8)
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2))
# egg
nblocks <- 200
#step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])<500 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
#step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])>500 & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
#step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])>250 & as.numeric(coord_all_clean[,8])<4000 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
cum.dist_egg <- 1 - cumsum(norm_freq_egg)/nblocks
#plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(cum.dist_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="Incubating birds")
egg.ls <- loess(log(cum.dist_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
#points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")
# linear regression (power law on log-log scale)
egg.lm.ls <- glm(log(cum.dist_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),family = gaussian)
#abline(egg.lm.ls,col="red")
summary(egg.lm.ls)
logLik(egg.lm.ls)
AIC(egg.lm.ls)
# chick
nblocks <- 200
#step_lens_chick <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])<500 & as.numeric(coord_all_clean[,8])<4000 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
step_lens_chick <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
norm_freq_chick <- step_lens_chick$counts[step_lens_chick$counts>1] / max(step_lens_chick$counts[step_lens_chick$counts>1])
cum.dist_chick <- 1 - cumsum(norm_freq_chick)/nblocks
#plot(log(step_lens_chick$mids[step_lens_chick$counts>1],10), log(cum.dist_chick,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="Chick rearing birds")
chick.ls <- loess(log(cum.dist_chick,10) ~ log(step_lens_chick$mids[step_lens_chick$counts>1],10),span = 0.2)
pred_chick <- predict(chick.ls, seq(0,5,.005))
#points(seq(0,5,.005), pred_chick, pch=19, cex=.25,col="orange")
# linear regression (power law on log-log scale)
chick.lm.ls <- glm(log(cum.dist_chick,10) ~ log(step_lens_chick$mids[step_lens_chick$counts>1],10),family = gaussian)
#abline(chick.lm.ls,col="red")
summary(chick.lm.ls)
logLik(chick.lm.ls)
AIC(chick.lm.ls)
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
#dev.off()


pdf("LevyPow_flight_2000.pdf",width=8,height=5.5)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
######################### egg #########################
x <- log(step_lens_egg$mids[step_lens_egg$counts>1],10)
y <- log(cum.dist_egg,10)
plot(y ~ x,type="l",,col="red",xlim=c(2.5,3.8),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="")
#a <- min(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),na.rm=T)
#a <- min(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])>500 & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),na.rm=T)
a <- min(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])>500 & as.numeric(coord_all_clean[,8])>0]),na.rm=T)
egg_ls <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(mu = 3),algorithm = "plinear")
summary(egg_ls)
pred_egg <- predict(egg_ls, seq(.5,3,.005))
lines(x, pred_egg,col="orange")
mu <- 1/(sum(norm_freq_egg)*13) # 1/(sum(norm_freq_egg)/length(norm_freq_egg))
exp.trace <- 10^(-mu*10^(seq(0,3,.005)))
#lines(log(10^(seq(0,3,.005)),10),log(exp.trace,10),lty=3,lwd=2,col="orange")
lp1.lm <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(a = 1, mu = -10))
summary(lp1.lm)
lp2.lm <- nls(y ~ lambda*exp(-lambda*(x-a)),start = list(a = .01, lambda = -1))
summary(lp2.lm)
#xv1 <- seq(0,3,.005); yv1 <- (predict(lp1.lm,list(x=xv))); lines(xv1, yv1,,lty=3,lwd=2,col="orange")
#xv2 <- seq(0,3,.005); yv2 <- (predict(lp2.lm,list(x=xv))); lines(xv2, yv2,,lty=3,lwd=2,col="orange")
#egg.lm.ls <- glm(log(cum.dist_egg,10)[log(step_lens_egg$mids[step_lens_egg$counts>1],10)<2] ~ 0+log(step_lens_egg$mids[step_lens_egg$counts>1],10)[log(step_lens_egg$mids[step_lens_egg$counts>1],10)<2],family = gaussian)
#abline(egg.lm.ls,col="yellow")
#summary(egg.lm.ls)
nb <- 13 + 10^x[1] #length(norm_freq_egg[norm_freq_egg > 0])
AIC.exp_egg <- -2 * (nb * log(mu) + nb*mu*10^x[1] - mu * sum(norm_freq_egg)) + 2
gamma <- 1 - nb/(nb*10^x[1] - sum(x))
AIC.pow_egg <- -2 * (nb*log(gamma-1)+nb*(gamma-1)*x[1]-gamma*sum(norm_freq_egg)) + 2
######################### chick #########################
x <- log(step_lens_chick$mids[step_lens_chick$counts>1],10)
y <- log(cum.dist_chick,10)
lines(y ~ x,type="l",,col="blue",xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]")
a <- min(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])>500 & as.numeric(coord_all_clean[,8])>0]),na.rm=T)
chick_ls <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(mu = 3),algorithm = "plinear")
summary(chick_ls)
pred_chick <- predict(chick_ls, seq(.5,3,.005))
lines(x, pred_chick,col="purple")
mu <- 1/(sum(norm_freq_chick)*50) # 1/(sum(norm_freq_chick)/length(norm_freq_chick))
exp.trace <- 10^(-mu*10^(seq(0,3,.005)))
#lines(log(10^(seq(0,3,.005)),10),log(exp.trace,10),lty=3,lwd=2,col="purple")
lp1.lm <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(a = 1, mu = -1))
summary(lp1.lm)
lp2.lm <- nls(y ~ lambda*exp(-lambda*(x-a)),start = list(a = .01, lambda = -1))
summary(lp2.lm)
xv1 <- seq(0,3,.005); yv1 <- (predict(lp1.lm,list(x=xv))); lines(xv1, yv1,,lty=3,lwd=2,col="purple")
#xv2 <- seq(0,3,.005); yv2 <- (predict(lp2.lm,list(x=xv))); lines(xv2, yv2,,lty=3,lwd=2,col="purple")
#chick.lm.ls <- glm(log(cum.dist_chick,10)[log(step_lens_chick$mids[step_lens_chick$counts>1],10)<2] ~ 0+log(step_lens_chick$mids[step_lens_chick$counts>1],10)[log(step_lens_chick$mids[step_lens_chick$counts>1],10)<2],family = gaussian)
#abline(chick.lm.ls,col="violet")
#summary(chick.lm.ls)
nb <- 50 + 10^x[1] #length(norm_freq_chick[norm_freq_chick > 0])
AIC.exp_chick <- -2 * (nb * log(mu) + nb*mu*10^x[1] - mu * sum(norm_freq_chick)) + 2
gamma <- 1 - nb/(nb*10^x[1] - sum(x))
AIC.pow_chick <- -2 * (nb*log(gamma-1)+nb*(gamma-1)*x[1]-gamma*sum(norm_freq_chick)) + 2
#legend(.5,-.04,c("Incubation: raw",paste("Incubation: exp (AIC=", format(AIC.exp_egg,digits=5),")",sep=""),paste("Incubation: pow (AIC=", format(AIC.pow_egg,digits=5),")",sep=""),
#legend(.5,-.1,c("Incubation: raw",paste("Incubation: exp (AIC=", format(AIC.exp_egg,digits=5),")",sep=""),paste("Incubation: pow (AIC=", format(AIC.pow_egg,digits=5),")",sep=""),
#                "Chick rearing: raw",paste("Chick rearing: exp (AIC=", format(AIC.exp_chick,digits=5),")",sep=""),paste("Chick rearing: pow (AIC=", format(AIC.pow_chick,digits=5),")",sep="")),col=c("red","orange","orange","blue","purple","purple"),lty=c(1,3,3,1,3,3))
legend(3.1,-.04,c("Incubation: raw data",paste("Incubation: fit (mu=", format(-summary(egg_ls)$coefficients[1,1],digits=3),")",sep=""),
#legend(.5,-.04,c("Incubation: raw data",paste("Incubation: fit (mu=", format(-summary(egg_ls)$coefficients[1,1],digits=3),")",sep=""),
                "Chick rearing: raw data",paste("Chick rearing: fit (mu=", format(-summary(chick_ls)$coefficients[1,1],digits=3),")",sep="")),col=c("red","orange","blue","purple"),lty=c(1,3,1,3))

format(-summary(egg_ls)$coefficients[1,1],digits=3)
format(-summary(chick_ls)$coefficients[1,1],digits=3)

dev.off()

format(AIC.exp_egg,digits=5)
format(AIC.pow_egg,digits=5)
format(AIC.exp_chick,digits=5)
format(AIC.pow_chick,digits=5)

AIC.exp_chick
AIC.pow_chick

# flight steps: "series of contiguous in-flight intervals; all single interval (i.e., *5* min) flight steps were ignored." (Humphries et al. PNAS 2012)
#pdf("flight_steps.pdf",width=12,height=8)
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2))
# egg
nblocks <- 200
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])>500 & as.numeric(coord_all_clean[,8])<7000 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
#step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
cum.dist_egg <- 1 - cumsum(norm_freq_egg)/nblocks
#plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(cum.dist_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="Incubating birds")
egg.ls <- loess(log(cum.dist_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
#points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")
# linear regression (power law on log-log scale)
egg.lm.ls <- glm(log(cum.dist_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),family = gaussian)
#abline(egg.lm.ls,col="red")
summary(egg.lm.ls)
logLik(egg.lm.ls)
AIC(egg.lm.ls)
# chick
nblocks <- 200
step_lens_chick <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])>500 & as.numeric(coord_all_clean[,8])<7000 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
#step_lens_chick <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])<2000 & as.numeric(coord_all_clean[,8])>0]),nblocks,xlab="Move-step length (m)",main="Incubating birds",xlim=c(0,10000))
norm_freq_chick <- step_lens_chick$counts[step_lens_chick$counts>1] / max(step_lens_chick$counts[step_lens_chick$counts>1])
cum.dist_chick <- 1 - cumsum(norm_freq_chick)/nblocks
#plot(log(step_lens_chick$mids[step_lens_chick$counts>1],10), log(cum.dist_chick,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="Chick rearing birds")
chick.ls <- loess(log(cum.dist_chick,10) ~ log(step_lens_chick$mids[step_lens_chick$counts>1],10),span = 0.2)
pred_chick <- predict(chick.ls, seq(0,5,.005))
#points(seq(0,5,.005), pred_chick, pch=19, cex=.25,col="orange")
# linear regression (power law on log-log scale)
chick.lm.ls <- glm(log(cum.dist_chick,10) ~ log(step_lens_chick$mids[step_lens_chick$counts>1],10),family = gaussian)
#abline(chick.lm.ls,col="red")
summary(chick.lm.ls)
logLik(chick.lm.ls)
AIC(chick.lm.ls)
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
#dev.off()


pdf("LevyPow_flight_7000-500.pdf",width=8,height=5.5)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
######################### egg #########################
x <- log(step_lens_egg$mids[step_lens_egg$counts>1],10)
y <- log(cum.dist_egg,10)
plot(y ~ x,type="l",,col="red",ylim=c(-.026,-.002),xlim=c(2.5,4),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="")
a <- min(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="egg" & as.numeric(coord_all_clean[,8])>500 & as.numeric(coord_all_clean[,8])<7000 & as.numeric(coord_all_clean[,8])>0]),na.rm=T)
egg_ls <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(mu = 3),algorithm = "plinear")
summary(egg_ls)
pred_egg <- predict(egg_ls, seq(.5,3,.005))
lines(x, pred_egg,col="orange")
mu <- 1/(sum(norm_freq_egg)*13) # 1/(sum(norm_freq_egg)/length(norm_freq_egg))
exp.trace <- 10^(-mu*10^(seq(0,3,.005)))
#lines(log(10^(seq(0,3,.005)),10),log(exp.trace,10),lty=3,lwd=2,col="orange")
lp1.lm <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(a = 1, mu = -10))
summary(lp1.lm)
lp2.lm <- nls(y ~ lambda*exp(-lambda*(x-a)),start = list(a = .01, lambda = -1))
summary(lp2.lm)
#xv1 <- seq(0,3,.005); yv1 <- (predict(lp1.lm,list(x=xv))); lines(xv1, yv1,,lty=3,lwd=2,col="orange")
#xv2 <- seq(0,3,.005); yv2 <- (predict(lp2.lm,list(x=xv))); lines(xv2, yv2,,lty=3,lwd=2,col="orange")
#egg.lm.ls <- glm(log(cum.dist_egg,10)[log(step_lens_egg$mids[step_lens_egg$counts>1],10)<2] ~ 0+log(step_lens_egg$mids[step_lens_egg$counts>1],10)[log(step_lens_egg$mids[step_lens_egg$counts>1],10)<2],family = gaussian)
#abline(egg.lm.ls,col="yellow")
#summary(egg.lm.ls)
nb <- 13 + 10^x[1] #length(norm_freq_egg[norm_freq_egg > 0])
AIC.exp_egg <- -2 * (nb * log(mu) + nb*mu*10^x[1] - mu * sum(norm_freq_egg)) + 2
gamma <- 1 - nb/(nb*10^x[1] - sum(x))
AIC.pow_egg <- -2 * (nb*log(gamma-1)+nb*(gamma-1)*x[1]-gamma*sum(norm_freq_egg)) + 2
######################### chick #########################
x <- log(step_lens_chick$mids[step_lens_chick$counts>1],10)
y <- log(cum.dist_chick,10)
lines(y ~ x,type="l",,col="blue",xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]")
a <- min(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick" & as.numeric(coord_all_clean[,8])>500 & as.numeric(coord_all_clean[,8])<7000 & as.numeric(coord_all_clean[,8])>0]),na.rm=T)
chick_ls <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(mu = 3),algorithm = "plinear")
summary(chick_ls)
pred_chick <- predict(chick_ls, seq(.5,3,.005))
lines(x, pred_chick,col="purple")
mu <- 1/(sum(norm_freq_chick)*50) # 1/(sum(norm_freq_chick)/length(norm_freq_chick))
exp.trace <- 10^(-mu*10^(seq(0,3,.005)))
#lines(log(10^(seq(0,3,.005)),10),log(exp.trace,10),lty=3,lwd=2,col="purple")
lp1.lm <- nls(y ~ ((mu-1)*a^(mu-1))*x^(-mu),start = list(a = 1, mu = -1))
summary(lp1.lm)
lp2.lm <- nls(y ~ lambda*exp(-lambda*(x-a)),start = list(a = .01, lambda = -1))
summary(lp2.lm)
xv1 <- seq(0,3,.005); yv1 <- (predict(lp1.lm,list(x=xv))); lines(xv1, yv1,,lty=3,lwd=2,col="purple")
#xv2 <- seq(0,3,.005); yv2 <- (predict(lp2.lm,list(x=xv))); lines(xv2, yv2,,lty=3,lwd=2,col="purple")
#chick.lm.ls <- glm(log(cum.dist_chick,10)[log(step_lens_chick$mids[step_lens_chick$counts>1],10)<2] ~ 0+log(step_lens_chick$mids[step_lens_chick$counts>1],10)[log(step_lens_chick$mids[step_lens_chick$counts>1],10)<2],family = gaussian)
#abline(chick.lm.ls,col="violet")
#summary(chick.lm.ls)
nb <- 50 + 10^x[1] #length(norm_freq_chick[norm_freq_chick > 0])
AIC.exp_chick <- -2 * (nb * log(mu) + nb*mu*10^x[1] - mu * sum(norm_freq_chick)) + 2
gamma <- 1 - nb/(nb*10^x[1] - sum(x))
AIC.pow_chick <- -2 * (nb*log(gamma-1)+nb*(gamma-1)*x[1]-gamma*sum(norm_freq_chick)) + 2
#legend(.5,-.04,c("Incubation: raw",paste("Incubation: exp (AIC=", format(AIC.exp_egg,digits=5),")",sep=""),paste("Incubation: pow (AIC=", format(AIC.pow_egg,digits=5),")",sep=""),
#legend(.5,-.1,c("Incubation: raw",paste("Incubation: exp (AIC=", format(AIC.exp_egg,digits=5),")",sep=""),paste("Incubation: pow (AIC=", format(AIC.pow_egg,digits=5),")",sep=""),
#                "Chick rearing: raw",paste("Chick rearing: exp (AIC=", format(AIC.exp_chick,digits=5),")",sep=""),paste("Chick rearing: pow (AIC=", format(AIC.pow_chick,digits=5),")",sep="")),col=c("red","orange","orange","blue","purple","purple"),lty=c(1,3,3,1,3,3))
legend(2.5,-.015,c("Incubation: raw data",paste("Incubation: fit (mu=", format(-summary(egg_ls)$coefficients[1,1],digits=3),")",sep=""),
#legend(.5,-.04,c("Incubation: raw data",paste("Incubation: fit (mu=", format(-summary(egg_ls)$coefficients[1,1],digits=3),")",sep=""),
                "Chick rearing: raw data",paste("Chick rearing: fit (mu=", format(-summary(chick_ls)$coefficients[1,1],digits=3),")",sep="")),col=c("red","orange","blue","purple"),lty=c(1,3,1,3))
dev.off()

format(AIC.exp_egg,digits=5)
format(AIC.pow_egg,digits=5)
format(AIC.exp_chick,digits=5)
format(AIC.pow_chick,digits=5)

AIC.pow_egg - AIC.exp_egg
AIC.pow_chick - AIC.exp_chick

AIC.exp_chick
AIC.pow_chick

#######################################################################
#######################################################################
#######################################################################



# dive steps
# egg
nblocks <- 200
step_lens_egg <- hist(res_all$maxdep[res_all$stage=="egg" & res_all$maxdep<15],60,xlab="Dive depth (m)",main="Incubating birds")
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
cum.dist_egg <- 1 - cumsum(norm_freq_egg)/nblocks
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(cum.dist_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="Incubating birds")
egg.ls <- loess(log(cum.dist_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")
# chick
nblocks <- 200
step_lens_chick <- hist(res_all$maxdep[res_all$stage=="chick" & res_all$maxdep<15],60,xlab="Dive depth (m)",main="Chick rearing birds")
norm_freq_chick <- step_lens_chick$counts[step_lens_chick$counts>1] / max(step_lens_chick$counts[step_lens_chick$counts>1])
cum.dist_chick <- 1 - cumsum(norm_freq_chick)/nblocks
plot(log(step_lens_chick$mids[step_lens_chick$counts>1],10), log(cum.dist_chick,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]",main="Chick rearing birds")
chick.ls <- loess(log(cum.dist_chick,10) ~ log(step_lens_chick$mids[step_lens_chick$counts>1],10),span = 0.2)
pred_chick <- predict(chick.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_chick, pch=19, cex=.25,col="orange")

pdf("LevyPow_maxdive.pdf",width=8,height=5.5)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# egg
x <- log(step_lens_egg$mids[step_lens_egg$counts>1],10)
y <- log(cum.dist_egg,10)
plot(y ~ x,type="l",,col="red",xlim=c(.5,1.2),xlab="log10[Max dive, x (m)]",ylab="log10[Normalized frequency > x]",main="")
mu <- 1/(sum(norm_freq_egg)*14) # 1/(sum(norm_freq_egg)/length(norm_freq_egg))
exp.trace <- 10^(-mu*10^(seq(0,3,.005)))
##lines(log(10^(seq(0,3,.005)),10),log(exp.trace,10),lty=3,lwd=2,col="orange")
egg.lm.ls <- glm(y ~ x,family = gaussian)
abline(egg.lm.ls,col="orange",lty=3,lwd=2)
#summary(egg.lm.ls)
nb <- 14 + 10^x[1] #length(norm_freq_egg[norm_freq_egg > 0])
AIC.exp_egg <- -2 * (nb * log(mu) + nb*mu*10^x[1] - mu * sum(norm_freq_egg)) + 2
gamma <- 1 - nb/(nb*10^x[1] - sum(x))
AIC.pow_egg <- -2 * (nb*log(gamma-1)+nb*(gamma-1)*x[1]-gamma*sum(norm_freq_egg)) + 2
# chick
x <- log(step_lens_chick$mids[step_lens_chick$counts>1],10)
y <- log(cum.dist_chick,10)
lines(y ~ x,type="l",,col="blue",xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency > x]")
mu <- 1/(sum(norm_freq_chick)*50) # 1/(sum(norm_freq_chick)/length(norm_freq_chick))
exp.trace <- 10^(-mu*10^(seq(0,3,.005)))
#lines(log(10^(seq(0,3,.005)),10),log(exp.trace,10),lty=3,lwd=2,col="purple")
chick.lm.ls <- glm(y ~ x,family = gaussian)
abline(chick.lm.ls,col="purple",lty=3,lwd=2)
#summary(chick.lm.ls)
nb <- 14 + 10^x[1] #length(norm_freq_chick[norm_freq_chick > 0])
AIC.exp_chick <- -2 * (nb * log(mu) + nb*mu*10^x[1] - mu * sum(norm_freq_chick)) + 2
gamma <- 1 - nb/(nb*10^x[1] - sum(x))
AIC.pow_chick <- -2 * (nb*log(gamma-1)+nb*(gamma-1)*x[1]-gamma*sum(norm_freq_chick)) + 2
legend(.9,0,c("Incubation: raw",paste("Incubation: exp (AIC=", format(AIC.exp_egg,digits=5),")",sep=""),paste("Incubation: pow (AIC=", format(AIC.pow_egg,digits=5),")",sep=""),
                "Chick rearing: raw",paste("Chick rearing: exp (AIC=", format(AIC.exp_chick,digits=5),")",sep=""),paste("Chick rearing: pow (AIC=", format(AIC.pow_chick,digits=5),")",sep="")),
                col=c("red","orange","orange","blue","purple","purple"),lty=c(1,3,3,1,3,3))
dev.off()








# theta -- edited 131112
# flight steps: "series of contiguous in-flight intervals; all single interval (i.e., *5* min) flight steps were ignored." (Humphries et al. PNAS 2012)
pdf("theta_steps.pdf",width=12,height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2))
# egg
step_lens_egg <- hist(as.numeric(ipq_time$theta),200,xlab="Move-step length (m)",main="Theta-all",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="Incubating birds")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")
# chick
step_lens_chick <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,10]=="chick"]),200,xlab="Move-step length (m)",main="Chick rearing birds",xlim=c(0,10000))
norm_freq_chick <- step_lens_chick$counts[step_lens_chick$counts>1] / max(step_lens_chick$counts[step_lens_chick$counts>1])
plot(log(step_lens_chick$mids[step_lens_chick$counts>1],10), log(norm_freq_chick,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="Chick rearing birds")
chick.ls <- loess(log(norm_freq_chick,10) ~ log(step_lens_chick$mids[step_lens_chick$counts>1],10),span = 0.2)
pred_chick <- predict(chick.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_chick, pch=19, cex=.25,col="orange")
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()


atan2(45,0) * 180/pi

# wrt coordinate (0,0)
angle <- atan2(as.numeric(coord_all[,4]), as.numeric(coord_all[,3])) * 180/pi
hist(angle,100)
pdf("angle_stage_GPS_00.pdf",width=6,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
egg.angle <- angle[coord_all[,10]== "egg"]
chick.angle <- angle[coord_all[,10]== "chick"]
plot(density(as.numeric(chick.angle),na.rm=T,bw=.035), xlab = expression(paste("Inclination ",theta," during flights",sep="")), col = "blue", main="")
lines(density(as.numeric(egg.angle),na.rm=T,bw=.035), col="red")
legend(95.55,3,c("Incubation","Chick rearing"), col=c("red","blue"), lty=1, cex=.7)
dev.off()

angle <- atan2(as.numeric(res_all[,7]), as.numeric(res_all[,6])) * 180/pi #trackAngle(as.matrix(x))
hist(angle,100)
pdf("angle_stage_TDR_00.pdf",width=6,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
egg.angle <- angle[res_all$stage== "egg"]
chick.angle <- angle[res_all$stage== "chick"]
plot(density(as.numeric(chick.angle),na.rm=T,bw=.035), xlab = expression(paste("Inclination ",theta," during dives",sep="")), col = "blue", main="")
lines(density(as.numeric(egg.angle),na.rm=T,bw=.035), col="red")
legend(95.55,3.5,c("Incubation","Chick rearing"), col=c("red","blue"), lty=1, cex=.7)
dev.off()

# wrt colony: col.loc1
col.loc <- rep(col.loc1, length(coord_all[,4]))
col.loc[as.numeric(coord_all[,5])==2012] <- col.loc2

angle <- atan2(as.numeric(res_all[,7])-col.loc[2], as.numeric(res_all[,6])-col.loc[1]) * 180/pi #trackAngle(as.matrix(x))
#hist(angle,100)
pdf("angle_stage_TDR-GPS_col.pdf",width=6,height=8)
#pdf("angle_stage_TDR_col.pdf",width=6,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,1))
egg.angle <- angle[res_all$stage== "egg"]
chick.angle <- angle[res_all$stage== "chick"]
plot(density(as.numeric(chick.angle),na.rm=T,bw=6.2), xlab = expression(paste("Inclination ",theta," during dives wrt colony location",sep="")), col = "blue", main="")
lines(density(as.numeric(egg.angle),na.rm=T,bw=6.2), col="red")
legend(-100,.0225,c("Incubation","Chick rearing"), col=c("red","blue"), lty=1, cex=.99)
abline(v=-180,col="gray",lty=2)
abline(v=180,col="gray",lty=2)
#dev.off()
ks.test(chick.angle,egg.angle)

angle <- atan2(as.numeric(coord_all[,4])-col.loc[2], as.numeric(coord_all[,3])-col.loc[1]) * 180/pi
#hist(angle,100)
#pdf("angle_stage_GPS_col.pdf",width=6,height=4)
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
egg.angle <- angle[coord_all[,10]== "egg"]
chick.angle <- angle[coord_all[,10]== "chick"]
plot(density(as.numeric(chick.angle),na.rm=T,bw=6.2), xlab = expression(paste("Inclination ",theta," during flights wrt colony location",sep="")), col = "blue", main="")
lines(density(as.numeric(egg.angle),na.rm=T,bw=6.2), col="red")
abline(v=-180,col="gray",lty=2)
abline(v=180,col="gray",lty=2)
legend(-100,.0225,c("Incubation","Chick rearing"), col=c("red","blue"), lty=1, cex=.99)
dev.off()
ks.test(chick.angle,egg.angle)





pdf("angle_individual.pdf")
angle <- angle
id <- coord_all[,11]
boxplot(as.numeric(angle)~id, xlab = "Angle by individuals")
dev.off()

pdf("angle_stage.pdf", width=15.5, height=11)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,1))
egg.angle <- angle[coord_all[,10]== "egg"]
egg.id <- coord_all[,11][coord_all[,10]== "egg"]
chick.angle <- angle[coord_all[,10]== "chick"]
chick.id <- coord_all[,11][coord_all[,10]== "chick"]
boxplot(as.numeric(egg.angle)~egg.id, xlab = "Angles during incubation")
boxplot(as.numeric(chick.angle)~chick.id, xlab = "Angle during chick-rearing")
dev.off()











# by individual
pdf("flight_steps_individual.pdf",width=4,height=6)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(4,2))
# k24793
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="k24793"]),200,xlab="Move-step length (m)",main="k24793",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="k24793_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

# k93793
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="k93793"]),200,xlab="Move-step length (m)",main="k93793",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="k93793_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

#m27931
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="m27931"]),200,xlab="Move-step length (m)",main="m27931",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="m27931_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

#m93896
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="m93896"]),200,xlab="Move-step length (m)",main="m93896",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="m93896_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

#m93924
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="m93924"]),200,xlab="Move-step length (m)",main="m93924",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="m93924_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

#m93998
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="m93998"]),200,xlab="Move-step length (m)",main="m93998",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="m93998_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

#m93925
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="m93925"]),200,xlab="Move-step length (m)",main="m93925",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="m93925_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

#k24827
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="k24827"]),200,xlab="Move-step length (m)",main="k24827",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="k24827_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

#k24828
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="k24828"]),200,xlab="Move-step length (m)",main="k24828",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="k24828_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

#k24830
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="k24830"]),200,xlab="Move-step length (m)",main="k24830",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="k24830_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

#k24835
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="k24835"]),200,xlab="Move-step length (m)",main="k24835",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="k24835_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

#m27788
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="m27788"]),200,xlab="Move-step length (m)",main="m27788",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="m27788_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

#m93992
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="m93992"]),200,xlab="Move-step length (m)",main="m93992",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="m93992_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

#m93994
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="m93994"]),200,xlab="Move-step length (m)",main="m93994",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="m93994_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")

#m93995
step_lens_egg <- hist(as.numeric(coord_all_clean[,8][coord_all_clean[,11]=="m93995"]),200,xlab="Move-step length (m)",main="m93995",xlim=c(0,10000))
norm_freq_egg <- step_lens_egg$counts[step_lens_egg$counts>1] / max(step_lens_egg$counts[step_lens_egg$counts>1])
plot(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]",main="m93995_flight")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(step_lens_egg$mids[step_lens_egg$counts>1],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.005))
points(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="orange")
dev.off()

pdf("flight_steps_overlayed.pdf",width=6,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# egg
plot(seq(0,5,.005), pred_egg, pch=19, cex=.25,col="blue",xlim=c(min(log(step_lens_egg$mids[step_lens_egg$counts>1],10),log(step_lens_chick$mids[step_lens_chick$counts>1],10)),max(log(step_lens_egg$mids[step_lens_egg$counts>1],10),log(step_lens_chick$mids[step_lens_chick$counts>1],10))),ylim=c(min(log(norm_freq_egg,10),log(norm_freq_chick,10)),max(log(norm_freq_egg,10),log(norm_freq_chick,10))),xlab="log10[Move step, x (m)]",ylab="log10[Normalized frequency]")
points(log(step_lens_egg$mids[step_lens_egg$counts>1],10), log(norm_freq_egg,10),col="blue")
# chick
points(seq(0,5,.005), pred_chick, pch=19, cex=.25,col="red")
points(log(step_lens_chick$mids[step_lens_chick$counts>1],10), log(norm_freq_chick,10),col="red")
legend(1.8,-1.5,c("Incubating","Chick rearing"),lty=1,col=c("blue","red"),lwd=2)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()

# dive steps
pdf("dive_depths.pdf",width=12,height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2))
# egg
dive_lens_egg <- hist(res_all$maxdep[res_all$stage=="egg"],60,xlab="Dive depth (m)",main="Incubating birds")
norm_freq_egg <- dive_lens_egg$counts[dive_lens_egg$counts>2.5] / max(dive_lens_egg$counts[dive_lens_egg$counts>2.5])
plot(log(dive_lens_egg$mids[dive_lens_egg$counts>2.5],10), log(norm_freq_egg,10),xlab="log10[Dive depth, x (m)]",ylab="log10[Frequency]",main="Incubating birds")
egg.ls <- loess(log(norm_freq_egg,10) ~ log(dive_lens_egg$mids[dive_lens_egg$counts>2.5],10),span = 0.2)
pred_egg <- predict(egg.ls, seq(0,5,.0005))
points(seq(0,5,.0005), pred_egg, pch=19, cex=.25,col="orange")
# chick
dive_lens_chick <- hist(res_all$maxdep[res_all$stage=="chick"],60,xlab="Dive depth (m)",main="Chick rearing birds")
norm_freq_chick <- dive_lens_chick$counts[dive_lens_chick$counts>2.5] / max(dive_lens_chick$counts[dive_lens_chick$counts>2.5])
plot(log(dive_lens_chick$mids[dive_lens_chick$counts>2.5],10), log(norm_freq_chick,10),xlab="log10[Dive depth, x (m)]",ylab="log10[Frequency]",main="Chick rearing birds")
chick.ls <- loess(log(norm_freq_chick,10) ~ log(dive_lens_chick$mids[dive_lens_chick$counts>2.5],10),span = 0.2)
pred_chick <- predict(chick.ls, seq(0,5,.0005))
points(seq(0,5,.0005), pred_chick, pch=19, cex=.25,col="orange")
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()


pdf("dive_dives_overlayed.pdf",width=6,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# egg
plot(seq(0,5,.0005), pred_egg, pch=19, cex=.25,col="blue",xlim=c(min(log(dive_lens_egg$mids[dive_lens_egg$counts>2.5],10),log(dive_lens_chick$mids[dive_lens_chick$counts>2.5],10)),max(log(dive_lens_egg$mids[dive_lens_egg$counts>2.5],10),log(dive_lens_chick$mids[dive_lens_chick$counts>2.5],10))),ylim=c(min(log(norm_freq_egg,10),log(norm_freq_chick,10)),max(log(norm_freq_egg,10),log(norm_freq_chick,10))),xlab="log10[Dive depth, x (m)]",ylab="log10[Normalized frequency]")
points(log(dive_lens_egg$mids[dive_lens_egg$counts>2.5],10), log(norm_freq_egg,10),col="blue")
# chick
points(seq(0,5,.0005), pred_chick, pch=19, cex=.25,col="red")
points(log(dive_lens_chick$mids[dive_lens_chick$counts>2.5],10), log(norm_freq_chick,10),col="red")
legend(.6,-1.5,c("Incubating","Chick rearing"),lty=1,col=c("blue","red"),lwd=2)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()



####################################################################
save.image(paste("tdr_new", mytdr_d.thresh,".", mytdr_t.thresh, ".RData", sep=""))
####################################################################




##################################
# IPQ as a function of trip time #
##################################
# 2011
ipq_time_2011 <- data.frame(ipq_time[,1][ipq_time$year == 2011],ipq_time[,3][ipq_time$year == 2011],ipq_time[,5][ipq_time$year == 2011],ipq_time[,6][ipq_time$year == 2011])
colnames(ipq_time_2011) <- c("BirdID","tripID","tripdur","ipq")
# line plots
p_ipq_time_2011_trips <- ggplot(ipq_time_2011, aes(x = tripdur, y = ipq, linetype=factor(BirdID), col = factor(tripID))) + 
geom_line() +
scale_x_continuous(name="Trip duration") + 
scale_y_continuous(name="log IPQ") + 
ggtitle("log IPQ vs. trip duration in 2011")
# regressions
p_ipq_time_2011_reg <- ggplot(ipq_time_2011, aes(x = tripdur, y = ipq, linetype=factor(BirdID), col = factor(tripID))) + 
geom_smooth(method="rlm",se=F) +
#theme_bw() +
scale_x_continuous(name="Trip duration") + 
scale_y_continuous(name="log IPQ") + 
ggtitle("log IPQ vs. trip duration in 2011")
# 2012
ipq_time_2012 <- data.frame(ipq_time[,1][ipq_time$year == 2012],ipq_time[,3][ipq_time$year == 2012],ipq_time[,5][ipq_time$year == 2012],ipq_time[,6][ipq_time$year == 2012])
colnames(ipq_time_2012) <- c("BirdID","tripID","tripdur","ipq")
# line plots
p_ipq_time_2012_trips <- ggplot(ipq_time_2012, aes(x = tripdur, y = ipq, linetype=factor(BirdID), col = factor(tripID))) + 
geom_line() +
scale_x_continuous(name="Trip duration") + 
scale_y_continuous(name="log IPQ") + 
ggtitle("log IPQ vs. trip duration in 2012")
# regressions
p_ipq_time_2012_reg <- ggplot(ipq_time_2012, aes(x = tripdur, y = ipq, linetype=factor(BirdID), col = factor(tripID))) + 
geom_smooth(method="rlm",se=F) +
#theme_bw() +
scale_x_continuous(name="Trip duration") + 
scale_y_continuous(name="log IPQ") + 
ggtitle("log IPQ vs. trip duration in 2012")
# plotting 2011 and 2012 together
p <- arrangeGrob(p_ipq_time_2011_trips, p_ipq_time_2011_reg, p_ipq_time_2012_trips, p_ipq_time_2012_reg, ncol=2)
ggsave("ipq_time.pdf", plot=p, width=16, height=8)

# regression slopes per trip per year
slopes_2011 <- rep(NA,max(ipq_time_2011$tripID))
for(i in 1:max(ipq_time_2011$tripID)){
	if(sum(ipq_time_2011$tripID == i)>2){
		slopes_2011[i] <- summary(lmRob(ipq_time_2011$ipq[ipq_time_2011$tripID == i] ~ ipq_time_2011$tripdur[ipq_time_2011$tripID == i]))$coefficients[2,1]
	}
	
}
mean(slopes_2011,na.rm=T)
sum(slopes_2011>0,na.rm=T)/sum(!is.na(slopes_2011))
slopes_2012 <- rep(NA,max(ipq_time_2012$tripID))
for(i in 1:max(ipq_time_2012$tripID)){
	if(sum(ipq_time_2012$tripID == i)>2){
		slopes_2012[i] <- summary(lm(ipq_time_2012$ipq[ipq_time_2012$tripID == i] ~ ipq_time_2012$tripdur[ipq_time_2012$tripID == i]))$coefficients[2,1]
	}
	
}
mean(slopes_2012,na.rm=T)
sum(slopes_2012>0,na.rm=T)/sum(!is.na(slopes_2012))

# regression slopes per bird per year
slopes_2011 <- rep(NA,max(as.numeric(ipq_time_2011$BirdID)))
for(i in 1:max(as.numeric(ipq_time_2011$BirdID))){
	if(sum(as.numeric(ipq_time_2011$BirdID) == i)>2){
		slopes_2011[i] <- summary(lm(ipq_time_2011$ipq[as.numeric(ipq_time_2011$BirdID) == i] ~ ipq_time_2011$tripdur[as.numeric(ipq_time_2011$BirdID) == i]))$coefficients[2,1]
	}
	
}
mean(slopes_2011,na.rm=T)
sum(slopes_2011>0,na.rm=T)/sum(!is.na(slopes_2011))
slopes_2012 <- rep(NA,max(as.numeric(ipq_time_2012$BirdID)))
for(i in 1:max(as.numeric(ipq_time_2012$BirdID))){
	if(sum(as.numeric(ipq_time_2012$BirdID) == i)>2){
		slopes_2012[i] <- summary(lm(ipq_time_2012$ipq[as.numeric(ipq_time_2012$BirdID) == i] ~ ipq_time_2012$tripdur[as.numeric(ipq_time_2012$BirdID) == i]))$coefficients[2,1]
	}
	
}
mean(slopes_2012,na.rm=T)
sum(slopes_2012>0,na.rm=T)/sum(!is.na(slopes_2012))


#################
# IPQ and speed #
#################

pdf("ipq.speed.pdf")
plot(in.out.times[,19], in.out.times[,11],pch=16, xlab="Mean speed", ylab="Mean IPQ", main="From in.out.times")
plot(coord_all[,6], coord_all[,9],pch=16, xlab="Mean speed", ylab="Mean IPQ", main="From coord_all")
summary(lmRob(as.numeric(coord_all[,9])~as.numeric(coord_all[,6])))
dev.off()


################
# time budgets #
################

pdf("time.budget_linear.pdf", width=12, height=8)
mycols <- in.out.times[,4]
mycols[mycols=="2011"] <- "red"
mycols[mycols=="2012"] <- "blue"
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2))
plot(in.out.times[,11]~ in.out.times[,13],pch=16, col=mycols, xlab="Time total (s)", ylab="IPQ")
lm_tot1 <- lmRob(as.numeric(in.out.times[,11][in.out.times[,4]==2011])~ as.numeric(in.out.times[,13][in.out.times[,4] == 2011]) )
lm_tot2 <- lmRob(as.numeric(in.out.times[,11][in.out.times[,4]==2012])~ as.numeric(in.out.times[,13][in.out.times[,4] == 2012]) )
abline(lm_tot1, col=mycols[in.out.times[,4]==2011])
abline(lm_tot2, col=mycols[in.out.times[,4]==2012])
legend("topright",col=c("red","blue"),c("2011","2012"),pch=16)
text(40000,-.18,paste("2011: P=",format(summary(lm_tot1)$coefficients[2,4], scientific = F, digits=4),"; 2012: P=",format(summary(lm_tot2)$coefficients[2,4], scientific = F, digits=4),sep=""))
plot(in.out.times[,11]~ in.out.times[,14],pch=16, col=mycols, xlab="Time diving (s)", ylab="IPQ")
lm_div1 <- lmRob(as.numeric(in.out.times[,11][in.out.times[,4]==2011])~ as.numeric(in.out.times[,14][in.out.times[,4]==2011]) )
lm_div2 <- lmRob(as.numeric(in.out.times[,11][in.out.times[,4]==2012])~ as.numeric(in.out.times[,14][in.out.times[,4]==2012]) )
abline(lm_div1, col=mycols[in.out.times[,4]==2011])
abline(lm_div2, col=mycols[in.out.times[,4]==2012])
legend("topright",col=c("red","blue"),c("2011","2012"),pch=16)
text(5000,-.18,paste("2011: P=",format(summary(lm_div1)$coefficients[2,4], scientific = F, digits=4),"; 2012: P=",format(summary(lm_div2)$coefficients[2,4], scientific = F, digits=4),sep=""))
plot(in.out.times[,11]~ in.out.times[,15],pch=16, col=mycols, xlab="Time flying (s)", ylab="IPQ")
lm_fly1 <- lmRob(as.numeric(in.out.times[,11][in.out.times[,4]==2011])~ as.numeric(in.out.times[,15][in.out.times[,4]==2011]) )
lm_fly2 <- lmRob(as.numeric(in.out.times[,11][in.out.times[,4]==2012])~ as.numeric(in.out.times[,15][in.out.times[,4]==2012]) )
abline(lm_fly1, col=mycols[in.out.times[,4]==2011])
abline(lm_fly2, col=mycols[in.out.times[,4]==2012])
legend("topright",col=c("red","blue"),c("2011","2012"),pch=16)
text(2000,-.18,paste("2011: P=",format(summary(lm_fly1)$coefficients[2,4], scientific = F, digits=4),"; 2012: P=",format(summary(lm_fly2)$coefficients[2,4], scientific = F, digits=4),sep=""))
plot(in.out.times[,11]~ in.out.times[,16],pch=16, col=mycols, xlab="Time floating (s)", ylab="IPQ")
lm_flo1 <- lmRob(as.numeric(in.out.times[,11][in.out.times[,4]==2011])~ as.numeric(in.out.times[,16][in.out.times[,4]==2011]) )
lm_flo2 <- lmRob(as.numeric(in.out.times[,11][in.out.times[,4]==2012])~ as.numeric(in.out.times[,16][in.out.times[,4]==2012]) )
abline(lm_flo1, col=mycols[in.out.times[,4]==2011])
abline(lm_flo2, col=mycols[in.out.times[,4]==2012])
legend("topright",col=c("red","blue"),c("2011","2012"),pch=16)
text(20000,-.18,paste("2011: P=",format(summary(lm_flo1)$coefficients[2,4], scientific = F, digits=4),"; 2012: P=",format(summary(lm_flo1)$coefficients[2,4], scientific = F, digits=4),sep=""))
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()



pdf("time.budget_quadratic.pdf", width=12, height=8)
mycols <- in.out.times[,4]
mycols[mycols=="2011"] <- "red"
mycols[mycols=="2012"] <- "blue"
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2))
# total time
plot(in.out.times[,11]~ in.out.times[,13],pch=16, col=mycols, xlab="Time total (s)", ylab="IPQ")
x <- as.numeric(in.out.times[,13][in.out.times[,4] == 2011])
y <- as.numeric(in.out.times[,11][in.out.times[,4]==2011])
lm2_tot1 <- lm(y ~ poly(x,2)); summary(lm2_tot1)
xv <- seq(min(as.numeric(in.out.times[,13])),max(as.numeric(in.out.times[,13])),10)
yv <- predict(lm2_tot1,list(x=xv))
points(xv, yv,pch=".",col=mycols[in.out.times[,4]==2011])
x <- as.numeric(in.out.times[,13][in.out.times[,4] == 2012])
y <- as.numeric(in.out.times[,11][in.out.times[,4]==2012])
lm2_tot2 <- lm(y ~ poly(x,2)); summary(lm2_tot2)
xv <- seq(min(as.numeric(in.out.times[,13])),max(as.numeric(in.out.times[,13])),10)
yv <- predict(lm2_tot2,list(x=xv))
points(xv, yv,pch=".",col=mycols[in.out.times[,4]==2012])
legend("topright",col=c("red","blue"),c("2011","2012"),pch=16)
# time diving
plot(in.out.times[,11]~ in.out.times[,14],pch=16, col=mycols, xlab="Time diving (s)", ylab="IPQ")
x <- as.numeric(in.out.times[,14][in.out.times[,4]==2011])
y <- as.numeric(in.out.times[,11][in.out.times[,4]==2011])
lm2_tot1 <- lm(y ~ poly(x,2)); summary(lm2_tot1)
xv <- seq(min(as.numeric(in.out.times[,14])),max(as.numeric(in.out.times[,14])),10)
yv <- predict(lm2_tot1,list(x=xv))
points(xv, yv,pch=".",col=mycols[in.out.times[,4]==2011])
x <- as.numeric(in.out.times[,14][in.out.times[,4]==2012])
y <- as.numeric(in.out.times[,11][in.out.times[,4]==2012])
lm2_tot2 <- lm(y ~ poly(x,2)); summary(lm2_tot2)
xv <- seq(min(as.numeric(in.out.times[,14])),max(as.numeric(in.out.times[,14])),10)
yv <- predict(lm2_tot2,list(x=xv))
points(xv, yv,pch=".",col=mycols[in.out.times[,4]==2012])
legend("topright",col=c("red","blue"),c("2011","2012"),pch=16)
# time flying
plot(in.out.times[,11]~ in.out.times[,15],pch=16, col=mycols, xlab="Time flying (s)", ylab="IPQ")
x <- as.numeric(in.out.times[,15][in.out.times[,4]==2011])
y <- as.numeric(in.out.times[,11][in.out.times[,4]==2011])
lm2_tot1 <- lm(y ~ poly(x,2)); summary(lm2_tot1)
xv <- seq(min(as.numeric(in.out.times[,15])),max(as.numeric(in.out.times[,15])),10)
yv <- predict(lm2_tot1,list(x=xv))
points(xv, yv,pch=".",col=mycols[in.out.times[,4]==2011])
x <- as.numeric(in.out.times[,15][in.out.times[,4]==2012])
y <- as.numeric(in.out.times[,11][in.out.times[,4]==2012])
lm2_tot2 <- lm(y ~ poly(x,2)); summary(lm2_tot2)
xv <- seq(min(as.numeric(in.out.times[,15])),max(as.numeric(in.out.times[,15])),10)
yv <- predict(lm2_tot2,list(x=xv))
points(xv, yv,pch=".",col=mycols[in.out.times[,4]==2012])
legend("topright",col=c("red","blue"),c("2011","2012"),pch=16)
# time floating
plot(in.out.times[,11]~ in.out.times[,16],pch=16, col=mycols, xlab="Time floating (s)", ylab="IPQ")
x <- as.numeric(in.out.times[,16][in.out.times[,4]==2011])
y <- as.numeric(in.out.times[,11][in.out.times[,4]==2011])
lm2_tot1 <- lm(y ~ poly(x,2)); summary(lm2_tot1)
xv <- seq(min(as.numeric(in.out.times[,16])),max(as.numeric(in.out.times[,16])),10)
yv <- predict(lm2_tot1,list(x=xv))
points(xv, yv,pch=".",col=mycols[in.out.times[,4]==2011])
x <- as.numeric(in.out.times[,16][in.out.times[,4]==2012])
y <- as.numeric(in.out.times[,11][in.out.times[,4]==2012])
lm2_tot2 <- lm(y ~ poly(x,2)); summary(lm2_tot2)
xv <- seq(min(as.numeric(in.out.times[,16])),max(as.numeric(in.out.times[,16])),10)
yv <- predict(lm2_tot2,list(x=xv))
points(xv, yv,pch=".",col=mycols[in.out.times[,4]==2012])
legend("topright",col=c("red","blue"),c("2011","2012"),pch=16)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()





pdf("maxdist2col.hist.pdf", width=6, height=4)
hist(as.numeric(in.out.times[,7]),20, xlab=paste("maximum distance to colony (mode at",format(half.range.mode(as.numeric(in.out.times[,7]),B=100),digits=7),"m)"), main="")
dev.off()


plot(coord_all[,8], col=coord_all[,2],type="l", ylab="travel distance from GPS files")

goingouttime <- unix2POSIXct(in.out.times[,5])
comingbacktime <- unix2POSIXct(in.out.times[,6])

# just to check
goingouttime[1] == unix2POSIXct(coord_all[70,7]) 
comingbacktime[1] == unix2POSIXct(coord_all[82,7])

tripduration <- difftime(comingbacktime, goingouttime, units="secs")

hist(as.numeric(tripduration), nclass=20, xlab="Trip duration (sec)", main="")

pdf("tripdurations.pdf", height=5, width=6.5)
plot(density(as.numeric(tripduration)),xlab="Trip duration (sec)", main="", ylim=c(0,.000015))
lines(density(as.numeric(tripduration[in.out.times[,4] == 2011])), col="red")
lines(density(as.numeric(tripduration[in.out.times[,4] == 2012])), col="blue")
lines(density(as.numeric(tripduration[in.out.times[,4] == 2013])), col="orange")
legend(100000,.00003, col=c("black","red","blue", "orange"), c("all data","2011","2012", "2013"), lty=1)
dev.off()

t.test(as.numeric(tripduration[in.out.times[,4] == 2011]),as.numeric(tripduration[in.out.times[,4] == 2012]))

pdf("tripdurations_stage.pdf", height=5, width=6.5)
plot(density(as.numeric(tripduration)),xlab="Trip duration (sec)", main="", ylim=c(0,.000015))
lines(density(as.numeric(tripduration[in.out.times[,24] == "chick"])), col="red")
lines(density(as.numeric(tripduration[in.out.times[,24] == "egg"])), col="blue")
#lines(density(as.numeric(tripduration[in.out.times[,4] == 2012])), col="blue")
#lines(density(as.numeric(tripduration[in.out.times[,4] == 2013])), col="orange")
legend(150000,.00002, col=c("black","red","blue"), c("all data","chick","egg"), lty=1)
dev.off()

pdf("tripdurations_stage&year.pdf", height=5, width=6.5)
plot(density(as.numeric(tripduration)[in.out.times[,24] == "egg"]),xlab="Trip duration (sec)", main="", ylim=c(0,.000010))
lines(density(as.numeric(tripduration[in.out.times[,24] == "egg" & in.out.times[,4] == 2012])), col="red")
lines(density(as.numeric(tripduration[in.out.times[,24] == "egg" & in.out.times[,4] == 2013])), col="blue")
#lines(density(as.numeric(tripduration[in.out.times[,4] == 2012])), col="blue")
#lines(density(as.numeric(tripduration[in.out.times[,4] == 2013])), col="orange")
legend(150000,.000010, col=c("black","red","blue"), c("all data","egg_2012","egg_2013"), lty=1)
dev.off()


#	Welch Two Sample t-test
#
#data:  as.numeric(tripduration[in.out.times[, 4] == 2011]) and as.numeric(tripduration[in.out.times[, 4] == 2012]) 
#t = -4.3461, df = 22.02, p-value = 0.0002586
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# -48366.99 -17119.62 
#sample estimates:
#mean of x mean of y 
# 17984.30  50727.61 



################
# IPQ modeling #
################
aov1 <- aov(res_all$ipq ~ res_all$dist2col + res_all$year)
summary(aov1)

# for all data together
plot(res_all$dist2col ~ res_all$ipq, xlab="IPQ", ylab="Distance to colony (m)")
ipq.lm <- lm(res_all$dist2col ~ res_all$ipq)
abline(ipq.lm)
summary(ipq.lm)

ipq.dist2col <- data.frame(res_all$ipq,res_all$dist2col)
ipq.dist2col <- ipq.dist2col[!is.na(ipq.dist2col[,2]),]
m1 <- ggplot(data= ipq.dist2col, aes(x = res_all.ipq, y = res_all.dist2col)) + theme_bw() + ylim(05, 56) 
m1ga <- m1 + geom_density2d() + xlab("IPQ") + ylab("Distance to colony (km)") + ggtitle("All data") + geom_abline(intercept = coefficients(ipq.lm)[[1]], slope = coefficients(ipq.lm)[[2]])
#m1ga
ggsave("ipq.dist2col.pdf", plot= m1ga, width=7, height=4.5)


# year by year
ipq.dist2col.2011 <- data.frame(res_all$ipq[res_all$year == 2011],res_all$dist2col[res_all$year == 2011])
ipq.dist2col.2011 <- ipq.dist2col.2011[!is.na(ipq.dist2col.2011[,2]),]
colnames(ipq.dist2col.2011) <- c("res_all.ipq","res_all.dist2col")
ipq.lm.2011 <- lm(ipq.dist2col.2011[,2] ~ ipq.dist2col.2011[,1])
ipq.dist2col.2012 <- data.frame(res_all$ipq[res_all$year == 2012],res_all$dist2col[res_all$year == 2012])
ipq.dist2col.2012 <- ipq.dist2col.2012[!is.na(ipq.dist2col.2012[,2]),]
colnames(ipq.dist2col.2012) <- c("res_all.ipq","res_all.dist2col")
ipq.lm.2012 <- lm(ipq.dist2col.2012[,2] ~ ipq.dist2col.2012[,1])
ipq.dist2col.2013 <- data.frame(res_all$ipq[res_all$year == 2013],res_all$dist2col[res_all$year == 2013])
ipq.dist2col.2013 <- ipq.dist2col.2013[!is.na(ipq.dist2col.2013[,2]),]
colnames(ipq.dist2col.2013) <- c("res_all.ipq","res_all.dist2col")
ipq.lm.2013 <- lm(ipq.dist2col.2013[,2] ~ ipq.dist2col.2013[,1])
#
m2 <- ggplot(data= ipq.dist2col.2011, aes(x = res_all.ipq, y = res_all.dist2col)) + theme_bw() + ylim(05, 56)
m2ga <- m2 + geom_density2d() + xlab("IPQ") + ylab("Distance to colony (km)") + ggtitle("2011") + geom_abline(intercept = coefficients(ipq.lm.2011)[[1]], slope = coefficients(ipq.lm.2011)[[2]])
#m2ga
m3 <- ggplot(data= ipq.dist2col.2012, aes(x = res_all.ipq, y = res_all.dist2col)) + theme_bw() + ylim(05, 56)
m3ga <- m3 + geom_density2d() + xlab("IPQ") + ylab("Distance to colony (km)") + ggtitle("2012") + geom_abline(intercept = coefficients(ipq.lm.2012)[[1]], slope = coefficients(ipq.lm.2012)[[2]])
#m3ga
m4 <- ggplot(data= ipq.dist2col.2013, aes(x = res_all.ipq, y = res_all.dist2col)) + theme_bw() + ylim(05, 56)
m4ga <- m4 + geom_density2d() + xlab("IPQ") + ylab("Distance to colony (km)") + ggtitle("2013") + geom_abline(intercept = coefficients(ipq.lm.2013)[[1]], slope = coefficients(ipq.lm.2013)[[2]])
#m3ga
ggsave("ipq.dist2col.2011.pdf", plot= m2ga, width=7, height=4.5)
ggsave("ipq.dist2col.2012.pdf", plot= m3ga, width=7, height=4.5)
ggsave("ipq.dist2col.2013.pdf", plot= m4ga, width=7, height=4.5)

# by stages
ipq.dist2col.egg <- data.frame(res_all$ipq[res_all$stage == "egg"],res_all$dist2col[res_all$stage == "egg"])
ipq.dist2col.egg <- ipq.dist2col.egg[!is.na(ipq.dist2col.egg[,2]),]
colnames(ipq.dist2col.egg) <- c("res_all.ipq","res_all.dist2col")
ipq.lm.egg <- lm(ipq.dist2col.egg[,2] ~ ipq.dist2col.egg[,1])
ipq.dist2col.chick <- data.frame(res_all$ipq[res_all$stage == "chick"],res_all$dist2col[res_all$stage == "chick"])
ipq.dist2col.chick <- ipq.dist2col.chick[!is.na(ipq.dist2col.chick[,2]),]
colnames(ipq.dist2col.chick) <- c("res_all.ipq","res_all.dist2col")
ipq.lm.chick <- lm(ipq.dist2col.chick[,2] ~ ipq.dist2col.chick[,1])
#
m2 <- ggplot(data= ipq.dist2col.2011, aes(x = res_all.ipq, y = res_all.dist2col)) + theme_bw() + ylim(05, 56)
m2ga <- m2 + geom_density2d() + xlab("IPQ") + ylab("Distance to colony (km)") + ggtitle("incubation") + geom_abline(intercept = coefficients(ipq.lm.2011)[[1]], slope = coefficients(ipq.lm.2011)[[2]])
#m2ga
m3 <- ggplot(data= ipq.dist2col.2012, aes(x = res_all.ipq, y = res_all.dist2col)) + theme_bw() + ylim(05, 56)
m3ga <- m3 + geom_density2d() + xlab("IPQ") + ylab("Distance to colony (km)") + ggtitle("chick rearing") + geom_abline(intercept = coefficients(ipq.lm.2012)[[1]], slope = coefficients(ipq.lm.2012)[[2]])
#m3ga
ggsave("ipq.dist2col.egg.pdf", plot= m2ga, width=7, height=4.5)
ggsave("ipq.dist2col.chick.pdf", plot= m3ga, width=7, height=4.5)



#########
# bouts #
#########

nbirds <- max(res_all$ID,na.rm=T)
bec <- rep(NA,nbirds)

pdf("bouts.pdf", width= 8, height= 3*nbirds)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(nbirds,2))

# computes individual BECs (for each bird)
for(i in 1:nbirds){
	postdives <- res_all$postdive.dur[res_all$ID == i] 		# Postdiveduration based BEC
	if(length(postdives)>10){
		postdives.diff <- abs(diff(postdives))
		postdives.diff <- postdives.diff[postdives.diff < 2000] # remove isolated dives
		lnfreq <- boutfreqs(postdives.diff, bw=0.1, plot=F)
		startval <- boutinit(lnfreq, 50, plot=F)
		p <- startval[[1]]["a"] / (startval[[1]]["a"] + startval[[2]]["a"])
		init.parms <- list(p=as.vector(logit(p)),
		                   lambda1=as.vector(log(startval[[1]]["lambda"])),
		                   lambda2=as.vector(log(startval[[2]]["lambda"])))
		bout.fit1 <- bouts.mle(bouts2.LL, start=init.parms, x=postdives.diff,
		                       method="L-BFGS-B", lower=c(-2, -5, -10))
		coefs <- as.vector(coef(bout.fit1))
		
		## Un-transform and fit the original parameterization
		init.parms <- list(p=unLogit(coefs[1]), lambda1=exp(coefs[2]),
		                   lambda2=exp(coefs[3]))
		bout.fit2 <- bouts.mle(bouts2.ll, x=postdives.diff, start=init.parms,
		                       method="L-BFGS-B", lower=rep(1e-08, 3),
		                       control=list(parscale=c(1, 0.1, 0.01)))
		plotBouts(bout.fit2, postdives.diff, main=paste("Bird no. ",i," (",res_all$year[res_all$ID == i][1],")",sep=""))
		plotBouts2.cdf(bout.fit2, postdives.diff)
		bec[i] <- bec2(bout.fit2)
		labelBouts(postdives, rep(bec[i],length(postdives)), bec.method = "seq.diff")
	}
}

par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()

BEC <- mean(bec,na.rm=T)

#########################################
# Multiple bouts, dive activity and IPQ #
#########################################

# this is using time between consecutive dives -- all dives
ndives <- length(res_all$postdive.dur)
#mb.thres <- 30 # multiple bout interval threshold_plus (time unit: seconds)
mb.thres <- BEC
mb.status <- numeric(ndives)
lag <- numeric(ndives)

for(i in 1:(ndives)){
	lag[i] <- res_all$postdive.dur[i]
	if(lag[i] < mb.thres) mb.status[i] <- 1
}
mb.status
# getting patch IDs
counter <- 1
patch.id <- numeric(ndives)
for(i in 2:ndives){
	if((mb.status[i-1]==0) & (mb.status[i]==1)){ # entering a new patch
		patch.id[i] <- counter
	}else{
		if((mb.status[i-1]==1) & (mb.status[i]==1)){
			patch.id[i] <- counter
		}else{
			if((mb.status[i-1]==1) & (mb.status[i]==0)){
				counter <- counter + 1
			}
		}
	}
}
patch.id


############################
# number of dives per bout #
############################

nbouts <- max(patch.id)
mean.dist <- numeric(nbouts)
mean.bottime <- numeric(nbouts)
mean.spd <- numeric(nbouts)
n.dives <- numeric(nbouts)
for(i in 1:nbouts){
	mean.dist[i] <- mean(res_all$dist2col[patch.id==i])
	n.dives[i] <- sum(patch.id==i)
	mean.bottime[i] <- mean(res_all$botttim[patch.id==i])
	mean.spd[i] <- mean(res_all$postdive.dur[patch.id==i])
}

plot(mean.dist, n.dives, pch=3)
lm1 <- lm(n.dives~ mean.dist)
summary(lm1)
abline(lm1,col="red")

plot(mean.dist, mean.bottime, pch=3)
lm2 <- lm(mean.bottime ~ mean.dist)
summary(lm2)
abline(lm2,col="red")

plot(mean.dist, mean.spd, pch=3)
lm3 <- lm(mean.spd ~ mean.dist)
summary(lm3)
abline(lm3,col="red")



#######################################
# IPQ of last bout vs IPQ other bouts #
#######################################

# get passage IDs from passageID for use with multiple bouts
birdID.mb <- res_all[,2]
passageID.mb <- rep(NA,ndives)
for(j in 1:ndives){
	dive.beg  <- res_all$begdesc[j]
	dive.end  <- res_all$begdesc[j] + res_all$divetim[j]
	rightbird <- coord_all[,2]==birdID.mb[j]
	maxntrips <- max(passageID[coord_all[,2]==birdID.mb[j]],na.rm=T)
	for(i in 1:maxntrips){
		rightpassid         <- passageID==i
		righttimeinterval.1 <- min(as.numeric(coord_all[which(rightbird & rightpassid),7]),na.rm=T)
		righttimeinterval.2 <- max(as.numeric(coord_all[which(rightbird & rightpassid),7]),na.rm=T)
		#passageID.mb[res_all[,2]==birdID.mb[j] & res_all$begdesc>=righttimeinterval.1 & (res_all$begdesc+res_all$divetim)<=righttimeinterval.2] <- i#birdID.mb[j]#i
		passageID.mb[res_all$begdesc>=righttimeinterval.1 & (res_all$begdesc+res_all$divetim)<=righttimeinterval.2] <- i#birdID.mb[j]#i
	}
}
passageID.mb
# a quick plot just to check the match btw passageID and passageID.mb
# incomplete match as some trips have no dives
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,1))
plot(birdID.mb,passageID.mb,ylim=c(0,15))
plot(coord_all[,2],passageID,ylim=c(0,15))
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# this extracts IPQ values for the last trip of each log -- meaningless
last.ipq <- rep(NA,max(res_all[,2],na.rm=T))
other.ipq <- rep(NA,max(res_all[,2],na.rm=T))
for(i in 1:max(res_all[,2],na.rm=T)){
	maxtrip <- max(passageID[coord_all[,2]==i],na.rm=T)
	if(sum(birdID.mb==i & passageID.mb==maxtrip,na.rm=T)>0){ # this should not be here...
		last.ipq[i] <- mean(res_all$ipq[which(birdID.mb==i & passageID.mb==maxtrip)],na.rm=T)
	}
	other.ipq[i] <- mean(res_all$ipq[which(birdID.mb==i & passageID.mb==c(1:(maxtrip-1)))],na.rm=T)
}
t.test(other.ipq,last.ipq,paired=T)
wilcox.test(other.ipq,last.ipq, paired=T) 

# this extracts IPQ values of the last bout of each trip
last_ipq <- other_ipq <- NA
for(i in 1:max(res_all[,2],na.rm=T)){
	maxtrip <- max(passageID[coord_all[,2]==i],na.rm=T)
	if(sum(birdID.mb==i & passageID.mb== maxtrip,na.rm=T)>0){ # in case there's no dive during last trip
		for(j in 1:maxtrip){
			if(sum(birdID.mb==i & passageID.mb==j,na.rm=T)>0){ # in case there's no dive during one trip, but not the last one
				maxpatchid <- max(patch.id[birdID.mb==i & passageID.mb==j],na.rm=T)
				minpatchid <- min(patch.id[birdID.mb==i & passageID.mb==j & patch.id>0],na.rm=T)
				last_ipq <- c(last_ipq,res_all$ipq[which(patch.id==maxpatchid)])
				other_ipq <- c(other_ipq,res_all$ipq[which(patch.id==c(minpatchid:(maxpatchid-1)))])
			}
		}
	}else{
		for(j in 1:(maxtrip-1)){
			if(sum(birdID.mb==i & passageID.mb==j,na.rm=T)>0){
				maxpatchid <- max(patch.id[birdID.mb==i & passageID.mb==j],na.rm=T)
				minpatchid <- min(patch.id[birdID.mb==i & passageID.mb==j & patch.id>0],na.rm=T)
				last_ipq <- c(last_ipq,res_all$ipq[which(patch.id==maxpatchid)])
				other_ipq <- c(other_ipq,res_all$ipq[which(patch.id==c(minpatchid:(maxpatchid-1)))])
			}
		}
	}
}
t.test(other_ipq,last_ipq)

# this extracts IPQ values of the last bout of each trip
last_ipq <- other_ipq <- NA
for(i in 1:max(res_all[,2],na.rm=T)){
	maxtrip <- max(passageID[coord_all[,2]==i],na.rm=T)
	for(j in 1:maxtrip){
		if(sum(birdID.mb==i & passageID.mb==j,na.rm=T)>0){ # in case there's no dives during one trip [skip it]
			maxpatchid <- max(patch.id[birdID.mb==i & passageID.mb==j],na.rm=T)
			minpatchid <- min(patch.id[birdID.mb==i & passageID.mb==j & patch.id>0],na.rm=T)
			last_ipq <- c(last_ipq,res_all$ipq[which(patch.id==maxpatchid)])
			for(k in minpatchid:(maxpatchid-1)){
				other_ipq <- c(other_ipq,res_all$ipq[which(patch.id==k)])
			}
		}
	}
}
t.test(other_ipq,last_ipq)

# this extracts IPQ values of the last bout of each trip *2011*
last_ipq <- other_ipq <- NA
for(i in 1:max(res_all[,2],na.rm=T)){
	if(as.numeric(coord_all[,5][coord_all[,2]==i][1])==2011){
		maxtrip <- max(passageID[coord_all[,2]==i],na.rm=T)
		for(j in 1:maxtrip){
			if(sum(birdID.mb==i & passageID.mb==j,na.rm=T)>0){ # in case there's no dives during one trip [skip it]
				maxpatchid <- max(patch.id[birdID.mb==i & passageID.mb==j],na.rm=T)
				minpatchid <- min(patch.id[birdID.mb==i & passageID.mb==j & patch.id>0],na.rm=T)
				last_ipq <- c(last_ipq,res_all$ipq[which(patch.id==maxpatchid)])
				for(k in minpatchid:(maxpatchid-1)){
					other_ipq <- c(other_ipq,res_all$ipq[which(patch.id==k)])
				}
			}
		}
	}
}
t.test(other_ipq,last_ipq)

# this extracts IPQ values of the last bout of each trip *2012*
last_ipq <- other_ipq <- NA
for(i in 1:max(res_all[,2],na.rm=T)){
	if(as.numeric(coord_all[,5][coord_all[,2]==i][1])==2012){
		maxtrip <- max(passageID[coord_all[,2]==i],na.rm=T)
		for(j in 1:maxtrip){
			if(sum(birdID.mb==i & passageID.mb==j,na.rm=T)>0){ # in case there's no dives during one trip [skip it]
				maxpatchid <- max(patch.id[birdID.mb==i & passageID.mb==j],na.rm=T)
				minpatchid <- min(patch.id[birdID.mb==i & passageID.mb==j & patch.id>0],na.rm=T)
				last_ipq <- c(last_ipq,res_all$ipq[which(patch.id==maxpatchid)])
				for(k in minpatchid:(maxpatchid-1)){
					other_ipq <- c(other_ipq,res_all$ipq[which(patch.id==k)])
				}
			}
		}
	}
}
t.test(other_ipq,last_ipq)


# this extracts IPQ values of each bout of each trip and plots them
mean_bout_ipq <- numeric(0)
for(k in 1:max(patch.id)){
	if(patch.id[k]>0){
		mean_bout_ipq[k] <- mean(res_all$ipq[which(patch.id==k)],na.rm=T)
	}
}
mean_bout_ipq


plot(1000,1000,xlim=c(0,100),ylim=c(min(res_all$ipq,na.rm=T),max(res_all$ipq,na.rm=T)),xlab="Bouts",ylab="Mean IPQ")
last_bout <- 1
for(i in 1:max(res_all[,2],na.rm=T)){
	maxtrip <- max(passageID[coord_all[,2]==i],na.rm=T)
	for(j in 1:maxtrip){
		for(k in 1:max(patch.id)){
			points(k-last_bout,mean(res_all$ipq[which(birdID.mb==i & passageID.mb==j & patch.id==k)],na.rm=T),pch=j)
		}
		last_bout <- zzz # unfinished
	}
}


##################
# something else #
##################

# get mean IPQ within each patch
mean_patch.ipq <- numeric(nbouts)
mean_patch.ndives <- numeric(nbouts)
mean_patch.duration <- numeric(nbouts)
mean_patch.maxdepth <- numeric(nbouts)
# new stuff (13/04/07)
mean_patch.year <- numeric(nbouts)
mean_patch.birdID <- numeric(nbouts)
mean_patch.dist2col <- numeric(nbouts)
all.maxdepth <- numeric(nbouts)         # tracks all individual maxdepths within each patch
sd_patch.maxdepth <- numeric(nbouts)
patch.length <- 0
for(i in 2:ndives){
	if((patch.id[i-1]==0) & !(patch.id[i]==0)){
		patch.length <- patch.length + 1
		mean_patch.ipq[patch.id[i]] <- mean_patch.ipq[patch.id[i]] + res_all$ipq[i]
		mean_patch.duration[patch.id[i]] <- res_all$begdesc[i]
		mean_patch.maxdepth[patch.id[i]] <- res_all$maxdep[i]
		all.maxdepth <- rep(NA, nbouts)
		all.maxdepth[patch.length] <- res_all$maxdep[i]
		mean_patch.year[patch.id[i]] <- res_all$year[i]
		mean_patch.birdID[patch.id[i]] <- res_all$ID[i]
		mean_patch.dist2col[patch.id[i]] <- mean_patch.dist2col[patch.id[i]] + res_all$dist2col[i]
	}else{
		mean_patch.year[patch.id[i]] <- res_all$year[i]
		mean_patch.birdID[patch.id[i]] <- res_all$ID[i]
		if(!(patch.id[i-1]==0) & !(patch.id[i]==0)){
			patch.length <- patch.length + 1
			mean_patch.ipq[patch.id[i]] <- mean_patch.ipq[patch.id[i]] + res_all$ipq[i]
			mean_patch.maxdepth[patch.id[i]] <- mean_patch.maxdepth[patch.id[i]] + res_all$maxdep[i]
			all.maxdepth[patch.length] <- res_all$maxdep[i]
			mean_patch.dist2col[patch.id[i]] <- mean_patch.dist2col[patch.id[i]] + res_all$dist2col[i]
		}else{
			if(!(patch.id[i-1]==0) & (patch.id[i]==0)){
				mean_patch.ipq[patch.id[i-1]] <- mean_patch.ipq[patch.id[i-1]] / patch.length
				mean_patch.ndives[patch.id[i-1]] <- patch.length
				mean_patch.duration[patch.id[i-1]] <- (res_all$begdesc[i-1] + res_all$divetim[i-1]) - mean_patch.duration[patch.id[i-1]]
				mean_patch.maxdepth[patch.id[i-1]] <- mean_patch.maxdepth[patch.id[i-1]] / patch.length
				mean_patch.dist2col[patch.id[i-1]] <- mean_patch.dist2col[patch.id[i-1]] / patch.length
				patch.length <- 0
				#sd_patch.maxdepth[patch.id[i-1]] <- sd(all.maxdepth,na.rm=T)
				sd_patch.maxdepth[patch.id[i-1]] <- var(all.maxdepth,na.rm=T)
			}
		}
	}
}
mean_patch.ipq
mean_patch.ndives
mean_patch.duration
mean_patch.maxdepth
sd_patch.maxdepth # NA values when only one dive in a given patch
mean_patch.year
mean_patch.birdID
mean_patch.dist2col

cor.test(sd_patch.maxdepth, mean_patch.ipq)

pdf("ipq.patch.pdf", width=10, height=6)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2))
plot(mean_patch.ipq ~ mean_patch.ndives, pch=3, xlab="number of dives per patch", ylab="mean patch IPQ")
lm1 <- lm(mean_patch.ipq ~ mean_patch.ndives)
abline(lm1, lty=2, col="gray")
summary(lm1)

plot(mean_patch.ndives ~ mean_patch.duration, pch=3, ylab="Mean number of dives per patch", xlab="Sojourn time in each patch (s)")
lm3 <- lm(mean_patch.ndives ~ mean_patch.duration)
abline(lm3, lty=2, col="gray")
summary(lm3)

#plot(mean_patch.ipq ~ mean_patch.duration, pch=3, xlab="length of sojourn in patch (s)", ylab="mean patch IPQ", log="x")
plot(mean_patch.ipq ~ mean_patch.duration, pch=3, xlab="Sojourn time in each patch (s)", ylab="mean patch IPQ")
lm2.ls <- loess(mean_patch.ipq ~ mean_patch.duration)
xv <- seq(min(mean_patch.duration,na.rm=T),max(mean_patch.duration,na.rm=T),.01)
pred.lm2 <- predict(lm2.ls,xv)
points(xv, pred.lm2, pch=19, cex=.05, col="red")
lm2 <- lm(mean_patch.ipq ~ mean_patch.duration)
abline(lm2, lty=2, col="gray")
summary(lm2)

#plot(mean_patch.ipq ~ mean_patch.maxdepth, pch=3, xlab="maximum depth of dives per patch", ylab="mean patch IPQ")
#lm4 <- lm(mean_patch.ipq ~ mean_patch.maxdepth)
#abline(lm4, lty=2, col="gray")
#summary(lm4)

#plot(mean_patch.ipq ~ sd_patch.maxdepth, pch=3, xlab="SD of dive depths per patch", ylab="mean patch IPQ")
plot(mean_patch.ipq ~ sd_patch.maxdepth, pch=3, xlab="Variance of dive depths per patch (log)", ylab="mean patch IPQ",log="x")
lm4.ls <- loess(mean_patch.ipq ~ sd_patch.maxdepth)
xv <- seq(min(sd_patch.maxdepth,na.rm=T),max(sd_patch.maxdepth,na.rm=T),.01)
pred.lm4 <- predict(lm4.ls,xv)
points(xv, pred.lm4, pch=19, cex=.05, col="red")
#lm4 <- lm(mean_patch.ipq ~ sd_patch.maxdepth)
lm4 <- lm(mean_patch.ipq[sd_patch.maxdepth>0] ~ log(sd_patch.maxdepth[sd_patch.maxdepth>0]))
abline(lm4, lty=2, col="gray")
summary(lm4)

par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,1))
dev.off()


t.test(res_all$ipq[(mb.status == 1) & (res_all$dist2col> threshold_plus)], res_all$ipq[(mb.status == 0) & (res_all$dist2col> threshold_plus)])
# 2 sec: p-value = 2.144e-07	| 0.09283325 0.15889242  
# 5 sec: p-value = 1.188e-11	| 0.09655938 0.16141458
#10 sec: p-value = 1.005e-10	| 0.1000700 0.1636039
#30 sec: p-value = 0.10846		| 0.1523783 0.1619223
#60 sec: p-value = 0.9268		| 0.1559623 0.1549635 


# this is using space between consecutive dives
#mb.thres <- 100 # multiple bout interval threshold_plus (space unit: meter)
#mb.status <- numeric(ndives)
#lag <- numeric(ndives)

#for(i in 1:(ndives-1)){
#	lag[i] <- distMeeus(c(interpol.lon[i],interpol.lat[i]),c(interpol.lon[i+1],interpol.lat[i+1]))
#	if(lag[i] < mb.thres) mb.status[i] <- 1
#}
#mb.status
#plot(mb.status)
#t.test(ipq[mb.status == 1], ipq[mb.status == 0])


ncolor <- 7
mypalette<-brewer.pal(ncolor,"BuPu")
ipq <- as.numeric(coord_all[,9])
ipq.col <- rep(mypalette[1],length(ipq))
min.ipq <- summary(ipq)[[2]]
ipq.col[ipq < min.ipq] <- mypalette[1]
ipq.col[ipq > min.ipq] <- mypalette[1 + floor(ncolor * (-min(ipq[ipq > min.ipq]) + ipq[ipq > min.ipq])/max(-min(ipq[ipq > min.ipq]) + ipq[ipq > min.ipq]))]

pdf("ipq.locations.pdf", width=10, height=5)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,2))
plot(coord_all[,3], coord_all[,4], pch=3, col=ipq.col, xlab="longitude", ylab="latitude")
legend(min(as.numeric(coord_all[,3])),max(as.numeric(coord_all[,4])),col=mypalette, pch=3, c("low ipq","2","3","4","5","6","high ipq"))
plot(as.numeric(coord_all[,3])[ipq<summary(ipq)[[2]]], as.numeric(coord_all[,4])[ipq<summary(ipq)[[2]]], pch=3, col=mypalette[1], xlab="longitude", ylab="latitude")
points(as.numeric(coord_all[,3])[ipq>summary(ipq)[[5]]], as.numeric(coord_all[,4])[ipq>summary(ipq)[[5]]], pch=3, col=mypalette[ncolor])
legend(min(as.numeric(coord_all[,3])[ipq>summary(ipq)[[5]]]),max(as.numeric(coord_all[,4])[ipq>summary(ipq)[[5]]]),col=c(mypalette[1], mypalette[ncolor]), pch=3, c("ipq < 1st quartile of ipq","ipq > 3rd quartile of ipq"))
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()

#map('world', interior=T, col="gray80", fill=T, xlim=c(min(as.numeric(coord_all[,3]))-.5,max(as.numeric(coord_all[,3]))+.5), ylim=c(min(as.numeric(coord_all[,4]))-.5,max(as.numeric(coord_all[,4]))+.5))
#points(coord_all[,3], coord_all[,4], pch=3, col=ipq.col, xlab="longitude", ylab="latitude")
#legend(min(as.numeric(coord_all[,3]))-.5,max(as.numeric(coord_all[,4]))+.5,col=mypalette, pch=3, c("low ipq","2","3","4","5","6","high ipq"))


plot(unix2POSIXct(res_all$begdesc[res_all$ID==1]), res_all$ilon[res_all$ID==1])


################################################################
# these line are experimental and should not be run -- for now #
################################################################
#
#birdID <- 3
#df4pca <- data.frame(res_all$begdesc[res_all$ID== birdID],res_all$dist2col[res_all$ID== birdID],res_all$ilon[res_all$ID== birdID],res_all$ilat[res_all$ID== birdID])
#colnames(df4pca) <- c("desctimes","dist2col")
#pcadat <- princomp(df4pca, cor = TRUE)
#bestk <- silcheck(data.frame(pcadat$scores[,1],pcadat$scores[,2],pcadat$scores[,3]),kmax=100)[1]
#bestk
#bouts <- pam(data.frame(pcadat$scores[,1],pcadat$scores[,2],pcadat$scores[,3]), k=bestk)
#plot(bouts, which.plots = 1)




#########################################
# number of dives on out-/inbound trips #
#########################################
# this depends on threshold_plus as defined above (thru passageID)

# initiallize vector of indices for out- and inbound trips
ndives <- length(res_all_ori[,1])
outinbound_index <- numeric(ndives)

# place 1 each time distance to colony > threshold_plus: this defines individual trips
#outinbound_index[res_all_ori$dist2col > 80] <- 1
outinbound_index[res_all_ori$dist2col > 1500] <- 1
#outinbound_index[res_all_ori$dist2col > threshold_plus] <- 1

maxdist <- 0
diving_bird <- 0
for(i in 1:(ndives-1)){
	if((outinbound_index[i]==1) & (outinbound_index[i+1]==1) & (res_all_ori$dist2col[i] > maxdist)){
		maxdist <- res_all_ori$dist2col[i]
	}
	else{
		if((outinbound_index[i]==1) & (outinbound_index[i+1]==1)& (res_all_ori$dist2col[i] < (maxdist-2000))){
			outinbound_index[i] <- -1
		}
		else{
			if((outinbound_index[i]==0) & (outinbound_index[i+1]==1)){
				maxdist <- 0
			}
			else{
				if((outinbound_index[i]==1) & (outinbound_index[i+1]==0)){
					outinbound_index[i] <- -1
				}
			}
		}
	}
}
outinbound_index

pdf("outinbound.pdf",width=12,height=8)
mycolor <- rep("blue", ndives)
mycolor[outinbound_index<0] <- "red"
mycolor[outinbound_index==0] <- "gray"
plot(res_all_ori$dist2col, type="o", cex=.5, col=mycolor,xlim=c(0,1000))
abline(h= threshold_plus, col="gray")
par(new=TRUE)
mycolor[mycolor=="red"] <- "orange"
mycolor[mycolor=="blue"] <- "green"
plot(res_all_ori$ipq, type="o", cex=.5, col=mycolor,xlim=c(0,1000),ylim=c(-1,1))
axis(4)
mtext("ipq",side=4,line=3)
dev.off()

# difference in IPQ between outbound and inbound trips?
ipq.blue <- res_all_ori$ipq[outinbound_index>0]
ipq.red <- res_all_ori$ipq[outinbound_index<0]
t.test(ipq.blue, ipq.red)
#	Welch Two Sample t-test
#
#data:  ipq.blue and ipq.red 
#t = 1.5717, df = 5228.902, p-value = 0.1161
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# -0.002387728  0.021695883 
#sample estimates:
# mean of x  mean of y 
#0.09967007 0.09001599 

# - when red (inbound) is followed by blue (outbound), red IPQ increases?
# - when red (inbound) is last, IPQ decreases?
in.ipq <- rep(NA, ndives)
in1.ipq <- rep(NA, ndives) # ipq inbound values when bird reaches colony
in1res.ipq <- rep(NA, ndives) # ipq inbound values when bird reaches colony
in2.ipq <- rep(NA, ndives) # ipq inbound values when bird goes back out fishing
in2res.ipq <- rep(NA, ndives) # ipq inbound values when bird goes back out fishing
i_vect <- c(1:ndives)
for(i in 1:(ndives-1)){
	if( ((outinbound_index[i]>0) & (outinbound_index[i+1]<0)) || ((outinbound_index[i]<0) & (outinbound_index[i+1]<0)) ){
		in.ipq[i] <- res_all_ori$ipq[i]
	}else{
		if((outinbound_index[i]<0) & (outinbound_index[i+1]>0)){ # goes back out fishing
			in.ipq[i] <- res_all_ori$ipq[i]
			in2.ipq <- in.ipq
			in.ipq <- rep(NA, ndives)
			
			if(length(in2.ipq[!is.na(in2.ipq)])>4){ #60
				x <- summary(lmRob(in2.ipq~i_vect))$coefficients
				#if(x[2,4]<.05){
					in2res.ipq[i] <- x[2,1] # average slope of IPQ values when bird goes back out fishing
				#}
			}
		}else{
			if((outinbound_index[i]<0) & (outinbound_index[i+1]==0)){ # returned to colony
				in.ipq[i] <- res_all_ori$ipq[i]
				in1.ipq <- in.ipq
				in.ipq <- rep(NA, ndives)
				if(length(in1.ipq[!is.na(in1.ipq)])>4){
					x <- summary(lmRob(in1.ipq~i_vect))$coefficients
					in1res.ipq[i] <- x[2,1] # average slope of IPQ values when bird goes back to colony
				}
			}
		}
	}
}


#plot(density(in2res.ipq,na.rm=T),col="red",xlim=c(-.1,max(c(in1res.ipq,in2res.ipq),na.rm=T)))
#lines(density(in1res.ipq,na.rm=T))

bigdata <- data.frame(cbind(res_all_ori, in1res.ipq, in2res.ipq))
p1 <- qplot(factor(year), in1res.ipq, fill=factor(year), data=bigdata, geom="boxplot", position="dodge")+theme_bw()+scale_x_discrete(name="Year")+scale_y_continuous(name="mean IPQ slope of returning bird\n that goes back to colony")
p2 <- qplot(factor(year), in2res.ipq, fill=factor(year), data=bigdata, geom="boxplot", position="dodge")+theme_bw()+scale_x_discrete(name="Year")+scale_y_continuous(name="mean IPQ slope of returning bird\n that goes back out fishing again")
mp <- grid.arrange(p1, p2, nrow=1)
#pdf("ipq.dive.choice.pdf",width=11,height=4.25)
#print(grid.arrange(p1, p2, nrow=1))
#dev.off()
p <- arrangeGrob(p1, p2, ncol=2)
ggsave("ipq.dive.choice.pdf", plot=p, width=11, height=4.25)


t.test(in1res.ipq[res_all_ori$year == 2011],in1res.ipq[res_all_ori$year == 2012])
t.test(in2res.ipq[res_all_ori$year == 2011],in2res.ipq[res_all_ori$year == 2012])

#################################################
# Plots all tracks on same figure by ipq values #
#################################################

if(online){
	all_tracks <- data.frame(coord_all[,1],as.numeric(coord_all[,2]),as.numeric(coord_all[,3]),as.numeric(coord_all[,4]),as.numeric(coord_all[,9]),as.numeric(coord_all[,5]))
	colnames(all_tracks) <- c("name","ID","lon","lat","ipq","year")
	SHmap <- qmap(c(lon=mean(as.numeric(coord_all[,3]), na.rm=T), lat=mean(as.numeric(coord_all[,4]), na.rm=T)), maptype = 'satellite', zoom=9)
	p <- SHmap + 
		geom_path(data= all_tracks, aes(lon, lat, col=round(-min(all_tracks$ipq)+all_tracks$ipq)), lineend="round") #+ 
		#+ scale_colour_gradient(low="red",high="blue",mid="green",midpoint=.1)
		#opts(legend.position="left", legend.text = theme_text(colour="black", size = 6)) +
		#geom_point(data=data.frame(res_all$ilon,res_all$ilat), aes(res_all$ilon,res_all$ilat), col="black", pch=3, cex=.9)
	#p
	ggsave("all.google.trax.ipq.pdf", plot=p, width=11, height=11)
}


# add bathymetry data
 gebco_dir  = "/Users/stephane/Documents/data/akiko/130305_bathymetry/bakabathy/"
 gebco_file = paste(gebco_dir, "gebco_08_-10_48_-2_53.nc", sep="") 

 ex.nc = open.ncdf(gebco_file)

 # Retrieve chl_oc5 variable:
 x_range = get.var.ncdf( ex.nc, "x_range")
 y_range = get.var.ncdf( ex.nc, "y_range")
 z_range = get.var.ncdf( ex.nc, "z_range")
 spacing = get.var.ncdf( ex.nc, "spacing")
 z = get.var.ncdf( ex.nc, "z")

 # Convert to a 2D array
 bath_pa <- array(z,dim=c(960,610))

 # Reverse latitudes (order of columns)
 bath_pa <- bath_pa[, rev(seq_len(ncol(bath_pa)))]

 # Construct longitude and latitude grids
 # 
 # Longitude: -10 to -2 
 # Latitude: 48 to 53.0833
 #
 # Lon/Lat spacing: 0.0083333 0.0083333
 #
 # Boxes span 30 arc seconds (corresponding to 0.0083333 degrees). Data coverage down to 
 # a longitude of -10 is therefore provided by a point centred at -9deg,59',45" (-9.995833).
 # Similarly, data coverage up to -2 is provided by a point at -2deg,00',15" (-2.004167).
 #
 # The lon/lat grid therefore has the following range of elements:
 # Longitude range (elements) -9.995833 (-9 deg,59',45") to -2.004167 (-2 deg,00',15") 
 # Latitude range (elements) 48.00417 (48 deg,59',15") to 53.07917 (53 deg,59',45")

 # 960 longitude (x) elements
 # 610 latitude (y) elements

 # Construct x/y sequences with 0.0083333 spacing
 geb_x <- seq(x_range[1], x_range[2], by = spacing[1])
 geb_y <- seq(y_range[1], y_range[2], by = spacing[2])

 # Offset by (0.0083333/2)
 geb_x <- geb_x + (spacing[1]/2)
 geb_y <- geb_y + (spacing[2]/2)

 # Determine the correct length
 last_x = length(geb_x)-1
 last_y = length(geb_y)-1

 # Discard the last element
 geb_lon <- geb_x[1:last_x]
 geb_lat <- geb_y[1:last_y]

 # Replicate longitude over columns
 geb_lon <- cbind( do.call(cbind, rep(list(geb_lon), 610)) ) 

 # Replicate latitude over rows
 geb_lat <- cbind( do.call(cbind, rep(list(geb_lat), 960)) )
 geb_lat <- t(geb_lat)

 # Concatinate into one large array: 
 # Latitude, Longitude, Bathymetry

 bath_arr <- array(0,dim=c(960,610,3)) 

 bath_arr[,,1] <- geb_lat
 bath_arr[,,2] <- geb_lon
 bath_arr[,,3] <- bath_pa

 #---------------- PLOT BATHYMETRY -----------------------

  # Reverse columns for plotting the correct way up:
  #bath_pa <- bath_pa[, rev(seq_len(ncol(bath_pa)))]

  min_val <- 900 #min(bath_pa,na.rm = TRUE)        
  max_val <- -2000 #max(bath_pa,na.rm = TRUE)

	col.loc1 <- c(-5.304686, 51.732764) # 2011
	col.loc2 <- c(-5.300254, 51.745299) # 2012
	# Plot chl_oc5 data: 
	pdf("bathymetry_tracks1.pdf")
	filled.contour(
  				#x = seq(-10, -2, length.out = nrow(bath_pa)),
               	#y = seq(48, 53, length.out = ncol(bath_pa)),
  				x = seq(-6, -4.8, length.out = nrow(bath_pa)),
               	y = seq(51, 52.5, length.out = ncol(bath_pa)),
               	bath_pa,
               	zlim = range(min_val, max_val, finite = TRUE),
               	color.palette = rainbow, nlevels = 200,
               	plot.title = title(main = "Bathymetry"), 
               	key.title = title(main = "m")
               	)
    dev.off()
	pdf("bathymetry_tracks2.pdf")
	plot(NULL,xlim=c(-6, -4.8),ylim=c(51, 52.5),xlab="Longitude",ylab="Latitude")
 	abline(v=-6); abline(v=-4.8)
 	abline(h=51.2); abline(h=52.4)
	lines(all_tracks[,3][all_tracks[,6]==2011], all_tracks[,4][all_tracks[,6]==2011],col="yellow")
	lines(all_tracks[,3][all_tracks[,6]==2012], all_tracks[,4][all_tracks[,6]==2012],col="green")
	dev.off()


##############
# bathymetry #
##############

bathy_grid <- read.csv("/Users/stephane/Documents/data/akiko/130305_bathymetry/bathy_arrays/bathy_depth.csv")
lon_grid <- read.csv("/Users/stephane/Documents/data/akiko/130305_bathymetry/bathy_arrays/bathy_lon.csv")
lat_grid <- read.csv("/Users/stephane/Documents/data/akiko/130305_bathymetry/bathy_arrays/bathy_lat.csv")

# first for all GPS locations -- regression limited to bathy < 0
ncoords <- length(coord_all[,1])
bathy <- rep(NA, ncoords)
x1 <- y1 <- rep(NA, ncoords)
for(i in 1: ncoords){
	x1[i] <- which(lon_grid[,2]>=as.numeric(coord_all[i,3]))[1] #lon_grid[which(lon_grid[,2]>=as.numeric(coord_all[i,3]))[1],2]
	y1[i] <- which(lat_grid[2,]>=as.numeric(coord_all[i,4]))[1] #lat_grid[2,which(lat_grid[2,]>=as.numeric(coord_all[i,4]))[1]]
	bathy[i] <- (bathy_grid[x1[i],y1[i]] + bathy_grid[x1[i]+1,y1[i]] +
					bathy_grid[x1[i],y1[i]+1] + bathy_grid[x1[i]+1,y1[i]+1])/4
	print(paste("Now doing",i,": (",x1[i],",",y1[i],") ->",bathy[i],sep=" "))
}

pdf("ipq_bathy1.pdf",width=6,height=4)
smoothScatter(bathy[bathy<0],coord_all[,9][bathy<0],xlab="Bathymetry (m)",ylab="log IPQ")
lm.bath1 <- lm(as.numeric(coord_all[,9][bathy<0])~bathy[bathy<0])
summary(lm.bath1)
abline(lm.bath1,col="red")
dev.off()

pdf("ipq_dist2col.pdf",width=6,height=4)
smoothScatter(res_all[,4]~res_all[,5],xlab="Dist2col (m)",ylab="log IPQ")
#lmer.fit <- lmer(res_all[,4] ~ res_all[,5] + res_all[,3] + (1|res_all[,1]))
lmer.fit <- lmer(res_all[,4] ~ res_all[,5] + (1|res_all[,1]))
summary(lmer.fit)
abline(coef=c(lmer.fit@fixef[1],lmer.fit@fixef[2]), col = "red")

smoothScatter(res_all[,4][res_all[,5]>10000]~res_all[,5][res_all[,5]>10000],xlab="Dist2col (m)",ylab="log IPQ")
#lmer.fit <- lmer(res_all[,4] ~ res_all[,5] + res_all[,3] + (1|res_all[,1]))
lmer.fit <- lmer(res_all[,4][res_all[,5]>10000] ~ res_all[,5][res_all[,5]>10000] + (1|res_all[,1][res_all[,5]>10000]))
summary(lmer.fit)
abline(coef=c(lmer.fit@fixef[1],lmer.fit@fixef[2]), col = "red")

smoothScatter(res_all[,4][res_all[,3]==2011]~res_all[,5][res_all[,3]==2011],xlab="Dist2col (m)",ylab="log IPQ")
#lmer.fit <- lmer(res_all[,4] ~ res_all[,5] + res_all[,3] + (1|res_all[,1]))
lmer.fit <- lmer(res_all[,4][res_all[,3]==2011] ~ res_all[,5][res_all[,3]==2011] + (1|res_all[,1][res_all[,3]==2011]))
summary(lmer.fit)
abline(coef=c(lmer.fit@fixef[1],lmer.fit@fixef[2]), col = "red")

smoothScatter(res_all[,4][res_all[,3]==2012]~res_all[,5][res_all[,3]==2012],xlab="Dist2col (m)",ylab="log IPQ")
#lmer.fit <- lmer(res_all[,4] ~ res_all[,5] + res_all[,3] + (1|res_all[,1]))
lmer.fit <- lmer(res_all[,4][res_all[,3]==2012] ~ res_all[,5][res_all[,3]==2012] + (1|res_all[,1][res_all[,3]==2012]))
summary(lmer.fit)
abline(coef=c(lmer.fit@fixef[1],lmer.fit@fixef[2]), col = "red")
dev.off()

# only for dives
ndives <- length(res_all[,1])
bathy2 <- rep(NA, ndives)
x2 <- y2 <- rep(NA, ndives)
for(i in 1: ndives){
	x2[i] <- which(lon_grid[,2]>=res_all$ilon[i])[1] #lon_grid[which(lon_grid[,2]>=as.numeric(coord_all[i,3]))[1],2]
	y2[i] <- which(lat_grid[2,]>=res_all$ilat[i])[1] #lat_grid[2,which(lat_grid[2,]>=as.numeric(coord_all[i,4]))[1]]
	bathy2[i] <- (bathy_grid[x2[i],y2[i]] + bathy_grid[x2[i]+1,y2[i]] +
					bathy_grid[x2[i],y2[i]+1] + bathy_grid[x2[i]+1,y2[i]+1])/4
	print(paste("Now doing",i,": (",x2[i],",",y2[i],") ->",bathy2[i],sep=" "))
}

#pdf("ipq_bathy2.pdf",width=6,height=4)
#smoothScatter(bathy2[bathy2<0],res_all$ipq[bathy2<0],xlab="Bathymetry (m -- only for dives)",ylab="log IPQ")
#lm.bath2 <- lm(res_all$ipq[bathy2<0]~bathy2[bathy2<0])
#summary(lm.bath2)
#abline(lm.bath2,col="red")
#dev.off()

pdf("ipq_bathy3.pdf",width=6,height=4)
smoothScatter(bathy2[bathy2<0],-res_all$maxdep[bathy2<0],xlab="Bathymetry (m)",ylab="Dive depth (m)")
lm.bath <- lm(-res_all$maxdep[bathy2<0]~bathy2[bathy2<0])
summary(lm.bath)
abline(lm.bath,col="red")
dev.off()


##########
# 150810 #
##########


mean(as.numeric(coord_all[,8]), na.rm=T)
mean(as.numeric(coord_all[,8][coord_all[,10] == "egg"]), na.rm=T)
sd(as.numeric(coord_all[,8][coord_all[,10] == "egg"]), na.rm=T)
mean(as.numeric(coord_all[,8][coord_all[,10] == "chick"]), na.rm=T)
sd(as.numeric(coord_all[,8][coord_all[,10] == "chick"]), na.rm=T)

mean(as.numeric(in.out.times[,7]), na.rm=T)
sd(as.numeric(in.out.times[,7]), na.rm=T)
min(as.numeric(in.out.times[,7]), na.rm=T)
max(as.numeric(in.out.times[,7]), na.rm=T)
mean(as.numeric(in.out.times[,7][in.out.times[,24] == "chick"]), na.rm=T)
sd(as.numeric(in.out.times[,7][in.out.times[,24] == "chick"]), na.rm=T)
mean(as.numeric(in.out.times[,7][in.out.times[,24] == "egg"]), na.rm=T)
sd(as.numeric(in.out.times[,7][in.out.times[,24] == "egg"]), na.rm=T)


########################
# Primary productivity #
# 131109               #
########################

# approx location of Region Of Interest (ROI; contra Rest Of Area: ROA)
roi_x1 <- -5.75; roi_y1 <- 51.65
roi_x2 <- -5.3; roi_y2 <- 51.825

pdf("chl_zonetest.pdf")
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="fly"], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="fly"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: fly")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="fly"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="fly"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
abline(h= roi_y1,col="red")
abline(h= roi_y2,col="red")
abline(v= roi_x1,col="red")
abline(v= roi_x2,col="red")
dev.off()

# begin/end of each period
begin_egg <- unix2POSIXct(min(as.numeric(coord_all[,7][coord_all[,10]=="egg"])))
end_egg <- unix2POSIXct(max(as.numeric(coord_all[,7][coord_all[,10]=="egg"])))
begin_chick <- unix2POSIXct(min(as.numeric(coord_all[,7][coord_all[,10]=="chick"])))
end_chick <- unix2POSIXct(max(as.numeric(coord_all[,7][coord_all[,10]=="chick"])))
# and of course...
end_egg < begin_chick
# mid point around June 15
mid_point <- as.Date("06/15/2011", "%m/%d/%Y")
# take mid point? keep overlap?
begin_egg <- min(as.numeric(coord_all[,7][coord_all[,10]=="egg"]))
end_egg <- max(as.numeric(coord_all[,7][coord_all[,10]=="egg"]))
begin_chick <- min(as.numeric(coord_all[,7][coord_all[,10]=="chick"]))
end_chick <- max(as.numeric(coord_all[,7][coord_all[,10]=="chick"]))


############
# OC5 data #
############
#load("/Users/stephane/Documents/data/akiko/130305_bathymetry/whoisbaka.RData")
# Recall:
# ocean_arr [dimensions (512,512,8)]: date(ddmmyyyy) time(hhmm) lat lon chl_oc5 sst windspeed windangle 
# e.g. ocean_arr[100,100,]: 01012010 1415 49.0349349975586 -8.453125 -999 11.1999998092651 3.37868094444275 -69.492790222168
# NB. MODIS screwed up, lat is lon and lon is lat...
ndives <- length(res_all[,1])
lon_roi_egg <- lon_roi_chick <- lon_notroi_egg <- lon_notroi_chick <- c()
lat_roi_egg <- lat_roi_chick <- lat_notroi_egg <- lat_notroi_chick <- c()
chl_roi_egg <- chl_roi_chick <- chl_notroi_egg <- chl_notroi_chick <- c()
sst_roi_egg <- sst_roi_chick <- sst_notroi_egg <- sst_notroi_chick <- c()
x2 <- y2 <- rep(NA, ndives)


#library(foreach)
#library(doMC)
#registerDoMC(cores=4)

#foreach(k=1: ndives, .combine=c) %dopar%{
for(k in 1:ndives){
	print(paste("Now doing dive ",k," out of ",ndives," in ", res_all$name[k],sep=""))
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]]    
	mth <- as.numeric(substr(as.character(datetime),6,7)) #unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)
	
	# check stage first
	if(res_all$stage[k] == "egg"){
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     			   
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
			   #ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
			   lon_roi_egg <- c(lon_roi_egg,as.numeric(ocean_arr[1,mylon,4]))
			   lat_roi_egg <- c(lat_roi_egg,as.numeric(ocean_arr[mylat,1,3]))
			   chl_roi_egg <- c(chl_roi_egg,as.numeric(ocean_arr[mylat,mylon,5]))
			   sst_roi_egg <- c(sst_roi_egg,as.numeric(ocean_arr[mylat,mylon,6]))
			   
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
		   # MODIS script ends here
		
		   # sab script starts over from here
		   #k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}else{ # outside of ROI
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
				#ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
				lon_notroi_egg <- c(lon_notroi_egg,as.numeric(ocean_arr[1,mylon,4]))
				lat_notroi_egg <- c(lat_notroi_egg,as.numeric(ocean_arr[mylat,1,3]))
				chl_notroi_egg <- c(chl_notroi_egg,as.numeric(ocean_arr[mylat,mylon,5]))
				sst_notroi_egg <- c(sst_notroi_egg,as.numeric(ocean_arr[mylat,mylon,6]))
				
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
			# MODIS script ends here
		
			# sab script starts over from here
			#k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}
	}else{ # res_all$stage == "chick"
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
				#ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
				lon_roi_chick <- c(lon_roi_chick,as.numeric(ocean_arr[1,mylon,4]))
				lat_roi_chick <- c(lat_roi_chick,as.numeric(ocean_arr[mylat,1,3]))
				chl_roi_chick <- c(chl_roi_chick,as.numeric(ocean_arr[mylat,mylon,5]))
				sst_roi_chick <- c(sst_roi_chick,as.numeric(ocean_arr[mylat,mylon,6]))
				
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
			# MODIS script ends here
		
			# sab script starts over from here
			#k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}else{ # outside of ROI
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
				#ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
				lon_notroi_chick <- c(lon_notroi_chick,as.numeric(ocean_arr[1,mylon,4]))
				lat_notroi_chick <- c(lat_notroi_chick,as.numeric(ocean_arr[mylat,1,3]))
				chl_notroi_chick <- c(chl_notroi_chick,as.numeric(ocean_arr[mylat,mylon,5]))
				sst_notroi_chick <- c(sst_notroi_chick,as.numeric(ocean_arr[mylat,mylon,6]))
				
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
			# MODIS script ends here
		
			# sab script starts over from here
			#k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}
	}

}

sum(!is.na(chl_notroi_egg[chl_notroi_egg>0]))
sum(!is.na(chl_notroi_chick[chl_notroi_chick>0]))
sum(!is.na(chl_roi_egg[chl_roi_egg>0]))
sum(!is.na(chl_roi_chick[chl_roi_chick>0]))


plot3d(lat_roi_chick[chl_roi_chick>0], lon_roi_chick[chl_roi_chick>0],chl_roi_chick[chl_roi_chick>0])

plot(lat_roi_chick,lon_roi_chick)

pdf("chl_satsweeps.pdf")
smoothScatter(ocean_arr[,,4][ocean_arr[,,5]>0],ocean_arr[,,3][ocean_arr[,,5]>0],colramp = colorRampPalette(c("white", "green")),xlab="Longitude",ylab="Latitude")
points(ocean_arr[,,4][ocean_arr[,,5]>0],ocean_arr[,,3][ocean_arr[,,5]>0],pch=".")
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
abline(h= roi_y1,col="red")
abline(h= roi_y2,col="red")
abline(v= roi_x1,col="red")
abline(v= roi_x2,col="red")
abline(h=max(as.numeric(coord_all[,4]),na.rm=T),col="blue")
abline(h=min(as.numeric(coord_all[,4]),na.rm=T),col="blue")
abline(v=max(as.numeric(coord_all[,3]),na.rm=T),col="blue")
abline(v=min(as.numeric(coord_all[,3]),na.rm=T),col="blue")
dev.off()


unix2POSIXct(min(as.numeric(ocean_arr[,,1])))

ocean_arr[,,5][as.numeric(ocean_arr[,,1])>= begin_egg & as.numeric(ocean_arr[,,1])<= end_egg]

####################################################################
save.image(paste("tdr_new", mytdr_d.thresh,".", mytdr_t.thresh, ".RData", sep=""))
####################################################################

# 1
roi_x1 <- -5.75; roi_y1 <- 51.65
roi_x2 <- -5.3; roi_y2 <- 51.825
ndives <- length(res_all[,1])
lon_roi_egg <- lon_roi_chick <- lon_notroi_egg <- lon_notroi_chick <- c()
lat_roi_egg <- lat_roi_chick <- lat_notroi_egg <- lat_notroi_chick <- c()
chl_roi_egg <- chl_roi_chick <- chl_notroi_egg <- chl_notroi_chick <- c()
sst_roi_egg <- sst_roi_chick <- sst_notroi_egg <- sst_notroi_chick <- c()
x2 <- y2 <- rep(NA, ndives)
for(k in 1:ndives){
	print(paste("Now doing dive ",k," out of ",ndives," in ", res_all$name[k],sep=""))
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]]    
	mth <- as.numeric(substr(as.character(datetime),6,7)) #unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)
	
	# check stage first
	if(res_all$stage[k] == "egg"){
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     			   
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
			   #ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
			   lon_roi_egg <- c(lon_roi_egg,as.numeric(ocean_arr[1,mylon,4]))
			   lat_roi_egg <- c(lat_roi_egg,as.numeric(ocean_arr[mylat,1,3]))
			   chl_roi_egg <- c(chl_roi_egg,as.numeric(ocean_arr[mylat,mylon,5]))
			   sst_roi_egg <- c(sst_roi_egg,as.numeric(ocean_arr[mylat,mylon,6]))
			   
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
		   # MODIS script ends here
		
		   # sab script starts over from here
		   #k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}else{ # outside of ROI

		}
	}else{ # res_all$stage == "chick"
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI

		}
	}

}
save(lon_roi_egg, lat_roi_egg, chl_roi_egg, sst_roi_egg, file="roi_egg.RData")

# 2
roi_x1 <- -5.75; roi_y1 <- 51.65
roi_x2 <- -5.3; roi_y2 <- 51.825
ndives <- length(res_all[,1])
lon_roi_egg <- lon_roi_chick <- lon_notroi_egg <- lon_notroi_chick <- c()
lat_roi_egg <- lat_roi_chick <- lat_notroi_egg <- lat_notroi_chick <- c()
chl_roi_egg <- chl_roi_chick <- chl_notroi_egg <- chl_notroi_chick <- c()
sst_roi_egg <- sst_roi_chick <- sst_notroi_egg <- sst_notroi_chick <- c()
x2 <- y2 <- rep(NA, ndives)
for(k in 1:ndives){
	print(paste("Now doing dive ",k," out of ",ndives," in ", res_all$name[k],sep=""))
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]]    
	mth <- as.numeric(substr(as.character(datetime),6,7)) #unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)
	
	# check stage first
	if(res_all$stage[k] == "egg"){
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
				#ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
				lon_notroi_egg <- c(lon_notroi_egg,as.numeric(ocean_arr[1,mylon,4]))
				lat_notroi_egg <- c(lat_notroi_egg,as.numeric(ocean_arr[mylat,1,3]))
				chl_notroi_egg <- c(chl_notroi_egg,as.numeric(ocean_arr[mylat,mylon,5]))
				sst_notroi_egg <- c(sst_notroi_egg,as.numeric(ocean_arr[mylat,mylon,6]))
				
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
			# MODIS script ends here
		
			# sab script starts over from here
			#k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}
	}else{ # res_all$stage == "chick"
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI

		}
	}

}
save(lon_notroi_egg, lat_notroi_egg, chl_notroi_egg, sst_notroi_egg, file="notroi_egg.RData")

# 3
roi_x1 <- -5.75; roi_y1 <- 51.65
roi_x2 <- -5.3; roi_y2 <- 51.825
ndives <- length(res_all[,1])
lon_roi_egg <- lon_roi_chick <- lon_notroi_egg <- lon_notroi_chick <- c()
lat_roi_egg <- lat_roi_chick <- lat_notroi_egg <- lat_notroi_chick <- c()
chl_roi_egg <- chl_roi_chick <- chl_notroi_egg <- chl_notroi_chick <- c()
sst_roi_egg <- sst_roi_chick <- sst_notroi_egg <- sst_notroi_chick <- c()
x2 <- y2 <- rep(NA, ndives)
for(k in 1:ndives){
	print(paste("Now doing dive ",k," out of ",ndives," in ", res_all$name[k],sep=""))
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]]    
	mth <- as.numeric(substr(as.character(datetime),6,7)) #unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)
	
	# check stage first
	if(res_all$stage[k] == "egg"){
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI

		
		}
	}else{ # res_all$stage == "chick"
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
				#ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
				lon_roi_chick <- c(lon_roi_chick,as.numeric(ocean_arr[1,mylon,4]))
				lat_roi_chick <- c(lat_roi_chick,as.numeric(ocean_arr[mylat,1,3]))
				chl_roi_chick <- c(chl_roi_chick,as.numeric(ocean_arr[mylat,mylon,5]))
				sst_roi_chick <- c(sst_roi_chick,as.numeric(ocean_arr[mylat,mylon,6]))
				
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
			# MODIS script ends here
		
			# sab script starts over from here
			#k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}else{ # outside of ROI
		}
	}

}
save(lon_roi_chick, lat_roi_chick, chl_roi_chick, sst_roi_chick, file="roi_chick.RData")

# 4
roi_x1 <- -5.75; roi_y1 <- 51.65
roi_x2 <- -5.3; roi_y2 <- 51.825
ndives <- length(res_all[,1])
lon_roi_egg <- lon_roi_chick <- lon_notroi_egg <- lon_notroi_chick <- c()
lat_roi_egg <- lat_roi_chick <- lat_notroi_egg <- lat_notroi_chick <- c()
chl_roi_egg <- chl_roi_chick <- chl_notroi_egg <- chl_notroi_chick <- c()
sst_roi_egg <- sst_roi_chick <- sst_notroi_egg <- sst_notroi_chick <- c()
x2 <- y2 <- rep(NA, ndives)
for(k in 1:ndives){
	print(paste("Now doing dive ",k," out of ",ndives," in ", res_all$name[k],sep=""))
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]]    
	mth <- as.numeric(substr(as.character(datetime),6,7)) #unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)
	
	# check stage first
	if(res_all$stage[k] == "egg"){
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI
		}
	}else{ # res_all$stage == "chick"
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
				#ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
				lon_notroi_chick <- c(lon_notroi_chick,as.numeric(ocean_arr[1,mylon,4]))
				lat_notroi_chick <- c(lat_notroi_chick,as.numeric(ocean_arr[mylat,1,3]))
				chl_notroi_chick <- c(chl_notroi_chick,as.numeric(ocean_arr[mylat,mylon,5]))
				sst_notroi_chick <- c(sst_notroi_chick,as.numeric(ocean_arr[mylat,mylon,6]))
				
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
			# MODIS script ends here
		
			# sab script starts over from here
			#k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}
	}

}
save(lon_notroi_chick, lat_notroi_chick, chl_notroi_chick, sst_notroi_chick, file="notroi_chick.RData")



load("roi_egg.RData")
load("notroi_egg.RData")
load("roi_chick.RData")
load("notroi_chick.RData")

scale_chl=max(chl_notroi_egg, chl_roi_egg, chl_notroi_chick, chl_roi_chick,na.rm=T)
plot(lon_notroi_egg, lat_notroi_egg,col=rainbow(10))

xxx

pdf("chl_density.pdf",width=6,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(density(chl_notroi_egg,na.rm=T,bw=.125),col="red",lty=1,xlim=c(0,4),ylim=c(0,2.2),main="",xlab="Primary productivity (Chl OC5)")
lines(density(chl_roi_egg,na.rm=T,bw=.125),col="red",lty=2)
lines(density(chl_notroi_chick,na.rm=T,bw=.125),col="blue",lty=1)
lines(density(chl_roi_chick,na.rm=T,bw=.125),col="blue",lty=2)
legend(2,2,c("Incubation - not ROI","Incubation - ROI","Chick rearing - not ROI","Chick rearing - ROI"),col=c("red","red","blue","blue"),lty=c(1,2,1,2))
dev.off()

par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
cne_name <- rep("notroi_egg",length(chl_notroi_egg))
cre_name <- rep("roi_egg",length(chl_roi_egg))
cnc_name <- rep("notroi_chick",length(chl_notroi_chick))
crc_name <- rep("roi_chick",length(chl_roi_chick))
chl.df <- data.frame(cbind(c(chl_notroi_egg, chl_roi_egg, chl_notroi_chick, chl_roi_chick),c(cne_name, cre_name, cnc_name, crc_name)))
colnames(chl.df) <- c("Chl","Stage")
p1 <- qplot(factor(Stage),as.numeric(Chl), fill=factor(Stage), data= chl.df, geom="boxplot", position="dodge")+theme_bw()+scale_x_discrete(name="Stage")+scale_y_continuous(name="Chl") + theme(axis.text.x=element_text(angle=-90))
pdf("chl_boxplot.pdf")
print(grid.arrange(p1, ncol=1))
dev.off()


t.test(chl_notroi_egg, chl_roi_egg)
#data:  chl_notroi_egg and chl_roi_egg 
#t = -0.6346, df = 1591.58, p-value = 0.5258
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# -0.05127282  0.02620641 
#sample estimates:
#mean of x mean of y 
# 1.032088  1.044621 

t.test(chl_notroi_chick, chl_roi_chick)
#data:  chl_notroi_chick and chl_roi_chick 
#t = -10.4826, df = 889.258, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# -0.3199408 -0.2190306 
#sample estimates:
#mean of x mean of y 
#0.8688117 1.1382974 
t.test(chl_notroi_chick, chl_roi_chick,alternative="less")
#data:  chl_notroi_chick and chl_roi_chick 
#t = -10.4826, df = 889.258, p-value < 2.2e-16
#alternative hypothesis: true difference in means is less than 0 
#95 percent confidence interval:
#      -Inf -0.227156 
#sample estimates:
#mean of x mean of y 
#0.8688117 1.1382974 


t.test(chl_notroi_egg, chl_notroi_chick)
#data:  chl_notroi_egg and chl_notroi_chick 
#t = 7.6563, df = 1053.832, p-value = 4.316e-14
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# 0.1214302 0.2051218 
#sample estimates:
#mean of x mean of y 
#1.0320877 0.8688117 

t.test(chl_roi_egg, chl_roi_chick)
#data:  chl_roi_egg and chl_roi_chick 
#t = -3.8366, df = 910.182, p-value = 0.0001334
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# -0.14159638 -0.04575662 
#sample estimates:
#mean of x mean of y 
# 1.044621  1.138297 
t.test(chl_roi_egg, chl_roi_chick,alternative="less")


roi_x1 <- -5.75; roi_y1 <- 51.65
roi_x2 <- -5.3; roi_y2 <- 51.825
# defines sample space
min_lon <- min(c(lon_notroi_chick, lon_roi_chick, lon_notroi_egg, lon_roi_egg),na.rm=T)
max_lon <- max(c(lon_notroi_chick, lon_roi_chick, lon_notroi_egg, lon_roi_egg),na.rm=T)
min_lat <- min(c(lat_notroi_chick, lat_roi_chick, lat_notroi_egg, lat_roi_egg),na.rm=T)
max_lat <- max(c(lat_notroi_chick, lat_roi_chick, lat_notroi_egg, lat_roi_egg),na.rm=T)
# all data
all_lon <- c(lon_notroi_chick, lon_roi_chick, lon_notroi_egg, lon_roi_egg)
all_lat <- c(lat_notroi_chick, lat_roi_chick, lat_notroi_egg, lat_roi_egg)
all_chl <- c(chl_notroi_chick, chl_roi_chick, chl_notroi_egg, chl_roi_egg)
all_lon_chick <- c(lon_notroi_chick, lon_roi_chick)
all_lat_chick <- c(lat_notroi_chick, lat_roi_chick)
all_chl_chick <- c(chl_notroi_chick, chl_roi_chick)
all_lon_egg <- c(lon_notroi_egg, lon_roi_egg)
all_lat_egg <- c(lat_notroi_egg, lat_roi_egg)
all_chl_egg <- c(chl_notroi_egg, chl_roi_egg)


counter <- 1
Nreps <- 100000
t_stat_eggegg <- c()
while(counter < Nreps){
	# draw coordinate of top left corner of ROI uniformly thru positional fixes
	rand_lon <- runif(1, min_lon, max_lon)
	rand_lat <- runif(1, min_lat, max_lat)
	# contruct ROI
	rand_lon1 <- rand_lon
	rand_lat1 <- rand_lat
	rand_lon2 <- rand_lon1 + (roi_x2 - roi_x1)
	rand_lat2 <- rand_lat1 + (roi_y2 - roi_y1)
	# extract corresponding OC5 data
	rand_chl_notroi_chick <- all_chl_chick[((all_lon_chick < rand_lon1) | (all_lon_chick > rand_lon2)) & ((all_lat_chick < rand_lat1) | (all_lat_chick > rand_lat2))]
	rand_chl_roi_chick    <- all_chl_chick[(all_lon_chick > rand_lon1) & (all_lon_chick < rand_lon2) & (all_lat_chick > rand_lat1) & (all_lat_chick < rand_lat2)]
	rand_chl_notroi_egg   <- all_chl_egg[((all_lon_egg < rand_lon1) | (all_lon_egg > rand_lon2)) & ((all_lat_egg < rand_lat1) | (all_lat_egg > rand_lat2))]
	rand_chl_roi_egg      <- all_chl_egg[(all_lon_egg > rand_lon1) & (all_lon_egg < rand_lon2) & (all_lat_egg > rand_lat1) & (all_lat_egg < rand_lat2)]
	if( (sum(!is.na(rand_chl_notroi_egg)) > 5) & (sum(!is.na(rand_chl_roi_egg)) > 5)){
		t_stat_eggegg[counter] <- t.test(rand_chl_notroi_egg, rand_chl_roi_egg)$statistic[[1]]
		#t_stat_eggegg[counter] <- ks.test(rand_chl_notroi_egg, rand_chl_roi_egg)$statistic[[1]]
		counter <- counter + 1
	}
}
hist(t_stat_eggegg,20)
abline(v=t.test(chl_notroi_egg, chl_roi_egg)$statistic[[1]],col="red")

counter <- 1
Nreps <- 100000
t_stat_chickchick <- c()
while(counter < Nreps){
	# draw coordinate of top left corner of ROI uniformly thru positional fixes
	rand_lon <- runif(1, min_lon, max_lon)
	rand_lat <- runif(1, min_lat, max_lat)
	# contruct ROI
	rand_lon1 <- rand_lon
	rand_lat1 <- rand_lat
	rand_lon2 <- rand_lon1 + (roi_x2 - roi_x1)
	rand_lat2 <- rand_lat1 + (roi_y2 - roi_y1)
	# extract corresponding OC5 data
	rand_chl_notroi_chick <- all_chl_chick[((all_lon_chick < rand_lon1) | (all_lon_chick > rand_lon2)) & ((all_lat_chick < rand_lat1) | (all_lat_chick > rand_lat2))]
	rand_chl_roi_chick    <- all_chl_chick[(all_lon_chick > rand_lon1) & (all_lon_chick < rand_lon2) & (all_lat_chick > rand_lat1) & (all_lat_chick < rand_lat2)]
	rand_chl_notroi_chick   <- all_chl_chick[((all_lon_chick < rand_lon1) | (all_lon_chick > rand_lon2)) & ((all_lat_chick < rand_lat1) | (all_lat_chick > rand_lat2))]
	rand_chl_roi_chick      <- all_chl_chick[(all_lon_chick > rand_lon1) & (all_lon_chick < rand_lon2) & (all_lat_chick > rand_lat1) & (all_lat_chick < rand_lat2)]
	if( (sum(!is.na(rand_chl_notroi_chick)) > 5) & (sum(!is.na(rand_chl_roi_chick)) > 5)){
		t_stat_chickchick[counter] <- t.test(rand_chl_notroi_chick, rand_chl_roi_chick,alternative="less")$statistic[[1]]
		#t_stat_chickchick[counter] <- ks.test(rand_chl_notroi_chick, rand_chl_roi_chick,alternative="less")$statistic[[1]]
		counter <- counter + 1
	}
}
hist(t_stat_chickchick,20)
abline(v=t.test(chl_notroi_chick, chl_roi_chick,alternative="less")$statistic[[1]],col="red")

counter <- 1
Nreps <- 100000
t_stat_eggchick <- c()
while(counter < Nreps){
	# draw coordinate of top left corner of ROI uniformly thru positional fixes
	rand_lon <- runif(1, min_lon, max_lon)
	rand_lat <- runif(1, min_lat, max_lat)
	# contruct ROI
	rand_lon1 <- rand_lon
	rand_lat1 <- rand_lat
	rand_lon2 <- rand_lon1 + (roi_x2 - roi_x1)
	rand_lat2 <- rand_lat1 + (roi_y2 - roi_y1)
	# extract corresponding OC5 data
	rand_chl_notroi_chick <- all_chl_chick[((all_lon_chick < rand_lon1) | (all_lon_chick > rand_lon2)) & ((all_lat_chick < rand_lat1) | (all_lat_chick > rand_lat2))]
	rand_chl_roi_chick    <- all_chl_chick[(all_lon_chick > rand_lon1) & (all_lon_chick < rand_lon2) & (all_lat_chick > rand_lat1) & (all_lat_chick < rand_lat2)]
	rand_chl_notroi_chick   <- all_chl_chick[((all_lon_chick < rand_lon1) | (all_lon_chick > rand_lon2)) & ((all_lat_chick < rand_lat1) | (all_lat_chick > rand_lat2))]
	rand_chl_roi_chick      <- all_chl_chick[(all_lon_chick > rand_lon1) & (all_lon_chick < rand_lon2) & (all_lat_chick > rand_lat1) & (all_lat_chick < rand_lat2)]
	if( (sum(!is.na(rand_chl_roi_egg)) > 5) & (sum(!is.na(rand_chl_roi_chick)) > 5)){
		t_stat_eggchick[counter] <- t.test(rand_chl_roi_egg, rand_chl_roi_chick,alternative="less")$statistic[[1]]
		#t_stat_eggchick[counter] <- ks.test(rand_chl_roi_egg, rand_chl_roi_chick,alternative="less")$statistic[[1]]
		counter <- counter + 1
	}
}
hist(t_stat_eggchick,20)
abline(v=t.test(chl_roi_egg, chl_roi_chick,alternative="less")$statistic[[1]],col="red")

pdf("chl_boot_tstats.pdf",width=6,height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(3,1))
hist(t_stat_eggegg,50,xlab="Test statistic (t-test)",main="ROI vs. not ROI during incubation")
abline(v=t.test(chl_notroi_egg, chl_roi_egg)$statistic[[1]],col="red",lty=2)
#abline(v=ks.test(chl_notroi_egg, chl_roi_egg)$statistic[[1]],col="red",lty=2)
P1 <- sum(t_stat_eggegg < t.test(chl_notroi_egg, chl_roi_egg)$statistic[[1]])/Nreps
text(-60,14000,paste("P = ",P1,sep=""),adj = c(0,0))
hist(t_stat_chickchick,50,xlab="Test statistic (t-test)",main="ROI vs. not ROI during chick rearing")
abline(v=t.test(chl_notroi_chick, chl_roi_chick,alternative="less")$statistic[[1]],col="red",lty=2)
#abline(v=ks.test(chl_notroi_chick, chl_roi_chick,alternative="less")$statistic[[1]],col="red",lty=2)
P2 <- sum(t_stat_chickchick < t.test(chl_notroi_chick, chl_roi_chick,alternative="less")$statistic[[1]])/Nreps
text(-30,8000,paste("P = ",P2,sep=""),adj = c(0,0))
hist(t_stat_eggchick,50,xlab="Test statistic (t-test)",main="Incubation vs. chick rearing within ROI")
abline(v=t.test(chl_roi_egg, chl_roi_chick,alternative="less")$statistic[[1]],col="red",lty=2)
#abline(v=ks.test(chl_roi_egg, chl_roi_chick,alternative="less")$statistic[[1]],col="red",lty=2)
P3 <- sum(t_stat_eggchick < t.test(chl_roi_egg, chl_roi_chick,alternative="less")$statistic[[1]])/Nreps
text(-4,15000,paste("P = ",P3,sep=""),adj = c(0,0))
dev.off()

####################################################################
# alternative ROI #1
# 1

roi_x1 <- -5.75; roi_y1 <- 52.00
roi_x2 <- -5.3; roi_y2 <- 52.175

tmp_res_all <- res_all[complete.cases(res_all),]
tmp_coord_all <- coord_all[complete.cases(coord_all),]
pdf("activity_maps_colony_trimmed2_rect_alt1.pdf",width=10, height=10)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2))
# rest
bathy <- getNOAA.bathy(lon1 = min(tmp_res_all$ilon,na.rm=T), lon2 = max(tmp_res_all$ilon,na.rm=T), lat1 = min(tmp_res_all$ilat,na.rm=T), lat2 = max(tmp_res_all$ilat,na.rm=T),resolution = 1)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: fly/rest")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
rect(roi_x1, roi_y1, roi_x2, roi_y2, border ="red")
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: fly/rest")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
rect(roi_x1, roi_y1, roi_x2, roi_y2, border ="red")
# dive
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: dive")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
rect(roi_x1, roi_y1, roi_x2, roi_y2, border ="red")
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: dive")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
rect(roi_x1, roi_y1, roi_x2, roi_y2, border ="red")
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()
rm(z)


roi_x1 <- -5.75; roi_y1 <- 52.00
roi_x2 <- -5.3; roi_y2 <- 52.175
ndives <- length(res_all[,1])
lon_roi_egg <- lon_roi_chick <- lon_notroi_egg <- lon_notroi_chick <- c()
lat_roi_egg <- lat_roi_chick <- lat_notroi_egg <- lat_notroi_chick <- c()
chl_roi_egg <- chl_roi_chick <- chl_notroi_egg <- chl_notroi_chick <- c()
sst_roi_egg <- sst_roi_chick <- sst_notroi_egg <- sst_notroi_chick <- c()
x2 <- y2 <- rep(NA, ndives)
for(k in 1:ndives){
	print(paste("Now doing dive ",k," out of ",ndives," in ", res_all$name[k],sep=""))
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]]    
	mth <- as.numeric(substr(as.character(datetime),6,7)) #unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)
	
	# check stage first
	if(res_all$stage[k] == "egg"){
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     			   
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
			   #ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
			   lon_roi_egg <- c(lon_roi_egg,as.numeric(ocean_arr[1,mylon,4]))
			   lat_roi_egg <- c(lat_roi_egg,as.numeric(ocean_arr[mylat,1,3]))
			   chl_roi_egg <- c(chl_roi_egg,as.numeric(ocean_arr[mylat,mylon,5]))
			   sst_roi_egg <- c(sst_roi_egg,as.numeric(ocean_arr[mylat,mylon,6]))
			   
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
		   # MODIS script ends here
		
		   # sab script starts over from here
		   #k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}else{ # outside of ROI

		}
	}else{ # res_all$stage == "chick"
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI

		}
	}

}
save(lon_roi_egg, lat_roi_egg, chl_roi_egg, sst_roi_egg, file="roi_egg_alt1.RData")

# 2
roi_x1 <- -5.75; roi_y1 <- 52.00
roi_x2 <- -5.3; roi_y2 <- 52.175
ndives <- length(res_all[,1])
lon_roi_egg <- lon_roi_chick <- lon_notroi_egg <- lon_notroi_chick <- c()
lat_roi_egg <- lat_roi_chick <- lat_notroi_egg <- lat_notroi_chick <- c()
chl_roi_egg <- chl_roi_chick <- chl_notroi_egg <- chl_notroi_chick <- c()
sst_roi_egg <- sst_roi_chick <- sst_notroi_egg <- sst_notroi_chick <- c()
x2 <- y2 <- rep(NA, ndives)
for(k in 1:ndives){
	print(paste("Now doing dive ",k," out of ",ndives," in ", res_all$name[k],sep=""))
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]]    
	mth <- as.numeric(substr(as.character(datetime),6,7)) #unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)
	
	# check stage first
	if(res_all$stage[k] == "egg"){
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
				#ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
				lon_notroi_egg <- c(lon_notroi_egg,as.numeric(ocean_arr[1,mylon,4]))
				lat_notroi_egg <- c(lat_notroi_egg,as.numeric(ocean_arr[mylat,1,3]))
				chl_notroi_egg <- c(chl_notroi_egg,as.numeric(ocean_arr[mylat,mylon,5]))
				sst_notroi_egg <- c(sst_notroi_egg,as.numeric(ocean_arr[mylat,mylon,6]))
				
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
			# MODIS script ends here
		
			# sab script starts over from here
			#k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}
	}else{ # res_all$stage == "chick"
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI

		}
	}

}
save(lon_notroi_egg, lat_notroi_egg, chl_notroi_egg, sst_notroi_egg, file="notroi_egg_alt1.RData")

# 3
roi_x1 <- -5.75; roi_y1 <- 52.00
roi_x2 <- -5.3; roi_y2 <- 52.175
ndives <- length(res_all[,1])
lon_roi_egg <- lon_roi_chick <- lon_notroi_egg <- lon_notroi_chick <- c()
lat_roi_egg <- lat_roi_chick <- lat_notroi_egg <- lat_notroi_chick <- c()
chl_roi_egg <- chl_roi_chick <- chl_notroi_egg <- chl_notroi_chick <- c()
sst_roi_egg <- sst_roi_chick <- sst_notroi_egg <- sst_notroi_chick <- c()
x2 <- y2 <- rep(NA, ndives)
for(k in 1:ndives){
	print(paste("Now doing dive ",k," out of ",ndives," in ", res_all$name[k],sep=""))
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]]    
	mth <- as.numeric(substr(as.character(datetime),6,7)) #unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)
	
	# check stage first
	if(res_all$stage[k] == "egg"){
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI

		
		}
	}else{ # res_all$stage == "chick"
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
				#ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
				lon_roi_chick <- c(lon_roi_chick,as.numeric(ocean_arr[1,mylon,4]))
				lat_roi_chick <- c(lat_roi_chick,as.numeric(ocean_arr[mylat,1,3]))
				chl_roi_chick <- c(chl_roi_chick,as.numeric(ocean_arr[mylat,mylon,5]))
				sst_roi_chick <- c(sst_roi_chick,as.numeric(ocean_arr[mylat,mylon,6]))
				
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
			# MODIS script ends here
		
			# sab script starts over from here
			#k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}else{ # outside of ROI
		}
	}

}
save(lon_roi_chick, lat_roi_chick, chl_roi_chick, sst_roi_chick, file="roi_chick_alt1.RData")

# 4
roi_x1 <- -5.75; roi_y1 <- 52.00
roi_x2 <- -5.3; roi_y2 <- 52.175
ndives <- length(res_all[,1])
lon_roi_egg <- lon_roi_chick <- lon_notroi_egg <- lon_notroi_chick <- c()
lat_roi_egg <- lat_roi_chick <- lat_notroi_egg <- lat_notroi_chick <- c()
chl_roi_egg <- chl_roi_chick <- chl_notroi_egg <- chl_notroi_chick <- c()
sst_roi_egg <- sst_roi_chick <- sst_notroi_egg <- sst_notroi_chick <- c()
x2 <- y2 <- rep(NA, ndives)
for(k in 1:ndives){
	print(paste("Now doing dive ",k," out of ",ndives," in ", res_all$name[k],sep=""))
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]]    
	mth <- as.numeric(substr(as.character(datetime),6,7)) #unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)
	
	# check stage first
	if(res_all$stage[k] == "egg"){
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI
		}
	}else{ # res_all$stage == "chick"
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
				#ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
				lon_notroi_chick <- c(lon_notroi_chick,as.numeric(ocean_arr[1,mylon,4]))
				lat_notroi_chick <- c(lat_notroi_chick,as.numeric(ocean_arr[mylat,1,3]))
				chl_notroi_chick <- c(chl_notroi_chick,as.numeric(ocean_arr[mylat,mylon,5]))
				sst_notroi_chick <- c(sst_notroi_chick,as.numeric(ocean_arr[mylat,mylon,6]))
				
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
			# MODIS script ends here
		
			# sab script starts over from here
			#k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}
	}

}
save(lon_notroi_chick, lat_notroi_chick, chl_notroi_chick, sst_notroi_chick, file="notroi_chick_alt1.RData")


####################################################################
# alternative ROI #2
# 1

roi_x1 <- -5.75; roi_y1 <- 51.30
roi_x2 <- -5.3; roi_y2 <- 51.475

tmp_res_all <- res_all[complete.cases(res_all),]
tmp_coord_all <- coord_all[complete.cases(coord_all),]
pdf("activity_maps_colony_trimmed2_rect_alt2.pdf",width=10, height=10)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2))
# rest
bathy <- getNOAA.bathy(lon1 = min(tmp_res_all$ilon,na.rm=T), lon2 = max(tmp_res_all$ilon,na.rm=T), lat1 = min(tmp_res_all$ilat,na.rm=T), lat2 = max(tmp_res_all$ilat,na.rm=T),resolution = 1)
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: fly/rest")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
rect(roi_x1, roi_y1, roi_x2, roi_y2, border ="red")
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: fly/rest")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & (tmp_coord_all[,12]=="rest" | tmp_coord_all[,12]=="fly")]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
rect(roi_x1, roi_y1, roi_x2, roi_y2, border ="red")
# dive
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"], tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Incubation: dive")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="egg" & tmp_coord_all[,12]=="dive"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
rect(roi_x1, roi_y1, roi_x2, roi_y2, border ="red")
smoothScatter(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"], tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"],xlim=c(min(tmp_res_all$ilon,na.rm=T),max(tmp_res_all$ilon,na.rm=T)),ylim=c(min(tmp_res_all$ilat,na.rm=T),max(tmp_res_all$ilat,na.rm=T)),colramp = colorRampPalette(c("white", "blue")),xlab="Longitude",ylab="Latitude",main="Chick rearing: dive")
z <- kde2d(as.numeric(tmp_coord_all[,3][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), as.numeric(tmp_coord_all[,4][tmp_coord_all[,10]=="chick" & tmp_coord_all[,12]=="dive"]), n=50)
plot(bathy,add=T,col="gray")
contour(z, drawlabels=T, nlevels=8, add=TRUE)
map('worldHires',c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'),fill=T,col="gray",add=T)
points(col.loc1[1],col.loc1[2],pch=7); points(col.loc2[1],col.loc2[2],pch=9)
rect(roi_x1, roi_y1, roi_x2, roi_y2, border ="red")
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()
rm(z)


roi_x1 <- -5.75; roi_y1 <- 51.30
roi_x2 <- -5.3; roi_y2 <- 51.475
ndives <- length(res_all[,1])
lon_roi_egg <- lon_roi_chick <- lon_notroi_egg <- lon_notroi_chick <- c()
lat_roi_egg <- lat_roi_chick <- lat_notroi_egg <- lat_notroi_chick <- c()
chl_roi_egg <- chl_roi_chick <- chl_notroi_egg <- chl_notroi_chick <- c()
sst_roi_egg <- sst_roi_chick <- sst_notroi_egg <- sst_notroi_chick <- c()
x2 <- y2 <- rep(NA, ndives)
for(k in 1:ndives){
	print(paste("Now doing dive ",k," out of ",ndives," in ", res_all$name[k],sep=""))
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]]    
	mth <- as.numeric(substr(as.character(datetime),6,7)) #unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)
	
	# check stage first
	if(res_all$stage[k] == "egg"){
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     			   
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
			   #ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
			   lon_roi_egg <- c(lon_roi_egg,as.numeric(ocean_arr[1,mylon,4]))
			   lat_roi_egg <- c(lat_roi_egg,as.numeric(ocean_arr[mylat,1,3]))
			   chl_roi_egg <- c(chl_roi_egg,as.numeric(ocean_arr[mylat,mylon,5]))
			   sst_roi_egg <- c(sst_roi_egg,as.numeric(ocean_arr[mylat,mylon,6]))
			   
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
		   # MODIS script ends here
		
		   # sab script starts over from here
		   #k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}else{ # outside of ROI

		}
	}else{ # res_all$stage == "chick"
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI

		}
	}

}
save(lon_roi_egg, lat_roi_egg, chl_roi_egg, sst_roi_egg, file="roi_egg_alt2.RData")

# 2
roi_x1 <- -5.75; roi_y1 <- 51.30
roi_x2 <- -5.3; roi_y2 <- 51.475
ndives <- length(res_all[,1])
lon_roi_egg <- lon_roi_chick <- lon_notroi_egg <- lon_notroi_chick <- c()
lat_roi_egg <- lat_roi_chick <- lat_notroi_egg <- lat_notroi_chick <- c()
chl_roi_egg <- chl_roi_chick <- chl_notroi_egg <- chl_notroi_chick <- c()
sst_roi_egg <- sst_roi_chick <- sst_notroi_egg <- sst_notroi_chick <- c()
x2 <- y2 <- rep(NA, ndives)
for(k in 1:ndives){
	print(paste("Now doing dive ",k," out of ",ndives," in ", res_all$name[k],sep=""))
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]]    
	mth <- as.numeric(substr(as.character(datetime),6,7)) #unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)
	
	# check stage first
	if(res_all$stage[k] == "egg"){
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
				#ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
				lon_notroi_egg <- c(lon_notroi_egg,as.numeric(ocean_arr[1,mylon,4]))
				lat_notroi_egg <- c(lat_notroi_egg,as.numeric(ocean_arr[mylat,1,3]))
				chl_notroi_egg <- c(chl_notroi_egg,as.numeric(ocean_arr[mylat,mylon,5]))
				sst_notroi_egg <- c(sst_notroi_egg,as.numeric(ocean_arr[mylat,mylon,6]))
				
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
			# MODIS script ends here
		
			# sab script starts over from here
			#k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}
	}else{ # res_all$stage == "chick"
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI

		}
	}

}
save(lon_notroi_egg, lat_notroi_egg, chl_notroi_egg, sst_notroi_egg, file="notroi_egg_alt2.RData")

# 3
roi_x1 <- -5.75; roi_y1 <- 51.30
roi_x2 <- -5.3; roi_y2 <- 51.475
ndives <- length(res_all[,1])
lon_roi_egg <- lon_roi_chick <- lon_notroi_egg <- lon_notroi_chick <- c()
lat_roi_egg <- lat_roi_chick <- lat_notroi_egg <- lat_notroi_chick <- c()
chl_roi_egg <- chl_roi_chick <- chl_notroi_egg <- chl_notroi_chick <- c()
sst_roi_egg <- sst_roi_chick <- sst_notroi_egg <- sst_notroi_chick <- c()
x2 <- y2 <- rep(NA, ndives)
for(k in 1:ndives){
	print(paste("Now doing dive ",k," out of ",ndives," in ", res_all$name[k],sep=""))
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]]    
	mth <- as.numeric(substr(as.character(datetime),6,7)) #unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)
	
	# check stage first
	if(res_all$stage[k] == "egg"){
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI

		
		}
	}else{ # res_all$stage == "chick"
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
				#ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
				lon_roi_chick <- c(lon_roi_chick,as.numeric(ocean_arr[1,mylon,4]))
				lat_roi_chick <- c(lat_roi_chick,as.numeric(ocean_arr[mylat,1,3]))
				chl_roi_chick <- c(chl_roi_chick,as.numeric(ocean_arr[mylat,mylon,5]))
				sst_roi_chick <- c(sst_roi_chick,as.numeric(ocean_arr[mylat,mylon,6]))
				
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
			# MODIS script ends here
		
			# sab script starts over from here
			#k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}else{ # outside of ROI
		}
	}

}
save(lon_roi_chick, lat_roi_chick, chl_roi_chick, sst_roi_chick, file="roi_chick_alt2.RData")

# 4
roi_x1 <- -5.75; roi_y1 <- 51.30
roi_x2 <- -5.3; roi_y2 <- 51.475
ndives <- length(res_all[,1])
lon_roi_egg <- lon_roi_chick <- lon_notroi_egg <- lon_notroi_chick <- c()
lat_roi_egg <- lat_roi_chick <- lat_notroi_egg <- lat_notroi_chick <- c()
chl_roi_egg <- chl_roi_chick <- chl_notroi_egg <- chl_notroi_chick <- c()
sst_roi_egg <- sst_roi_chick <- sst_notroi_egg <- sst_notroi_chick <- c()
x2 <- y2 <- rep(NA, ndives)
for(k in 1:ndives){
	print(paste("Now doing dive ",k," out of ",ndives," in ", res_all$name[k],sep=""))
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]]    
	mth <- as.numeric(substr(as.character(datetime),6,7)) #unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)
	
	# check stage first
	if(res_all$stage[k] == "egg"){
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI
		}
	}else{ # res_all$stage == "chick"
		# then check location
		if((res_all$ilon[k] >= roi_x1) & (res_all$ilon[k] <= roi_x2) & (res_all$ilat[k] >= roi_y1) & (res_all$ilat[k] <= roi_y2)){

		}else{ # outside of ROI
		   # MODIS script	
		   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
		   sst_dir      = paste("/Users/stephane/Documents/data/akiko/131004/130305_bathymetry/AVHRR_1km_SST/",year, sep="")
		   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
		   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
		   if (length(file_list) == 0) next
		   for(i in 1:length(file_list)){
			    jul_date <- substr(file_list[i], 6, 8)
			    time <- substr(file_list[i], 10, 13)
			    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
			    sst_file  <- paste("M",year,jul_date,"-",year,jul_date,".*.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.*.data.nc", sep="")
			    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
			    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
			    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
			    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
			    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
			    file_paths <- c(chl_oc5_filepath, sst_filepath)
			    variable   <- c("chl_oc5","sstp")
			    # Only read in files if data is available for ALL variables
			    #if (length(file_paths)==3) {
			    #if (length(file_paths)==2) {
			    if (length(file_paths)>0) {
			     count = count + 1
			     ex.nc <- open.ncdf(file_paths[1])
			     numRows<-length(get.var.ncdf(ex.nc,"latitude"))
			     numColumns<-length(get.var.ncdf(ex.nc,"longitude"))
			     # Define empty matrix to hold all variables
			     variable_matrix  <- array(0,dim=c(numRows,numColumns,4))
			     #for(j in 1:3){
			     for(j in 1:2){
	                  # Open netcdf file
	                  ex.nc = open.ncdf(file_paths[j])
	                  # Extract latitude and longitude
	                  long = get.var.ncdf(ex.nc,"longitude")
	                  lat = get.var.ncdf(ex.nc,"latitude")
	                  # Retrieve variable:
	                  # sab
	                  #z  = get.var.ncdf(ex.nc, variable[j])
	                  if(j==1){
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }else{
	                  	z  = get.var.ncdf(ex.nc, variable[j])
	                  }
	                  # variables get read in by get.var.ncdf on their side for 
	                  # some reason, so transpose
	                  z <- t(z)
	                  # Save to variable matrix
	                  variable_matrix[,,j] <- z
	                  # If reading the modis wind netcdf, extract the second variable
	                  if (j==3) {
	                     z2 <- get.var.ncdf(ex.nc, variable[j+1])
	                     # transpose
	                     z2<-t(z2)
	                     # Save to variable matrix
	                     variable_matrix[,,j+1] <- z2
	                  } # end if
	               } # end for (j loop)
			    } else next 
               # Replicate longitude over rows, then turn on side to fit everything else
               long <- matrix(rep(long,numRows),nrow=numColumns) 
               long <-t (long)
               # Replicate latitude over columns
               lat <- matrix(rep(lat,numColumns),nrow=numRows) 
               # Create date array
               date <- paste(day, month, year,sep="")
               date_vec <- rep(date, numRows)
               date_arr <- cbind( do.call(cbind, rep(list(date_vec), numColumns)) )
               # Create time array
               time_vec <- rep(time, numRows)
               time_arr <- cbind( do.call(cbind, rep(list(time_vec), numColumns)) )
               # Concatinate into one large array: 
               # Date (ddmmyyy), time (hhmm), Latitude, Longitude, chl_oc5, sst, windspeed, windangle
               ocean_arr <- array(0,dim=c(numRows,numColumns,8))
               ocean_arr[,,1] <- date_arr
               ocean_arr[,,2] <- time_arr
               ocean_arr[,,3] <- lat
               ocean_arr[,,4] <- long
               ocean_arr[,,5:8] <- variable_matrix
			     
				mylat <- which(as.numeric(ocean_arr[,1,3])<=res_all$ilat[k])[1]
				mylon <- which(as.numeric(t(ocean_arr[1,,4]))>=res_all$ilon[k])[1]
				#ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
				lon_notroi_chick <- c(lon_notroi_chick,as.numeric(ocean_arr[1,mylon,4]))
				lat_notroi_chick <- c(lat_notroi_chick,as.numeric(ocean_arr[mylat,1,3]))
				chl_notroi_chick <- c(chl_notroi_chick,as.numeric(ocean_arr[mylat,mylon,5]))
				sst_notroi_chick <- c(sst_notroi_chick,as.numeric(ocean_arr[mylat,mylon,6]))
				
			   # Assign 'ocean_arr' for different times of the same day with a unique name
			   #assign(paste("ocean_arr_",count,sep=""), ocean_arr)
			   # Create list of arrays
			   print(c(paste(year,"/", month,"/", day," @ ",time," lon=",format(as.numeric(ocean_arr[1,mylon,4]),digits=4)," lat=",format(as.numeric(ocean_arr[mylat,1,3]),digits=4)," chl=",format(as.numeric(ocean_arr[mylat,mylon,5]),digits=4)," sst=",format(as.numeric(ocean_arr[mylat,mylon,6]),digits=4),sep=""))   )
			   rm(ocean_arr, time_arr, time_vec, lat, long, variable_matrix, z, ex.nc)
		   } # i for loop
		   # Reset counter
		   count = 0
			# MODIS script ends here
		
			# sab script starts over from here
			#k<-465;res_all$ilon[k];res_all$ilat[k]
		
		}
	}

}
save(lon_notroi_chick, lat_notroi_chick, chl_notroi_chick, sst_notroi_chick, file="notroi_chick_alt2.RData")


load("roi_egg_alt1.RData")
load("notroi_egg_alt1.RData")
load("roi_chick_alt1.RData")
load("notroi_chick_alt1.RData")

pdf("chl_density_alt1.pdf",width=6,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(density(chl_notroi_egg,na.rm=T,bw=.125),col="red",lty=1,xlim=c(0,4),ylim=c(0,2.2),main="",xlab="Primary productivity (Chl OC5)")
lines(density(chl_roi_egg,na.rm=T,bw=.125),col="red",lty=2)
lines(density(chl_notroi_chick,na.rm=T,bw=.125),col="blue",lty=1)
lines(density(chl_roi_chick,na.rm=T,bw=.125),col="blue",lty=2)
legend(2,2,c("Incubation - not ROI","Incubation - ROI","Chick rearing - not ROI","Chick rearing - ROI"),col=c("red","red","blue","blue"),lty=c(1,2,1,2))
dev.off()

t.test(chl_notroi_egg, chl_roi_egg)
#data:  chl_notroi_egg and chl_roi_egg 
#t = -3.6467, df = 928.519, p-value = 0.0002804
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# -0.10816993 -0.03247828 
#sample estimates:
#mean of x mean of y 
# 1.064113  1.134437 

t.test(chl_notroi_chick, chl_roi_chick)
#data:  chl_notroi_chick and chl_roi_chick 
#t = 15.9506, df = 253.485, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# 0.2517060 0.3226157 
#sample estimates:
#mean of x mean of y 
#1.1562226 0.8690618 


t.test(chl_notroi_egg, chl_notroi_chick)
#data:  chl_notroi_egg and chl_notroi_chick 
#t = -4.1762, df = 1491.759, p-value = 3.136e-05
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# -0.13537321 -0.04884533 
#sample estimates:
#mean of x mean of y 
# 1.064113  1.156223 

t.test(chl_roi_egg, chl_roi_chick)
#data:  chl_roi_egg and chl_roi_chick 
#t = 18.333, df = 108.478, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# 0.2366846 0.2940667 
#sample estimates:
#mean of x mean of y 
#1.1344374 0.8690618 


load("roi_egg_alt2.RData")
load("notroi_egg_alt2.RData")
load("roi_chick_alt2.RData")
load("notroi_chick_alt2.RData")

pdf("chl_density_alt2.pdf",width=6,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(density(chl_notroi_egg,na.rm=T,bw=.125),col="red",lty=1,xlim=c(0,4),ylim=c(0,2.2),main="",xlab="Primary productivity (Chl OC5)")
lines(density(chl_roi_egg,na.rm=T,bw=.125),col="red",lty=2)
lines(density(chl_notroi_chick,na.rm=T,bw=.125),col="blue",lty=1)
lines(density(chl_roi_chick,na.rm=T,bw=.125),col="blue",lty=2)
legend(2,2,c("Incubation - not ROI","Incubation - ROI","Chick rearing - not ROI","Chick rearing - ROI"),col=c("red","red","blue","blue"),lty=c(1,2,1,2))
dev.off()

t.test(chl_notroi_egg, chl_roi_egg)
#data:  chl_notroi_egg and chl_roi_egg 
#t = 0.5411, df = 157.585, p-value = 0.5892
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# -0.06710149  0.11773509 
#sample estimates:
#mean of x mean of y 
# 1.054614  1.029297 

t.test(chl_notroi_chick, chl_roi_chick)
#data:  chl_notroi_chick and chl_roi_chick 
#t = -5.9858, df = 608.268, p-value = 3.681e-09
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# -0.15810787 -0.07999024 
#sample estimates:
#mean of x mean of y 
# 1.145846  1.264895 


t.test(chl_notroi_egg, chl_notroi_chick)
#data:  chl_notroi_egg and chl_notroi_chick 
#t = -3.9434, df = 1246.581, p-value = 8.48e-05
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# -0.13662045 -0.04584353 
#sample estimates:
#mean of x mean of y 
# 1.054614  1.145846 

t.test(chl_roi_egg, chl_roi_chick)
#data:  chl_roi_egg and chl_roi_chick 
#t = -5.2039, df = 137.994, p-value = 6.916e-07
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# -0.3251177 -0.1460780 
#sample estimates:
#mean of x mean of y 
# 1.029297  1.264895 



####################################################################
####################################################################

# resolution
(as.numeric(ocean_arr[512,512,3]) - as.numeric(ocean_arr[1,1,3]))/512 # lon
(as.numeric(ocean_arr[512,512,4]) - as.numeric(ocean_arr[1,1,4]))/512 # lat

pdf("ipq_chl.pdf",width=6,height=4)
smoothScatter(chl[chl>0],res_all$ipq[chl>0],xlab="Primary productivity [chl_oc5]",ylab="log IPQ")
lm.chl <- lm(res_all$ipq[chl>0]~chl[chl>0])
summary(lm.chl)
abline(lm.chl,col="red")
dev.off()

pdf("ipq_sst.pdf",width=6,height=4)
smoothScatter(sst[sst>0],res_all$ipq[sst>0],xlab="Sea surface temperature",ylab="log IPQ")
lm.sst <- lm(res_all$ipq[sst>0]~sst[sst>0])
summary(lm.sst)
abline(lm.sst,col="red")
dev.off()

# sst prevvious day
sst_prev <- rep(NA, ndives)
x2 <- y2 <- rep(NA, ndives)
count <- 0
for(k in 1:ndives){
	datetime <- unix2POSIXct(res_all$begdesc[k]) # unix2POSIXct(1307593743)
	dy <- unlist(datetime)[[4]] - 1
	mth <- unlist(datetime)[[5]]
	year <- as.numeric(paste("2",unlist(datetime)[[6]],sep="")) - 100
	year_short <- substr(year,3,4)
	month <- sprintf("%02d",mth)
	day <- sprintf("%02d",dy)

	# MODIS script	
   chl_oc5_dir  = paste("/Users/stephane/Documents/data/akiko/130305_bathymetry/MODIS_1km_OC5/",year,"/",month,"/",day,"/", sep="")
   sst_dir      = paste("/Users/stephane/Documents/data/akiko/130305_bathymetry/AVHRR_1km_SST/",year,"/AVHRR_1km_",month,"_",year,".pa", sep="")
   #wind_dir     = paste("",year,"/",month,"/",day,"/", sep="")
   file_list <- list.files(path = chl_oc5_dir, pattern = ".nc", all.files = FALSE, full.names = FALSE, recursive = FALSE,  ignore.case = FALSE, include.dirs = FALSE)
   month <- sprintf("%02d",mth)
   mth_arr <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
   if (length(file_list) == 0) next
   for(i in 1:length(file_list)){
	    jul_date <- substr(file_list[i], 6, 8)
	    time <- substr(file_list[i], 10, 13)
	    chl_oc5_filepath = paste(chl_oc5_dir, file_list[i], sep="") 
	    sst_file      = paste("M",year,jul_date,"-",year,jul_date,".pa.sstp.AVH.L3_median.",day,mth_arr[mth],year_short,"-",day,mth_arr[mth],year_short,".v1.",year,jul_date,"*.data.nc", sep="")
	    sst_filepath  <- Sys.glob(paste(sst_dir,"/", sst_file, sep=""))
	    #wind_file     = paste("M",year,jul_date,".",time,".pa.all_products.MYD.",day,mth_arr[mth],year_short,time,".v1.*.data.nc", sep="")  
	    #wind_filepath <- Sys.glob(paste(wind_dir,wind_file,sep=""))
	    #file_paths <- c(chl_oc5_filepath, sst_filepath, wind_filepath)
	    #variable   <- c("chl_oc5","sstp","windspeed","windangle")
	    file_paths <- c(chl_oc5_filepath, sst_filepath)
	    variable   <- c("chl_oc5","sstp")
	    # Only read in files if data is available for ALL variables
	    #if (length(file_paths)==3) {
	    #if (length(file_paths)==2) {
	    if (length(file_paths)>0) {
	     count = count + 1
	     z_mat  <- array(0,dim=c(512,512,4))
	     #for(j in 1:3){
	     for(j in 1:2){
	     	if(!is.na(file_paths[j])){
		      ex.nc = open.ncdf(file_paths[j])
		      x = get.var.ncdf( ex.nc, "longitude")
		      y = get.var.ncdf( ex.nc, "latitude")
		      z  = get.var.ncdf( ex.nc, variable[j])
		      z <- z[, rev(seq_len(ncol(z)))]
		      z_mat[,,j] <- z
		      # If reading the modis wind netcdf, extract the second variable
		      #if (j==3) {
		      #  z2 = get.var.ncdf( ex.nc, variable[j+1])
		      #
		      #  # Reverse latitudes (order of columns)
		      #  z2 <- z2[, rev(seq_len(ncol(z2)))]
		      #
		      #  # Save to variable matrix
		      #  z_mat[,,j+1] <- z2
		      #}
		      close.ncdf(ex.nc)
	     	}
	     } # j for loop
	    } else next 
	     lon <- cbind( do.call(cbind, rep(list(x), 512)) ) 
	     y <- rev(y)
	     lat <- cbind( do.call(cbind, rep(list(y), 512)) )
	     lat <- t(lat)
	     date <- paste(day, month, year,sep="")
	     date_vec <- rep(date, 512)
	     date_arr <- cbind( do.call(cbind, rep(list(date_vec), 512)) )
	     time_vec <- rep(time, 512)
	     time_arr <- cbind( do.call(cbind, rep(list(time_vec), 512)) )
	     ocean_arr <- array(0,dim=c(512,512,8))
	     ocean_arr[,,1] <- date_arr
	     ocean_arr[,,2] <- time_arr
	     ocean_arr[,,3] <- lat
	     ocean_arr[,,4] <- lon
	     ocean_arr[,,5:8] <- z_mat
	   # Assign 'ocean_arr' for different times of the same day with a unique name
	   assign(paste("ocean_arr_",count,sep=""), ocean_arr)
	   # Create list of arrays
	   print(c(paste("ocean_arr_",count,sep=""), time, day, month, year))
   } # i for loop
   # Reset counter
   count = 0
	# MODIS script ends here

	# sab script starts over from here
	#k<-465;res_all$ilon[k];res_all$ilat[k]

	mylat <- which(ocean_arr[1,,3]>=res_all$ilat[k])[1]
	mylon <- which(as.numeric(t(ocean_arr[,1,4]))>=res_all$ilon[k])[1]
	#ocean_arr[mylon, mylat,3];ocean_arr[mylon, mylat,4]
	sst_prev[k] <- as.numeric(ocean_arr[mylon,mylat,6])
}

sst_diff <- sst - sst_prev

pdf("ipq_sst_diff.pdf",width=6,height=4)
smoothScatter(sst_diff[sst>0 & sst_prev>0],res_all$ipq[sst>0 & sst_prev>0],xlab="Sea surface temperature difference",ylab="log IPQ")
lm.sst_diff <- lm(res_all$ipq[sst>0 & sst_prev>0]~sst_diff[sst>0 & sst_prev>0])
summary(lm.sst_diff)
abline(lm.sst_diff,col="red")
dev.off()

####################################################################
####################################################################
# from Akiko to plot log IPQ vs. dist to colony with quadratic fits

par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,1))
#2011
plot(as.numeric(res_all[,4][res_all$year == "2011"])~as.numeric(res_all[,5][res_all$year == "2011"]))
x <- as.numeric(res_all[,5][res_all$year == "2011"])
y <- as.numeric(res_all[,4][res_all$year == "2011"])
lm2_tot1 <- lm(y ~ poly(x,2)); summary(lm2_tot1)
xv <- seq(min(x), max(x), 10)
yv <- predict(lm2_tot1, list(x=xv))
points(xv, yv, col = "gray")
#2012
plot(as.numeric(res_all[,4][res_all$year == "2012"])~as.numeric(res_all[,5][res_all$year == "2012"]))
x <- as.numeric(res_all[,5][res_all$year == "2012"])
y <- as.numeric(res_all[,4][res_all$year == "2012"])
lm2_tot1 <- lm(y ~ poly(x,2)); summary(lm2_tot1)
xv <- seq(min(x), max(x), 10)
yv <- predict(lm2_tot1, list(x=xv))
points(xv, yv)


####################################################################
####################################################################
# circadian IPQ variation
# http://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data

library(timeSeries)
library(Hmisc)

substr(unix2POSIXct(res_all$begdesc[1]),12,16) # hour of day + minutes
as.duration(hms(substr(unix2POSIXct(res_all$begdesc[1]),12,19)))
plot(as.duration(hms(substr(unix2POSIXct(res_all$begdesc),12,19))), res_all$ipq)

ipq_mean_hour <- rep(0,24)
ipq_sd_hour <- rep(0,24)
ipq_count_hour <- rep(0,24)
for(i in 1:length(res_all$ipq)){
	ipq_mean_hour[as.numeric(substr(unix2POSIXct(res_all$begdesc[i]),12,13))] <- ipq_mean_hour[as.numeric(substr(unix2POSIXct(res_all$begdesc[i]),12,13))] + res_all$ipq[i]
	ipq_sd_hour[as.numeric(substr(unix2POSIXct(res_all$begdesc[i]),12,13))] <- ipq_sd_hour[as.numeric(substr(unix2POSIXct(res_all$begdesc[i]),12,13))] + (res_all$ipq[i])^2
	ipq_count_hour[as.numeric(substr(unix2POSIXct(res_all$begdesc[i]),12,13))] <- ipq_count_hour[as.numeric(substr(unix2POSIXct(res_all$begdesc[i]),12,13))] + 1
}

ipq_mean_hour <- ipq_mean_hour / ipq_count_hour
ipq_sd_hour <- sqrt((ipq_sd_hour / ipq_count_hour) - (ipq_mean_hour)^2)
ipq_sd_hour <- ipq_sd_hour / sqrt(ipq_count_hour)

plot(ipq_mean_hour)
for(i in 1:24){
	segments(i, ipq_mean_hour[i]-ipq_sd_hour[i],i, ipq_mean_hour[i]+ipq_sd_hour[i])
}

ipq2days <- c(ipq_mean_hour, ipq_mean_hour, ipq_mean_hour)
ipq2days <- c(ipq_mean_hour)
plot(ipq2days, type="l",xlab="Hours of days",ylab="log IPQ")

ts_table <- ts(ipq2days, frequency=24, start=0)
ts_table[which(ts_table=="NaN")] <- NA
ts_table <- ts(interpNA(as.timeSeries(ts_table), method = "linear"), frequency=24, start=0) # interpolate missing values

y <- ts_table[2:69]
t <- 2:69
ssp <- spectrum(y)  
per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
reslm <- lm(y ~ sin(2*pi/per*t)+cos(2*pi/per*t))
summary(reslm)
#Call:
#lm(formula = y ~ sin(2 * pi/per * t) + cos(2 * pi/per * t))
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.15757 -0.02773  0.01691  0.04315  0.05931 
#
#Coefficients:
#                     Estimate Std. Error  t value Pr(>|t|)    
#(Intercept)         -1.987609   0.006672 -297.924   <2e-16 ***
#sin(2 * pi/per * t)  0.023858   0.009189    2.596   0.0116 *  
#cos(2 * pi/per * t) -0.118647   0.009668  -12.272   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
#
#Residual standard error: 0.05483 on 65 degrees of freedom
#Multiple R-squared: 0.7066,	Adjusted R-squared: 0.6975 
#F-statistic: 78.26 on 2 and 65 DF,  p-value: < 2.2e-16 



# 140625
# counting non 2013 birds
nbirds_ID <- unique(res_all$ID[res_all$year != 2013])
nbirds <- length(nbirds_ID)
# computing hourly averages per bird
ipq_mean_hour_pb <- rep(0,24* nbirds)
ipq_sd_hour_pb <- rep(0,24* nbirds)
ipq_count_hour_pb <- rep(0,24* nbirds)
hours <- rep(seq(1:24),nbirds)
ipq_birdID <- rep(0,24* nbirds)
for(j in 1:nbirds){
	for(i in 1:length(res_all$ipq)){
		if(res_all$ID[i] == nbirds_ID[j]){
			ipq_mean_hour_pb[j * as.numeric(substr(unix2POSIXct(res_all$begdesc[i]),12,13))] <- ipq_mean_hour_pb[j * as.numeric(substr(unix2POSIXct(res_all$begdesc[i]),12,13))] + res_all$ipq[i]
			ipq_sd_hour_pb[j * as.numeric(substr(unix2POSIXct(res_all$begdesc[i]),12,13))] <- ipq_sd_hour_pb[j * as.numeric(substr(unix2POSIXct(res_all$begdesc[i]),12,13))] + (res_all$ipq[i])^2
			ipq_count_hour_pb[j * as.numeric(substr(unix2POSIXct(res_all$begdesc[i]),12,13))] <- ipq_count_hour_pb[j * as.numeric(substr(unix2POSIXct(res_all$begdesc[i]),12,13))] + 1
			ipq_birdID[j * as.numeric(substr(unix2POSIXct(res_all$begdesc[i]),12,13))] <- nbirds_ID[j]
		}
	}
}
ipq_mean_hour_pb <- ipq_mean_hour_pb / ipq_count_hour_pb
ipq_sd_hour_pb <- sqrt((ipq_sd_hour / ipq_count_hour) - (ipq_mean_hour)^2)
ipq_sd_hour_pb <- ipq_sd_hour / sqrt(ipq_count_hour)
plot(ipq_mean_hour_pb ~ hours)

reslme <- lme(y ~ sin(2*pi/per*t)+cos(2*pi/per*t), random=~1|)



rg <- diff(range(y))
t <- t -24
pdf("ipq_circadian.pdf",width=6,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(y~t,ylim=c(min(y)-0.1*rg,max(y)+0.1*rg),pch=20,xlim=c(0,24),,xlab="Hours of day",ylab="log IPQ")
lines(fitted(reslm)~t,col="red",lty=2,lwd=2)   # dashed blue line is sin fit
for(i in 1:24){
	segments(i, ipq_mean_hour[i]-ipq_sd_hour[i],i, ipq_mean_hour[i]+ipq_sd_hour[i])
}
minor.tick(nx=5, tick.ratio=.5)
dev.off()



barplot(ipq_count_hour)

####################################################################
save.image(paste("tdr_new", mytdr_d.thresh,".", mytdr_t.thresh, ".RData", sep=""))
####################################################################










