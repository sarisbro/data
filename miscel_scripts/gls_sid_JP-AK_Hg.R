# This script was developed to analyze data in Shoji et al. (unpublished)
# comparing GLS data for rhinoceros auklets sampled on both sides of the Pacific
# (Japan [JP] and Alaska [AK]) in relation to environmental [THg] (total mercury)
# during wintering.
# GLS data are available from the first author.


# R CMD BATCH gls_sid_JP-AK_Hg.R &

library(chron)
library(maps)
library(lubridate)
library(geosphere)
library(fields)
library(MASS)
library(ggmap)
library(caTools)
library(colorspace)
library(TeachingDemos)
library(robust)
library(Hmisc)
library(raster)
library(RColorBrewer)
library(corrplot)
library(PerformanceAnalytics)
library(heatmaply)
library(gridExtra)
library(dunn.test)
library(imputeTS)
library(suncalc)
library(foreach)
library(plot3D)

library(doMC)
registerDoMC(cores=6)


colo_AK_location <- c(-146.33888889, 59.43277778) # c(146° 20′ 20″ O, 59° 25′ 58″ N, ), according to https://www.latlong.net/degrees-minutes-seconds-to-decimal-degrees
colo_JP_location <- c(141.30000000, 44.41666667) # 44°25′N 141°18′E


if(file.exists("salt_2017_ak_jp.RData")){
	load("salt_2017_ak_jp.RData")
}else{
												#########################
												# read in original data #
												#########################
	
	####################################################################################################################
	####################################################################################################################
	
														###########
														# AK data #
														###########
	
	#bird_hg_ak <- read.csv("../180717_Middelton_Hg/MDO_2018_Hg.csv", header=T)
	#colnames(bird_hg_ak) <- c("nb", "GLS_Data", "blood_Hg", "breast_Hg", "primary_Hg", "retrix_Hg", "feces_Hg", "ID", "Ring", "SID.file.name", "z")
	bird_hg_ak <- read.csv("190726_MDO_2018_Hg_new.csv", header=T)
	colnames(bird_hg_ak) <- c("nb", "GLS_Data", "blood_log10Hg", "breast_log10Hg", "primary_log10Hg", "retrix_log10Hg", "feces_log10Hg", "ID", "Ring", "SID.file.name")
	bird_hg_ak$blood_log10Hg <- log10(bird_hg_ak$blood_log10Hg)
	bird_hg_ak$breast_log10Hg <- log10(bird_hg_ak$breast_log10Hg)
	bird_hg_ak$primary_log10Hg <- log10(bird_hg_ak$primary_log10Hg)
	bird_hg_ak$retrix_log10Hg <- log10(bird_hg_ak$retrix_log10Hg)
	bird_hg_ak$feces_log10Hg <- log10(bird_hg_ak$feces_log10Hg)
	
	# eliminate row w/out GLS data
	na_pos <- which(is.na(bird_hg_ak$SID.file.name))
	if(length(na_pos) > 0){
		bird_hg_ak <- bird_hg_ak[-na_pos,]
	}
	
	# get bird_name
	bird_name <- bird_hg_ak$SID.file.name
	bird_name <- sub(".deg", "", bird_name)
	bird_name <- sub("md_rhau_BA", "", bird_name)
	bird_name <- substr(bird_name, 5, length(bird_name))
	bird_hg_ak <- cbind(bird_hg_ak, bird_name)
	
	bird_hg_ak_dup <- bird_hg_ak
	
	######### GLS ##############
	
	dat_AK_files <- list.files(path = "../180516_GLS_2017/middletown/2day_median/", pattern = ".csv")
	gls_AK_dat <- NULL
	
	# read GLS files
	for(i in dat_AK_files){
		bird_id <- substr(i, 9, 15)
		bird_sp <- substr(i, 4, 7)
		print(paste0("Now doing ", bird_id, " (", bird_sp, ")"))
		dat <- read.csv(paste0("../180516_GLS_2017/middletown/2day_median/", i), header=T)
		datetime <- as.character(dat$date)
		datetime1 <- substr(datetime, 1, 11)
		datetime2 <- substr(datetime, 13, 20)
		dat$date <- as.POSIXct(chron(dates = datetime1, times = datetime2, format=c("d-mon-y","h:m:s")))
		distances <- ave_speed <- c()
		for(j in 1:length(dat[,1])){
			if(j == 1){
				distances[j] <- distMeeus(c(colo_AK_location[1], colo_AK_location[2]), c(dat$med_lon[j],dat$med_lat[j])) / 1000
				ave_speed[j] <- NA
			}else{
				distances[j] <- distMeeus(c(dat$med_lon[j-1],dat$med_lat[j-1]), c(dat$med_lon[j],dat$med_lat[j])) / 1000
				ave_speed[j] <- distances[j] / (as.numeric(difftime(dat$date[j], dat$date[j-1], unit="hours")))
			}
		}
		gls_AK_dat <- rbind(gls_AK_dat, cbind(dat, distances, ave_speed))
	}
	
	gls_bird_AK_IDs <- as.character(unique(gls_AK_dat$bird))
	


	########## SID ###########
	salt_files_ak <- list.files(path = "../180516_GLS_2017/middletown/salt/", pattern=".deg")
	
	# read files
	all_salt_ak <- data.frame(
		status = character(),
		datetime = character(),
		col3 = numeric(),
		salt = numeric(),
		filename = character(),
		bird_name = character(),
		hg_blood = numeric(),
		hg_breast = numeric(),
		hg_rectrix = numeric(),
		# not in JP data
		hg_primary = numeric(),
		hg_feces = numeric()
	)
	
	for(i in salt_files_ak){
		print(paste0("Now doing: ", i))
		tmp_salt <- bird_name <- filename <- NULL
		tmp_salt <- read.table(paste0("../180516_GLS_2017/middletown/salt/", i), skip=19, header=T)
		tmp_salt[,1] <- as.POSIXct(chron(dates = as.character(tmp_salt[,1]), times = as.character(tmp_salt[,2]), format=c("d/m/y","h:m:s")))
		# get filename, bird_name, sp_name
		filename <- i #sub(".act", "", i)
		bird_name <- substr(filename, 15, nchar(filename))
		bird_name <- sub(".deg", "", bird_name)
		tmp_salt <- cbind(tmp_salt, filename, bird_name, bird_hg_ak[bird_hg_ak$SID.file.name == i,c(3:7)])
		all_salt_ak <- rbind(all_salt_ak, tmp_salt)
		
	}
	colnames(all_salt_ak) <- c("datetime", "datetime2", "salt", "filename", "bird_name", "hg_blood", "hg_breast", "hg_primary", "hg_rectrix", "hg_feces")
	bird_IDs_sid_ak <- unique(as.character(all_salt_ak$bird_name))
	
	
	
	
														###########
														# JP data #
														###########
	
	bird_hg_jp <- read.csv("../171213_salt2017/teuri_gls_salt_table.csv", header=T)
	names(bird_hg_jp) <- c("filename", "Nest", "bird_name", "blood_log10Hg", "breast_log10Hg", "retrix_log10Hg", "blood_d13C", "primary_d13C", "blood_d15N", "primary_d15N")
	bird_hg_jp$blood_log10Hg <- log10(bird_hg_jp$blood_log10Hg)
	bird_hg_jp$breast_log10Hg <- log10(bird_hg_jp$breast_log10Hg)
	bird_hg_jp$retrix_log10Hg <- log10(bird_hg_jp$retrix_log10Hg)
	
	# eliminate row w/out GLS data
	na_pos <- which(is.na(bird_hg_jp$filename))
	if(length(na_pos) > 0){
		bird_hg_jp <- bird_hg_jp[-na_pos,]
	}
	
	bird_hg_jp_dup <- bird_hg_jp
	
	######### GLS ##############
	dat_JP_files <- list.files(path = "../180516_GLS_2017/teuri/2day_median/", pattern = ".csv")
	gls_JP_dat <- NULL
	
	# read GLS files
	for(i in dat_JP_files){
		bird_id <- substr(i, 9, 16)
		bird_sp <- substr(i, 4, 7)
		print(paste0("Now doing ", bird_id, " (", bird_sp, ")"))
		dat <- read.csv(paste0("../180516_GLS_2017/teuri/2day_median/", i), header=T)
		datetime <- as.character(dat$date)
		dat$date <- as.POSIXct(strptime(dat$date, format = "%Y-%m-%d"))
		distances <- ave_speed <- c()
		for(j in 1:length(dat[,1])){
			if(j == 1){
				distances[j] <- distMeeus(c(colo_JP_location[1], colo_JP_location[2]), c(dat$med_lon[j],dat$med_lat[j])) / 1000
				ave_speed[j] <- NA
			}else{
				distances[j] <- distMeeus(c(dat$med_lon[j-1],dat$med_lat[j-1]), c(dat$med_lon[j],dat$med_lat[j])) / 1000
				ave_speed[j] <- distances[j] / (as.numeric(difftime(dat$date[j], dat$date[j-1], unit="hours")))
			}
		}
		gls_JP_dat <- rbind(gls_JP_dat, cbind(dat, distances, ave_speed))
	}
	
	gls_bird_JP_IDs <- as.character(unique(gls_JP_dat$bird))
	
	########## SID ###########
	salt_files_jp <- list.files(path="../180516_GLS_2017/teuri/salt", pattern=".act")
	
	# read files
	all_salt_jp <- data.frame(
		status = character(),
		datetime = character(),
		col3 = numeric(),
		salt = numeric(),
		filename = character(),
		bird_name = character(),
		hg_blood = numeric(),
		hg_breast = numeric(),
		hg_rectrix = numeric()
	)
	
	for(i in salt_files_jp){
		print(paste0("Now doing: ", i))
		tmp_salt <- bird_name <- filename <- NULL
		tmp_salt <- read.csv(paste0("../180516_GLS_2017/teuri/salt/", i), header=F)
		# some files have > 4 columns: ignore these
		tmp_salt <- tmp_salt[,1:4]
		tmp_salt[,2] <- as.POSIXct(strptime(tmp_salt[,2], format = "%d-%B-%Y"))
		# get filename, bird_name, sp_name
		filename <- sub(".act", "", i)
		bird_name <- substr(filename, 9, nchar(filename))
		tmp_salt <- cbind(tmp_salt, filename, bird_name, bird_hg_jp[bird_hg_jp$bird_name == bird_name,4:10])
		all_salt_jp <- rbind(all_salt_jp, tmp_salt)
	}
	colnames(all_salt_jp) <- c("status", "datetime", "xls_time", "salt", "filename", "bird_name", "hg_blood", "hg_breast", "hg_rectrix", "d13C_blood", "d13C_prime", "d15N_blood", "d15N_prime")
	bird_IDs_sid_jp <- unique(as.character(all_salt_jp$bird_name))
	#plot(as.numeric(all_salt_jp$hg_blood))
	#plot(all_salt_jp$datetime)
	
	
	
	
	save.image("salt_2017_ak_jp.RData")
	####################################################################################################################
	####################################################################################################################
}

####################################################################################################################
####################################################################################################################


												######################
												# filter AK GLS data #
												######################
colnames(gls_AK_dat)[2] <- "datetime"

# define migratory seasons; from end-Jul to mid-Apr
gls_years <- year(gls_AK_dat$datetime)
gls_years_m <- min(gls_years, na.rm=T)
gls_years_M <- max(gls_years, na.rm=T)
gls_months <- month(gls_AK_dat$datetime)
gls_days <- day(gls_AK_dat$datetime)
# NB. no entries are invalid...
which(is.na(gls_years))
# get seasons
gls_year_list <- gls_years_m:gls_years_M
gls_seasons <- rep("RS_colony", length(gls_AK_dat$datetime))
start_month <- 8#7
start_day   <- 15#15
end_month   <- 4#3
end_day     <- 15#15
for(i in 1:(length(unique(gls_years))-1)){
	gls_seasons[(gls_months >= start_month) & (gls_years == gls_year_list[i])] <- paste0(gls_year_list[i], "/", gls_year_list[i+1])
	gls_seasons[(gls_months == start_month) & (gls_days < start_day) & (gls_years == gls_year_list[i])] <- "RS_colony"
	gls_seasons[(gls_months <= end_month) & (gls_years == gls_year_list[i+1])] <- paste0(gls_year_list[i], "/", gls_year_list[i+1])
	gls_seasons[(gls_months == end_month) & (gls_days > end_day) & (gls_years == gls_year_list[i+1])] <- "RS_colony"
}
sort(unique(gls_seasons))
gls_AK_dat <- cbind(gls_AK_dat, gls_seasons)
	
# filter GLS data
gls_AK_dat_ori <- gls_AK_dat
length(gls_AK_dat[,1])
# keep only confidence >= 5
gls_AK_dat  <- gls_AK_dat[gls_AK_dat$conf >= 5,]
length(gls_AK_dat[,1])
# exclude breeding period from mid-March to mid-July
gls_AK_dat  <- gls_AK_dat[gls_AK_dat$gls_seasons != "RS_colony",]
length(gls_AK_dat[,1])
	
# re-computes distances (in km); checked at: http://www.nhc.noaa.gov/gccalc.shtml
distances <- c()
for(i in 2:length(gls_AK_dat[,1])){
	# for each bird / season
	if(gls_AK_dat$bird[i] == gls_AK_dat$bird[i-1]){
		distances[i] <- distMeeus(c(gls_AK_dat$med_lon[i-1],gls_AK_dat$med_lat[i-1]), c(gls_AK_dat$med_lon[i],gls_AK_dat$med_lat[i])) / 1000
	}else{
		distances[i] <- 0
	}
}

gls_AK_dat <- cbind(gls_AK_dat, distances)
gls_AK_dat <- gls_AK_dat[complete.cases(gls_AK_dat),]



												######################
												# filter AK SID data #
												######################



# define migratory seasons; from mid-July to mid-March
sid_years <- year(all_salt_ak$datetime)
sid_years_m <- min(sid_years, na.rm=T)
sid_years_M <- max(sid_years, na.rm=T)
sid_months <- month(all_salt_ak$datetime)
sid_days <- day(all_salt_ak$datetime)
# NB. no entries are invalid...
no_sid_dates <- which(is.na(sid_years))
no_sid_dates
# eliminate entries with invalid dates
all_sid <- all_salt_ak #[-no_sid_dates,]
sid_years <- year(all_sid$datetime)
sid_years_m <- min(sid_years, na.rm=T)
sid_years_M <- max(sid_years, na.rm=T)
sid_months <- month(all_sid$datetime)
sid_days <- day(all_sid$datetime)
# get seasons
sid_year_list <- sid_years_m:sid_years_M
sid_seasons <- rep("RS_colony", length(all_sid$datetime))
for(i in 1:(length(unique(sid_years))-1)){
	sid_seasons[(sid_months >= 8) & (sid_years == sid_year_list[i])] <- paste0(sid_year_list[i], "/", sid_year_list[i+1])
	sid_seasons[(sid_months == 8) & (sid_days < 15) & (sid_years == sid_year_list[i])] <- "RS_colony"
	sid_seasons[(sid_months <= 4) & (sid_years == sid_year_list[i+1])] <- paste0(sid_year_list[i], "/", sid_year_list[i+1])
	sid_seasons[(sid_months == 4) & (sid_days > 15) & (sid_years == sid_year_list[i+1])] <- "RS_colony"
}
sort(unique(sid_seasons))
all_sid_ak <- cbind(all_sid, sid_seasons)

# eliminate SID at colony
at_colony <- which(all_sid_ak$sid_seasons == "RS_colony")
all_sid_ak <- all_sid_ak[-at_colony,]
dim(all_sid_ak)


												######################
												# filter JP GLS data #
												######################

colnames(gls_JP_dat)[2] <- "datetime"

# define migratory seasons; from end-Jul to mid-Apr
gls_years <- year(gls_JP_dat$datetime)
gls_years_m <- min(gls_years, na.rm=T)
gls_years_M <- max(gls_years, na.rm=T)
gls_months <- month(gls_JP_dat$datetime)
gls_days <- day(gls_JP_dat$datetime)
# NB. no entries are invalid...
which(is.na(gls_years))
# get seasons
gls_year_list <- gls_years_m:gls_years_M
gls_seasons <- rep("RS_colony", length(gls_JP_dat$datetime))
start_month <- 7#7
start_day   <- 31#15
end_month   <- 3#3
end_day     <- 31#15
for(i in 1:(length(unique(gls_years))-1)){
	gls_seasons[(gls_months >= start_month) & (gls_years == gls_year_list[i])] <- paste0(gls_year_list[i], "/", gls_year_list[i+1])
	gls_seasons[(gls_months == start_month) & (gls_days < start_day) & (gls_years == gls_year_list[i])] <- "RS_colony"
	gls_seasons[(gls_months <= end_month) & (gls_years == gls_year_list[i+1])] <- paste0(gls_year_list[i], "/", gls_year_list[i+1])
	gls_seasons[(gls_months == end_month) & (gls_days > end_day) & (gls_years == gls_year_list[i+1])] <- "RS_colony"
}
sort(unique(gls_seasons))
gls_JP_dat <- cbind(gls_JP_dat, gls_seasons)
	
# filter GLS data
gls_JP_dat_ori <- gls_JP_dat
length(gls_JP_dat[,1])
# keep only confidence >= 5
gls_JP_dat  <- gls_JP_dat[gls_JP_dat$conf >= 5,]
length(gls_JP_dat[,1])
# exclude breeding period from mid-March to mid-July
gls_JP_dat  <- gls_JP_dat[gls_JP_dat$gls_seasons != "RS_colony",]
length(gls_JP_dat[,1])
	
# re-computes distances (in km); checked at: http://www.nhc.noaa.gov/gccalc.shtml
distances <- c()
for(i in 2:length(gls_JP_dat[,1])){
	# for each bird / season
	if(gls_JP_dat$bird[i] == gls_JP_dat$bird[i-1]){
		distances[i] <- distMeeus(c(gls_JP_dat$med_lon[i-1],gls_JP_dat$med_lat[i-1]), c(gls_JP_dat$med_lon[i],gls_JP_dat$med_lat[i])) / 1000
	}else{
		distances[i] <- 0
	}
}

gls_JP_dat <- cbind(gls_JP_dat, distances)
gls_JP_dat <- gls_JP_dat[complete.cases(gls_JP_dat),]



												######################
												# filter JP SID data #
												######################



# define migratory seasons; from mid-July to mid-March
sid_years <- year(all_salt_jp$datetime)
sid_years_m <- min(sid_years, na.rm=T)
sid_years_M <- max(sid_years, na.rm=T)
sid_months <- month(all_salt_jp$datetime)
sid_days <- day(all_salt_jp$datetime)
# NB. no entries are invalid...
no_sid_dates <- which(is.na(sid_years))
no_sid_dates
# eliminate entries with invalid dates
all_sid <- all_salt_jp #[-no_sid_dates,]
sid_years <- year(all_sid$datetime)
sid_years_m <- min(sid_years, na.rm=T)
sid_years_M <- max(sid_years, na.rm=T)
sid_months <- month(all_sid$datetime)
sid_days <- day(all_sid$datetime)
# get seasons
sid_year_list <- sid_years_m:sid_years_M
sid_seasons <- rep("RS_colony", length(all_sid$datetime))
for(i in 1:(length(unique(sid_years))-1)){
	sid_seasons[(sid_months >= 7) & (sid_years == sid_year_list[i])] <- paste0(sid_year_list[i], "/", sid_year_list[i+1])
	sid_seasons[(sid_months == 7) & (sid_days < 31) & (sid_years == sid_year_list[i])] <- "RS_colony"
	sid_seasons[(sid_months <= 3) & (sid_years == sid_year_list[i+1])] <- paste0(sid_year_list[i], "/", sid_year_list[i+1])
	sid_seasons[(sid_months == 3) & (sid_days > 31) & (sid_years == sid_year_list[i+1])] <- "RS_colony"
}
sort(unique(sid_seasons))
all_sid_jp <- cbind(all_sid, sid_seasons)

# eliminate SID at colony
at_colony <- which(all_sid_jp$sid_seasons == "RS_colony")
all_sid_jp <- all_sid_jp[-at_colony,]
dim(all_sid_jp)


save.image("salt_2017_ak_jp_filtered.RData")

####################################################################################################################
####################################################################################################################
#i=151858
											#############################
											# remove nighttime from SID #
											#############################

## AK SID 
if(!file.exists("salt_2017_ak_jp_filtered_NoNite_AK.done")){
	if(file.exists("salt_2017_ak_jp_filtered_NoNite_AK.RData") ){
		load("salt_2017_ak_jp_filtered_NoNite_AK.RData")
		start_i <- i
	}else{
		all_sid_ak_with_night <- all_sid_ak
		all_sid_tmp <- NULL
		start_i <- 1
	}
	sunrise_AK <- sunset_AK <- character()
	for(i in start_i:length(all_sid_ak[,1])){
		sunTimes <- cur_lat <- cur_lon <- cur_date <- cur_bird <- cur_gls <- date_diffs <- pos <- NULL
		if(!(i %% 1000)){
			print(paste0("Now doing AK sunTimes: ", i))
			save(all_sid_tmp, i, file="salt_2017_ak_jp_filtered_NoNite_AK.RData")
		}
		cur_bird <- all_sid_ak$bird_name[i]
		cur_date <- lubridate::date(all_sid_ak$datetime[i])
		cur_gls <- gls_AK_dat[gls_AK_dat$bird == cur_bird,]
		date_diffs <- difftime(cur_gls$datetime, cur_date)
		pos <- which.min(date_diffs)
		cur_lat <- cur_gls$med_lat[pos]
		cur_lon <- cur_gls$med_lon[pos]
	
		sunTimes <- try(getSunlightTimes(date = cur_date, lat = cur_lat, lon = cur_lon, tz="America/Anchorage"))
		sunrise_AK[i] <- as.character(sunTimes$sunrise)
		sunset_AK[i] <- as.character(sunTimes$sunset)
		
		if((all_sid_ak$datetime[i] > force_tz(sunTimes$sunrise, tzone = "")) & (all_sid_ak$datetime[i] < force_tz(sunTimes$sunset, tzone = ""))){
			all_sid_tmp <- rbind(all_sid_tmp, all_sid_ak[i,])
		}
	}
	
	all_sid_ak <- all_sid_tmp
	system("touch salt_2017_ak_jp_filtered_NoNite_AK.done")
	save.image("salt_2017_ak_jp_filtered.RData")
}


## JP SID 

time2 <- as.POSIXct(all_sid_jp$xls_time * (60*60*24), origin="1899-12-29", tz="Japan")
all_sid_jp <- cbind(all_sid_jp, time2)

if(!file.exists("salt_2017_ak_jp_filtered_NoNite_JP.done")){
	if(file.exists("salt_2017_ak_jp_filtered_NoNite_JP.RData")){
		load("salt_2017_ak_jp_filtered_NoNite_JP.RData")
		start_i <- i
	}else{
		all_sid_jp_with_night <- all_sid_jp
		all_sid_tmp <- NULL
		start_i <- 1
	}
	sunrise_JP <- sunset_JP <- character()
	for(i in start_i:length(all_sid_jp[,1])){
		sunTimes <- cur_lat <- cur_lon <- cur_date <- cur_bird <- cur_gls <- date_diffs <- pos <- NULL
		if(!(i %% 1000)){
			print(paste0("Now doing JP sunTimes: ", i))
			save(all_sid_tmp, i, file="salt_2017_ak_jp_filtered_NoNite_JP.RData")
		}
		cur_bird <- all_sid_jp$bird_name[i]
		cur_date <- lubridate::date(all_sid_jp$time2[i])
		cur_gls <- gls_JP_dat[gls_JP_dat$bird == cur_bird,]
		date_diffs <- difftime(cur_gls$datetime, cur_date)
		pos <- which.min(date_diffs)
		cur_lat <- cur_gls$med_lat[pos]
		cur_lon <- cur_gls$med_lon[pos]
	
		sunTimes <- try(getSunlightTimes(date = cur_date, lat = cur_lat, lon = cur_lon, tz="Japan"))
		sunrise_JP[i] <- as.character(sunTimes$sunrise)
		sunset_JP[i] <- as.character(sunTimes$sunset)
		
		if((all_sid_jp$time2[i] > (force_tz(sunTimes$sunrise, tzone = "Japan") - (60*60*24))) & (all_sid_jp$time2[i] < (force_tz(sunTimes$sunset, tzone = "Japan") - (60*60*24)))){
			all_sid_tmp <- rbind(all_sid_tmp, all_sid_jp[i,])
		}
	}
	
	all_sid_jp <- all_sid_tmp
	system("touch salt_2017_ak_jp_filtered_NoNite_JP.done")
	save.image("salt_2017_ak_jp_filtered.RData")
}

load("salt_2017_ak_jp_filtered.RData")

####################################################################################################################
####################################################################################################################

										#################################
										# filtering by dist from colony #
										#################################

pdf("distcolony.pdf", width=6, height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,1))
# AK
hist(gls_AK_dat$distcolony, 100, main="", xlab="Distance from colony -- all AK birds (km)", xlim=c(0,max(c(gls_AK_dat$distcolony, gls_JP_dat$distcolony))), col="#56B4E9")
distcolo_AK_cutoff <- quantile(gls_AK_dat$distcolony, probs=.95, na.rm=T)
abline(v = distcolo_AK_cutoff, col="orange", lty=2)
legend("topright", legend=paste0("95% cutoff at ", format(distcolo_AK_cutoff, nsmall=2, digits=2), "km"), pch=NULL)
# JP
hist(gls_JP_dat$distcolony, 55, main="", xlab="Distance from colony -- all JP birds (km)", xlim=c(0,max(c(gls_AK_dat$distcolony, gls_JP_dat$distcolony))), col="orange")
distcolo_JP_cutoff <- quantile(gls_JP_dat$distcolony, probs=.95, na.rm=T)
abline(v = distcolo_JP_cutoff, col="orange", lty=2)
legend("topright", legend=paste0("95% cutoff at ", format(distcolo_JP_cutoff, nsmall=2, digits=2), "km"), pch=NULL)
dev.off()

gls_AK_dat_distColCutoff <- gls_AK_dat[gls_AK_dat$distcolony <= distcolo_AK_cutoff, ]
min_AK_lon2 <- min(gls_AK_dat_distColCutoff$med_lon)
max_AK_lon2 <- max(gls_AK_dat_distColCutoff$med_lon)
min_AK_lat2 <- min(gls_AK_dat_distColCutoff$med_lat)
max_AK_lat2 <- max(gls_AK_dat_distColCutoff$med_lat)

gls_JP_dat_distColCutoff <- gls_JP_dat[gls_JP_dat$distcolony <= distcolo_JP_cutoff, ]
min_JP_lon2 <- min(gls_JP_dat_distColCutoff$med_lon)
max_JP_lon2 <- max(gls_JP_dat_distColCutoff$med_lon)
min_JP_lat2 <- min(gls_JP_dat_distColCutoff$med_lat)
max_JP_lat2 <- max(gls_JP_dat_distColCutoff$med_lat)

min_lon2 <- min(min_AK_lon2, min_JP_lon2)
max_lon2 <- max(max_AK_lon2, max_JP_lon2)
min_lat2 <- min(min_AK_lat2, min_JP_lat2)
max_lat2 <- max(max_AK_lat2, max_JP_lat2)

pdf("gls_global_map_all0_zoom_distColCutoff.pdf", width=12, height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,2))
# AK
for(i in gls_bird_AK_IDs){
	map('world', interior=F, fill=F, col="gray80", xlim=c(min_AK_lon2, max_AK_lon2), ylim=c(min_AK_lat2, max_AK_lat2))
	points(gls_AK_dat_distColCutoff$med_lon[gls_AK_dat_distColCutoff$bird == i], gls_AK_dat_distColCutoff$med_lat[gls_AK_dat_distColCutoff$bird == i], pch=20, cex=.6, col="#56B4E9")
	title(paste0("GLS map for ", i))
}
# JP
for(i in gls_bird_JP_IDs){
	map('world', interior=F, fill=F, col="gray80", xlim=c(min_JP_lon2, max_JP_lon2), ylim=c(min_JP_lat2, max_JP_lat2))
	points(gls_JP_dat_distColCutoff$med_lon[gls_JP_dat_distColCutoff$bird == i], gls_JP_dat_distColCutoff$med_lat[gls_JP_dat_distColCutoff$bird == i], pch=20, cex=.6, col="orange")
	title(paste0("GLS map for ", i))
}
dev.off()

###########################################
# https://stackoverflow.com/questions/5353184/fixing-maps-library-data-for-pacific-centred-0-360-longitude-display
plot.map<- function(database,center,...){
    Obj <- map(database,...,plot=F)
    coord <- cbind(Obj[[1]],Obj[[2]])

    # split up the coordinates
    id <- rle(!is.na(coord[,1]))
    id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
    polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})

    # split up polygons that differ too much
    polygons <- lapply(polygons,function(x){
        x[,1] <- x[,1] + center
        x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
        if(sum(diff(x[,1])>300,na.rm=T) >0){
          id <- x[,1] < 0
          x <- rbind(x[id,],c(NA,NA),x[!id,])
       }
       x
    })
    # reconstruct the object
    polygons <- do.call(rbind,polygons)
    Obj[[1]] <- polygons[,1]
    Obj[[2]] <- polygons[,2]

    map(Obj,...)
}
###########################################



pdf("gls_global_map_all0_zoom_distColCutoff_density.pdf", width=12, height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
lon_offset <- 200
plot.map("world", center= lon_offset, col="black", bg="white", fill=F, xlim=c(min_lon2, max_lon2), ylim=c(min_lat2, max_lat2), mar=c(0,0,0,0))
# AK
z_AK <- kde2d(gls_AK_dat_distColCutoff$med_lon, gls_AK_dat_distColCutoff$med_lat, n=150)
contour((z_AK$x + lon_offset), z_AK$y, z_AK$z, drawlabels=F, nlevels=13, add=T, col="#56B4E9")
# JP
z_JP <- kde2d(gls_JP_dat_distColCutoff$med_lon, gls_JP_dat_distColCutoff$med_lat, n=150)
contour((z_JP$x + lon_offset -360), z_JP$y, z_JP$z, drawlabels=F, nlevels=13, add=T, col="orange")
plot.map("world", center= lon_offset, col="gray95", bg="white", fill=T, xlim=c(min_lon2, max_lon2), ylim=c(min_lat2, max_lat2), mar=c(0,0,0,0), add=T)
dev.off()

####################################################################################################################
####################################################################################################################

								#########################################
								# activity data: flight / rest / forage #
								#########################################


### AK ###
daily_activity_ak <- data.frame()
max_salt_ak <- max(all_sid_ak$salt)
upper_salt_ak <- .98 * max_salt_ak
lower_salt_ak <- .02 * max_salt_ak

for(i in bird_IDs_sid_ak){
	print(paste0("Bird: ", i))
	cur_bird <- all_sid_ak[all_sid_ak$bird_name == i,]
	yeardays <- unique(yday(cur_bird$datetime))
	fulldates <- as.character(cur_bird$datetime)
	fulldates <- substr(fulldates, 1, 10)
	fulldates <- unique(fulldates)
	for(j in yeardays){
		cur_day_salt <- cur_bird$salt[yeardays == j]
		sum(cur_day_salt >= upper_salt_ak)
		daily_activity_ak <- rbind(daily_activity_ak, cbind(bird=i, yearday=j, Date= fulldates[which(yeardays == j)],
					dive=sum(cur_day_salt >= upper_salt_ak),
					fly=sum(cur_day_salt <= lower_salt_ak), 
					rest=sum((cur_day_salt > lower_salt_ak) & (cur_day_salt < upper_salt_ak)), 
					propdive=sum(cur_day_salt >= upper_salt_ak)/length(cur_day_salt),
					propfly=sum(cur_day_salt <= lower_salt_ak)/length(cur_day_salt),
					proprest=sum((cur_day_salt > lower_salt_ak) & (cur_day_salt < upper_salt_ak))/length(cur_day_salt)
		))
	}
}
# convert to dates
daily_activity_ak$Date <- base::as.Date(as.character(daily_activity_ak$Date))
# convert to numeric
daily_activity_ak[, c(2,4:9)] <- sapply(daily_activity_ak[, c(2,4:9)], as.character)
daily_activity_ak[, c(2,4:9)] <- sapply(daily_activity_ak[, c(2,4:9)], as.numeric)



### JP ###
daily_activity_jp <- data.frame()
max_salt_jp <- max(all_sid_jp$salt)
upper_salt_jp <- .98 * max_salt_jp
lower_salt_jp <- .02 * max_salt_jp

for(i in bird_IDs_sid_jp){
	print(paste0("Bird: ", i))
	cur_bird <- all_sid_jp[all_sid_jp$bird_name == i,]
	yeardays <- unique(yday(cur_bird$datetime))
	fulldates <- as.character(unique(cur_bird$datetime))
	for(j in yeardays){
		cur_day_salt <- cur_bird$salt[yeardays == j]
		sum(cur_day_salt >= upper_salt_jp)
		daily_activity_jp <- rbind(daily_activity_jp, cbind(bird=i, yearday=j, Date= fulldates[which(yeardays == j)],
					dive=sum(cur_day_salt >= upper_salt_jp),
					fly=sum(cur_day_salt <= lower_salt_jp), 
					rest=sum((cur_day_salt > lower_salt_jp) & (cur_day_salt < upper_salt_jp)), 
					propdive=sum(cur_day_salt >= upper_salt_jp)/length(cur_day_salt),
					propfly=sum(cur_day_salt <= lower_salt_jp)/length(cur_day_salt),
					proprest=sum((cur_day_salt > lower_salt_jp) & (cur_day_salt < upper_salt_jp))/length(cur_day_salt)
		))
	}
}
# convert to dates
daily_activity_jp$Date <- base::as.Date(as.character(daily_activity_jp$Date))
# convert to numeric
daily_activity_jp[, c(2,4:9)] <- sapply(daily_activity_jp[, c(2,4:9)], as.character)
daily_activity_jp[, c(2,4:9)] <- sapply(daily_activity_jp[, c(2,4:9)], as.numeric)

# compare activity ak vs. jp
a1 <- dunn.test(list(daily_activity_ak$propfly, daily_activity_jp$propfly))
a2 <- dunn.test(list(daily_activity_ak$proprest, daily_activity_jp$proprest))
a3 <- dunn.test(list(daily_activity_ak$propdive, daily_activity_jp$propdive))

# plot activity
dailyFly_ak <- ggplot(daily_activity_ak, aes(x= bird, y= propfly)) +
	     geom_boxplot(outlier.colour=NA, fill="#56B4E9", colour="grey40") +
	     theme_bw() + theme(legend.position="none") +
	     labs(title="AK", x = NULL, y="Daily fly %") + ylim(0., 1.)
dailyRest_ak <- ggplot(daily_activity_ak, aes(x= bird, y= proprest)) +
	     geom_boxplot(outlier.colour=NA, fill="#56B4E9", colour="grey40") +
	     theme_bw() + theme(legend.position="none") +
	     labs(x = NULL, y="Daily rest %") + ylim(0., 1.)
dailyDive_ak <- ggplot(daily_activity_ak, aes(x= bird, y= propdive)) +
	     geom_boxplot(outlier.colour=NA, fill="#56B4E9", colour="grey40") +
	     theme_bw() + theme(legend.position="none") +
	     labs(x = "Birds", y="Daily dive %") + ylim(0., 1.)

dailyFly_jp <- ggplot(daily_activity_jp, aes(x= bird, y= propfly)) +
	     geom_boxplot(outlier.colour=NA, fill="orange", colour="grey40") +
	     theme_bw() + theme(legend.position="none") +
	     annotate("text", label = paste0("AK vs. JP: Z = ", format(a1$Z, nsmall=2, digits=2), "; P = ", format(a1$P, nsmall=4, digits=2)), x = 8, y = 1, size = 4, colour = "black") +
	     labs(title="JP", x = NULL, y= NULL) + ylim(0., 1.)
dailyRest_jp <- ggplot(daily_activity_jp, aes(x= bird, y= proprest)) +
	     geom_boxplot(outlier.colour=NA, fill="orange", colour="grey40") +
	     theme_bw() + theme(legend.position="none") +
	     annotate("text", label = paste0("AK vs. JP: Z = ", format(a2$Z, nsmall=2, digits=2), "; P = ", format(a2$P, nsmall=4, digits=2)), x = 8, y = 1, size = 4, colour = "black") +
	     labs(x = NULL, y= NULL) + ylim(0., 1.)
dailyDive_jp <- ggplot(daily_activity_jp, aes(x= bird, y= propdive)) +
	     geom_boxplot(outlier.colour=NA, fill="orange", colour="grey40") +
	     theme_bw() + theme(legend.position="none") +
	     annotate("text", label = paste0("AK vs. JP: Z = ", format(a3$Z, nsmall=2, digits=2), "; P = ", format(a3$P, nsmall=4, digits=2)), x = 8, y = 1, size = 4, colour = "black") +
	     labs(x = "Birds", y= NULL) + ylim(0., 1.)



pdf("ActivityBirdsBoxplots.pdf", height=8, width=12)
grid.arrange(dailyFly_ak, dailyFly_jp,
			dailyRest_ak, dailyRest_jp,
			dailyDive_ak, dailyDive_jp,
			nrow=3)
dev.off()






####################################################################################################################
####################################################################################################################

												##############
												# in-bird Hg #
												##############

t1 <- dunn.test(list(bird_hg_ak$blood_log10Hg[complete.cases(bird_hg_ak$blood_log10Hg)], bird_hg_jp$blood_log10Hg[complete.cases(bird_hg_jp$blood_log10Hg)]))
t2 <- dunn.test(list(bird_hg_ak$breast_log10Hg[complete.cases(bird_hg_ak$breast_log10Hg)], bird_hg_jp$breast_log10Hg[complete.cases(bird_hg_jp$breast_log10Hg)]))
t3 <- dunn.test(list(bird_hg_ak$retrix_log10Hg[complete.cases(bird_hg_ak$retrix_log10Hg)], bird_hg_jp$retrix_log10Hg[complete.cases(bird_hg_jp$retrix_log10Hg)]))

dat_hg_ak <- bird_hg_ak[,c(3,4,6)]
dat_hg_jp <- bird_hg_jp[,c(4:6)]
dat_hg_df <- data.frame(rbind(cbind(Region = "AK", dat_hg_ak), cbind(Region = "JP", dat_hg_jp)))

p1 <- ggplot(dat_hg_df, aes(x= Region, y= blood_log10Hg)) +
	     geom_boxplot(outlier.colour=NA, fill="grey90", colour="grey40") +
	     geom_point(aes(fill= Region), size=4, shape=21, colour="grey20",
	                position=position_jitter(width=0.2, height=0.1)) +
	     scale_fill_manual(values=c("#56B4E9", "orange")) +
	     theme_bw() + theme(legend.position="none") +
	     labs(title=paste0("blood [Hg]: P = ", format(t1$P, nsmall=4, digits=2)), x = NULL)
p2 <- ggplot(dat_hg_df, aes(x= Region, y= breast_log10Hg)) +
	     geom_boxplot(outlier.colour=NA, fill="grey90", colour="grey40") +
	     geom_point(aes(fill= Region), size=4, shape=21, colour="grey20",
	                position=position_jitter(width=0.2, height=0.1)) +
	     scale_fill_manual(values=c("#56B4E9", "orange")) +
	     theme_bw() + theme(legend.position="none") +
	     labs(title=paste0("breast [Hg]: P = ", format(t2$P, nsmall=4, digits=2)), x = NULL)
p3 <- ggplot(dat_hg_df, aes(x= Region, y= retrix_log10Hg)) +
	     geom_boxplot(outlier.colour=NA, fill="grey90", colour="grey40") +
	     geom_point(aes(fill= Region), size=4, shape=21, colour="grey20",
	                position=position_jitter(width=0.2, height=0.1)) +
	     scale_fill_manual(values=c("#56B4E9", "orange")) +
	     theme_bw() + theme(legend.position="none") +
	     labs(title=paste0("rectrix [Hg]: P = ", format(t3$P, nsmall=4, digits=2)), x = NULL)

#pdf("Hg_jittering.pdf", height=4, width=12)
#grid.arrange(p1, p2, p3, nrow=1)
#dev.off()

											####################
											# environmental Hg #
											####################


# from: https://www.amap.no/mercury-emissions/datasets
#system("wget https://www.amap.no/documents/download/1447")
#system("unzip 1447")
cell_ID <- read.table("../171213_salt2017/cellSurfaceArea05_14feb2013_FS.txt", header=T)
cell_Hg <- read.csv("../171213_salt2017/globalHgEmission2010ALL_14feb2013_FS_SW.txt", header=T)

# AK
# restrict cell IDs to those of interest here
min_lon <- min_AK_lon2 - 5
max_lon <- max_AK_lon2 + 5
min_lat <- min_AK_lat2 - 5
max_lat <- max_AK_lat2 + 5

my_cell_ID <- cell_ID[(cell_ID$centerLat >= min_lat) & (cell_ID$centerLat <= max_lat) & (cell_ID$centerLon >= min_lon) & (cell_ID$centerLon <= max_lon),]

# extract corresponding [THg]
matching_AK_Hg <- matching_AK_lon <- matching_AK_lat <- c()
for(i in 1:length(my_cell_ID[,1])){
	pos <- which(cell_Hg$cell_code == my_cell_ID$cellCode[i])
	if(length(pos) > 0){
		matching_AK_Hg[i] <- cell_Hg$HgT[pos]
		matching_AK_lon[i] <- my_cell_ID$centerLon[i]
		matching_AK_lat[i] <- my_cell_ID$centerLat[i]
	}
}

Hg_AK_df <- data.frame(lon = matching_AK_lon, lat = matching_AK_lat, THg = matching_AK_Hg, Region="AK")
Hg_AK_df <- Hg_AK_df[complete.cases(Hg_AK_df),]


# JP
# restrict cell IDs to those of interest here
min_lon <- min_JP_lon2 - 5
max_lon <- max_JP_lon2 + 5
min_lat <- min_JP_lat2 - 5
max_lat <- max_JP_lat2 + 5

my_cell_ID <- cell_ID[(cell_ID$centerLat >= min_lat) & (cell_ID$centerLat <= max_lat) & (cell_ID$centerLon >= min_lon) & (cell_ID$centerLon <= max_lon),]

# extract corresponding [THg]
matching_JP_Hg <- matching_JP_lon <- matching_JP_lat <- c()
for(i in 1:length(my_cell_ID[,1])){
	pos <- which(cell_Hg$cell_code == my_cell_ID$cellCode[i])
	if(length(pos) > 0){
		matching_JP_Hg[i] <- cell_Hg$HgT[pos]
		matching_JP_lon[i] <- my_cell_ID$centerLon[i]
		matching_JP_lat[i] <- my_cell_ID$centerLat[i]
	}
}

Hg_JP_df <- data.frame(lon = matching_JP_lon, lat = matching_JP_lat, THg = matching_JP_Hg, Region="JP")
Hg_JP_df <- Hg_JP_df[complete.cases(Hg_JP_df),]

t4 <- dunn.test(list(log10(Hg_AK_df$THg[complete.cases(Hg_AK_df$THg)]), log10(Hg_JP_df$THg[complete.cases(Hg_JP_df$THg)])))

Hg_AK_JP_df <- data.frame(rbind(cbind(Region = "AK", Hg_AK_df), cbind(Region = "JP", Hg_JP_df)))
names(Hg_AK_JP_df)[4] <- "log10THg"
Hg_AK_JP_df$log10THg <- log10(Hg_AK_JP_df$log10THg)



p4 <- ggplot(Hg_AK_JP_df, aes(x= Region, y= log10THg)) +
	     geom_boxplot(outlier.colour=NA, fill="grey90", colour="grey40") +
	     geom_point(aes(fill= Region), size=4, shape=21, colour="grey20",
	                position=position_jitter(width=0.2, height=0.1)) +
	     scale_fill_manual(values=c("#56B4E9", "orange")) +
	     theme_bw() + theme(legend.position=c(.85, .15)) +
	     labs(title=paste0("Environmental [Hg]: P = ", format(t4$P, nsmall=4, digits=2)), x = NULL)



pdf("Hg_jittering_all.pdf", height=4, width=16)
grid.arrange(p2, p3, p1, p4, nrow=1)
dev.off()


# added 200524
plant2017 <- read.csv("200524_2017_AKJP_plant.csv", header=T)
plant2017_df <- data.frame(plant2017)
plant2017_df2 <- cbind(plant2017_df, Site_TRT = paste0(plant2017_df$Site, "_", plant2017_df$Plot))
pos2rm <- which(plant2017_df2$Part != "roots")
plant2017_df3 <- plant2017_df2[-pos2rm,]
names(plant2017_df3)[which(names(plant2017_df3) == "Hg_microgram")] <- "log_10Hg_microgram"

t5 <- dunn.test(plant2017_df3$log_10Hg_microgram, plant2017_df3$Site_TRT, method="BH")

p5 <- ggplot(plant2017_df3, aes(x= Site_TRT, y= log_10Hg_microgram)) +
	     geom_boxplot(outlier.colour=NA, fill="grey90", colour="grey40") +
	     geom_point(aes(fill= Site_TRT), size=4, shape=21, colour="grey20",
	                position=position_jitter(width=0.2, height=0.1)) +
	     scale_fill_manual(values=c("darkblue", "#56B4E9", "darkred", "orange")) +
	     theme_bw() + theme(legend.position="none") + xlab("") +
	     labs(title=paste0("In plantae (root) [Hg]", x = NULL)) +
	     annotate("text", label = "(a)", x = 1, y = .22, size = 4, colour = "black") +
	     annotate("text", label = "(b)", x = 2, y = .22, size = 4, colour = "black") +
	     annotate("text", label = "(c)", x = 3, y = .22, size = 4, colour = "black") +
	     annotate("text", label = "(b)", x = 4, y = .22, size = 4, colour = "black") 

pdf("Hg_jittering_all_200524.pdf", height=4, width=20)
grid.arrange(p4, p3, p2, p1, p5, nrow=1)
dev.off()

####################################################################################################################
####################################################################################################################

										###################################
										# average GLS track decomposition #
										#     aka "Mean trajectories"     #
										###################################

moving_ave <- 20

# AK
# compute raw averages
yearday <- yday(gls_AK_dat$datetime)
gls_AK_dat_yday <- cbind(gls_AK_dat, yearday)
ave_ak_lat <- tapply(gls_AK_dat_yday$med_lat, gls_AK_dat_yday$yearday, mean)
ave_ak_lon <- tapply(gls_AK_dat_yday$med_lon, gls_AK_dat_yday$yearday, mean)
sd_ak_lat <- tapply(gls_AK_dat_yday$sd_lat, gls_AK_dat_yday$yearday, mean)
sd_ak_lon <- tapply(gls_AK_dat_yday$sd_lon, gls_AK_dat_yday$yearday, mean)
# re-order by the 1st day of nonbreeding season
first_day_NBS <- yearday[1]
pos_NBS <- which(names(ave_ak_lat) == first_day_NBS)
ave_ak_lat <- c(ave_ak_lat[pos_NBS:length(ave_ak_lat)], ave_ak_lat[1:(pos_NBS-1)])
ave_ak_lon <- c(ave_ak_lon[pos_NBS:length(ave_ak_lon)], ave_ak_lon[1:(pos_NBS-1)])
# compute running averages
ma_ak_lon <- runmean(as.numeric(ave_ak_lon), moving_ave)
ma_ak_lat <- runmean(as.numeric(ave_ak_lat), moving_ave)

plot(ma_ak_lon, ma_ak_lat, pch=20, type="l", xlab="Longitude", ylab="Latitude")
map('world', interior=F, fill=T, col="gray99", add=T)
points(ma_ak_lon, ma_ak_lat, pch=20, type="l")
points(ma_ak_lon, ma_ak_lat, pch=20, col=diverge_hcl(length(ma_ak_lon), h=c(130,43), c=100, l=c(70,90)))
text(colo_AK_location[1], colo_AK_location[2], "+", col="red", cex=2)


# JP
# compute raw averages
yearday <- yday(gls_JP_dat$datetime)
gls_JP_dat_yday <- cbind(gls_JP_dat, yearday)
ave_jp_lat <- tapply(gls_JP_dat_yday$med_lat, gls_JP_dat_yday$yearday, mean)
ave_jp_lon <- tapply(gls_JP_dat_yday$med_lon, gls_JP_dat_yday$yearday, mean)
sd_jp_lat <- tapply(gls_JP_dat_yday$sd_lat, gls_JP_dat_yday$yearday, mean)
sd_jp_lon <- tapply(gls_JP_dat_yday$sd_lon, gls_JP_dat_yday$yearday, mean)
# re-order by the 1st day of nonbreeding season
first_day_NBS <- yearday[1]
pos_NBS <- which(names(ave_jp_lat) == first_day_NBS)
ave_jp_lat <- c(ave_jp_lat[pos_NBS:length(ave_jp_lat)], ave_jp_lat[1:(pos_NBS-1)])
ave_jp_lon <- c(ave_jp_lon[pos_NBS:length(ave_jp_lon)], ave_jp_lon[1:(pos_NBS-1)])
# compute running averages
ma_jp_lon <- runmean(as.numeric(ave_jp_lon), moving_ave)
ma_jp_lat <- runmean(as.numeric(ave_jp_lat), moving_ave)

plot(ma_jp_lon, ma_jp_lat, pch=20, type="l", xlab="Longitude", ylab="Latitude")
map('world', interior=F, fill=T, col="gray99", add=T)
points(ma_jp_lon, ma_jp_lat, pch=20, type="l")
points(ma_jp_lon, ma_jp_lat, pch=20, col=diverge_hcl(length(ma_jp_lon), h=c(130,43), c=100, l=c(70,90)))
text(colo_JP_location[1], colo_JP_location[2], "+", col="red", cex=2)


pdf("MeanTrajectories.pdf",width=4,height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,1))
# AK
plot(ma_ak_lon, ma_ak_lat, type="l", xlab="Longitude", ylab="Latitude")
points(ma_ak_lon, ma_ak_lat, pch=20, col=diverge_hcl(length(ma_ak_lon), h=c(130,43), c=100, l=c(70,90)))
# JP
plot(ma_jp_lon, ma_jp_lat, type="l", xlab="Longitude", ylab="Latitude")
points(ma_jp_lon, ma_jp_lat, pch=20, col=diverge_hcl(length(ma_jp_lon), h=c(130,43), c=100, l=c(70,90)))
dev.off()

# defining 3-part trips (outbound, inbound and wintering grounds) in terms of yearday dates
# AK
yearday <- c(yday(min(gls_AK_dat$datetime)):365,1:yday(max(gls_AK_dat$datetime)))
maxsplit <- ceiling(length(yearday)/3)
x <- seq_along(yearday)
#ak_splits <- split(yearday, ceiling(x/maxsplit))
ak_splits <- split(yearday, c(rep(1,137), rep(2, 70), rep(3, 35)))
# JP
yearday <- c(yday(min(gls_JP_dat$datetime)):365,1:yday(max(gls_JP_dat$datetime)))
maxsplit <- ceiling(length(yearday)/3)
x <- seq_along(yearday)
#jp_splits <- split(yearday, ceiling(x/maxsplit))
jp_splits <- split(yearday, c(rep(1, 125), rep(2, 75), rep(3, 15)))

####################################################################################################################
####################################################################################################################

								############################################
								# matching env Hg for each mean trajectory #
								############################################

# AK
chunk_1_ak <- gls_AK_dat_yday[gls_AK_dat_yday$yearday %in% ak_splits[[1]],]
chunk_2_ak <- gls_AK_dat_yday[gls_AK_dat_yday$yearday %in% ak_splits[[2]],]
chunk_3_ak <- gls_AK_dat_yday[gls_AK_dat_yday$yearday %in% ak_splits[[3]],]

# JP
chunk_1_jp <- gls_JP_dat_yday[gls_JP_dat_yday$yearday %in% jp_splits[[1]],]
chunk_2_jp <- gls_JP_dat_yday[gls_JP_dat_yday$yearday %in% jp_splits[[2]],]
chunk_3_jp <- gls_JP_dat_yday[gls_JP_dat_yday$yearday %in% jp_splits[[3]],]


closestEnvHgReadings <- function(cur_chunk, Hg_neighborhood){
	# find closest [Hg] readings to each GLS data point
	#cur_chunk <- chunk_1_ak
	cur_chunk <- cur_chunk[complete.cases(cur_chunk),]
	#Hg_neighborhood <- 500000  #1000 #100000 # size of neighborhood to compute mean [Hg], in meters
	Hg_env_mean <- Hg_env_sd <- log10Hg_env_mean <- log10Hg_env_sd <- NULL
	for(i in 1:length(cur_chunk[,1])){
		print(paste0("Doing ", i, " out of ", length(cur_chunk[,1]), " hood size ", Hg_neighborhood))
		# find closest Hg reading for focal GLS point
		cur_dist_vect <- NULL
		for(j in 1:length(Hg_AK_JP_df[,1])){
			cur_dist_vect[j] <- distMeeus(c(cur_chunk$med_lon[i],cur_chunk$med_lat[i]), c(Hg_AK_JP_df$lon[j], Hg_AK_JP_df$lat[j]))
		}
		closest_pos <- which.min(cur_dist_vect)
		# find neighborhood readings
		tmp_Hg_env <- NULL
		cur_dist_vect_neighbor <- NULL
		for(j in 1:length(Hg_AK_JP_df[,1])){
			cur_dist_vect_neighbor[j] <- distMeeus(c(Hg_AK_JP_df$lon[closest_pos], Hg_AK_JP_df$lat[closest_pos]), c(Hg_AK_JP_df$lon[j], Hg_AK_JP_df$lat[j]))
		}
		neighborhood_pos <- which(cur_dist_vect_neighbor <= Hg_neighborhood)
		Hg_env_mean[i] <- mean(10^(Hg_AK_JP_df$log10THg[neighborhood_pos]), na.rm=T)
		Hg_env_sd[i] <- sd(10^(Hg_AK_JP_df$log10THg[neighborhood_pos]), na.rm=T)
		log10Hg_env_mean[i] <- mean(Hg_AK_JP_df$log10THg[neighborhood_pos], na.rm=T)
		log10Hg_env_sd[i] <- sd(Hg_AK_JP_df$log10THg[neighborhood_pos], na.rm=T)
	}
	return(cbind(cur_chunk, Hg_env_mean, Hg_env_sd, log10Hg_env_mean, log10Hg_env_sd))
}

closestEnvHgReadingsPar <- function(cur_chunk, Hg_neighborhood){
	# find closest [Hg] readings to each GLS data point
	cur_chunk <- cur_chunk[complete.cases(cur_chunk),]
	#Hg_neighborhood <- 500000  #1000 #100000 # size of neighborhood to compute mean [Hg], in meters
	env_features <- NULL
	env_features <- foreach(i = 1:length(cur_chunk[,1]), .combine=rbind)%dopar%{
	#for(i in 1:length(cur_chunk[,1])){
		Hg_env_mean <- Hg_env_sd <- log10Hg_env_mean <- log10Hg_env_sd <- NULL
		print(paste0("Doing ", i, " out of ", length(cur_chunk[,1]), " hood size ", Hg_neighborhood))
		# find closest Hg reading for focal GLS point
		cur_dist_vect <- NULL
		for(j in 1:length(Hg_AK_JP_df[,1])){
			cur_dist_vect[j] <- distMeeus(c(cur_chunk$med_lon[i],cur_chunk$med_lat[i]), c(Hg_AK_JP_df$lon[j], Hg_AK_JP_df$lat[j]))
		}
		closest_pos <- which.min(cur_dist_vect)
		# find neighborhood readings
		tmp_Hg_env <- NULL
		cur_dist_vect_neighbor <- NULL
		for(j in 1:length(Hg_AK_JP_df[,1])){
			cur_dist_vect_neighbor[j] <- distMeeus(c(Hg_AK_JP_df$lon[closest_pos], Hg_AK_JP_df$lat[closest_pos]), c(Hg_AK_JP_df$lon[j], Hg_AK_JP_df$lat[j]))
		}
		neighborhood_pos <- which(cur_dist_vect_neighbor <= Hg_neighborhood)
		Hg_env_mean <- mean(10^(Hg_AK_JP_df$log10THg[neighborhood_pos]), na.rm=T)
		Hg_env_sd <- sd(10^(Hg_AK_JP_df$log10THg[neighborhood_pos]), na.rm=T)
		log10Hg_env_mean <- mean(Hg_AK_JP_df$log10THg[neighborhood_pos], na.rm=T)
		log10Hg_env_sd <- sd(Hg_AK_JP_df$log10THg[neighborhood_pos], na.rm=T)
		return(cbind(Hg_env_mean, Hg_env_sd, log10Hg_env_mean, log10Hg_env_sd))
	}
	return(cbind(cur_chunk, env_features))
}

chunk_1_ak_envHg <- chunk_2_ak_envHg <- chunk_3_ak_envHg <- list()
chunk_1_jp_envHg <- chunk_2_jp_envHg <- chunk_3_jp_envHg <- list()

# # this one takes *a while*
# #for(neigh_size in seq(500, 30000, 10000)){
# for(neigh_size in seq(500, 3000000, 10000)){
# 	# ak
# 	print(paste0("Now doing chunk_1_ak_envHg ", neigh_size))
# 	chunk_1_ak_envHg[[neigh_size]] <- closestEnvHgReadings(chunk_1_ak, neigh_size)
# 	print(paste0("Now doing chunk_2_ak_envHg ", neigh_size))
# 	chunk_2_ak_envHg[[neigh_size]] <- closestEnvHgReadings(chunk_2_ak, neigh_size)
# 	print(paste0("Now doing chunk_3_ak_envHg ", neigh_size))
# 	chunk_3_ak_envHg[[neigh_size]] <- closestEnvHgReadings(chunk_3_ak, neigh_size)
# 	# jp
# 	print(paste0("Now doing chunk_1_jp_envHg ", neigh_size))
# 	chunk_1_jp_envHg[[neigh_size]] <- closestEnvHgReadings(chunk_1_jp, neigh_size)
# 	print(paste0("Now doing chunk_2_jp_envHg ", neigh_size))
# 	chunk_2_jp_envHg[[neigh_size]] <- closestEnvHgReadings(chunk_2_jp, neigh_size)
# 	print(paste0("Now doing chunk_3_jp_envHg ", neigh_size))
# 	chunk_3_jp_envHg[[neigh_size]] <- closestEnvHgReadings(chunk_3_jp, neigh_size)
# }


# this one is much faster!
#for(neigh_size in seq(500, 30000, 10000)){
for(neigh_size in seq(500, 3000000, 10000)){
	# ak
	print(paste0("Now doing chunk_1_ak_envHg ", neigh_size))
	chunk_1_ak_envHg[[neigh_size]] <- closestEnvHgReadingsPar(chunk_1_ak, neigh_size)
	print(paste0("Now doing chunk_2_ak_envHg ", neigh_size))
	chunk_2_ak_envHg[[neigh_size]] <- closestEnvHgReadingsPar(chunk_2_ak, neigh_size)
	print(paste0("Now doing chunk_3_ak_envHg ", neigh_size))
	chunk_3_ak_envHg[[neigh_size]] <- closestEnvHgReadingsPar(chunk_3_ak, neigh_size)
	# jp
	print(paste0("Now doing chunk_1_jp_envHg ", neigh_size))
	chunk_1_jp_envHg[[neigh_size]] <- closestEnvHgReadingsPar(chunk_1_jp, neigh_size)
	print(paste0("Now doing chunk_2_jp_envHg ", neigh_size))
	chunk_2_jp_envHg[[neigh_size]] <- closestEnvHgReadingsPar(chunk_2_jp, neigh_size)
	print(paste0("Now doing chunk_3_jp_envHg ", neigh_size))
	chunk_3_jp_envHg[[neigh_size]] <- closestEnvHgReadingsPar(chunk_3_jp, neigh_size)
}




save.image("210215_salt_2017_ak_jp_all.RData")



####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
bird_hg_ak <- bird_hg_ak_dup
bird_hg_jp <- bird_hg_jp_dup

chunk_Hg_data <- function(cur_chunk, bird_hg){
	#cur_chunk <- chunk_1_ak_envHg_500000
	#bird_hg <- bird_hg_ak
	chunk_bird_blood_log10Hg <- chunk_bird_breast_log10Hg <- chunk_bird_retrix_log10Hg <- NULL
	for(i in 1:length(cur_chunk[,1])){
		bird_chunk <- as.character(cur_chunk$bird[i])
		bird_chunk_pos <- which(bird_hg$bird_name == bird_chunk)
		chunk_bird_blood_log10Hg[i] <- bird_hg$blood_log10Hg[bird_chunk_pos]
		chunk_bird_breast_log10Hg[i] <- bird_hg$breast_log10Hg[bird_chunk_pos]
		chunk_bird_retrix_log10Hg[i] <- bird_hg$retrix_log10Hg[bird_chunk_pos]
	}
	return(cbind(cur_chunk, chunk_bird_blood_log10Hg, chunk_bird_breast_log10Hg, chunk_bird_retrix_log10Hg))
}

chunk_1_ak_envHg_birdHg <- chunk_2_ak_envHg_birdHg <- chunk_3_ak_envHg_birdHg <- list()
chunk_1_jp_envHg_birdHg <- chunk_2_jp_envHg_birdHg <- chunk_3_jp_envHg_birdHg <- list()

# super quick to run
for(neigh_size in seq(500, 210500, 10000)){
#for(neigh_size in seq(500, 3000000, 10000)){
	# ak
	print(paste0("Now doing chunk_1_ak_envHg ", neigh_size))
	chunk_1_ak_envHg_birdHg[[neigh_size]] <- chunk_Hg_data(chunk_1_ak_envHg[[neigh_size]], bird_hg_ak)
	print(paste0("Now doing chunk_2_ak_envHg ", neigh_size))
	chunk_2_ak_envHg_birdHg[[neigh_size]] <- chunk_Hg_data(chunk_2_ak_envHg[[neigh_size]], bird_hg_ak)
	print(paste0("Now doing chunk_3_ak_envHg ", neigh_size))
	chunk_3_ak_envHg_birdHg[[neigh_size]] <- chunk_Hg_data(chunk_3_ak_envHg[[neigh_size]], bird_hg_ak)
	# jp
	print(paste0("Now doing chunk_1_jp_envHg ", neigh_size))
	chunk_1_jp_envHg_birdHg[[neigh_size]] <- chunk_Hg_data(chunk_1_jp_envHg[[neigh_size]], bird_hg_jp)
	print(paste0("Now doing chunk_2_jp_envHg ", neigh_size))
	chunk_2_jp_envHg_birdHg[[neigh_size]] <- chunk_Hg_data(chunk_2_jp_envHg[[neigh_size]], bird_hg_jp)
	print(paste0("Now doing chunk_3_jp_envHg ", neigh_size))
	chunk_3_jp_envHg_birdHg[[neigh_size]] <- chunk_Hg_data(chunk_3_jp_envHg[[neigh_size]], bird_hg_jp)

}


# ploting results

BirdBlood_vs_envHg1 <- BirdBreast_vs_envHg1 <- BirdRectrix_vs_envHg1 <- list()
BirdBlood_vs_envHg2 <- BirdBreast_vs_envHg2 <- BirdRectrix_vs_envHg2 <- list()
BirdBlood_vs_envHg3 <- BirdBreast_vs_envHg3 <- BirdRectrix_vs_envHg3 <- list()

for(neigh_size in seq(500, 210500, 10000)){
#for(neigh_size in seq(500, 3000000, 10000)){
	# stage_1
	dat_Hg_env_bird <- data.frame(rbind(
		cbind(chunk_1_ak_envHg_birdHg[[neigh_size]], Region="AK"), 
		cbind(chunk_1_jp_envHg_birdHg[[neigh_size]], Region="JP")
	))
	BirdBlood_vs_envHg1[[neigh_size]] <- ggplot(dat_Hg_env_bird, aes(x=log10Hg_env_mean, y= chunk_bird_blood_log10Hg, col=Region)) +
	  geom_point() + scale_color_manual(values=c("#56B4E9", "orange")) +
	  stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
	  #xlim(0,160) +
	  theme_bw() + theme(legend.position="none") +
	  labs(title=paste0("stage_1:", neigh_size,"m; P_AK=",format(coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4], nsmall=4, digits=2),
	                    "; P_JP=",format(coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4], nsmall=4, digits=2)
	  ))
	BirdBreast_vs_envHg1[[neigh_size]] <- ggplot(dat_Hg_env_bird, aes(x=log10Hg_env_mean, y= chunk_bird_breast_log10Hg, col=Region)) +
	  geom_point() + scale_color_manual(values=c("#56B4E9", "orange")) +
	  stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
	  #xlim(0,160) +
	  theme_bw() + theme(legend.position="none") +
	  labs(title=paste0("stage_1:", neigh_size,"m; P_AK=",format(coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4], nsmall=4, digits=2),
	                    "; P_JP=",format(coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4], nsmall=4, digits=2)
	  ))
	BirdRectrix_vs_envHg1[[neigh_size]] <- ggplot(dat_Hg_env_bird, aes(x=log10Hg_env_mean, y= chunk_bird_retrix_log10Hg, col=Region)) +
	  geom_point() + scale_color_manual(values=c("#56B4E9", "orange")) +
	  stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
	  #xlim(0,160) +
	  theme_bw() + theme(legend.position=c(.1, .2)) +
	  labs(title=paste0("stage_1:", neigh_size,"m; P_AK=",format(coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4], nsmall=4, digits=2),
	                    "; P_JP=",format(coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4], nsmall=4, digits=2)
	  ))
	# stage_2
	dat_Hg_env_bird <- data.frame(rbind(
		cbind(chunk_2_ak_envHg_birdHg[[neigh_size]], Region="AK"), 
		cbind(chunk_2_jp_envHg_birdHg[[neigh_size]], Region="JP")
	))
	BirdBlood_vs_envHg2[[neigh_size]] <- ggplot(dat_Hg_env_bird, aes(x=log10Hg_env_mean, y= chunk_bird_blood_log10Hg, col=Region)) +
	  geom_point() + scale_color_manual(values=c("#56B4E9", "orange")) +
	  stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
	  #xlim(0,160) +
	  theme_bw() + theme(legend.position="none") +
	  labs(title=paste0("stage_2:", neigh_size,"m; P_AK=",format(coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4], nsmall=4, digits=2),
	                    "; P_JP=",format(coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4], nsmall=4, digits=2)
	  ))
	BirdBreast_vs_envHg2[[neigh_size]] <- ggplot(dat_Hg_env_bird, aes(x=log10Hg_env_mean, y= chunk_bird_breast_log10Hg, col=Region)) +
	  geom_point() + scale_color_manual(values=c("#56B4E9", "orange")) +
	  stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
	  #xlim(0,160) +
	  theme_bw() + theme(legend.position="none") +
	  labs(title=paste0("stage_2:", neigh_size,"m; P_AK=",format(coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4], nsmall=4, digits=2),
	                    "; P_JP=",format(coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4], nsmall=4, digits=2)
	  ))
	BirdRectrix_vs_envHg2[[neigh_size]] <- ggplot(dat_Hg_env_bird, aes(x=log10Hg_env_mean, y= chunk_bird_retrix_log10Hg, col=Region)) +
	  geom_point() + scale_color_manual(values=c("#56B4E9", "orange")) +
	  stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
	  #xlim(0,160) +
	  theme_bw() + theme(legend.position="none") +
	  labs(title=paste0("stage_2:", neigh_size,"m; P_AK=",format(coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4], nsmall=4, digits=2),
	                    "; P_JP=",format(coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4], nsmall=4, digits=2)
	  ))
	# stage_3
	dat_Hg_env_bird <- data.frame(rbind(
		cbind(chunk_3_ak_envHg_birdHg[[neigh_size]], Region="AK"), 
		cbind(chunk_3_jp_envHg_birdHg[[neigh_size]], Region="JP")
	))
	BirdBlood_vs_envHg3[[neigh_size]] <- ggplot(dat_Hg_env_bird, aes(x=log10Hg_env_mean, y= chunk_bird_blood_log10Hg, col=Region)) +
	  geom_point() + scale_color_manual(values=c("#56B4E9", "orange")) +
	  stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
	  #xlim(0,160) +
	  theme_bw() + theme(legend.position="none") +
	  labs(title=paste0("stage_3:", neigh_size,"m; P_AK=",format(coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4], nsmall=4, digits=2),
	                    "; P_JP=",format(coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4], nsmall=4, digits=2)
	  ))
	BirdBreast_vs_envHg3[[neigh_size]] <- ggplot(dat_Hg_env_bird, aes(x=log10Hg_env_mean, y= chunk_bird_breast_log10Hg, col=Region)) +
	  geom_point() + scale_color_manual(values=c("#56B4E9", "orange")) +
	  stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
	  #xlim(0,160) +
	  theme_bw() + theme(legend.position="none") +
	  labs(title=paste0("stage_3:", neigh_size,"m; P_AK=",format(coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4], nsmall=4, digits=2),
	                    "; P_JP=",format(coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4], nsmall=4, digits=2)
	  ))
	BirdRectrix_vs_envHg3[[neigh_size]] <- ggplot(dat_Hg_env_bird, aes(x=log10Hg_env_mean, y= chunk_bird_retrix_log10Hg, col=Region)) +
	  geom_point() + scale_color_manual(values=c("#56B4E9", "orange")) +
	  stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
	  #xlim(0,160) +
	  theme_bw() + theme(legend.position="none") +
	  labs(title=paste0("stage_3:", neigh_size,"m; P_AK=",format(coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4], nsmall=4, digits=2),
	                    "; P_JP=",format(coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4], nsmall=4, digits=2)
	  ))

	pdf(paste0("210215_bird_envHg_all_", neigh_size,".pdf"), height=10, width=15)
	grid.arrange(
				 BirdRectrix_vs_envHg1[[neigh_size]], BirdRectrix_vs_envHg2[[neigh_size]], BirdRectrix_vs_envHg3[[neigh_size]], 
				 BirdBreast_vs_envHg1[[neigh_size]], BirdBreast_vs_envHg2[[neigh_size]], BirdBreast_vs_envHg3[[neigh_size]], 
				 BirdBlood_vs_envHg1[[neigh_size]], BirdBlood_vs_envHg2[[neigh_size]], BirdBlood_vs_envHg3[[neigh_size]], 
				 nrow=3)
	dev.off()

	pdf(paste0("210215_bird_envHg_focused_", neigh_size,".pdf"), height=4, width=12)
	grid.arrange(
				 BirdRectrix_vs_envHg1[[neigh_size]], 
				 BirdBreast_vs_envHg2[[neigh_size]], 
				 nrow=1)
	dev.off()


}




####################################################################################################################
####################################################################################################################
save.image("210215_salt_2017_ak_jp_all.RData")


signif_exposure <- data.frame()

for(neigh_size in seq(500, 210500, 10000)){
#for(neigh_size in seq(500, 3000000, 10000)){
	size <- neigh_size
	chunk <- 1
	dat_Hg_env_bird <- data.frame(rbind(cbind(chunk_1_ak_envHg_birdHg[[neigh_size]], Region="AK"), cbind(chunk_1_jp_envHg_birdHg[[neigh_size]], Region="JP")))
	region <- "AK"
	breast <- coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4]
	retrix <- coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4]
	blood <- coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4]
	breast_sgn <- sign(coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,1])
	retrix_sgn <- sign(coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,1])
	blood_sgn <- sign(coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,1])
	signif_exposure <- rbind(signif_exposure, cbind(region, size, chunk, breast, retrix, blood, breast_sgn, retrix_sgn, blood_sgn))
	region <- "JP"
	breast <- coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4]
	retrix <- coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4]
	blood <- coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4]
	breast_sgn <- sign(coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,1])
	retrix_sgn <- sign(coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,1])
	blood_sgn <- sign(coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,1])
	signif_exposure <- rbind(signif_exposure, cbind(region, size, chunk, breast, retrix, blood, breast_sgn, retrix_sgn, blood_sgn))
	chunk <- 2
	dat_Hg_env_bird <- data.frame(rbind(cbind(chunk_2_ak_envHg_birdHg[[neigh_size]], Region="AK"), cbind(chunk_2_jp_envHg_birdHg[[neigh_size]], Region="JP")))
	region <- "AK"
	breast <- coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4]
	retrix <- coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4]
	blood <- coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4]
	breast_sgn <- sign(coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,1])
	retrix_sgn <- sign(coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,1])
	blood_sgn <- sign(coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,1])
	signif_exposure <- rbind(signif_exposure, cbind(region, size, chunk, breast, retrix, blood, breast_sgn, retrix_sgn, blood_sgn))
	region <- "JP"
	breast <- coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4]
	retrix <- coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4]
	blood <- coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4]
	breast_sgn <- sign(coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,1])
	retrix_sgn <- sign(coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,1])
	blood_sgn <- sign(coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,1])
	signif_exposure <- rbind(signif_exposure, cbind(region, size, chunk, breast, retrix, blood, breast_sgn, retrix_sgn, blood_sgn))
	chunk <- 3
	dat_Hg_env_bird <- data.frame(rbind(cbind(chunk_3_ak_envHg_birdHg[[neigh_size]], Region="AK"), cbind(chunk_3_jp_envHg_birdHg[[neigh_size]], Region="JP")))
	region <- "AK"
	breast <- coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4]
	retrix <- coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4]
	blood <- coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,4]
	breast_sgn <- sign(coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,1])
	retrix_sgn <- sign(coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,1])
	blood_sgn <- sign(coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",])))[2,1])
	signif_exposure <- rbind(signif_exposure, cbind(region, size, chunk, breast, retrix, blood, breast_sgn, retrix_sgn, blood_sgn))
	region <- "JP"
	breast <- coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4]
	retrix <- coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4]
	blood <- coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,4]
	breast_sgn <- sign(coef(summary(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,1])
	retrix_sgn <- sign(coef(summary(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,1])
	blood_sgn <- sign(coef(summary(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",])))[2,1])
	signif_exposure <- rbind(signif_exposure, cbind(region, size, chunk, breast, retrix, blood, breast_sgn, retrix_sgn, blood_sgn))
}



signif_exposure[, c(2, 4:9)] <- sapply(signif_exposure[, c(2, 4:9)], as.character)
signif_exposure[, c(2, 4:9)] <- sapply(signif_exposure[, c(2, 4:9)], as.numeric)

signif_exposure_log10 <- signif_exposure
signif_exposure_log10[, c(2)] <- sapply(signif_exposure_log10[, c(2)], function(x){x/1000})
signif_exposure_log10[, c(4:6)] <- sapply(signif_exposure_log10[, c(4:6)], log10)
signif_exposure_log10[, 4] <- -1 * signif_exposure_log10[, 4] * signif_exposure_log10[, 7]
signif_exposure_log10[, 5] <- -1 * signif_exposure_log10[, 5] * signif_exposure_log10[, 8]
signif_exposure_log10[, 6] <- -1 * signif_exposure_log10[, 6] * signif_exposure_log10[, 9]


SI_log10 <- list()

# chunk 1
SI_log10[[1]] <- ggplot(signif_exposure_log10[signif_exposure_log10$chunk == 1,], aes(x=size, y= retrix, col=region)) +
  geom_point(shape=18, col="black") + scale_color_manual(values=c("#56B4E9", "orange")) +
  geom_line(size=2, alpha = 0.7) + geom_hline(yintercept = -log10(0.01), col="black", linetype=2) +
  theme_bw() + theme(legend.position="none") +
  xlab("") + ylab("Exposure vs rectrix [-log10(P)]") + labs(title="(A) -- stage 1 (early winter)")
SI_log10[[2]] <- ggplot(signif_exposure_log10[signif_exposure_log10$chunk == 1,], aes(x=size, y= breast, col=region)) +
  geom_point(shape=18, col="black") + scale_color_manual(values=c("#56B4E9", "orange")) +
  geom_line(size=2, alpha = 0.7) + geom_hline(yintercept = -log10(0.01), col="black", linetype=2) +
  theme_bw() + theme(legend.position="none") +
  xlab("") + ylab("Exposure vs breast [-log10(P)]") + labs(title="(B)")
SI_log10[[3]] <- ggplot(signif_exposure_log10[signif_exposure_log10$chunk == 1,], aes(x=size, y= blood, col=region)) +
  geom_point(shape=18, col="black") + scale_color_manual(values=c("#56B4E9", "orange")) +
  geom_line(size=2, alpha = 0.7) + geom_hline(yintercept = -log10(0.01), col="black", linetype=2) +
  theme_bw() + theme(legend.position="none") +
  xlab("Neigborhood size (km)") + ylab("Exposure vs blood [-log10(P)]") + labs(title="(C)")
# chunk 2
SI_log10[[4]] <- ggplot(signif_exposure_log10[signif_exposure_log10$chunk == 2,], aes(x=size, y= retrix, col=region)) +
  geom_point(shape=18, col="black") + scale_color_manual(values=c("#56B4E9", "orange")) +
  geom_line(size=2, alpha = 0.7) + geom_hline(yintercept = -log10(0.01), col="black", linetype=2) +
  theme_bw() + theme(legend.position="none") +
  xlab("") + ylab("Exposure vs rectrix [-log10(P)]") + labs(title="(D) -- stage 2 (mid-winter)")
SI_log10[[5]] <- ggplot(signif_exposure_log10[signif_exposure_log10$chunk == 2,], aes(x=size, y= breast, col=region)) +
  geom_point(shape=18, col="black") + scale_color_manual(values=c("#56B4E9", "orange")) +
  geom_line(size=2, alpha = 0.7) + geom_hline(yintercept = -log10(0.01), col="black", linetype=2) +
  theme_bw() + theme(legend.position="none") +
  xlab("") + ylab("Exposure vs breast [-log10(P)]") + labs(title="(E)")
SI_log10[[6]] <- ggplot(signif_exposure_log10[signif_exposure_log10$chunk == 2,], aes(x=size, y= blood, col=region)) +
  geom_point(shape=18, col="black") + scale_color_manual(values=c("#56B4E9", "orange")) +
  geom_line(size=2, alpha = 0.7) + geom_hline(yintercept = -log10(0.01), col="black", linetype=2) +
  theme_bw() + theme(legend.position="none") +
  xlab("Neigborhood size (km)") + ylab("Exposure vs blood [-log10(P)]") + labs(title="(F)")
# chunk 3
SI_log10[[7]] <- ggplot(signif_exposure_log10[signif_exposure_log10$chunk == 3,], aes(x=size, y= retrix, col=region)) +
  geom_point(shape=18, col="black") + scale_color_manual(values=c("#56B4E9", "orange")) +
  geom_line(size=2, alpha = 0.7) + geom_hline(yintercept = -log10(0.01), col="black", linetype=2) +
  theme_bw() + theme(legend.position="none") +
  xlab("") + ylab("Exposure vs rectrix [-log10(P)]") + labs(title="(G)  -- stage 3 (late winter)")
SI_log10[[8]] <- ggplot(signif_exposure_log10[signif_exposure_log10$chunk == 3,], aes(x=size, y= breast, col=region)) +
  geom_point(shape=18, col="black") + scale_color_manual(values=c("#56B4E9", "orange")) +
  geom_line(size=2, alpha = 0.7) + geom_hline(yintercept = -log10(0.01), col="black", linetype=2) +
  theme_bw() + theme(legend.position="none") +
  xlab("") + ylab("Exposure vs breast [-log10(P)]") + labs(title="(H)")
SI_log10[[9]] <- ggplot(signif_exposure_log10[signif_exposure_log10$chunk == 3,], aes(x=size, y= blood, col=region)) +
  geom_point(shape=18, col="black") + scale_color_manual(values=c("#56B4E9", "orange")) +
  geom_line(size=2, alpha = 0.7) + geom_hline(yintercept = -log10(0.01), col="black", linetype=2) +
  theme_bw() + theme(legend.position="none") +
  xlab("Neigborhood size (km)") + ylab("Exposure vs blood [-log10(P)]") + labs(title="(I)")


SI_log10[[10]] <- ggplot(signif_exposure_log10[signif_exposure_log10$chunk == 1,], aes(x=size, y= retrix, col=region)) +
  geom_point(shape=18, col="black") + scale_color_manual(values=c("#56B4E9", "orange")) +
  geom_line(size=2, alpha = 0.7) + geom_hline(yintercept = -log10(0.01), col="black", linetype=2) +
  theme_bw() + theme(legend.position="none") +
  xlab("") + ylab("Exposure vs rectrix [-log10(P)]") + labs(title="(A) -- stage 1 (early winter)")
SI_log10[[11]] <- ggplot(signif_exposure_log10[signif_exposure_log10$chunk == 2,], aes(x=size, y= breast, col=region)) +
  geom_point(shape=18, col="black") + scale_color_manual(values=c("#56B4E9", "orange")) +
  geom_line(size=2, alpha = 0.7) + geom_hline(yintercept = -log10(0.01), col="black", linetype=2) +
  theme_bw() + theme(legend.position="none") +
  xlab("Neigborhood size (km)") + ylab("Exposure vs breast [-log10(P)]") + labs(title="(B) -- stage 2 (mid-winter)")





pdf("210215_bird_envHg_significance_all.pdf", height=10, width=15)
grid.arrange(SI_log10[[1]], SI_log10[[4]], SI_log10[[7]],
			 SI_log10[[2]], SI_log10[[5]], SI_log10[[8]],
			 SI_log10[[3]], SI_log10[[6]], SI_log10[[9]],
			 nrow=3)
dev.off()

pdf("210215_bird_envHg_significance_focused.pdf", height=8, width=6)
grid.arrange(SI_log10[[10]], 
			 SI_log10[[11]], 
			 nrow=2)
dev.off()




####################################################################################################################
####################################################################################################################
save.image("210215_salt_2017_ak_jp_all.RData")
####################################################################################################################
####################################################################################################################


####################################################################################################################
####################################################################################################################

											############################
											# back to SID and activity #
											############################

### extracting env Hg data and activity data ###

# AK
yearday <- yday(all_salt_ak$datetime)
all_salt_ak_yday <- cbind(all_salt_ak, yearday)
chunk_sid_1_ak <- all_salt_ak_yday[all_salt_ak_yday$yearday %in% ak_splits[[1]],]
chunk_sid_2_ak <- all_salt_ak_yday[all_salt_ak_yday$yearday %in% ak_splits[[2]],]
chunk_sid_3_ak <- all_salt_ak_yday[all_salt_ak_yday$yearday %in% ak_splits[[3]],]

# JP
yearday <- yday(all_salt_jp$datetime)
all_salt_jp_yday <- cbind(all_salt_jp, yearday)
chunk_sid_1_jp <- all_salt_jp_yday[all_salt_jp_yday$yearday %in% jp_splits[[1]],]
chunk_sid_2_jp <- all_salt_jp_yday[all_salt_jp_yday$yearday %in% jp_splits[[2]],]
chunk_sid_3_jp <- all_salt_jp_yday[all_salt_jp_yday$yearday %in% jp_splits[[3]],]



all_salt_ak
daily_activity_jp

jp_sid_envHg_500000 <- data.frame()
jp_envHg_birdHg_500000 <- rbind(chunk_1_jp_envHg_birdHg_500000, chunk_2_jp_envHg_birdHg_500000, chunk_3_jp_envHg_birdHg_500000)

for(i in 1:length(daily_activity_jp[,1])){
	Bird_ID <- as.character(daily_activity_jp$bird[i])
	Date <- as.character(daily_activity_jp$Date[i])
	propdive <- daily_activity_jp$propdive[i]
	propfly <- daily_activity_jp$propfly[i]
	proprest <- daily_activity_jp$proprest[i]
	env_cur_bird <- jp_envHg_birdHg_500000[jp_envHg_birdHg_500000$bird == Bird_ID,]
	env_pos <- which(env_cur_bird$yearday == yday(Date))
	if(length(env_pos) > 0){
		EnvHg <- env_cur_bird$log10Hg_env_mean[env_pos]
	}else{
		EnvHg <- NA
	}
	jp_sid_envHg_500000 <- rbind(jp_sid_envHg_500000, data.frame(Region="JP", Bird_ID, Date, EnvHg, propdive, propfly, proprest))
}
jp_sid_envHg_500000$Date <- as.Date(jp_sid_envHg_500000$Date)

all_jp_birds <- unique(as.character(jp_sid_envHg_500000$Bird_ID))
all_jp_birds_Hg <- all_jp_birds_fly <- all_jp_birds_rest <- all_jp_birds_dive <- list()
for(i in 1:length(all_jp_birds)){
	cur_df <- jp_sid_envHg_500000[jp_sid_envHg_500000$Bird_ID == all_jp_birds[i],]
	all_jp_birds_Hg[[i]] <- ggplot(cur_df, aes(x=Date, y=EnvHg)) +
							scale_x_date(date_labels = "%b") +
							geom_point(col="orange") + geom_smooth(span = 0.3, col="orange") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_jp_birds[i])
	all_jp_birds_fly[[i]] <- ggplot(cur_df, aes(x=Date, y=propfly)) +
							scale_x_date(date_labels = "%b") +
							geom_point(col="orange") + geom_smooth(span = 0.3, col="orange") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_jp_birds[i])
	all_jp_birds_rest[[i]] <- ggplot(cur_df, aes(x=Date, y=proprest)) +
							scale_x_date(date_labels = "%b") +
							geom_point(col="orange") + geom_smooth(span = 0.3, col="orange") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_jp_birds[i])
	all_jp_birds_dive[[i]] <- ggplot(cur_df, aes(x=Date, y=propdive)) +
							scale_x_date(date_labels = "%b") +
							geom_point(col="orange") + geom_smooth(span = 0.3, col="orange") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_jp_birds[i])
}

pdf("jp_sid_envHg_500000.pdf", height=25, width=20)
grid.arrange(all_jp_birds_Hg[[1]], all_jp_birds_fly[[1]], all_jp_birds_rest[[1]], all_jp_birds_dive[[1]], 
			 all_jp_birds_Hg[[2]], all_jp_birds_fly[[2]], all_jp_birds_rest[[2]], all_jp_birds_dive[[2]], 
			 all_jp_birds_Hg[[3]], all_jp_birds_fly[[3]], all_jp_birds_rest[[3]], all_jp_birds_dive[[3]], 
			 all_jp_birds_Hg[[4]], all_jp_birds_fly[[4]], all_jp_birds_rest[[4]], all_jp_birds_dive[[4]], 
			 all_jp_birds_Hg[[5]], all_jp_birds_fly[[5]], all_jp_birds_rest[[5]], all_jp_birds_dive[[5]], 
			 all_jp_birds_Hg[[6]], all_jp_birds_fly[[6]], all_jp_birds_rest[[6]], all_jp_birds_dive[[6]], 
			 all_jp_birds_Hg[[7]], all_jp_birds_fly[[7]], all_jp_birds_rest[[7]], all_jp_birds_dive[[7]], 
			 all_jp_birds_Hg[[8]], all_jp_birds_fly[[8]], all_jp_birds_rest[[8]], all_jp_birds_dive[[8]], 
			 all_jp_birds_Hg[[9]], all_jp_birds_fly[[9]], all_jp_birds_rest[[9]], all_jp_birds_dive[[9]], 
			 all_jp_birds_Hg[[10]], all_jp_birds_fly[[10]], all_jp_birds_rest[[10]], all_jp_birds_dive[[10]], 
			 nrow=10)
dev.off()

ak_sid_envHg_500000 <- data.frame()
ak_envHg_birdHg_500000 <- rbind(chunk_1_ak_envHg_birdHg_500000, chunk_2_ak_envHg_birdHg_500000, chunk_3_ak_envHg_birdHg_500000)

for(i in 1:length(daily_activity_ak[,1])){
	Bird_ID <- as.character(daily_activity_ak$bird[i])
	Date <- as.character(daily_activity_ak$Date[i])
	propdive <- daily_activity_ak$propdive[i]
	propfly <- daily_activity_ak$propfly[i]
	proprest <- daily_activity_ak$proprest[i]
	env_cur_bird <- ak_envHg_birdHg_500000[ak_envHg_birdHg_500000$bird == Bird_ID,]
	env_pos <- which(env_cur_bird$yearday == yday(Date))
	if(length(env_pos) > 0){
		EnvHg <- env_cur_bird$log10Hg_env_mean[env_pos]
	}else{
		EnvHg <- NA
	}
	ak_sid_envHg_500000 <- rbind(ak_sid_envHg_500000, data.frame(Region="JP", Bird_ID, Date, EnvHg, propdive, propfly, proprest))
}
ak_sid_envHg_500000$Date <- as.Date(ak_sid_envHg_500000$Date)

all_ak_birds <- unique(as.character(ak_sid_envHg_500000$Bird_ID))
all_ak_birds_Hg <- all_ak_birds_fly <- all_ak_birds_rest <- all_ak_birds_dive <- list()
for(i in 1:length(all_ak_birds)){
	cur_df <- ak_sid_envHg_500000[ak_sid_envHg_500000$Bird_ID == all_ak_birds[i],]
	all_ak_birds_Hg[[i]] <- ggplot(cur_df, aes(x=Date, y=EnvHg)) +
							scale_x_date(date_labels = "%b") + 
							geom_point(col="#56B4E9") + geom_smooth(span = 0.3, col="#56B4E9") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_ak_birds[i])
	all_ak_birds_fly[[i]] <- ggplot(cur_df, aes(x=Date, y=propfly)) +
							scale_x_date(date_labels = "%b") + 
							geom_point(col="#56B4E9") + geom_smooth(span = 0.3, col="#56B4E9") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_ak_birds[i])
	all_ak_birds_rest[[i]] <- ggplot(cur_df, aes(x=Date, y=proprest)) +
							scale_x_date(date_labels = "%b") + 
							geom_point(col="#56B4E9") + geom_smooth(span = 0.3, col="#56B4E9") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_ak_birds[i])
	all_ak_birds_dive[[i]] <- ggplot(cur_df, aes(x=Date, y=propdive)) +
							scale_x_date(date_labels = "%b") + 
							geom_point(col="#56B4E9") + geom_smooth(span = 0.3, col="#56B4E9") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_ak_birds[i])
}

pdf("ak_sid_envHg_500000.pdf", height=25, width=20)
grid.arrange(all_ak_birds_Hg[[1]], all_ak_birds_fly[[1]], all_ak_birds_rest[[1]], all_ak_birds_dive[[1]], 
			 all_ak_birds_Hg[[2]], all_ak_birds_fly[[2]], all_ak_birds_rest[[2]], all_ak_birds_dive[[2]], 
			 all_ak_birds_Hg[[3]], all_ak_birds_fly[[3]], all_ak_birds_rest[[3]], all_ak_birds_dive[[3]], 
			 all_ak_birds_Hg[[4]], all_ak_birds_fly[[4]], all_ak_birds_rest[[4]], all_ak_birds_dive[[4]], 
			 all_ak_birds_Hg[[5]], all_ak_birds_fly[[5]], all_ak_birds_rest[[5]], all_ak_birds_dive[[5]], 
			 all_ak_birds_Hg[[6]], all_ak_birds_fly[[6]], all_ak_birds_rest[[6]], all_ak_birds_dive[[6]], 
			 all_ak_birds_Hg[[7]], all_ak_birds_fly[[7]], all_ak_birds_rest[[7]], all_ak_birds_dive[[7]], 
			 all_ak_birds_Hg[[8]], all_ak_birds_fly[[8]], all_ak_birds_rest[[8]], all_ak_birds_dive[[8]], 
			 all_ak_birds_Hg[[9]], all_ak_birds_fly[[9]], all_ak_birds_rest[[9]], all_ak_birds_dive[[9]], 
			 all_ak_birds_Hg[[10]], all_ak_birds_fly[[10]], all_ak_birds_rest[[10]], all_ak_birds_dive[[10]], 
			 nrow=10)
dev.off()




### CCF analysis ###

# AK
dates <- seq(as.Date(min(ak_sid_envHg_500000$Date)), as.Date(max(ak_sid_envHg_500000$Date)), by=1)
ave_ak_sid_envHg_500000 <- data.frame()
for(i in dates){
	print(paste0("Doing ", as.Date(i), " / ", max(ak_sid_envHg_500000$Date)))
	cur_dat <- ak_sid_envHg_500000[ak_sid_envHg_500000$Date == i,]
	cur_dat <- cur_dat[complete.cases(cur_dat), ]
	ave_ak_sid_envHg_500000 <- rbind(ave_ak_sid_envHg_500000, cbind(Region= cur_dat$Region[1], Date=as.Date(i), 
	EnvHg=mean(cur_dat$EnvHg, na.rm=T), 
	propdive=mean(cur_dat$propdive, na.rm=T), propfly=mean(cur_dat$propfly, na.rm=T), proprest=mean(cur_dat$proprest, na.rm=T)))
}
# impute missing means
ave_ak_sid_envHg_500000$EnvHg <- na.kalman(ave_ak_sid_envHg_500000$EnvHg)
ave_ak_sid_envHg_500000$propdive <- na.kalman(ave_ak_sid_envHg_500000$propdive)
ave_ak_sid_envHg_500000$proprest <- na.kalman(ave_ak_sid_envHg_500000$proprest)
ave_ak_sid_envHg_500000$propfly <- na.kalman(ave_ak_sid_envHg_500000$propfly)



# JP
dates <- seq(as.Date(min(jp_sid_envHg_500000$Date)), as.Date(max(jp_sid_envHg_500000$Date)), by=1)
ave_jp_sid_envHg_500000 <- data.frame()
for(i in dates){
	print(paste0("Doing ", as.Date(i), " / ", max(jp_sid_envHg_500000$Date)))
	cur_dat <- jp_sid_envHg_500000[jp_sid_envHg_500000$Date == i,]
	cur_dat <- cur_dat[complete.cases(cur_dat), ]
	ave_jp_sid_envHg_500000 <- rbind(ave_jp_sid_envHg_500000, cbind(Region= cur_dat$Region[1], Date=as.Date(i), 
	EnvHg=mean(cur_dat$EnvHg, na.rm=T), 
	propdive=mean(cur_dat$propdive, na.rm=T), propfly=mean(cur_dat$propfly, na.rm=T), proprest=mean(cur_dat$proprest, na.rm=T)))
}
# impute missing means
ave_jp_sid_envHg_500000$EnvHg <- na.kalman(ave_jp_sid_envHg_500000$EnvHg)
ave_jp_sid_envHg_500000$propdive <- na.kalman(ave_jp_sid_envHg_500000$propdive)
ave_jp_sid_envHg_500000$proprest <- na.kalman(ave_jp_sid_envHg_500000$proprest)
ave_jp_sid_envHg_500000$propfly <- na.kalman(ave_jp_sid_envHg_500000$propfly)


pdf("ccf.pdf", width=15, height=6)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 4),mfrow=c(2,3))
zz <- ccf(ts(ave_ak_sid_envHg_500000$EnvHg), ts(ave_ak_sid_envHg_500000$propfly), type = "correlation", ci = 0.95, cex.axis=0.4, col="#56B4E9", lwd=3)
text(-18, max(zz$acf), "AK propfly")
zz <- ccf(ts(ave_ak_sid_envHg_500000$EnvHg), ts(ave_ak_sid_envHg_500000$proprest), type = "correlation", ci = 0.95, cex.axis=0.4, col="#56B4E9", lwd=3)
text(-18, max(zz$acf), "AK proprest")
zz <- ccf(ts(ave_ak_sid_envHg_500000$EnvHg), ts(ave_ak_sid_envHg_500000$propdive), type = "correlation", ci = 0.95, cex.axis=0.4, col="#56B4E9", lwd=3)
text(-18, max(zz$acf), "AK propdive")

zz <- ccf(ts(ave_jp_sid_envHg_500000$EnvHg), ts(ave_jp_sid_envHg_500000$propfly), type = "correlation", ci = 0.95, cex.axis=0.4, col="orange", lwd=3)
text(-18, max(zz$acf), "JP propfly")
zz <- ccf(ts(ave_jp_sid_envHg_500000$EnvHg), ts(ave_jp_sid_envHg_500000$proprest), type = "correlation", ci = 0.95, cex.axis=0.4, col="orange", lwd=3)
text(-18, .1, "JP proprest")
zz <- ccf(ts(ave_jp_sid_envHg_500000$EnvHg), ts(ave_jp_sid_envHg_500000$propdive), type = "correlation", ci = 0.95, cex.axis=0.4, col="orange", lwd=3)
text(-18, max(zz$acf), "JP propdive")
dev.off()



envHg_activity <- list()
envHg_activity[[1]] <- ggplot(ave_ak_sid_envHg_500000, aes(x= EnvHg, y= propfly)) +
   geom_point(col="#56B4E9") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(propfly ~ EnvHg, data= ave_ak_sid_envHg_500000)))[2,4], nsmall=4, digits=2) ))
envHg_activity[[2]] <- ggplot(ave_ak_sid_envHg_500000, aes(x= EnvHg, y= proprest)) +
   geom_point(col="#56B4E9") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(proprest ~ EnvHg, data= ave_ak_sid_envHg_500000)))[2,4], nsmall=4, digits=2) ))
envHg_activity[[3]] <- ggplot(ave_ak_sid_envHg_500000, aes(x= EnvHg, y= propdive)) +
   geom_point(col="#56B4E9") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(propdive ~ EnvHg, data= ave_ak_sid_envHg_500000)))[2,4], nsmall=4, digits=2) ))

envHg_activity[[4]] <- ggplot(ave_jp_sid_envHg_500000, aes(x= EnvHg, y= propfly)) +
   geom_point(col="orange") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(propfly ~ EnvHg, data= ave_jp_sid_envHg_500000)))[2,4], nsmall=4, digits=2) ))
envHg_activity[[5]] <- ggplot(ave_jp_sid_envHg_500000, aes(x= EnvHg, y= proprest)) +
   geom_point(col="orange") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(proprest ~ EnvHg, data= ave_jp_sid_envHg_500000)))[2,4], nsmall=4, digits=2) ))
envHg_activity[[6]] <- ggplot(ave_jp_sid_envHg_500000, aes(x= EnvHg, y= propdive)) +
   geom_point(col="orange") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(propdive ~ EnvHg, data= ave_jp_sid_envHg_500000)))[2,4], nsmall=4, digits=2) ))

pdf("envHg_acvtivity_500000.pdf", height=8, width=15)
grid.arrange(envHg_activity[[1]], envHg_activity[[2]], envHg_activity[[3]], 
			 envHg_activity[[4]], envHg_activity[[5]], envHg_activity[[6]],  
			 nrow=2)
dev.off()



####################################################################################################################
####################################################################################################################
# 200527 update
											############################
											# back to SID and activity #
											############################

yearday <- c(yday(min(gls_AK_dat$datetime)):365,1:yday(max(gls_AK_dat$datetime)))
ak_splits <- split(yearday, c(rep(1,93), rep(2, 90), rep(3, 59))) # according to Fig.3: mid-Nov--mid-Feb

yearday <- c(yday(min(gls_JP_dat$datetime)):365,1:yday(max(gls_JP_dat$datetime)))
jp_splits <- split(yearday, c(rep(1, 154), rep(2, 46), rep(3, 15))) # according to Fig.3: Jan-Feb

### extracting env Hg data and activity data ###

# AK
yearday <- yday(all_salt_ak$datetime)
all_salt_ak_yday <- cbind(all_salt_ak, yearday)
chunk_sid_1_ak <- all_salt_ak_yday[all_salt_ak_yday$yearday %in% ak_splits[[1]],]
chunk_sid_2_ak <- all_salt_ak_yday[all_salt_ak_yday$yearday %in% ak_splits[[2]],]
chunk_sid_3_ak <- all_salt_ak_yday[all_salt_ak_yday$yearday %in% ak_splits[[3]],]

# JP
yearday <- yday(all_salt_jp$datetime)
all_salt_jp_yday <- cbind(all_salt_jp, yearday)
chunk_sid_1_jp <- all_salt_jp_yday[all_salt_jp_yday$yearday %in% jp_splits[[1]],]
chunk_sid_2_jp <- all_salt_jp_yday[all_salt_jp_yday$yearday %in% jp_splits[[2]],]
chunk_sid_3_jp <- all_salt_jp_yday[all_salt_jp_yday$yearday %in% jp_splits[[3]],]


# 500m
jp_sid_envHg_500000 <- data.frame()
jp_envHg_birdHg_500000 <- rbind(chunk_1_jp_envHg_birdHg_500, chunk_2_jp_envHg_birdHg_500, chunk_3_jp_envHg_birdHg_500)

for(i in 1:length(daily_activity_jp[,1])){
	Bird_ID <- as.character(daily_activity_jp$bird[i])
	Date <- as.character(daily_activity_jp$Date[i])
	propdive <- daily_activity_jp$propdive[i]
	propfly <- daily_activity_jp$propfly[i]
	proprest <- daily_activity_jp$proprest[i]
	env_cur_bird <- jp_envHg_birdHg_500000[jp_envHg_birdHg_500000$bird == Bird_ID,]
	env_pos <- which(env_cur_bird$yearday == yday(Date))
	if(length(env_pos) > 0){
		EnvHg <- env_cur_bird$log10Hg_env_mean[env_pos]
	}else{
		EnvHg <- NA
	}
	jp_sid_envHg_500000 <- rbind(jp_sid_envHg_500000, data.frame(Region="JP", Bird_ID, Date, EnvHg, propdive, propfly, proprest))
}
jp_sid_envHg_500000$Date <- as.Date(jp_sid_envHg_500000$Date)

all_jp_birds <- unique(as.character(jp_sid_envHg_500000$Bird_ID))
all_jp_birds_Hg <- all_jp_birds_fly <- all_jp_birds_rest <- all_jp_birds_dive <- list()
for(i in 1:length(all_jp_birds)){
	cur_df <- jp_sid_envHg_500000[jp_sid_envHg_500000$Bird_ID == all_jp_birds[i],]
	all_jp_birds_Hg[[i]] <- ggplot(cur_df, aes(x=Date, y=EnvHg)) +
							scale_x_date(date_labels = "%b") +
							geom_point(col="orange") + geom_smooth(span = 0.3, col="orange") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_jp_birds[i])
	all_jp_birds_fly[[i]] <- ggplot(cur_df, aes(x=Date, y=propfly)) +
							scale_x_date(date_labels = "%b") +
							geom_point(col="orange") + geom_smooth(span = 0.3, col="orange") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_jp_birds[i])
	all_jp_birds_rest[[i]] <- ggplot(cur_df, aes(x=Date, y=proprest)) +
							scale_x_date(date_labels = "%b") +
							geom_point(col="orange") + geom_smooth(span = 0.3, col="orange") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_jp_birds[i])
	all_jp_birds_dive[[i]] <- ggplot(cur_df, aes(x=Date, y=propdive)) +
							scale_x_date(date_labels = "%b") +
							geom_point(col="orange") + geom_smooth(span = 0.3, col="orange") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_jp_birds[i])
}

pdf("200527_jp_sid_envHg_500.pdf", height=25, width=20)
grid.arrange(all_jp_birds_Hg[[1]], all_jp_birds_fly[[1]], all_jp_birds_rest[[1]], all_jp_birds_dive[[1]], 
			 all_jp_birds_Hg[[2]], all_jp_birds_fly[[2]], all_jp_birds_rest[[2]], all_jp_birds_dive[[2]], 
			 all_jp_birds_Hg[[3]], all_jp_birds_fly[[3]], all_jp_birds_rest[[3]], all_jp_birds_dive[[3]], 
			 all_jp_birds_Hg[[4]], all_jp_birds_fly[[4]], all_jp_birds_rest[[4]], all_jp_birds_dive[[4]], 
			 all_jp_birds_Hg[[5]], all_jp_birds_fly[[5]], all_jp_birds_rest[[5]], all_jp_birds_dive[[5]], 
			 all_jp_birds_Hg[[6]], all_jp_birds_fly[[6]], all_jp_birds_rest[[6]], all_jp_birds_dive[[6]], 
			 all_jp_birds_Hg[[7]], all_jp_birds_fly[[7]], all_jp_birds_rest[[7]], all_jp_birds_dive[[7]], 
			 all_jp_birds_Hg[[8]], all_jp_birds_fly[[8]], all_jp_birds_rest[[8]], all_jp_birds_dive[[8]], 
			 all_jp_birds_Hg[[9]], all_jp_birds_fly[[9]], all_jp_birds_rest[[9]], all_jp_birds_dive[[9]], 
			 all_jp_birds_Hg[[10]], all_jp_birds_fly[[10]], all_jp_birds_rest[[10]], all_jp_birds_dive[[10]], 
			 nrow=10)
dev.off()

ak_sid_envHg_500000 <- data.frame()
ak_envHg_birdHg_500000 <- rbind(chunk_1_ak_envHg_birdHg_500, chunk_2_ak_envHg_birdHg_500, chunk_3_ak_envHg_birdHg_500)

for(i in 1:length(daily_activity_ak[,1])){
	Bird_ID <- as.character(daily_activity_ak$bird[i])
	Date <- as.character(daily_activity_ak$Date[i])
	propdive <- daily_activity_ak$propdive[i]
	propfly <- daily_activity_ak$propfly[i]
	proprest <- daily_activity_ak$proprest[i]
	env_cur_bird <- ak_envHg_birdHg_500000[ak_envHg_birdHg_500000$bird == Bird_ID,]
	env_pos <- which(env_cur_bird$yearday == yday(Date))
	if(length(env_pos) > 0){
		EnvHg <- env_cur_bird$log10Hg_env_mean[env_pos]
	}else{
		EnvHg <- NA
	}
	ak_sid_envHg_500000 <- rbind(ak_sid_envHg_500000, data.frame(Region="JP", Bird_ID, Date, EnvHg, propdive, propfly, proprest))
}
ak_sid_envHg_500000$Date <- as.Date(ak_sid_envHg_500000$Date)

all_ak_birds <- unique(as.character(ak_sid_envHg_500000$Bird_ID))
all_ak_birds_Hg <- all_ak_birds_fly <- all_ak_birds_rest <- all_ak_birds_dive <- list()
for(i in 1:length(all_ak_birds)){
	cur_df <- ak_sid_envHg_500000[ak_sid_envHg_500000$Bird_ID == all_ak_birds[i],]
	all_ak_birds_Hg[[i]] <- ggplot(cur_df, aes(x=Date, y=EnvHg)) +
							scale_x_date(date_labels = "%b") + 
							geom_point(col="#56B4E9") + geom_smooth(span = 0.3, col="#56B4E9") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_ak_birds[i])
	all_ak_birds_fly[[i]] <- ggplot(cur_df, aes(x=Date, y=propfly)) +
							scale_x_date(date_labels = "%b") + 
							geom_point(col="#56B4E9") + geom_smooth(span = 0.3, col="#56B4E9") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_ak_birds[i])
	all_ak_birds_rest[[i]] <- ggplot(cur_df, aes(x=Date, y=proprest)) +
							scale_x_date(date_labels = "%b") + 
							geom_point(col="#56B4E9") + geom_smooth(span = 0.3, col="#56B4E9") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_ak_birds[i])
	all_ak_birds_dive[[i]] <- ggplot(cur_df, aes(x=Date, y=propdive)) +
							scale_x_date(date_labels = "%b") + 
							geom_point(col="#56B4E9") + geom_smooth(span = 0.3, col="#56B4E9") +
							theme_bw() + theme(legend.position="none") +
							labs(title= all_ak_birds[i])
}

pdf("200527_ak_sid_envHg_500.pdf", height=25, width=20)
grid.arrange(all_ak_birds_Hg[[1]], all_ak_birds_fly[[1]], all_ak_birds_rest[[1]], all_ak_birds_dive[[1]], 
			 all_ak_birds_Hg[[2]], all_ak_birds_fly[[2]], all_ak_birds_rest[[2]], all_ak_birds_dive[[2]], 
			 all_ak_birds_Hg[[3]], all_ak_birds_fly[[3]], all_ak_birds_rest[[3]], all_ak_birds_dive[[3]], 
			 all_ak_birds_Hg[[4]], all_ak_birds_fly[[4]], all_ak_birds_rest[[4]], all_ak_birds_dive[[4]], 
			 all_ak_birds_Hg[[5]], all_ak_birds_fly[[5]], all_ak_birds_rest[[5]], all_ak_birds_dive[[5]], 
			 all_ak_birds_Hg[[6]], all_ak_birds_fly[[6]], all_ak_birds_rest[[6]], all_ak_birds_dive[[6]], 
			 all_ak_birds_Hg[[7]], all_ak_birds_fly[[7]], all_ak_birds_rest[[7]], all_ak_birds_dive[[7]], 
			 all_ak_birds_Hg[[8]], all_ak_birds_fly[[8]], all_ak_birds_rest[[8]], all_ak_birds_dive[[8]], 
			 all_ak_birds_Hg[[9]], all_ak_birds_fly[[9]], all_ak_birds_rest[[9]], all_ak_birds_dive[[9]], 
			 all_ak_birds_Hg[[10]], all_ak_birds_fly[[10]], all_ak_birds_rest[[10]], all_ak_birds_dive[[10]], 
			 nrow=10)
dev.off()

# Hg along mig trajectories
mean(ave_jp_sid_envHg_500000$EnvHg) - mean(ave_ak_sid_envHg_500000$EnvHg)
ak_weedays <- as.numeric(unlist(ak_splits))
jp_weedays <- as.numeric(unlist(jp_splits))
tseries <- c(min(ak_weedays[1], jp_weedays[1]):365, 1:max(ak_weedays[length(ak_weedays)], jp_weedays[length(jp_weedays)]))

envHg_500_df <- data.frame(tm = 1:length(ave_ak_sid_envHg_500000$EnvHg),ave_ak_sid_envHg_500000$EnvHg, ave_jp_sid_envHg_500000$EnvHg)

p1 <- ggplot() + theme_bw() + xlab("Time along migration trajectory") +
  ylab("Environmental Hg (log10)") + #scale_x_continuous(limits = c(0, 1)) +
  geom_smooth(data= envHg_500_df, aes(x=tm, y= ave_ak_sid_envHg_500000.EnvHg), method="loess", fill="#56B4E9", colour="#56B4E9", size=1) + 
  geom_point(data= envHg_500_df, aes(x= tm, y= ave_ak_sid_envHg_500000.EnvHg), colour="#56B4E9", size=1, alpha=.2, shape=20) + 
  geom_smooth(data= envHg_500_df, aes(x= tm, y= ave_jp_sid_envHg_500000.EnvHg), method="loess", fill="orange", colour="orange", size=1) + 
  geom_point(data= envHg_500_df, aes(x= tm, y= ave_jp_sid_envHg_500000.EnvHg), colour="orange", size=1, alpha=.2, shape=20) + 
  theme(plot.title = element_text(size=12, hjust=0))

pdf("figs/Hg_along_traj.pdf", height=4, width=6)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1))
grid.arrange(p1,  
			 nrow=1)
dev.off()


### CCF analysis ###

# AK
dates <- seq(as.Date(min(ak_sid_envHg_500000$Date)), as.Date(max(ak_sid_envHg_500000$Date)), by=1)
ave_ak_sid_envHg_500000 <- data.frame()
for(i in dates){
	print(paste0("Doing ", as.Date(i), " / ", max(ak_sid_envHg_500000$Date)))
	cur_dat <- ak_sid_envHg_500000[ak_sid_envHg_500000$Date == i,]
	cur_dat <- cur_dat[complete.cases(cur_dat), ]
	ave_ak_sid_envHg_500000 <- rbind(ave_ak_sid_envHg_500000, cbind(Region= cur_dat$Region[1], Date=as.Date(i), 
	EnvHg=mean(cur_dat$EnvHg, na.rm=T), 
	propdive=mean(cur_dat$propdive, na.rm=T), propfly=mean(cur_dat$propfly, na.rm=T), proprest=mean(cur_dat$proprest, na.rm=T)))
}
# impute missing means
ave_ak_sid_envHg_500000$EnvHg <- na.kalman(ave_ak_sid_envHg_500000$EnvHg)
ave_ak_sid_envHg_500000$propdive <- na.kalman(ave_ak_sid_envHg_500000$propdive)
ave_ak_sid_envHg_500000$proprest <- na.kalman(ave_ak_sid_envHg_500000$proprest)
ave_ak_sid_envHg_500000$propfly <- na.kalman(ave_ak_sid_envHg_500000$propfly)



# JP
dates <- seq(as.Date(min(jp_sid_envHg_500000$Date)), as.Date(max(jp_sid_envHg_500000$Date)), by=1)
ave_jp_sid_envHg_500000 <- data.frame()
for(i in dates){
	print(paste0("Doing ", as.Date(i), " / ", max(jp_sid_envHg_500000$Date)))
	cur_dat <- jp_sid_envHg_500000[jp_sid_envHg_500000$Date == i,]
	cur_dat <- cur_dat[complete.cases(cur_dat), ]
	ave_jp_sid_envHg_500000 <- rbind(ave_jp_sid_envHg_500000, cbind(Region= cur_dat$Region[1], Date=as.Date(i), 
	EnvHg=mean(cur_dat$EnvHg, na.rm=T), 
	propdive=mean(cur_dat$propdive, na.rm=T), propfly=mean(cur_dat$propfly, na.rm=T), proprest=mean(cur_dat$proprest, na.rm=T)))
}
# impute missing means
ave_jp_sid_envHg_500000$EnvHg <- na.kalman(ave_jp_sid_envHg_500000$EnvHg)
ave_jp_sid_envHg_500000$propdive <- na.kalman(ave_jp_sid_envHg_500000$propdive)
ave_jp_sid_envHg_500000$proprest <- na.kalman(ave_jp_sid_envHg_500000$proprest)
ave_jp_sid_envHg_500000$propfly <- na.kalman(ave_jp_sid_envHg_500000$propfly)


pdf("200527_ccf.pdf", width=15, height=6)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 4),mfrow=c(2,3))
zz <- ccf(ts(ave_ak_sid_envHg_500000$EnvHg), ts(ave_ak_sid_envHg_500000$propfly), type = "correlation", ci = 0.95, cex.axis=0.4, col="#56B4E9", lwd=3)
text(-18, max(zz$acf), "AK propfly")
zz <- ccf(ts(ave_ak_sid_envHg_500000$EnvHg), ts(ave_ak_sid_envHg_500000$proprest), type = "correlation", ci = 0.95, cex.axis=0.4, col="#56B4E9", lwd=3)
text(-18, max(zz$acf), "AK proprest")
zz <- ccf(ts(ave_ak_sid_envHg_500000$EnvHg), ts(ave_ak_sid_envHg_500000$propdive), type = "correlation", ci = 0.95, cex.axis=0.4, col="#56B4E9", lwd=3)
text(-18, max(zz$acf), "AK propdive")

zz <- ccf(ts(ave_jp_sid_envHg_500000$EnvHg), ts(ave_jp_sid_envHg_500000$propfly), type = "correlation", ci = 0.95, cex.axis=0.4, col="orange", lwd=3)
text(-18, max(zz$acf), "JP propfly")
zz <- ccf(ts(ave_jp_sid_envHg_500000$EnvHg), ts(ave_jp_sid_envHg_500000$proprest), type = "correlation", ci = 0.95, cex.axis=0.4, col="orange", lwd=3)
text(-18, .1, "JP proprest")
zz <- ccf(ts(ave_jp_sid_envHg_500000$EnvHg), ts(ave_jp_sid_envHg_500000$propdive), type = "correlation", ci = 0.95, cex.axis=0.4, col="orange", lwd=3)
text(-18, max(zz$acf), "JP propdive")
dev.off()



envHg_activity <- list()
envHg_activity[[1]] <- ggplot(ave_ak_sid_envHg_500000, aes(x= EnvHg, y= propfly)) +
   geom_point(col="#56B4E9") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(propfly ~ EnvHg, data= ave_ak_sid_envHg_500000)))[2,4], nsmall=4, digits=2) ))
envHg_activity[[2]] <- ggplot(ave_ak_sid_envHg_500000, aes(x= EnvHg, y= proprest)) +
   geom_point(col="#56B4E9") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(proprest ~ EnvHg, data= ave_ak_sid_envHg_500000)))[2,4], nsmall=4, digits=2) ))
envHg_activity[[3]] <- ggplot(ave_ak_sid_envHg_500000, aes(x= EnvHg, y= propdive)) +
   geom_point(col="#56B4E9") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(propdive ~ EnvHg, data= ave_ak_sid_envHg_500000)))[2,4], nsmall=4, digits=2) ))

envHg_activity[[4]] <- ggplot(ave_jp_sid_envHg_500000, aes(x= EnvHg, y= propfly)) +
   geom_point(col="orange") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(propfly ~ EnvHg, data= ave_jp_sid_envHg_500000)))[2,4], nsmall=4, digits=2) ))
envHg_activity[[5]] <- ggplot(ave_jp_sid_envHg_500000, aes(x= EnvHg, y= proprest)) +
   geom_point(col="orange") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(proprest ~ EnvHg, data= ave_jp_sid_envHg_500000)))[2,4], nsmall=4, digits=2) ))
envHg_activity[[6]] <- ggplot(ave_jp_sid_envHg_500000, aes(x= EnvHg, y= propdive)) +
   geom_point(col="orange") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(propdive ~ EnvHg, data= ave_jp_sid_envHg_500000)))[2,4], nsmall=4, digits=2) ))

pdf("200527_envHg_acvtivity_500.pdf", height=8, width=15)
grid.arrange(envHg_activity[[1]], envHg_activity[[2]], envHg_activity[[3]], 
			 envHg_activity[[4]], envHg_activity[[5]], envHg_activity[[6]],  
			 nrow=2)
dev.off()



####################################################################################################################
####################################################################################################################
####################################################################################################################
# new 200121
# data from:
# https://earthobservatory.nasa.gov/global-maps/MY1DMM_CHLORA/MOD17A2_M_PSN
# and more specifically
# https://neo.sci.gsfc.nasa.gov/view.php?datasetId=MY1DMM_CHLORA&date=2017-05-01
# see also: http://geog.uoregon.edu/bartlein/courses/geog490/week06-geospatial.html
library(plot3D)


PP_Hg_AK_df <- Hg_AK_df
PP_Hg_JP_df <- Hg_JP_df

list_earthData <- list.files("EarthData/", pattern=".CSV")
for(i in list_earthData){
	tmp <- tmp_mat <- cur_date <- curLon <- curLat <- NULL
	tmp <- read.csv(paste0("EarthData/", i), header=T, row.names=1)
	#
	cur_date <- substr(i, 15, 21)
	cur_date <- sub("-", "_", cur_date)
	print(paste0("Now doing ", cur_date))
	curLat <- as.numeric(rownames(tmp))
	curLon <- colnames(tmp)
	curLon <- gsub("X\\.", "-", curLon)
	curLon <- gsub("X", "", curLon)
	curLon <- as.numeric(curLon)
	#
	tmp[tmp == 99999] <- NA
	tmp_mat <- as.matrix(tmp)
	pdf(paste0("EarthData/map_", cur_date, ".pdf"), height=6, width=4)
	image2D(z=tmp_mat, main= cur_date, xlab="", ylab="", resfac=3)
	dev.off()
	# AK
	cur_PP <- NULL
	for(j in 1:length(Hg_AK_df[,1])){
		cur_lat_dist <- cur_lon_dist <- NULL
		cur_lat_dist <- abs(curLat - Hg_AK_df$lat[j])
		closest_lat <- base::which.min(cur_lat_dist)
		cur_lon_dist <- abs(curLon - Hg_AK_df$lon[j])
		closest_lon <- base::which.min(cur_lon_dist)
		cur_PP <- c(cur_PP, tmp_mat[closest_lat, closest_lon])
	}
	PP_Hg_AK_df <- cbind(PP_Hg_AK_df, cur_PP)
	names(PP_Hg_AK_df) <- sub("cur_PP", paste0("PP_", cur_date), names(PP_Hg_AK_df))
	# JP
	cur_PP <- NULL
	for(j in 1:length(Hg_JP_df[,1])){
		cur_lat_dist <- cur_lon_dist <- NULL
		cur_lat_dist <- abs(curLat - Hg_JP_df$lat[j])
		closest_lat <- base::which.min(cur_lat_dist)
		cur_lon_dist <- abs(curLon - Hg_JP_df$lon[j])
		closest_lon <- base::which.min(cur_lon_dist)
		cur_PP <- c(cur_PP, tmp_mat[closest_lat, closest_lon])
	}
	PP_Hg_JP_df <- cbind(PP_Hg_JP_df, cur_PP)
	names(PP_Hg_JP_df) <- sub("cur_PP", paste0("PP_", cur_date), names(PP_Hg_JP_df))
}

PP_columns <- grep("PP_", names(PP_Hg_AK_df))
PPmean_Hg_AK_df <- cbind(PP_Hg_AK_df, PPmean = rowMeans(PP_Hg_AK_df[, PP_columns], na.rm=T))  
PPmean_Hg_AK_df_complete <- PPmean_Hg_AK_df[complete.cases(PPmean),]
plot(log10(PPmean_Hg_AK_df_complete$THg) ~ log10(PPmean_Hg_AK_df_complete$PPmean))
mylm <- lm(log10(PPmean_Hg_AK_df_complete$THg) ~ log10(PPmean_Hg_AK_df_complete$PPmean))
summary(mylm)

PP_columns <- grep("PP_", names(PP_Hg_JP_df))
PPmean_Hg_JP_df <- cbind(PP_Hg_JP_df, PPmean = rowMeans(PP_Hg_JP_df[, PP_columns], na.rm=T))  
PPmean_Hg_JP_df_complete <- PPmean_Hg_JP_df[complete.cases(PPmean),]
plot(log10(PPmean_Hg_JP_df_complete$THg) ~ log10(PPmean_Hg_JP_df_complete$PPmean))
mylm <- lm(log10(PPmean_Hg_JP_df_complete$THg) ~ log10(PPmean_Hg_JP_df_complete$PPmean))
summary(mylm)

envHg_PPmean <- list()
envHg_PPmean[[1]] <- ggplot(PPmean_Hg_AK_df_complete, aes(x= log10(PPmean), y= log10(THg))) +
   geom_point(col="#56B4E9") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(log10(THg) ~ log10(PPmean), data= PPmean_Hg_AK_df_complete)))[2,4], nsmall=4, digits=2) ))
envHg_PPmean[[2]] <- ggplot(PPmean_Hg_JP_df_complete, aes(x= log10(PPmean), y= log10(THg))) +
   geom_point(col="orange") + 
   stat_smooth(method=function(formula,data,weights=weight) rlm(formula,data,weights=weight,method="MM"),fullrange=T) +
   #xlim(0,160) +
   theme_bw() + theme(legend.position="none") +
   labs(title=paste0("P=",format(coef(summary(lmRob(log10(THg) ~ log10(PPmean), data= PPmean_Hg_JP_df_complete)))[2,4], nsmall=4, digits=2) ))

pdf("envHg_PPmean.pdf", height=8, width=5)
grid.arrange(envHg_PPmean[[1]], envHg_PPmean[[2]],  
			 nrow=2)
dev.off()

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

save.image("salt_2017_ak_jp_all.RData")
q(save="no")

####################################################################################################################
####################################################################################################################
####################################################################################################################
# new 200511 -- grid search for potential molting dates

load("salt_2017_ak_jp_all.RData")

# lat-linear trajectories
pdf("MeanTrajectories_LatOnly.pdf",width=6,height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,1))
# AK
plot(ma_ak_lat, type="l", xlab="Time (# of days)", ylab="Latitude", main="AK birds")
points(ma_ak_lat, pch=20, col=diverge_hcl(length(ma_ak_lon), h=c(130,43), c=100, l=c(70,90)))
abline(v=25, col="gray", lwd=2, lty=2)
abline(v=110, col="gray", lwd=2, lty=2)
# JP
plot(ma_jp_lat, type="l", xlab="Time (# of days)", ylab="Latitude", main="JP birds")
points(ma_jp_lat, pch=20, col=diverge_hcl(length(ma_jp_lon), h=c(130,43), c=100, l=c(70,90)))
abline(v=45, col="gray", lwd=2, lty=2)
abline(v=123, col="gray", lwd=2, lty=2)
dev.off()


# let's dive!
all_chunks_ak_envHg_500 <- closestEnvHgReadings(gls_AK_dat_yday, 50)
all_chunks_ak_envHg_birdHg_500 <- chunk_Hg_data(all_chunks_ak_envHg_500, bird_hg_ak)

all_chunks_jp_envHg_500 <- closestEnvHgReadings(gls_JP_dat_yday, 50)
all_chunks_jp_envHg_birdHg_500 <- chunk_Hg_data(all_chunks_jp_envHg_500, bird_hg_jp)

dir.create("gridsearch")

yearday <- c(yday(min(gls_AK_dat$datetime)):365,1:yday(max(gls_AK_dat$datetime)))
ak_splits <- split(yearday, c(rep(1,137), rep(2, 70), rep(3, 35)))
len_ak_splits <- length(unlist(ak_splits))
yearday <- c(yday(min(gls_JP_dat$datetime)):365,1:yday(max(gls_JP_dat$datetime)))
jp_splits <- split(yearday, c(rep(1, 125), rep(2, 75), rep(3, 15)))
len_jp_splits <- length(unlist(jp_splits))




# AK
grid_res_AK <- NULL
grid_res_AK <- foreach(break1 = 1:len_ak_splits, .combine='rbind')%dopar%{
#for(break1 in 1:len_ak_splits){
	for(break2 in break1:len_ak_splits){
		print(paste0("break1: ", break1, " -- break2: ", break2))

		grid_P_AK_blood1 <- grid_AIC_AK_blood1 <- grid_P_AK_breast1 <- grid_AIC_AK_breast1 <- grid_P_AK_retrix1 <- grid_AIC_AK_retrix1 <- NULL
		grid_P_AK_blood2 <- grid_AIC_AK_blood2 <- grid_P_AK_breast2 <- grid_AIC_AK_breast2 <- grid_P_AK_retrix2 <- grid_AIC_AK_retrix2 <- NULL
		grid_P_AK_blood3 <- grid_AIC_AK_blood3 <- grid_P_AK_breast3 <- grid_AIC_AK_breast3 <- grid_P_AK_retrix3 <- grid_AIC_AK_retrix3 <- NULL

		# defining 3-part trips (chunk1, chunk2 and chunk3) in terms of yearday dates
		yearday <- c(yday(min(gls_AK_dat$datetime)):365,1:yday(max(gls_AK_dat$datetime)))
		ak_splits <- split(yearday, c(rep(1, break1), rep(2, (break2-break1)), rep(3, (length(yearday)-break2))))
		
		if(length(ak_splits) == 3){
			#if((length(ak_splits[[1]]) > 12) & (length(ak_splits[[2]]) > 12) & (length(ak_splits[[3]]) > 12)){
				# AK
				chunk_1_ak <- gls_AK_dat_yday[gls_AK_dat_yday$yearday %in% ak_splits[[1]],]
				chunk_2_ak <- gls_AK_dat_yday[gls_AK_dat_yday$yearday %in% ak_splits[[2]],]
				chunk_3_ak <- gls_AK_dat_yday[gls_AK_dat_yday$yearday %in% ak_splits[[3]],]
		
				# 500
				chunk_1_ak_envHg_500 <- all_chunks_ak_envHg_500[all_chunks_ak_envHg_500$yearday %in% ak_splits[[1]],] #closestEnvHgReadings(chunk_1_ak, 500)
				chunk_2_ak_envHg_500 <- all_chunks_ak_envHg_500[all_chunks_ak_envHg_500$yearday %in% ak_splits[[2]],] #closestEnvHgReadings(chunk_2_ak, 500)
				chunk_3_ak_envHg_500 <- all_chunks_ak_envHg_500[all_chunks_ak_envHg_500$yearday %in% ak_splits[[3]],] #closestEnvHgReadings(chunk_3_ak, 500)
		
		
				# 500
				chunk_1_ak_envHg_birdHg_500 <- chunk_Hg_data(chunk_1_ak_envHg_500, bird_hg_ak)
				chunk_2_ak_envHg_birdHg_500 <- chunk_Hg_data(chunk_2_ak_envHg_500, bird_hg_ak)
				chunk_3_ak_envHg_birdHg_500 <- chunk_Hg_data(chunk_3_ak_envHg_500, bird_hg_ak)
		
		
				# chunk_1:3
				dat_Hg_env_bird <- data.frame(rbind(
													cbind(chunk_1_ak_envHg_birdHg_500, Region="AK"), 
													cbind(chunk_2_ak_envHg_birdHg_500, Region="AK"), 
													cbind(chunk_3_ak_envHg_birdHg_500, Region="AK")
													))
				
				# chunk_1
				dat_Hg_env_bird <- data.frame(rbind(cbind(chunk_1_ak_envHg_birdHg_500, Region="AK"), cbind(chunk_1_jp_envHg_birdHg_500, Region="JP")))
				#
				lm_blood <- lm_breast <- lm_retrix <- NULL
				lm_blood <- try(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				lm_breast <- try(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				lm_retrix <- try(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				grid_P_AK_blood1 <- -log10(coef(summary(lm_blood))[2,4])
				grid_AIC_AK_blood1 <- AIC(lm(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				grid_P_AK_breast1 <- -log10(coef(summary(lm_breast))[2,4])
				grid_AIC_AK_breast1 <- AIC(lm(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				grid_P_AK_retrix1 <- -log10(coef(summary(lm_retrix))[2,4])
				grid_AIC_AK_retrix1 <- AIC(lm(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
		
				# chunk_2
				dat_Hg_env_bird <- data.frame(rbind(cbind(chunk_2_ak_envHg_birdHg_500, Region="AK"), cbind(chunk_2_jp_envHg_birdHg_500, Region="JP")))
				#
				lm_blood <- lm_breast <- lm_retrix <- NULL
				lm_blood <- try(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				lm_breast <- try(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				lm_retrix <- try(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				grid_P_AK_blood2 <- -log10(coef(summary(lm_blood))[2,4])
				grid_AIC_AK_blood2 <- AIC(lm(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				grid_P_AK_breast2 <- -log10(coef(summary(lm_breast))[2,4])
				grid_AIC_AK_breast2 <- AIC(lm(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				grid_P_AK_retrix2 <- -log10(coef(summary(lm_retrix))[2,4])
				grid_AIC_AK_retrix2 <- AIC(lm(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				
				# chunk_3
				dat_Hg_env_bird <- data.frame(rbind(cbind(chunk_3_ak_envHg_birdHg_500, Region="AK"), cbind(chunk_3_jp_envHg_birdHg_500, Region="JP")))
				#
				lm_blood <- lm_breast <- lm_retrix <- NULL
				lm_blood <- try(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				lm_breast <- try(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				lm_retrix <- try(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				grid_P_AK_blood3 <- -log10(coef(summary(lm_blood))[2,4])
				grid_AIC_AK_blood3 <- AIC(lm(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				grid_P_AK_breast3 <- -log10(coef(summary(lm_breast))[2,4])
				grid_AIC_AK_breast3 <- AIC(lm(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				grid_P_AK_retrix3 <- -log10(coef(summary(lm_retrix))[2,4])
				grid_AIC_AK_retrix3 <- AIC(lm(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="AK",]))
				
				save(
					#
					grid_P_AK_blood1,
					grid_AIC_AK_blood1,
					grid_P_AK_breast1,
					grid_AIC_AK_breast1,
					grid_P_AK_retrix1,
					grid_AIC_AK_retrix1,
					#
					grid_P_AK_blood2,
					grid_AIC_AK_blood2,
					grid_P_AK_breast2,
					grid_AIC_AK_breast2,
					grid_P_AK_retrix2,
					grid_AIC_AK_retrix2,
					#
					grid_P_AK_blood3,
					grid_AIC_AK_blood3,
					grid_P_AK_breast3,
					grid_AIC_AK_breast3,
					grid_P_AK_retrix3,
					grid_AIC_AK_retrix3,
		
					file=paste0("gridsearch/AK_x", break1, "_y", break2,".RData"))
		
			#}
		}

	}
}

# now putting it all back together
rdat_files <- list.files(path="gridsearch/", pattern="^AK")
map_P_AK_blood1 <- map_AIC_AK_blood1 <- matrix(data=NA, nrow = len_ak_splits, ncol = len_ak_splits)
map_P_AK_blood2 <- map_AIC_AK_blood2 <- matrix(data=NA, nrow = len_ak_splits, ncol = len_ak_splits)
map_P_AK_blood3 <- map_AIC_AK_blood3 <- matrix(data=NA, nrow = len_ak_splits, ncol = len_ak_splits)
map_P_AK_breast1 <- map_AIC_AK_breast1 <- matrix(data=NA, nrow = len_ak_splits, ncol = len_ak_splits)
map_P_AK_breast2 <- map_AIC_AK_breast2 <- matrix(data=NA, nrow = len_ak_splits, ncol = len_ak_splits)
map_P_AK_breast3 <- map_AIC_AK_breast3 <- matrix(data=NA, nrow = len_ak_splits, ncol = len_ak_splits)
map_P_AK_retrix1 <- map_AIC_AK_retrix1 <- matrix(data=NA, nrow = len_ak_splits, ncol = len_ak_splits)
map_P_AK_retrix2 <- map_AIC_AK_retrix2 <- matrix(data=NA, nrow = len_ak_splits, ncol = len_ak_splits)
map_P_AK_retrix3 <- map_AIC_AK_retrix3 <- matrix(data=NA, nrow = len_ak_splits, ncol = len_ak_splits)

for(i in rdat_files){
	load(paste0("gridsearch/", i))
	filename <- sub(".RData", "", i)
	filename <- sub("AK_", "", filename)
	coord <- unlist(base::strsplit(filename, "_"))
	x <- as.numeric(sub("x", "", coord[1]))
	y <- as.numeric(sub("y", "", coord[2]))
	
	# AK (P: upper triangular; AIC: lower triangular)
	map_P_AK_blood1[x,y] <- grid_P_AK_blood1
	map_AIC_AK_blood1[y,x] <- grid_AIC_AK_blood1
	map_P_AK_blood2[x,y] <- grid_P_AK_blood2
	map_AIC_AK_blood2[y,x] <- grid_AIC_AK_blood2
	map_P_AK_blood3[x,y] <- grid_P_AK_blood3
	map_AIC_AK_blood3[y,x] <- grid_AIC_AK_blood3

	map_P_AK_breast1[x,y] <- grid_P_AK_breast1
	map_AIC_AK_breast1[y,x] <- grid_AIC_AK_breast1
	map_P_AK_breast2[x,y] <- grid_P_AK_breast2
	map_AIC_AK_breast2[y,x] <- grid_AIC_AK_breast2
	map_P_AK_breast3[x,y] <- grid_P_AK_breast3
	map_AIC_AK_breast3[y,x] <- grid_AIC_AK_breast3

	map_P_AK_retrix1[x,y] <- grid_P_AK_retrix1
	map_AIC_AK_retrix1[y,x] <- grid_AIC_AK_retrix1
	map_P_AK_retrix2[x,y] <- grid_P_AK_retrix2
	map_AIC_AK_retrix2[y,x] <- grid_AIC_AK_retrix2
	map_P_AK_retrix3[x,y] <- grid_P_AK_retrix3
	map_AIC_AK_retrix3[y,x] <- grid_AIC_AK_retrix3

}

pdf("gridsearch/AK_AIC.pdf", width=18, height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(3,3))
#
map_AIC_AK_blood1_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_AIC_AK_blood1)
image2D(map_AIC_AK_blood1_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: blood chunk 1", clab = "AIC")
map_AIC_AK_blood2_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_AIC_AK_blood2)
image2D(map_AIC_AK_blood2_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: blood chunk 2", clab = "AIC")
map_AIC_AK_blood3_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_AIC_AK_blood3)
image2D(map_AIC_AK_blood3_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: blood chunk 3", clab = "AIC")
#
map_AIC_AK_breast1_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_AIC_AK_breast1)
image2D(map_AIC_AK_breast1_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: breast chunk 1", clab = "AIC")
map_AIC_AK_breast2_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_AIC_AK_breast2)
image2D(map_AIC_AK_breast2_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: breast chunk 2", clab = "AIC")
map_AIC_AK_breast3_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_AIC_AK_breast3)
image2D(map_AIC_AK_breast3_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: breast chunk 3", clab = "AIC")
#
map_AIC_AK_retrix1_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_AIC_AK_retrix1)
image2D(map_AIC_AK_retrix1_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: retrix chunk 1", clab = "AIC")
map_AIC_AK_retrix2_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_AIC_AK_retrix2)
image2D(map_AIC_AK_retrix2_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: retrix chunk 2", clab = "AIC")
map_AIC_AK_retrix3_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_AIC_AK_retrix3)
image2D(map_AIC_AK_retrix3_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: retrix chunk 3", clab = "AIC")
#
dev.off()


pdf("gridsearch/AK_P.pdf", width=14.5, height=13)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(3,3))
#
map_P_AK_blood1_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_P_AK_blood1)
image2D(map_P_AK_blood1_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: blood chunk 1", clab = "-log10 P")
map_P_AK_blood2_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_P_AK_blood2)
image2D(map_P_AK_blood2_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: blood chunk 2", clab = "-log10 P")
map_P_AK_blood3_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_P_AK_blood3)
image2D(map_P_AK_blood3_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: blood chunk 3", clab = "-log10 P")
#
map_P_AK_breast1_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_P_AK_breast1)
image2D(map_P_AK_breast1_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: breast chunk 1", clab = "-log10 P")
map_P_AK_breast2_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_P_AK_breast2)
image2D(map_P_AK_breast2_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: breast chunk 2", clab = "-log10 P")
map_P_AK_breast3_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_P_AK_breast3)
image2D(map_P_AK_breast3_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: breast chunk 3", clab = "-log10 P")
#
map_P_AK_retrix1_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_P_AK_retrix1)
image2D(map_P_AK_retrix1_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: retrix chunk 1", clab = "-log10 P")
map_P_AK_retrix2_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_P_AK_retrix2)
image2D(map_P_AK_retrix2_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: retrix chunk 2", clab = "-log10 P")
map_P_AK_retrix3_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_P_AK_retrix3)
image2D(map_P_AK_retrix3_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: retrix chunk 3", clab = "-log10 P")
#
dev.off()

pdf("gridsearch/AK_P_chunk2.pdf", width=14.5, height=4.33)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,3))
#
map_P_AK_breast2_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_P_AK_breast2)
image2D(map_P_AK_breast2_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: breast chunk 2", clab = "-log10 P")
#
map_P_AK_retrix2_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_P_AK_retrix2)
image2D(map_P_AK_retrix2_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: retrix chunk 2", clab = "-log10 P")
#
map_P_AK_blood2_lst <- list(x= 1:len_ak_splits, y= 1:len_ak_splits, z= map_P_AK_blood2)
image2D(map_P_AK_blood2_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "AK: blood chunk 2", clab = "-log10 P")
#
dev.off()





# JP
grid_res_JP <- NULL
grid_res_JP <- foreach(break1 = 1:len_jp_splits, .combine='rbind')%dopar%{
#for(break1 in 1:len_jp_splits){
	for(break2 in break1:len_jp_splits){
		print(paste0("break1: ", break1, " -- break2: ", break2))

		grid_P_JP_blood1 <- grid_AIC_JP_blood1 <- grid_P_JP_breast1 <- grid_AIC_JP_breast1 <- grid_P_JP_retrix1 <- grid_AIC_JP_retrix1 <- NA
		grid_P_JP_blood2 <- grid_AIC_JP_blood2 <- grid_P_JP_breast2 <- grid_AIC_JP_breast2 <- grid_P_JP_retrix2 <- grid_AIC_JP_retrix2 <- NA
		grid_P_JP_blood3 <- grid_AIC_JP_blood3 <- grid_P_JP_breast3 <- grid_AIC_JP_breast3 <- grid_P_JP_retrix3 <- grid_AIC_JP_retrix3 <- NA

		# defining 3-part trips (chunk1, chunk2 and chunk3) in terms of yearday dates
		yearday <- c(yday(min(gls_JP_dat$datetime)):365,1:yday(max(gls_JP_dat$datetime)))
		jp_splits <- split(yearday, c(rep(1, break1), rep(2, (break2-break1)), rep(3, (length(yearday)-break2))))
		
		if(length(jp_splits) == 3){
			#if((length(jp_splits[[1]]) > 12) & (length(jp_splits[[2]]) > 12) & (length(jp_splits[[3]]) > 12)){
				# JP
				chunk_1_jp <- gls_JP_dat_yday[gls_JP_dat_yday$yearday %in% jp_splits[[1]],]
				chunk_2_jp <- gls_JP_dat_yday[gls_JP_dat_yday$yearday %in% jp_splits[[2]],]
				chunk_3_jp <- gls_JP_dat_yday[gls_JP_dat_yday$yearday %in% jp_splits[[3]],]
		
				# 500
				chunk_1_jp_envHg_500 <- all_chunks_jp_envHg_500[all_chunks_jp_envHg_500$yearday %in% jp_splits[[1]],] #closestEnvHgReadings(chunk_1_jp, 500)
				chunk_2_jp_envHg_500 <- all_chunks_jp_envHg_500[all_chunks_jp_envHg_500$yearday %in% jp_splits[[2]],] #closestEnvHgReadings(chunk_2_jp, 500)
				chunk_3_jp_envHg_500 <- all_chunks_jp_envHg_500[all_chunks_jp_envHg_500$yearday %in% jp_splits[[3]],] #closestEnvHgReadings(chunk_3_jp, 500)
		
		
				# 500
				chunk_1_jp_envHg_birdHg_500 <- chunk_Hg_data(chunk_1_jp_envHg_500, bird_hg_jp)
				chunk_2_jp_envHg_birdHg_500 <- chunk_Hg_data(chunk_2_jp_envHg_500, bird_hg_jp)
				chunk_3_jp_envHg_birdHg_500 <- chunk_Hg_data(chunk_3_jp_envHg_500, bird_hg_jp)
		
		
				# chunk_1:3
				if(dim(chunk_2_jp_envHg_birdHg_500)[1] > 0){
					dat_Hg_env_bird <- data.frame(rbind(
														cbind(chunk_1_jp_envHg_birdHg_500, Region="JP"), 
														cbind(chunk_2_jp_envHg_birdHg_500, Region="JP"), 
														cbind(chunk_3_jp_envHg_birdHg_500, Region="JP")
														))
					
					# chunk_1
					dat_Hg_env_bird <- data.frame(rbind(cbind(chunk_1_jp_envHg_birdHg_500, Region="JP"), cbind(chunk_1_jp_envHg_birdHg_500, Region="JP")))
					#
					lm_blood <- lm_breast <- lm_retrix <- NULL
					if(length(dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]$chunk_bird_blood_log10Hg) > 6){
						lm_blood <- try(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
						if(coef(summary(lm_blood))[2,1]){
							grid_P_JP_blood1 <- -log10(coef(summary(lm_blood))[2,4])
						}
					}
					if(length(dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]$chunk_bird_breast_log10Hg) > 6){
						lm_breast <- try(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
						if(coef(summary(lm_breast))[2,1]){
							grid_P_JP_breast1 <- -log10(coef(summary(lm_breast))[2,4])
						}
					}
					if(length(dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]$chunk_bird_retrix_log10Hg) > 6){
						lm_retrix <- try(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
						if(coef(summary(lm_retrix))[2,1]){
							grid_P_JP_retrix1 <- -log10(coef(summary(lm_retrix))[2,4])
						}
					}
					grid_AIC_JP_blood1 <- AIC(lm(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
					grid_AIC_JP_breast1 <- AIC(lm(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
					grid_AIC_JP_retrix1 <- AIC(lm(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
			
					# chunk_2
					dat_Hg_env_bird <- data.frame(rbind(cbind(chunk_2_jp_envHg_birdHg_500, Region="JP"), cbind(chunk_2_jp_envHg_birdHg_500, Region="JP")))
					#
					lm_blood <- lm_breast <- lm_retrix <- NULL
					if(length(dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]$chunk_bird_blood_log10Hg) > 6){
						lm_blood <- try(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
						if(coef(summary(lm_blood))[2,1]){
							grid_P_JP_blood2 <- -log10(coef(summary(lm_blood))[2,4])
						}
					}
					if(length(dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]$chunk_bird_breast_log10Hg) > 6){
						lm_breast <- try(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
						if(coef(summary(lm_breast))[2,1]){
							grid_P_JP_breast2 <- -log10(coef(summary(lm_breast))[2,4])
						}
					}
					if(length(dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]$chunk_bird_retrix_log10Hg) > 6){
						lm_retrix <- try(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
						if(coef(summary(lm_retrix))[2,1]){
							grid_P_JP_retrix2 <- -log10(coef(summary(lm_retrix))[2,4])
						}
					}
					grid_AIC_JP_blood2 <- AIC(lm(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
					grid_AIC_JP_breast2 <- AIC(lm(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
					grid_AIC_JP_retrix2 <- AIC(lm(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
					
					# chunk_3
					dat_Hg_env_bird <- data.frame(rbind(cbind(chunk_3_jp_envHg_birdHg_500, Region="JP"), cbind(chunk_3_jp_envHg_birdHg_500, Region="JP")))
					#
					lm_blood <- lm_breast <- lm_retrix <- NULL
					if(length(dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]$chunk_bird_blood_log10Hg) > 6){
						lm_blood <- try(lmRob(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
						if(coef(summary(lm_blood))[2,1]){
							grid_P_JP_blood3 <- -log10(coef(summary(lm_blood))[2,4])
						}
					}
					if(length(dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]$chunk_bird_breast_log10Hg) > 6){
						lm_breast <- try(lmRob(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
						if(coef(summary(lm_breast))[2,1]){
							grid_P_JP_breast3 <- -log10(coef(summary(lm_breast))[2,4])
						}
					}
					if(length(dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]$chunk_bird_retrix_log10Hg) > 6){
						lm_retrix <- try(lmRob(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
						if(coef(summary(lm_retrix))[2,1]){
							grid_P_JP_retrix3 <- -log10(coef(summary(lm_retrix))[2,4])
						}
					}
					grid_AIC_JP_blood3 <- AIC(lm(chunk_bird_blood_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
					grid_AIC_JP_breast3 <- AIC(lm(chunk_bird_breast_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
					grid_AIC_JP_retrix3 <- AIC(lm(chunk_bird_retrix_log10Hg ~ log10Hg_env_mean, data= dat_Hg_env_bird[dat_Hg_env_bird$Region=="JP",]))
					
					save(
						#
						grid_P_JP_blood1,
						grid_AIC_JP_blood1,
						grid_P_JP_breast1,
						grid_AIC_JP_breast1,
						grid_P_JP_retrix1,
						grid_AIC_JP_retrix1,
						#
						grid_P_JP_blood2,
						grid_AIC_JP_blood2,
						grid_P_JP_breast2,
						grid_AIC_JP_breast2,
						grid_P_JP_retrix2,
						grid_AIC_JP_retrix2,
						#
						grid_P_JP_blood3,
						grid_AIC_JP_blood3,
						grid_P_JP_breast3,
						grid_AIC_JP_breast3,
						grid_P_JP_retrix3,
						grid_AIC_JP_retrix3,
			
						file=paste0("gridsearch/JP_x", break1, "_y", break2,".RData"))
				}
		
			#}
		}

	}
}

# now putting it all back together
rdat_files <- list.files(path="gridsearch/", pattern="^JP")
map_P_JP_blood1 <- map_AIC_JP_blood1 <- matrix(data=NA, nrow = len_jp_splits, ncol = len_jp_splits)
map_P_JP_blood2 <- map_AIC_JP_blood2 <- matrix(data=NA, nrow = len_jp_splits, ncol = len_jp_splits)
map_P_JP_blood3 <- map_AIC_JP_blood3 <- matrix(data=NA, nrow = len_jp_splits, ncol = len_jp_splits)
map_P_JP_breast1 <- map_AIC_JP_breast1 <- matrix(data=NA, nrow = len_jp_splits, ncol = len_jp_splits)
map_P_JP_breast2 <- map_AIC_JP_breast2 <- matrix(data=NA, nrow = len_jp_splits, ncol = len_jp_splits)
map_P_JP_breast3 <- map_AIC_JP_breast3 <- matrix(data=NA, nrow = len_jp_splits, ncol = len_jp_splits)
map_P_JP_retrix1 <- map_AIC_JP_retrix1 <- matrix(data=NA, nrow = len_jp_splits, ncol = len_jp_splits)
map_P_JP_retrix2 <- map_AIC_JP_retrix2 <- matrix(data=NA, nrow = len_jp_splits, ncol = len_jp_splits)
map_P_JP_retrix3 <- map_AIC_JP_retrix3 <- matrix(data=NA, nrow = len_jp_splits, ncol = len_jp_splits)

for(i in rdat_files){
	load(paste0("gridsearch/", i))
	filename <- sub(".RData", "", i)
	filename <- sub("JP_", "", filename)
	coord <- unlist(base::strsplit(filename, "_"))
	x <- as.numeric(sub("x", "", coord[1]))
	y <- as.numeric(sub("y", "", coord[2]))
	
	# JP (P: upper triangular; AIC: lower triangular)
	map_P_JP_blood1[x,y] <- grid_P_JP_blood1
	map_AIC_JP_blood1[y,x] <- grid_AIC_JP_blood1
	map_P_JP_blood2[x,y] <- grid_P_JP_blood2
	map_AIC_JP_blood2[y,x] <- grid_AIC_JP_blood2
	map_P_JP_blood3[x,y] <- grid_P_JP_blood3
	map_AIC_JP_blood3[y,x] <- grid_AIC_JP_blood3

	map_P_JP_breast1[x,y] <- grid_P_JP_breast1
	map_AIC_JP_breast1[y,x] <- grid_AIC_JP_breast1
	map_P_JP_breast2[x,y] <- grid_P_JP_breast2
	map_AIC_JP_breast2[y,x] <- grid_AIC_JP_breast2
	map_P_JP_breast3[x,y] <- grid_P_JP_breast3
	map_AIC_JP_breast3[y,x] <- grid_AIC_JP_breast3

	map_P_JP_retrix1[x,y] <- grid_P_JP_retrix1
	map_AIC_JP_retrix1[y,x] <- grid_AIC_JP_retrix1
	map_P_JP_retrix2[x,y] <- grid_P_JP_retrix2
	map_AIC_JP_retrix2[y,x] <- grid_AIC_JP_retrix2
	map_P_JP_retrix3[x,y] <- grid_P_JP_retrix3
	map_AIC_JP_retrix3[y,x] <- grid_AIC_JP_retrix3

}

pdf("gridsearch/JP_AIC.pdf", width=18, height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(3,3))
#
map_AIC_JP_blood1_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_AIC_JP_blood1)
image2D(map_AIC_JP_blood1_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: blood chunk 1", clab = "AIC")
map_AIC_JP_blood2_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_AIC_JP_blood2)
image2D(map_AIC_JP_blood2_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: blood chunk 2", clab = "AIC")
map_AIC_JP_blood3_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_AIC_JP_blood3)
image2D(map_AIC_JP_blood3_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: blood chunk 3", clab = "AIC")
#
map_AIC_JP_breast1_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_AIC_JP_breast1)
image2D(map_AIC_JP_breast1_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: breast chunk 1", clab = "AIC")
map_AIC_JP_breast2_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_AIC_JP_breast2)
image2D(map_AIC_JP_breast2_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: breast chunk 2", clab = "AIC")
map_AIC_JP_breast3_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_AIC_JP_breast3)
image2D(map_AIC_JP_breast3_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: breast chunk 3", clab = "AIC")
#
map_AIC_JP_retrix1_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_AIC_JP_retrix1)
image2D(map_AIC_JP_retrix1_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: retrix chunk 1", clab = "AIC")
map_AIC_JP_retrix2_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_AIC_JP_retrix2)
image2D(map_AIC_JP_retrix2_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: retrix chunk 2", clab = "AIC")
map_AIC_JP_retrix3_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_AIC_JP_retrix3)
image2D(map_AIC_JP_retrix3_lst, xlab = "Breakpoint 2 along wintering tracks", ylab = "Breakpoint 1 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: retrix chunk 3", clab = "AIC")
#
dev.off()


pdf("gridsearch/JP_P.pdf", width=14.5, height=13)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(3,3))
#
map_P_JP_blood1_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_P_JP_blood1)
image2D(map_P_JP_blood1_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: blood chunk 1", clab = "-log10 P")
map_P_JP_blood2_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_P_JP_blood2)
image2D(map_P_JP_blood2_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: blood chunk 2", clab = "-log10 P")
map_P_JP_blood3_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_P_JP_blood3)
image2D(map_P_JP_blood3_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: blood chunk 3", clab = "-log10 P")
#
map_P_JP_breast1_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_P_JP_breast1)
image2D(map_P_JP_breast1_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: breast chunk 1", clab = "-log10 P")
map_P_JP_breast2_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_P_JP_breast2)
image2D(map_P_JP_breast2_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: breast chunk 2", clab = "-log10 P")
map_P_JP_breast3_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_P_JP_breast3)
image2D(map_P_JP_breast3_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: breast chunk 3", clab = "-log10 P")
#
map_P_JP_retrix1_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_P_JP_retrix1)
image2D(map_P_JP_retrix1_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: retrix chunk 1", clab = "-log10 P")
map_P_JP_retrix2_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_P_JP_retrix2)
image2D(map_P_JP_retrix2_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: retrix chunk 2", clab = "-log10 P")
map_P_JP_retrix3_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_P_JP_retrix3)
image2D(map_P_JP_retrix3_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: retrix chunk 3", clab = "-log10 P")
#
dev.off()

pdf("gridsearch/JP_P2.pdf", width=14.5, height=4.33)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,3))
#
map_P_JP_breast2_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_P_JP_breast2)
image2D(map_P_JP_breast2_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: breast chunk 2", clab = "-log10 P")
#
map_P_JP_retrix2_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_P_JP_retrix2)
image2D(map_P_JP_retrix2_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: retrix chunk 2", clab = "-log10 P")
#
map_P_JP_blood2_lst <- list(x= 1:len_jp_splits, y= 1:len_jp_splits, z= map_P_JP_blood2)
image2D(map_P_JP_blood2_lst, xlab = "Breakpoint 1 along wintering tracks", ylab = "Breakpoint 2 along wintering tracks", contour = list(levels = 0, col = "black", lwd = 1), shade = 0.1, main = "JP: blood chunk 2", clab = "-log10 P")
#
dev.off()






####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

save.image("salt_2017_ak_jp_all_grid.RData")
q(save="no")

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################


