library(locfit)
library(car)
library(animation)
library(xts)
library(circlize)
library(segmented)
library(modes)
library(datamart)
library(tree)
library(adabag)
library(corrplot)


# NB: cleanup sample names in econn dat file
# replace: "" with _
# replace: '' with _




												###############
												# E-conn data #
												###############

assay <- "Assay_A"
datafile <- "../data/Assay_A_0316_0614.csv"


# no edits required from here
system(paste0("rm -fr ", assay))
system(paste0("mkdir ", assay))
assay_dir <- paste0(assay, "/")
econn_dat <- read.csv(datafile)
head(econn_dat)
length(econn_dat[,1])
length(econn_dat[,1][econn_dat$Assay == assay])
save.image(paste0(assay_dir, "econn_cust_", assay, ".RData"))

# sort data by datetime
econn_dat_ori <- econn_dat
econn_dat <- econn_dat_ori[order(as.POSIXct(strptime(econn_dat_ori$Time.Metering, format = "%Y-%m-%d %H:%M:%S"), tz = "EST")),]


# timing
datetime <- as.POSIXct(strptime(econn_dat$Time.Metering, format = "%Y-%m-%d %H:%M:%S"), tz = "EST")
day_of_month <- as.numeric(format(datetime, "%d"))
hour_of_day <- as.numeric(format(datetime, "%I"))
time_of_day <- as.numeric(format(datetime, "%H"))
weekday <- weekdays(datetime)

which(is.na(datetime))
datetime[which(is.na(datetime))[1]]
econn_dat$Time.Metering[which(is.na(datetime))[1]]


min(datetime)
max(datetime)
max(datetime) - min(datetime) # 90.98508 days

pdf(paste0(assay_dir, assay, "_timing.pdf"), width=18, height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,3))
plot(datetime, pch=".", xlab="Ordered distribution of dates/times", ylab="Dates/times")
plot(day_of_month, datetime, pch=".")
plot(time_of_day, datetime, pch=".")
hist(time_of_day,24, xlab="Time of day (0-24h)", main="")
hist(day_of_month)
plot(table(weekday), type="h")
dev.off()

save.image(paste0(assay_dir, "econn_cust_", assay, ".RData"))

###########################
# split data at mid point #
###########################
midpoint_date <- min(datetime) + (max(datetime) - min(datetime))/2
econn_dat <- econn_dat[datetime > midpoint_date, ]

# sort data by datetime
econn_dat_ori <- econn_dat
econn_dat <- econn_dat_ori[order(as.POSIXct(strptime(econn_dat_ori$Time.Metering, format = "%Y-%m-%d %H:%M:%S"), tz = "EST")),]


# timing
datetime <- as.POSIXct(strptime(econn_dat$Time.Metering, format = "%Y-%m-%d %H:%M:%S"), tz = "EST")
day_of_month <- as.numeric(format(datetime, "%d"))
hour_of_day <- as.numeric(format(datetime, "%I"))
time_of_day <- as.numeric(format(datetime, "%H"))
weekday <- weekdays(datetime)

which(is.na(datetime))
datetime[which(is.na(datetime))[1]]
econn_dat$Time.Metering[which(is.na(datetime))[1]]


min(datetime)
max(datetime)
max(datetime) - min(datetime) # 90.98508 days

pdf(paste0(assay_dir, assay, "_timing_split.pdf"), width=18, height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,3))
plot(datetime, pch=".", xlab="Ordered distribution of dates/times", ylab="Dates/times")
plot(day_of_month, datetime, pch=".")
plot(time_of_day, datetime, pch=".")
hist(time_of_day,24, xlab="Time of day (0-24h)", main="")
hist(day_of_month)
plot(table(weekday), type="h")
dev.off()

save.image(paste0(assay_dir, "econn_cust_", assay, ".RData"))

												#################
												# Customer data #
												#################

cust_data <- read.csv("../reagent_complaints_012016_062016.csv")
head(cust_data)

# prune all Z* related call areas
cust_data <- cust_data[-grep("^Z", cust_data$Call_Area),]

# groups call areas
cust_data$Call_Area[cust_data$Call_Area == "QCDH"] <- as.factor("QCH")
cust_data$Call_Area[cust_data$Call_Area == "QCH"] <- as.factor("QCH")
cust_data$Call_Area[cust_data$Call_Area == "QCSH"] <- as.factor("QCH")
#
cust_data$Call_Area[cust_data$Call_Area == "QCDL"] <- as.factor("QCL")
cust_data$Call_Area[cust_data$Call_Area == "QCL"] <- as.factor("QCL")
cust_data$Call_Area[cust_data$Call_Area == "QCSL"] <- as.factor("QCL")

datetime_cust <- as.POSIXct(strptime(cust_data$Create_Audit_DT, format = "%Y-%m-%d %H:%M:%S"), tz = "EST")
min(datetime_cust)
max(datetime_cust)
max(datetime_cust) - min(datetime_cust) # 165.4125 days

# sort data by datetime
cust_data_ori <- cust_data
cust_data <- cust_data_ori[order(as.POSIXct(strptime(cust_data_ori$Create_Audit_DT, format = "%Y-%m-%d %H:%M:%S"), tz = "EST")),]
datetime_cust <- as.POSIXct(strptime(cust_data$Create_Audit_DT, format = "%Y-%m-%d %H:%M:%S"), tz = "EST")

# match customer data with E-conn data
cust_data_match <- cust_data[(datetime_cust > min(datetime)),]
datetime_cust <- as.POSIXct(strptime(cust_data_match$Create_Audit_DT, format = "%Y-%m-%d %H:%M:%S"), tz = "EST")
pdf(paste0(assay_dir, assay, "_timing_cust.pdf"), width=24, height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,2))
plot(datetime_cust, pch=20, col="red", main="Matching dates -- Customer data")
abline(h=datetime_cust[weekdays(datetime_cust) == "Sunday"], col="gray", lty=2)
plot(datetime, pch=20, col="blue", main="Matching dates -- Econn data")
dev.off()


save.image(paste0(assay_dir, "econn_cust_", assay, ".RData"))


												#################################
												# Error/resolution distribution #
												#################################


# list of error codes for DGXN:
sort(unique(cust_data_match$Call_Area))
err <- sort(unique(cust_data_match$Call_Area[cust_data_match$Call_Subject == assay]))
err

err_df <- list()
err_df_raw <- list()
err_name <- c()
pdf(paste0(assay_dir, assay, "_dist_errors.pdf"), width=24, height=16)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(4,4))
for(i in 1:length(err)){
	plot(datetime_cust, cumsum(cust_data_match$Call_Area == err[i]), main=err[i], xlab="Date/time", ylab="Cummulative sum", pch=20)
	abline(lm(cumsum(cust_data_match$Call_Area == err[i]) ~ datetime_cust), col="red")
	err_df[[i]] <- cumsum(cust_data_match$Call_Area == err[i])
	err_df_raw[[i]] <- (cust_data_match$Call_Area == err[i])
	print(length(cumsum(cust_data_match$Call_Area == err[i])))
	err_name[i] <- as.character(err[i])
}
dev.off()

err_df <- matrix(unlist(err_df), ncol = length(err), byrow = F)
colnames(err_df) <- err_name
err_df_raw <- matrix(unlist(err_df_raw), ncol = length(err), byrow = F)
colnames(err_df_raw) <- err_name

pdf(paste0(assay_dir, assay, "_dist_corr.pdf"), width=12, height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
err_M <- cor(err_df)
corrplot.mixed(err_M)
dev.off()


#pdf(paste0(assay_dir, assay, "_ccf_errors.pdf"), width=6*length(err_name), height=4*length(err_name))
#par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(length(err_name),length(err_name)))
#for(i in 1:length(err_name)){
#	for(j in 1:length(err_name)){
#		# bin data by day
#		time_interval <- 24*60*60 # in seconds
#		ts1 <- err_df_raw[,i]
#		ts1z <- zoo(ts1, datetime_cust)
#		ts1z_aggr <- aggregate(ts1z, time(ts1z) - as.numeric(time(ts1z)) %% time_interval, mean)
#		
#		ts2 <- err_df_raw[,j]
#		ts2z <- zoo(ts2, datetime_cust)
#		ts2z_aggr <- aggregate(ts2z, time(ts2z) - as.numeric(time(ts2z)) %% time_interval, mean)
#		
#		ccf_obj <- ccf(ts1z_aggr, ts2z_aggr, xlab="Lag (in seconds)", ylim=c(-1,1), lwd=2)
#		#text(-1e6, max(ccf_obj$acf), paste0(err_name[i], " vs. ", err_name[j]), font=2)
#		text(-1e6, 1, paste0(err_name[i], " vs. ", err_name[j]), font=2)
#	}
#}
#dev.off()


# pdf(paste0(assay_dir, assay, "_ccf_errors_ex.pdf"), width=12, height=4)
# par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,2))
# i <- 1
# j <- 15
# ts1 <- err_df_raw[,i]
# ts1z <- zoo(ts1, datetime_cust)
# ts1z_aggr <- aggregate(ts1z, time(ts1z) - as.numeric(time(ts1z)) %% time_interval, mean)
# ts2 <- err_df_raw[,j]
# ts2z <- zoo(ts2, datetime_cust)
# ts2z_aggr <- aggregate(ts2z, time(ts2z) - as.numeric(time(ts2z)) %% time_interval, mean)
# plot(ts1z_aggr, xlab="Date / time", ylab="Average number of daily calls")
# text(min(time(ts1z))+1e6, max(ts1z_aggr), err_name[i], font=2)
# plot(ts2z_aggr, xlab="Date / time", ylab="Average number of daily calls")
# text(min(time(ts2z))+1e6, max(ts2z_aggr), err_name[j], font=2)
# dev.off()




# list of resolutions for DGXN:
res <- sort(unique(cust_data_match$Resolution[cust_data_match$Call_Subject == assay]))
res

pdf(paste0(assay_dir, assay, "_dist_resolutions.pdf"), width=24, height=16)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(4,4))
for(i in 1:length(res)){
	plot(datetime_cust, cumsum(cust_data_match$Resolution == res[i]), main=res[i], xlab="Date/time", ylab="Cummulative sum", pch=20)
	abline(lm(cumsum(cust_data_match$Resolution == res[i]) ~ datetime_cust), col="red")
}
dev.off()

# error codes vs. resolutions
pdf(paste0(assay_dir, assay, "_err_vs_res1.pdf"), width=12, height=18)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,1))
# err codes
distErr <- table(sort(cust_data_match$Call_Area[cust_data_match$Call_Subject == assay]))[table(sort(cust_data_match$Call_Area[cust_data_match$Call_Subject == assay])) > 0]
plot(distErr, main=paste0("Distribution of error codes for ", assay," (aplhabetical order)"), xlab="", ylab="Counts", xaxt = "n", type="h", lwd=5)
axis(1, at=1:length(distErr), labels=names(distErr), las=3, cex=.6)
# resolution codes
distRes <- table(sort(cust_data_match$Resolution[cust_data_match$Call_Subject == assay]))[table(sort(cust_data_match$Resolution[cust_data_match$Call_Subject == assay])) > 0]
plot(distRes, main=paste0("Distribution of resolution codes for", assay," (aplhabetical order)"), xlab="", ylab="Counts", xaxt = "n", type="h", lwd=5)
axis(1, at=1:length(distRes), labels=names(distRes), las=3, cex=.6)
dev.off()

# pdf(paste0(assay_dir, assay, "_err_vs_res2.pdf"), width=12, height=12)
# par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
# # heatmap
# err2 <- cust_data_match$Call_Area[cust_data_match$Call_Subject == assay][!is.na(cust_data_match$Call_Area[cust_data_match$Call_Subject == assay])]
# res <- cust_data_match$Resolution[cust_data_match$Call_Subject == assay][!is.na(cust_data_match$Call_Area[cust_data_match$Call_Subject == assay])]
# DF <- data.frame(err, res)
# mat <- xtabs(~ err2 + res)
# clean_mat <- mat[rowSums(mat != 0) != 0, ]
# heatmap(log10(1/(clean_mat+.9)))
# dev.off()




pdf(paste0(assay_dir, assay, "_err_vs_res3.pdf"), width=24, height=16)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(4,4))
for(i in 1:length(err)){
	plot(table(cust_data_match$Resolution[cust_data_match$Call_Area[cust_data_match$Call_Subject == assay] == err[i] ]), ylab="Counts", las=2, cex=.7)
	text(10, (max(table(cust_data_match$Resolution[cust_data_match$Call_Area[cust_data_match$Call_Subject == assay] == err[i] ]), na.rm=T)*.75), as.character(err[i]), font=2)
}
dev.off()




												########################
												# Sample data matching #
												########################

uniq_SampName <- unique(econn_dat$Sample.Name)
pdf(paste0(assay_dir, assay, "_sample_dist.pdf"), width=(length(uniq_SampName)/100), height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
plot(table(econn_dat$Sample.Name))
dev.off()

pdf(paste0(assay_dir, assay, "_sample_dist2.pdf"), width=6, height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
hist(table(econn_dat$Sample.Name), 50, xlab="Number of complaints", main="")
dev.off()





# find "normal" testing concentrations
means <- NULL
for(i in 1:length(uniq_SampName)){
	means <- c(means, mean(econn_dat$F.Concentration[econn_dat$Sample.Name == uniq_SampName[i]], na.rm=T))
}
ddmeans <- density(means)
peak_locations <- which(diff(sign(diff(ddmeans$y))) == -2) + 1
pdf(paste0(assay_dir, assay, "_modes.pdf"), width=10, height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
plot(ddmeans, type="h", xlab="Concentration means by sample", main="")
empirical_QC <- ddmeans$x[peak_locations]
abline(v= empirical_QC, col="red", lty=2)
title(paste0("Modes at: ", paste(format(empirical_QC, digits=4), collapse=", ")), cex=.5  )
dev.off()


# pdf(paste0(assay_dir, assay, "_all_samples.pdf"), width=6, height=24)
# par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(6,1))
# #for(i in 1:20){
# for(i in 1:length(uniq_SampName)){
# 	# JNs should be shared; call subject is DGXN; shared Lot numbers
# 	cust_list <- unique(cust_data_match$Customer_Number[(unique(econn_dat$J.Number[econn_dat$Sample.Name == uniq_SampName[i]]) %in% cust_data_match$J_Number) & (cust_data_match$Call_Subject == "DGXN") & (sub("-XXXX-", "", as.character(cust_data_match$Lot_Serial_Number)) %in% as.character(econn_dat$Reagent.Lot.Number))])
# 	plot(econn_dat$F.Concentration[econn_dat$Sample.Name == uniq_SampName[i]] ~ datetime[econn_dat$Sample.Name == uniq_SampName[i]], pch=20, ylim=c(0,4), main=paste0("Sample: ", uniq_SampName[i]), xlab="Date / time", ylab="[DGXN]")
# 	abline(h= empirical_QC, col="blue")
# 	abline(h= empirical_QC + 1.96*QC_SDwithinday, col="blue", lty=2)
# 	abline(h= empirical_QC - 1.96*QC_SDwithinday, col="blue", lty=2)
# 	abline(h= empirical_QC + 1.96*QC_SDwithinlab, col="blue", lty=3)
# 	abline(h= empirical_QC - 1.96*QC_SDwithinlab, col="blue", lty=3)
# 	abline(v= datetime_cust[(cust_data_match$Customer_Number %in% cust_list) & (cust_data_match$Call_Subject == assay)], col="orange" )
# 	#text(datetime_cust[(cust_data_match$Customer_Number %in% cust_list[j]) & (cust_data_match$Call_Subject == "DGXN")], 4, paste0(cust_data_match$Call_Subject[(cust_data_match$Customer_Number %in% cust_list[j]) & (cust_data_match$Call_Subject == "DGXN")], "_",  cust_data_match$Call_Area[(cust_data_match$Customer_Number %in% cust_list[j]) & (cust_data_match$Call_Subject == "DGXN")]), cex=.75, adj = c(1,0), srt = 90, col="orange")
# }
# dev.off()






system(paste0("rm -fr ", assay_dir, "PotentialMatchings"))
system(paste0("mkdir ", assay_dir, "PotentialMatchings"))
matchings <- data.frame(
					customerID=character(),
					JN=character(),
					sampleID=character(),
					complained=character(),
					ave_testing_freq=numeric(),
					stringsAsFactors=F
				)
for(i in 1:length(uniq_SampName)){
#for(i in 1:10){
	JN_test_all <- JN_test <- customer <- testing_freq <- freq_units <- NULL
	JN_test_all <- unique(econn_dat$J.Number[econn_dat$Sample.Name == uniq_SampName[i]])
	for(JN_test in JN_test_all){
		customer <- unique(cust_data_match$Customer_Number[cust_data_match$J_Number == JN_test])
		if(length(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)]) > 1){
			testing_freq <- mean(diff(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)], 1))
			freq_units <- units(testing_freq)
			if(freq_units == "secs"){ # all time results in hours
				testing_freq <- as.numeric(testing_freq) / 3600
			}else{
				if(freq_units == "mins"){
					testing_freq <- as.numeric(testing_freq) / 60
				}else{
					if(freq_units == "hours"){
						testing_freq <- as.numeric(testing_freq)
					}else{
						if(freq_units == "days"){
							testing_freq <- as.numeric(testing_freq) * 24.
						}else{
							print("TIME!!!")
							break
						}
					}
				}
			}
		}else{
			testing_freq <- NA
		}
		if(length(customer) > 0){
			matchings[i,] <- list(as.character(customer), as.character(JN_test), as.character(uniq_SampName[i]), testing_freq)
			# positive samples: customer complained
			if(sum((cust_data_match$Customer_Number  == customer) & (cust_data_match$Call_Subject == assay)) > 0){
				matchings[i,] <- list(as.character(customer), as.character(JN_test), as.character(uniq_SampName[i]), "yes", testing_freq)
				pdf(paste0(assay_dir, "PotentialMatchings/pos_",assay, "_all_samples_PotentialMatchings",i, "_", JN_test,".pdf"), width=6, height=4)
				par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
				min_date <- min(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)], datetime_cust[(cust_data_match$Customer_Number  == customer) & (cust_data_match$Call_Subject == assay)])
				max_date <- max(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)], datetime_cust[(cust_data_match$Customer_Number  == customer) & (cust_data_match$Call_Subject == assay)])
				plot(econn_dat$F.Concentration[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)] ~ datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)], pch=20, xlim=c(min_date, max_date), ylim=c(0,max(means)), main=paste0("Sample: ", uniq_SampName[i], " for cust ", customer, " JN:", JN_test), xlab="Date / time", ylab="[DGXN]")
				abline(h= empirical_QC, col=c("blue"))
	#			abline(h= QC + 1.96*QC_SDwithinday, col=c("blue", "gray", "red"), lty=2)
	#			abline(h= QC - 1.96*QC_SDwithinday, col=c("blue", "gray", "red"), lty=2)
	#			abline(h= QC + 1.96*QC_SDwithinlab, col=c("blue", "gray", "red"), lty=3)
	#			abline(h= QC - 1.96*QC_SDwithinlab, col=c("blue", "gray", "red"), lty=3)
				abline(v= datetime_cust[(cust_data_match$Customer_Number == customer) & (cust_data_match$Call_Subject == assay)], col="orange" )
				text(datetime_cust[(cust_data_match$Customer_Number  == customer) & (cust_data_match$Call_Subject == assay)], max(means), paste0(cust_data_match$Call_Subject[(cust_data_match$Customer_Number  == customer) & (cust_data_match$Call_Subject == assay)], "_",cust_data_match$Complaint_Nbr___CH[(cust_data_match$Customer_Number  == customer) & (cust_data_match$Call_Subject == assay)], "_",  cust_data_match$Call_Area[(cust_data_match$Customer_Number  == customer) & (cust_data_match$Call_Subject == assay)]), cex=.75, adj = c(1,0), srt = 90, col="orange")
				# ID tining of changes of S.Gen
				running <- timing <- NULL
				running <- econn_dat$S.Gen[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				abline(v = datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], col="gray")
				text(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], max(means), paste0("S.Gen ",  econn_dat$S.Gen[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing]), cex=.75, adj = c(1,0), srt = 90, col="gray")
				# ID tining of changes of S.Lot
				running <- timing <- NULL
				running <- econn_dat$S.Lot[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				abline(v = datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], col="gray")
				text(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], max(means), paste0("S.Lot ",  econn_dat$S.Lot[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing]), cex=.75, adj = c(1,0), srt = 90, col="gray")
				# ID tining of changes of ERF.Lot
				running <- timing <- NULL
				running <- econn_dat$ERF.Lot[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				abline(v = datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], col="gray")
				text(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], max(means), paste0("ERF.Lot ",  econn_dat$ERF.Lot[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing]), cex=.75, adj = c(1,0), srt = 90, col="gray")
				# ID tining of changes of IWF.Lot
				running <- timing <- NULL
				running <- econn_dat$IWF.Lot[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				abline(v = datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], col="gray")
				text(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], max(means), paste0("IWF.Lot ",  econn_dat$IWF.Lot[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing]), cex=.75, adj = c(1,0), srt = 90, col="gray")
				# ID tining of changes of Control.Lot.Number
				running <- timing <- NULL
				running <- econn_dat$Control.Lot.Number[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				abline(v = datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], col="gray")
				text(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], max(means), paste0("Control.Lot.Number ",  econn_dat$Control.Lot.Number[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing]), cex=.75, adj = c(1,0), srt = 90, col="gray")
				# ID tining of changes of Cal.Curve.ID
				running <- timing <- NULL
				running <- econn_dat$Cal.Curve.ID[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				abline(v = datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], col="gray")
				text(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], max(means), paste0("Cal.Curve.ID ",  econn_dat$Cal.Curve.ID[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing]), cex=.75, adj = c(1,0), srt = 90, col="gray")
				dev.off()
			# negative samples: customer didn't complain
			}else{
				matchings[i,] <- list(as.character(customer), as.character(JN_test), as.character(uniq_SampName[i]), "no", testing_freq)
				pdf(paste0(assay_dir, "PotentialMatchings/neg_",assay, "_all_samples_PotentialMatchings",i, "_", JN_test,".pdf"), width=6, height=4)
				par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
				min_date <- min(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)], datetime_cust[(cust_data_match$Customer_Number  == customer) & (cust_data_match$Call_Subject == assay)])
				max_date <- max(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)], datetime_cust[(cust_data_match$Customer_Number  == customer) & (cust_data_match$Call_Subject == assay)])
				plot(econn_dat$F.Concentration[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)] ~ datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)], pch=20, xlim=c(min_date, max_date), ylim=c(0,max(means)), main=paste0("Sample: ", uniq_SampName[i], " for cust ", customer, " JN:", JN_test), xlab="Date / time", ylab="[DGXN]")
				abline(h= empirical_QC, col=c("blue"))
	#			abline(h= QC + 1.96*QC_SDwithinday, col=c("blue", "gray", "red"), lty=2)
	#			abline(h= QC - 1.96*QC_SDwithinday, col=c("blue", "gray", "red"), lty=2)
	#			abline(h= QC + 1.96*QC_SDwithinlab, col=c("blue", "gray", "red"), lty=3)
	#			abline(h= QC - 1.96*QC_SDwithinlab, col=c("blue", "gray", "red"), lty=3)
				abline(v= datetime_cust[(cust_data_match$Customer_Number == customer) & (cust_data_match$Call_Subject == assay)], col="orange" )
				text(datetime_cust[(cust_data_match$Customer_Number  == customer) & (cust_data_match$Call_Subject == assay)], 4, paste0(cust_data_match$Call_Subject[(cust_data_match$Customer_Number  == customer) & (cust_data_match$Call_Subject == assay)], "_",cust_data_match$Complaint_Nbr___CH[(cust_data_match$Customer_Number  == customer) & (cust_data_match$Call_Subject == assay)], "_",  cust_data_match$Call_Area[(cust_data_match$Customer_Number  == customer) & (cust_data_match$Call_Subject == assay)]), cex=.75, adj = c(1,0), srt = 90, col="orange")
				# ID tining of changes of S.Gen
				running <- timing <- NULL
				running <- econn_dat$S.Gen[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				abline(v = datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], col="gray")
				text(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], max(means), paste0("S.Gen ",  econn_dat$S.Gen[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing]), cex=.75, adj = c(1,0), srt = 90, col="gray")
				# ID tining of changes of S.Lot
				running <- timing <- NULL
				running <- econn_dat$S.Lot[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				abline(v = datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], col="gray")
				text(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], max(means), paste0("S.Lot ",  econn_dat$S.Lot[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing]), cex=.75, adj = c(1,0), srt = 90, col="gray")
				# ID tining of changes of ERF.Lot
				running <- timing <- NULL
				running <- econn_dat$ERF.Lot[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				abline(v = datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], col="gray")
				text(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], max(means), paste0("ERF.Lot ",  econn_dat$ERF.Lot[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing]), cex=.75, adj = c(1,0), srt = 90, col="gray")
				# ID tining of changes of IWF.Lot
				running <- timing <- NULL
				running <- econn_dat$IWF.Lot[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				abline(v = datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], col="gray")
				text(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], max(means), paste0("IWF.Lot ",  econn_dat$IWF.Lot[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing]), cex=.75, adj = c(1,0), srt = 90, col="gray")
				# ID tining of changes of Control.Lot.Number
				running <- timing <- NULL
				running <- econn_dat$Control.Lot.Number[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				abline(v = datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], col="gray")
				text(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], max(means), paste0("Control.Lot.Number ",  econn_dat$Control.Lot.Number[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing]), cex=.75, adj = c(1,0), srt = 90, col="gray")
				# ID tining of changes of Cal.Curve.ID
				running <- timing <- NULL
				running <- econn_dat$Cal.Curve.ID[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				abline(v = datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], col="gray")
				text(datetime[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing], max(means), paste0("Cal.Curve.ID ",  econn_dat$Cal.Curve.ID[(econn_dat$Sample.Name == uniq_SampName[i]) & (econn_dat$J.Number == JN_test)][timing]), cex=.75, adj = c(1,0), srt = 90, col="gray")
				dev.off()
			}
	
		}else{
			matchings[i,] <- list(NA, as.character(JN_test), as.character(uniq_SampName[i]), NA, testing_freq)
		}
	}
}
matchings


# no significant difference wrt log frequency between those who complained and who did not
pdf(paste0(assay_dir, assay, "_ave_testing_freq.pdf"), width=10, height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
hist(as.numeric(matchings$ave_testing_freq),100, xlab="E-conn log frequency (hours)", main="")
hist(as.numeric(matchings$ave_testing_freq)/24,100, xlab="E-conn log frequency (days)", main="")
logFreq_all <- locfit( ~ lp(as.numeric(matchings$ave_testing_freq)), maxk=500)
logFreq_nocall <- locfit( ~ lp(as.numeric(matchings$ave_testing_freq[matchings$complained == "no"])), maxk=500)
logFreq_call <- locfit( ~ lp(as.numeric(matchings$ave_testing_freq[matchings$complained == "yes"])), maxk=500)
plot(logFreq_all, xlab="E-conn log frequency (hours)", ylab="Density")
lines(logFreq_nocall, col="blue")
lines(logFreq_call, col="red")
t.test(as.numeric(matchings$ave_testing_freq[matchings$complained == "no"]), as.numeric(matchings$ave_testing_freq[matchings$complained == "yes"]))
dev.off()
#	Welch Two Sample t-test
#
#data:  as.numeric(matchings$ave_testing_freq[matchings$complained ==  and as.numeric(matchings$ave_testing_freq[matchings$complained ==     "no"]) and     "yes"])
#t = 2.5892, df = 1704.8, p-value = 0.009701
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  4.872734 35.315219
#sample estimates:
#mean of x mean of y 
#111.80394  91.70996 


################################################################################
################################################################################
################################## MODEL 1 #####################################
################################################################################
################################################################################
################################################################################

convert2hours <- function(mytime){
	freq_units <- units(mytime)
	if(freq_units == "secs"){ # all time results in hours
		mytime <- as.numeric(mytime) / 3600
	}else{
		if(freq_units == "mins"){
			mytime <- as.numeric(mytime) / 60
		}else{
			if(freq_units == "hours"){
				mytime <- as.numeric(mytime)
			}else{
				if(freq_units == "days"){
					mytime <- as.numeric(mytime) * 24.
				}else{
					print("TIME!!!")
					break
				}
			}
		}
	}
	return(mytime)
}


# create data for positive samples (those that generated a customer call)
# ecust_dat_neg is generated below
AssayEntries <- which(cust_data_match$Call_Subject == assay)
ecust_dat <- ecust_dat_NA <- ecust_dat_neg <- data.frame(
					i=numeric(), 
					j=numeric(), 
					# customer data
					CallArea =character(),
					CustomerNumber =character(),
					JN=character(),
					ComplaintID=character(),
					AuditDate =numeric(),
					# e-conn data
					SampleName =character(),
					DateSinceSampleStarted =numeric(),
					NbPriorSGenChange =numeric(),
					NbPriorSLotChange =numeric(),
					NbPriorERFLotChange =numeric(),
					NbPriorIWFLotChange =numeric(),
					NbPriorContLotNumChange =numeric(),
					NbPriorCalCurveChange =numeric(),
					TimeSinceLastSGenChange =numeric(),
					TimeSinceLastSLotChange =numeric(),
					TimeSinceLastERFLotChange =numeric(),
					TimeSinceLastIWFLotChange =numeric(),
					TimeSinceLastContLotNumChange =numeric(),
					TimeSinceLastCalCurveChange =numeric(),
					# e-conn concentrations
					TimeToComplain =numeric(),
					MostRecentConcentration =numeric(),
					FiveMostRecentConcentrationMean =numeric(),
					FiveMostRecentConcentrationSD =numeric(),
					TenMostRecentConcentrationMean =numeric(),
					TenMostRecentConcentrationSD =numeric(),
					TwoMostRecentConcentrationMean =numeric(),
					TwoMostRecentConcentrationSD =numeric(),
					# miscell
					stringsAsFactors=F
				)

for(i in 1:length(AssayEntries)){
	# customer data
	cur_i <- AssayEntries[i]
	CallArea <- CustomerNumber <- JN <- SampleNames <- NULL
	CallArea <- cust_data_match$Call_Area[cur_i]
	CustomerNumber <- cust_data_match$Customer_Number[cur_i]
	JN <- cust_data_match$J_Number[cur_i]
	ComplaintID <- cust_data_match$Complaint_Nbr___CH[cur_i]
	AuditDate <- datetime_cust[cur_i]
	# which(matchings$JN == JN) # to check w/ figures in PotentialMatchings
	# e-conn data
	#SampleNames <- matchings$sampleID[which(matchings$customerID == CustomerNumber)]
	SampleNames <- unique(econn_dat$Sample.Name[econn_dat$J.Number == JN])
	if(length(SampleNames) == 0){
		# not sure why this JN is not in econn...
				ecust_dat_NA[(length(ecust_dat_NA[,1])+1),] <- list(
						i,"JN_not_in_econn",
						# customer data
						as.character(CallArea),as.character(CustomerNumber), as.character(JN), as.character(ComplaintID), AuditDate,
						# e-conn data
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA,
						NA
					)
	}else{
		# this JN is in econn...
		for(j in 1:length(SampleNames)){
			cur_samp_name <- SampleNames[j]
			# we want the call to be placed after any e-conn logs
			time_constraint <- datetime_cust[cur_i] >= datetime[econn_dat$Sample.Name == cur_samp_name]
			DateSinceSampleStarted <- NULL
			if((time_constraint[1]) & (length(time_constraint) > 1)){
				####################################################
				# all samples with only one e-conn log are dropped #
				####################################################
				print(paste0(cur_i , " CurSamp:",cur_samp_name, ": ", sum(time_constraint), " logs"))
				DateSinceSampleStarted <- datetime[econn_dat$Sample.Name == cur_samp_name][1]
				# S.Gen
				running <- timing <- dateofchange <- NULL
				running <- econn_dat$S.Gen[econn_dat$Sample.Name == cur_samp_name]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][timing]
				NbPriorSGenChange <- sum(dateofchange <= datetime_cust[cur_i])
				timediff <- NULL
				if(!length(dateofchange)){
					dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][1]
				}
				timediff <- convert2hours(AuditDate - dateofchange)
				if(sum(timediff) > 0){
					# there was at least one change before the call
					TimeSinceLastSGenChange <- min(timediff[timediff > 0])
				}else{
					# all changes were after the call
					TimeSinceLastSGenChange <- timediff[which.min(abs(timediff))]
				}
				# S.Lot
				running <- timing <- dateofchange <- NULL
				running <- econn_dat$S.Lot[econn_dat$Sample.Name == cur_samp_name]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][timing]
				NbPriorSLotChange <- sum(dateofchange <= datetime_cust[cur_i])
				timediff <- NULL
				if(!length(dateofchange)){
					dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][1]
				}
				timediff <- convert2hours(AuditDate - dateofchange)
				if(sum(timediff) > 0){
					# there was at least one change before the call
					TimeSinceLastSLotChange <- min(timediff[timediff > 0])
				}else{
					# all changes were after the call
					TimeSinceLastSLotChange <- timediff[which.min(abs(timediff))]
				}
				# ERF.Lot
				running <- timing <- dateofchange <- NULL
				running <- econn_dat$ERF.Lot[econn_dat$Sample.Name == cur_samp_name]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][timing]
				NbPriorERFLotChange <- sum(dateofchange <= datetime_cust[cur_i])
				timediff <- NULL
				if(!length(dateofchange)){
					dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][1]
				}
				timediff <- convert2hours(AuditDate - dateofchange)
				if(sum(timediff) > 0){
					# there was at least one change before the call
					TimeSinceLastERFLotChange <- min(timediff[timediff > 0])
				}else{
					# all changes were after the call
					TimeSinceLastERFLotChange <- timediff[which.min(abs(timediff))]
				}
				# IWF.Lot
				running <- timing <- dateofchange <- NULL
				running <- econn_dat$IWF.Lot[econn_dat$Sample.Name == cur_samp_name]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][timing]
				NbPriorIWFLotChange <- sum(dateofchange <= datetime_cust[cur_i])
				timediff <- NULL
				if(!length(dateofchange)){
					dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][1]
				}
				timediff <- convert2hours(AuditDate - dateofchange)
				if(sum(timediff) > 0){
					# there was at least one change before the call
					TimeSinceLastIWFLotChange <- min(timediff[timediff > 0])
				}else{
					# all changes were after the call
					TimeSinceLastIWFLotChange <- timediff[which.min(abs(timediff))]
				}
				# Control.Lot.Number
				running <- timing <- dateofchange <- NULL
				running <- econn_dat$Control.Lot.Number[econn_dat$Sample.Name == cur_samp_name]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][timing]
				NbPriorContLotNumChange <- sum(dateofchange <= datetime_cust[cur_i])
				timediff <- NULL
				if(!length(dateofchange)){
					dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][1]
				}
				timediff <- convert2hours(AuditDate - dateofchange)
				if(sum(timediff) > 0){
					# there was at least one change before the call
					TimeSinceLastContLotNumChange <- min(timediff[timediff > 0])
				}else{
					# all changes were after the call
					TimeSinceLastContLotNumChange <- timediff[which.min(abs(timediff))]
				}
				# Cal.Curve.ID
				running <- timing <- dateofchange <- NULL
				running <- econn_dat$Cal.Curve.ID[econn_dat$Sample.Name == cur_samp_name]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][timing]
				NbPriorCalCurveChange <- sum(dateofchange <= datetime_cust[cur_i])
				timediff <- NULL
				if(!length(dateofchange)){
					dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][1]
				}
				timediff <- convert2hours(AuditDate - dateofchange)
				if(sum(timediff) > 0){
					# there was at least one change before the call
					TimeSinceLastCalCurveChange <- min(timediff[timediff > 0])
				}else{
					# all changes were after the call
					TimeSinceLastCalCurveChange <- timediff[which.min(abs(timediff))]
				}
				
				# concentrations
				MostRecentConcentration <- FiveMostRecentConcentrationMean <- FiveMostRecentConcentrationSD <- TenMostRecentConcentrationMean <- TenMostRecentConcentrationSD <- TwoMostRecentConcentrationMean <- TwoMostRecentConcentrationSD <- NULL
				running <- timing <- dateofchange <- NULL
				running <- time_constraint
				if((length(running) > 1) & (!all(is.na(running)))){
					if(length(time_constraint) == sum(time_constraint)){
						# no more e-conn logs after this call
						timing <- sum(time_constraint) + 1
					}else{
						for(k in 1:(length(running)-1)){
							if(!(running[k+1] == running[k])){
								timing <- c(timing, k+1)
							}
						}
					}
					MostRecentConcentration <- econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][timing - 1]
					FiveMostRecentConcentrationMean <- mean(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][max(1,(timing - 6)):(timing - 1)])
					FiveMostRecentConcentrationSD <- sd(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][max(1,(timing - 6)):(timing - 1)])
					if(is.na(FiveMostRecentConcentrationSD)){
						FiveMostRecentConcentrationSD <- 0
					}
					TenMostRecentConcentrationMean <- mean(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][max(1,(timing - 11)):(timing - 1)])
					TenMostRecentConcentrationSD <- sd(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][max(1,(timing - 11)):(timing - 1)])
					if(is.na(TenMostRecentConcentrationSD)){
						TenMostRecentConcentrationSD <- 0
					}
					TwoMostRecentConcentrationMean <- mean(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][max(1,(timing - 3)):(timing - 1)])
					TwoMostRecentConcentrationSD <- sd(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][max(1,(timing - 3)):(timing - 1)])
					if(is.na(TwoMostRecentConcentrationSD)){
						TwoMostRecentConcentrationSD <- 0
					}
				}else{
					timing <- 1
					MostRecentConcentration <- econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][timing]
					FiveMostRecentConcentrationMean <- NA
					FiveMostRecentConcentrationSD <- NA
					TenMostRecentConcentrationMean <- NA
					TenMostRecentConcentrationSD <- NA
					TwoMostRecentConcentrationMean <- NA
					TwoMostRecentConcentrationSD <- NA
				}
				TimeToComplain <- savelog <- NULL
				TimeToComplain <- datetime_cust[cur_i] - datetime[econn_dat$Sample.Name == cur_samp_name][timing - 1]
				TimeToComplain <- convert2hours(TimeToComplain)
				
				
				# fill in data frame for modeling
				#if(savelog){
				if(1){
					ConcentrationDist <- NULL
					ConcentrationDist <- hist(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name], plot=F)
					MostProbableQC <- ConcentrationDist$breaks[which.max(ConcentrationDist$density)]
					MostProbableQC <- empirical_QC[which.min(abs(MostProbableQC - empirical_QC))]
					ecust_dat[(length(ecust_dat[,1])+1),] <- list(
							i,j,
							# customer data
							as.character(CallArea),as.character(CustomerNumber), as.character(JN), as.character(ComplaintID), AuditDate,
							# e-conn data
							as.character(cur_samp_name),
							DateSinceSampleStarted,
							NbPriorSGenChange,
							NbPriorSLotChange,
							NbPriorERFLotChange,
							NbPriorIWFLotChange,
							NbPriorContLotNumChange,
							NbPriorCalCurveChange,
							TimeSinceLastSGenChange,
							TimeSinceLastSLotChange,
							TimeSinceLastERFLotChange,
							TimeSinceLastIWFLotChange,
							TimeSinceLastContLotNumChange,
							TimeSinceLastCalCurveChange,
							as.numeric(TimeToComplain),
							(MostRecentConcentration - MostProbableQC)/MostProbableQC, # smallest relative distance to assay specific MostProbableQC
							(FiveMostRecentConcentrationMean - MostProbableQC)/MostProbableQC,
							FiveMostRecentConcentrationSD,
							(TenMostRecentConcentrationMean - MostProbableQC)/MostProbableQC,
							TenMostRecentConcentrationSD,
							(TwoMostRecentConcentrationMean - MostProbableQC)/MostProbableQC,
							TwoMostRecentConcentrationSD
						)
				}
	
			}
		}
	}
}

head(ecust_dat_NA)
write.csv(ecust_dat_NA, paste0(assay_dir, "ecust_dat_NA_", assay, ".csv"), quote=F, row.names=F)
head(ecust_dat)
write.csv(ecust_dat, paste0(assay_dir, "ecust_dat_", assay, ".csv"), quote=F, row.names=F)


pdf(paste0(assay_dir, assay, "_RecentConcentrationMean.pdf"), width=10, height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
plot(ecust_dat$TwoMostRecentConcentrationMean, ecust_dat$FiveMostRecentConcentrationMean, pch=20)
plot(ecust_dat$TwoMostRecentConcentrationMean, ecust_dat$TenMostRecentConcentrationMean, pch=20)
plot(ecust_dat$FiveMostRecentConcentrationMean, ecust_dat$TenMostRecentConcentrationMean, pch=20)
dev.off()


# what is the most common timelag between issue and call?
pdf(paste0(assay_dir, assay, "_TTC_mode.pdf"), width=6, height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
TTC <- locfit( ~ lp(ecust_dat$TimeToComplain), maxk=500)
mode_TTC <- ecust_dat$TimeToComplain[which.max(predict(TTC, newdata=ecust_dat$TimeToComplain))]
plot(TTC, xlab="Time to complain", ylab="Density", main=paste0("Mode at ", mode_TTC, " hours"))
abline(v= mode_TTC, lty=2, col="red")
dev.off()


# create ecust_dat_neg, containing the negative data (no calls placed)
n_pos <- length(AssayEntries)
#econn_dat$Assay[econn_dat$Assay == "Na+"] <- assay
all_JNs <- unique(econn_dat$J.Number[econn_dat$Assay == assay])
pos_JNs <- unique(cust_data$J_Number[cust_data$Call_Subject == assay])
neg_JNs <- setdiff(all_JNs, pos_JNs)
neg_JNs_samp <- sample(neg_JNs, 1*n_pos, replace=F) # increasing the # of negative samples increases classification accuracy, but just b/c more OKs are well classified...
for(i in 1:length(neg_JNs_samp)){
	CallArea <- CustomerNumber <- JN <- SampleNames <- NULL
	cur_samp_name <- AuditDate <- time_constraint <- NULL
	# ID a sample
	counter <- 0
	repeat{
		cur_samp_name <- sample(unique(econn_dat$Sample.Name[econn_dat$J.Number == neg_JNs_samp[i]]), 1)
		if(length(datetime[econn_dat$Sample.Name == cur_samp_name]) > 10){
			break
		}
		counter <- counter + 1
		if(counter >= length(unique(econn_dat$Sample.Name[econn_dat$J.Number == neg_JNs_samp[i]]))){
			neg_JNs_samp[i] <- sample(neg_JNs, 1, replace=F)
		}
	}
	# negative sample	
	CallArea <- "OK"
	CustomerNumber <- "none"
	JN <- neg_JNs_samp[i]
	ComplaintID <- "none"
	# draw a random timepoint
	startDate <- endDate <- NULL
	startDate <- min(datetime[econn_dat$Sample.Name == cur_samp_name])
	endDate <- max(datetime[econn_dat$Sample.Name == cur_samp_name])
	AuditDate <- sample(seq(startDate, endDate, 3600/4), 1) # every 15 min, or 3600/4 sec
	time_constraint <- AuditDate >= datetime[econn_dat$Sample.Name == cur_samp_name]
	# get all details
				# S.Gen
				running <- timing <- dateofchange <- NULL
				running <- econn_dat$S.Gen[econn_dat$Sample.Name == cur_samp_name]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][timing]
				NbPriorSGenChange <- sum(dateofchange <= datetime_cust[cur_i])
				timediff <- NULL
				if(!length(dateofchange)){
					dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][1]
				}
				timediff <- convert2hours(AuditDate - dateofchange)
				if(sum(timediff) > 0){
					# there was at least one change before the call
					TimeSinceLastSGenChange <- min(timediff[timediff > 0])
				}else{
					# all changes were after the call
					TimeSinceLastSGenChange <- timediff[which.min(abs(timediff))]
				}
				# S.Lot
				running <- timing <- dateofchange <- NULL
				running <- econn_dat$S.Lot[econn_dat$Sample.Name == cur_samp_name]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][timing]
				NbPriorSLotChange <- sum(dateofchange <= datetime_cust[cur_i])
				timediff <- NULL
				if(!length(dateofchange)){
					dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][1]
				}
				timediff <- convert2hours(AuditDate - dateofchange)
				if(sum(timediff) > 0){
					# there was at least one change before the call
					TimeSinceLastSLotChange <- min(timediff[timediff > 0])
				}else{
					# all changes were after the call
					TimeSinceLastSLotChange <- timediff[which.min(abs(timediff))]
				}
				# ERF.Lot
				running <- timing <- dateofchange <- NULL
				running <- econn_dat$ERF.Lot[econn_dat$Sample.Name == cur_samp_name]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][timing]
				NbPriorERFLotChange <- sum(dateofchange <= datetime_cust[cur_i])
				timediff <- NULL
				if(!length(dateofchange)){
					dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][1]
				}
				timediff <- convert2hours(AuditDate - dateofchange)
				if(sum(timediff) > 0){
					# there was at least one change before the call
					TimeSinceLastERFLotChange <- min(timediff[timediff > 0])
				}else{
					# all changes were after the call
					TimeSinceLastERFLotChange <- timediff[which.min(abs(timediff))]
				}
				# IWF.Lot
				running <- timing <- dateofchange <- NULL
				running <- econn_dat$IWF.Lot[econn_dat$Sample.Name == cur_samp_name]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][timing]
				NbPriorIWFLotChange <- sum(dateofchange <= datetime_cust[cur_i])
				timediff <- NULL
				if(!length(dateofchange)){
					dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][1]
				}
				timediff <- convert2hours(AuditDate - dateofchange)
				if(sum(timediff) > 0){
					# there was at least one change before the call
					TimeSinceLastIWFLotChange <- min(timediff[timediff > 0])
				}else{
					# all changes were after the call
					TimeSinceLastIWFLotChange <- timediff[which.min(abs(timediff))]
				}
				# Control.Lot.Number
				running <- timing <- dateofchange <- NULL
				running <- econn_dat$Control.Lot.Number[econn_dat$Sample.Name == cur_samp_name]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][timing]
				NbPriorContLotNumChange <- sum(dateofchange <= datetime_cust[cur_i])
				timediff <- NULL
				if(!length(dateofchange)){
					dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][1]
				}
				timediff <- convert2hours(AuditDate - dateofchange)
				if(sum(timediff) > 0){
					# there was at least one change before the call
					TimeSinceLastContLotNumChange <- min(timediff[timediff > 0])
				}else{
					# all changes were after the call
					TimeSinceLastContLotNumChange <- timediff[which.min(abs(timediff))]
				}
				# Cal.Curve.ID
				running <- timing <- dateofchange <- NULL
				running <- econn_dat$Cal.Curve.ID[econn_dat$Sample.Name == cur_samp_name]
				# impute missing values
				if(sum(is.na(running)) > 0){
					na_pos <- which(is.na(running))
					if(length(na_pos) == length(running)){
						running <- rep(econn_dat$ERF.Lot[1], length(running))
					}else{
						for(mis in na_pos){
							# impute 1st position
							if(mis == 1){
								running[mis] <- running[which(!is.na(running))[1]] # the 1st one that is not NA
							}else{
								# impute last position
								if(mis == length(running)){
									running[mis] <- running[which(!is.na(running))[length(which(!is.na(running)))]] # the last one that was not NA
								}else{
									# impute other positions
									candidates <- which(!is.na(running))
									candidates <- candidates[which(candidates < mis)] # prev positions not NA
									running[mis] <- running[max(candidates)] # copy state at nearest not NA
								}
							}
						}
					}
				}
				if((length(running) > 1) & (!all(is.na(running)))){
					for(k in 1:(length(running)-1)){
						if(!(running[k+1] == running[k])){
							timing <- c(timing, k+1)
						}
					}
				}
				dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][timing]
				NbPriorCalCurveChange <- sum(dateofchange <= datetime_cust[cur_i])
				timediff <- NULL
				if(!length(dateofchange)){
					dateofchange <- datetime[econn_dat$Sample.Name == cur_samp_name][1]
				}
				timediff <- convert2hours(AuditDate - dateofchange)
				if(sum(timediff) > 0){
					# there was at least one change before the call
					TimeSinceLastCalCurveChange <- min(timediff[timediff > 0])
				}else{
					# all changes were after the call
					TimeSinceLastCalCurveChange <- timediff[which.min(abs(timediff))]
				}
				
				# concentrations
				MostRecentConcentration <- FiveMostRecentConcentrationMean <- FiveMostRecentConcentrationSD <- TenMostRecentConcentrationMean <- TenMostRecentConcentrationSD <- TwoMostRecentConcentrationMean <- TwoMostRecentConcentrationSD <- NULL
				running <- timing <- dateofchange <- NULL
				running <- time_constraint
				if((length(running) > 1) & (!all(is.na(running)))){
					if(length(time_constraint) == sum(time_constraint)){
						# no more e-conn logs after this call
						timing <- sum(time_constraint) + 1
					}else{
						for(k in 1:(length(running)-1)){
							if(!(running[k+1] == running[k])){
								timing <- c(timing, k+1)
							}
						}
					}
					MostRecentConcentration <- econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][timing - 1]
					FiveMostRecentConcentrationMean <- mean(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][max(1,(timing - 6)):(timing - 1)])
					FiveMostRecentConcentrationSD <- sd(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][max(1,(timing - 6)):(timing - 1)])
					if(is.na(FiveMostRecentConcentrationSD)){
						FiveMostRecentConcentrationSD <- 0
					}
					TenMostRecentConcentrationMean <- mean(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][max(1,(timing - 11)):(timing - 1)])
					TenMostRecentConcentrationSD <- sd(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][max(1,(timing - 11)):(timing - 1)])
					if(is.na(TenMostRecentConcentrationSD)){
						TenMostRecentConcentrationSD <- 0
					}
					TwoMostRecentConcentrationMean <- mean(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][max(1,(timing - 3)):(timing - 1)])
					TwoMostRecentConcentrationSD <- sd(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][max(1,(timing - 3)):(timing - 1)])
					if(is.na(TwoMostRecentConcentrationSD)){
						TwoMostRecentConcentrationSD <- 0
					}
				}else{
					timing <- 1
					MostRecentConcentration <- econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name][timing]
					FiveMostRecentConcentrationMean <- NA
					FiveMostRecentConcentrationSD <- NA
					TenMostRecentConcentrationMean <- NA
					TenMostRecentConcentrationSD <- NA
					TwoMostRecentConcentrationMean <- NA
					TwoMostRecentConcentrationSD <- NA
				}
				TimeToComplain <- savelog <- NULL
				TimeToComplain <- AuditDate - startDate
				TimeToComplain <- convert2hours(TimeToComplain)

	# fill out data
	ConcentrationDist <- NULL
	ConcentrationDist <- hist(econn_dat$F.Concentration[econn_dat$Sample.Name == cur_samp_name], plot=F)
	MostProbableQC <- ConcentrationDist$breaks[which.max(ConcentrationDist$density)]
	MostProbableQC <- empirical_QC[which.min(abs(MostProbableQC - empirical_QC))]
	ecust_dat_neg[(length(ecust_dat_neg[,1])+1),] <- list(
			i,0,
			# customer data
			as.character(CallArea),as.character(CustomerNumber), as.character(JN), as.character(ComplaintID), AuditDate,
			# e-conn data
			as.character(cur_samp_name),
			startDate,
			NbPriorSGenChange,
			NbPriorSLotChange,
			NbPriorERFLotChange,
			NbPriorIWFLotChange,
			NbPriorContLotNumChange,
			NbPriorCalCurveChange,
			TimeSinceLastSGenChange,
			TimeSinceLastSLotChange,
			TimeSinceLastERFLotChange,
			TimeSinceLastIWFLotChange,
			TimeSinceLastContLotNumChange,
			TimeSinceLastCalCurveChange,
			as.numeric(TimeToComplain),
			(MostRecentConcentration - MostProbableQC)/MostProbableQC, # smallest relative distance to assay specific MostProbableQC
			(FiveMostRecentConcentrationMean - MostProbableQC)/MostProbableQC,
			FiveMostRecentConcentrationSD,
			(TenMostRecentConcentrationMean - MostProbableQC)/MostProbableQC,
			TenMostRecentConcentrationSD,
			(TwoMostRecentConcentrationMean - MostProbableQC)/MostProbableQC,
			TwoMostRecentConcentrationSD
		)
	print(paste0("Done: ", i, " with sample ", as.character(cur_samp_name)))
}

ecust_dat_all <- rbind(ecust_dat, ecust_dat_neg)

					#############################
					# CART: with TimeToComplain #
					#############################


ecust_dat_noNA <- ecust_dat_all[complete.cases(ecust_dat_all), ]
ecust_dat_noNA$CallArea <- as.factor(ecust_dat_noNA$CallArea)

Nrep <- 2500
error.cart1 <- c()
cart_all1 <- tb1 <- list()
for(i in 1:Nrep){
	if(!(i %% 100)){
		print(paste0("Now doing: ", i, "/", Nrep))
	}
	l <- length(ecust_dat_noNA[,1])
	sub <- sample(1:l,2*l/3)
	mfinal <- 10 
	maxdepth <- 5
	
	cart_all1[[i]] <- tree(CallArea ~ NbPriorSGenChange+ NbPriorSLotChange+ NbPriorERFLotChange+ NbPriorIWFLotChange+ NbPriorContLotNumChange+ NbPriorCalCurveChange+ TimeToComplain+ MostRecentConcentration+ FiveMostRecentConcentrationMean+ FiveMostRecentConcentrationSD+ TenMostRecentConcentrationMean+ TenMostRecentConcentrationSD+ TwoMostRecentConcentrationMean+ TwoMostRecentConcentrationSD + TimeSinceLastSGenChange + TimeSinceLastSLotChange + TimeSinceLastERFLotChange + TimeSinceLastIWFLotChange + TimeSinceLastContLotNumChange + TimeSinceLastCalCurveChange, data= ecust_dat_noNA[sub,])
	cart_all1_pred <- predict(cart_all1[[i]], newdata= ecust_dat_noNA[-sub, ],type="class")
	tb1[[i]] <- table(cart_all1_pred, ecust_dat_noNA$CallArea[-sub])
	error.cart1[i] <- 1-(sum(diag(tb1[[i]]))/sum(tb1[[i]]))
}

# best model
best_tb1 <- tb1[[which.min(error.cart1)]]
best_tb1
best_cart_all1 <- cart_all1[[which.min(error.cart1)]]
summary(best_cart_all1)
pdf(paste0(assay_dir, assay, "_ML1_CART.pdf"), width=15, height=18)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,1))
hist(error.cart1, 20, xlab="CART error rate", xlim=c(0,.5), main="")
plot(best_cart_all1)
text(best_cart_all1, pretty=0)
dev.off()



					#################################
					# adaboost: with TimeToComplain #
					#################################

Nrep <- 2500
error.adaboost1 <- c()
adaboost1_all <- adaboost1_all_conf <- adaboost1_all_err <- adaboost1_all_pred <- importance_sorted1 <- list()
for(i in 1:Nrep){
	if(!(i %% 10)){
		print(paste0("Now doing: ", i, "/", Nrep))
	}
	l <- length(ecust_dat_noNA[,1])
	sub <- sample(1:l,2*l/3)
	mfinal <- 10 
	maxdepth <- 5
	
	adaboost1_all[[i]] <- boosting(CallArea ~ NbPriorSGenChange+ NbPriorSLotChange+ NbPriorERFLotChange+ NbPriorIWFLotChange+ NbPriorContLotNumChange+ NbPriorCalCurveChange+ TimeToComplain+ MostRecentConcentration+ FiveMostRecentConcentrationMean+ FiveMostRecentConcentrationSD+ TenMostRecentConcentrationMean+ TenMostRecentConcentrationSD+ TwoMostRecentConcentrationMean+ TwoMostRecentConcentrationSD + TimeSinceLastSGenChange + TimeSinceLastSLotChange + TimeSinceLastERFLotChange + TimeSinceLastIWFLotChange + TimeSinceLastContLotNumChange + TimeSinceLastCalCurveChange, data = ecust_dat_noNA[sub, ], 
					mfinal = 100, coeflearn = "Freund", control=rpart.control(minsplit=5))
	importance_sorted1[[i]] <- sort(adaboost1_all[[i]]$importance, decreasing=T)
	adaboost1_all_pred[[i]] <- predict.boosting(adaboost1_all[[i]], newdata= ecust_dat_noNA[-sub, ])
	adaboost1_all_conf[[i]] <- adaboost1_all_pred[[i]]$confusion
	adaboost1_all_err[[i]] <- adaboost1_all_pred[[i]]$error	
}


# comparing error evolution in training and test set
evol.train1 <- errorevol(adaboost1_all[[which.min(adaboost1_all_err)]],newdata= ecust_dat_noNA[sub, ])
evol.test1 <- errorevol(adaboost1_all[[which.min(adaboost1_all_err)]],newdata= ecust_dat_noNA[-sub, ])
pdf(paste0(assay_dir, assay, "_ML2_ADABOOST.pdf"), width=15, height=22)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(3,1))
hist(as.numeric(adaboost1_all_err), 20, xlab="AdaBoost error rate", xlim=c(0,.5), main="")
plot.errorevol(evol.test1,evol.train1)
plot(importance_sorted1[[which.min(adaboost1_all_err)]], type="h", xlab="Factors ranking", ylab="AdaBoost importance")
text(seq(1:length(importance_sorted1[[which.min(adaboost1_all_err)]])), importance_sorted1[[which.min(adaboost1_all_err)]], names(importance_sorted1[[which.min(adaboost1_all_err)]]), cex=.5, srt=90, pos=2)
dev.off()



# std 10-fold cross-validation
# adaboost.cv <- boosting.cv(CallArea ~ NbPriorSGenChange+ NbPriorSLotChange+ NbPriorERFLotChange+ NbPriorIWFLotChange+ NbPriorContLotNumChange+ NbPriorCalCurveChange+ TimeToComplain+ MostRecentConcentration+ FiveMostRecentConcentrationMean+ FiveMostRecentConcentrationSD, data = ecust_dat_noNA,
# 				 v=10, mfinal = 100, coeflearn = "Freund", control=rpart.control(minsplit=5))# adaboost.cv[-1] #   0.08196721
# 
# succRate <- 100 - 100* adaboost.cv[-1]$error
# succRate





save.image(paste0(assay_dir, "econn_cust_", assay, ".RData"))


					################################
					# CART: without TimeToComplain #
					################################


Nrep <- 2500
error.cart2 <- c()
cart_all2 <- tb2 <- list()
for(i in 1:Nrep){
	if(!(i %% 100)){
		print(paste0("Now doing: ", i, "/", Nrep))
	}
	l <- length(ecust_dat_noNA[,1])
	sub <- sample(1:l,2*l/3)
	mfinal <- 10 
	maxdepth <- 5
	
	cart_all2[[i]] <- tree(CallArea ~ NbPriorSGenChange+ NbPriorSLotChange+ NbPriorERFLotChange+ NbPriorIWFLotChange+ NbPriorContLotNumChange+ NbPriorCalCurveChange+ MostRecentConcentration+ FiveMostRecentConcentrationMean+ FiveMostRecentConcentrationSD+ TenMostRecentConcentrationMean+ TenMostRecentConcentrationSD+ TwoMostRecentConcentrationMean+ TwoMostRecentConcentrationSD + TimeSinceLastSGenChange + TimeSinceLastSLotChange + TimeSinceLastERFLotChange + TimeSinceLastIWFLotChange + TimeSinceLastContLotNumChange + TimeSinceLastCalCurveChange, data= ecust_dat_noNA[sub,])
	cart_all2_pred <- predict(cart_all2[[i]], newdata= ecust_dat_noNA[-sub, ],type="class")
	tb2[[i]] <- table(cart_all2_pred, ecust_dat_noNA$CallArea[-sub])
	error.cart2[i] <- 1-(sum(diag(tb2[[i]]))/sum(tb2[[i]]))
}

# best model
best_tb2 <- tb2[[which.min(error.cart2)]]
best_tb2
best_cart_all2 <- cart_all2[[which.min(error.cart2)]]
summary(best_cart_all2)
pdf(paste0(assay_dir, assay, "_ML3_CART.pdf"), width=15, height=18)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,1))
hist(error.cart2, 20, xlab="CART error rate", xlim=c(0,.5), main="")
plot(best_cart_all2)
text(best_cart_all2, pretty=0)
dev.off()



					####################################
					# adaboost: without TimeToComplain #
					####################################

Nrep <- 2500
error.adaboost2 <- c()
adaboost2_all <- adaboost2_all_conf <- adaboost2_all_err <- adaboost2_all_pred <- importance_sorted2 <- list()
for(i in 1:Nrep){
	if(!(i %% 10)){
		print(paste0("Now doing: ", i, "/", Nrep))
	}
	l <- length(ecust_dat_noNA[,1])
	sub <- sample(1:l,2*l/3)
	mfinal <- 10 
	maxdepth <- 5
	
	adaboost2_all[[i]] <- boosting(CallArea ~ NbPriorSGenChange+ NbPriorSLotChange+ NbPriorERFLotChange+ NbPriorIWFLotChange+ NbPriorContLotNumChange+ NbPriorCalCurveChange+ MostRecentConcentration+ FiveMostRecentConcentrationMean+ FiveMostRecentConcentrationSD+ TenMostRecentConcentrationMean+ TenMostRecentConcentrationSD+ TwoMostRecentConcentrationMean+ TwoMostRecentConcentrationSD + TimeSinceLastSGenChange + TimeSinceLastSLotChange + TimeSinceLastERFLotChange + TimeSinceLastIWFLotChange + TimeSinceLastContLotNumChange + TimeSinceLastCalCurveChange, data = ecust_dat_noNA[sub, ], 
					mfinal = 100, coeflearn = "Freund", control=rpart.control(minsplit=5))
	importance_sorted2[[i]] <- sort(adaboost2_all[[i]]$importance, decreasing=T)
	adaboost2_all_pred[[i]] <- predict.boosting(adaboost2_all[[i]], newdata= ecust_dat_noNA[-sub, ])
	adaboost2_all_conf[[i]] <- adaboost2_all_pred[[i]]$confusion
	adaboost2_all_err[[i]] <- adaboost2_all_pred[[i]]$error	
}


# comparing error evolution in training and test set
evol.train2 <- errorevol(adaboost2_all[[which.min(adaboost2_all_err)]],newdata= ecust_dat_noNA[sub, ])
evol.test2 <- errorevol(adaboost2_all[[which.min(adaboost2_all_err)]],newdata= ecust_dat_noNA[-sub, ])
pdf(paste0(assay_dir, assay, "_ML4_ADABOOST.pdf"), width=15, height=22)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(3,1))
hist(as.numeric(adaboost2_all_err), 20, xlab="AdaBoost error rate", xlim=c(0,.5), main="")
plot.errorevol(evol.test2,evol.train2)
plot(importance_sorted2[[which.min(adaboost2_all_err)]], type="h", xlab="Factors ranking", ylab="AdaBoost importance")
text(seq(1:length(importance_sorted2[[which.min(adaboost2_all_err)]])), importance_sorted2[[which.min(adaboost2_all_err)]], names(importance_sorted2[[which.min(adaboost2_all_err)]]), cex=.5, srt=90, pos=2)
dev.off()



# std 10-fold cross-validation
# adaboost.cv <- boosting.cv(CallArea ~ NbPriorSGenChange+ NbPriorSLotChange+ NbPriorERFLotChange+ NbPriorIWFLotChange+ NbPriorContLotNumChange+ NbPriorCalCurveChange+ TimeToComplain+ MostRecentConcentration+ FiveMostRecentConcentrationMean+ FiveMostRecentConcentrationSD, data = ecust_dat_noNA,
# 				 v=10, mfinal = 100, coeflearn = "Freund", control=rpart.control(minsplit=5))# adaboost.cv[-1] #   0.08196721
# 
# succRate <- 100 - 100* adaboost.cv[-1]$error
# succRate




################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################






####################################
# save and quit
save.image(paste0(assay_dir, "econn_cust_", assay, ".RData"))
q(save="no")
####################################
