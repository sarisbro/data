# This script was used in Shoji A., Aris-Brosou S., Culina A., Fayet A., Kirk H., Padget O., Juarez-Martinez I., Boyle D., Nakata T., Perrins C.M. and Guilford T. 2015. 
# Breeding phenology and winter activity predict subsequent breeding success in a trans-global migratory seabird. 
# Biology Letters 11:20150671.
# https://doi.org/10.1098/rsbl.2015.0671

library(adabag)
library(rpart)
library(mlbench)
library(lubridate)
library(MASS)
library(xts)
library(locfit)
library(Hmisc)  # for minor.tick()
library(doMC)
library(foreach)

rm(list = ls())


# load("coe_13.RData")

# read traits
#coe_trait <- read.csv("coe_trait.csv")
#coe_trait <- read.csv("coe_trait2.csv")
#coe_trait <- read.csv("coe_trait3.csv")
#coe_trait <- read.csv("coe_trait4.csv")
#coe_trait <- read.csv("coe_trait5_clean.csv")
coe_trait <- read.csv("coe_trait6.csv")


prior_lay.d <- prior_lay.m <- prior_lay.y <- numeric(length(coe_trait[,1]))
prior_hatch.d <- prior_hatch.m <- prior_hatch.y <- numeric(length(coe_trait[,1]))
prior_fledge.d <- prior_fledge.m <- prior_fledge.y <- numeric(length(coe_trait[,1]))
Colony_departure.d <- Colony_departure.m <- Colony_departure.y <- numeric(length(coe_trait[,1]))
Arrival_at_WG.d <- Arrival_at_WG.m <- Arrival_at_WG.y <- numeric(length(coe_trait[,1]))
WG_departure.d <- WG_departure.m <- WG_departure.y <- numeric(length(coe_trait[,1]))
Colony_arrival.d <- Colony_arrival.m <- Colony_arrival.y <- numeric(length(coe_trait[,1]))
post_hatch.d <- post_hatch.m <- post_hatch.y <- numeric(length(coe_trait[,1]))
post_fledge.d <- post_fledge.m <- post_fledge.y <- numeric(length(coe_trait[,1]))
for(i in 1:length(coe_trait[,1])){
	# prior_lay
	prior_lay.d[i] <- as.numeric(strsplit(as.character(coe_trait$prior_lay[i]), "/")[[1]][1])
	prior_lay.m[i] <- as.numeric(strsplit(as.character(coe_trait$prior_lay[i]), "/")[[1]][2])
	prior_lay.y[i] <- as.numeric(strsplit(as.character(coe_trait$prior_lay[i]), "/")[[1]][3])
	yday(dmy(coe_trait$prior_lay[i]))
	# prior_hatch
	prior_hatch.d[i] <- as.numeric(strsplit(as.character(coe_trait$prior_hatch[i]), "/")[[1]][1])
	prior_hatch.m[i] <- as.numeric(strsplit(as.character(coe_trait$prior_hatch[i]), "/")[[1]][2])
	prior_hatch.y[i] <- as.numeric(strsplit(as.character(coe_trait$prior_hatch[i]), "/")[[1]][3])
	# prior_fledge
	prior_fledge.d[i] <- as.numeric(strsplit(as.character(coe_trait$prior_fledge[i]), "/")[[1]][1])
	prior_fledge.m[i] <- as.numeric(strsplit(as.character(coe_trait$prior_fledge[i]), "/")[[1]][2])
	prior_fledge.y[i] <- as.numeric(strsplit(as.character(coe_trait$prior_fledge[i]), "/")[[1]][3])
	# Colony_departure
	Colony_departure.d[i] <- as.numeric(strsplit(as.character(coe_trait$Colony_departure[i]), "/")[[1]][1])
	Colony_departure.m[i] <- as.numeric(strsplit(as.character(coe_trait$Colony_departure[i]), "/")[[1]][2])
	Colony_departure.y[i] <- as.numeric(strsplit(as.character(coe_trait$Colony_departure[i]), "/")[[1]][3])
	# Arrival_at_WG
	Arrival_at_WG.d[i] <- as.numeric(strsplit(as.character(coe_trait$Arrival_of_WG[i]), "/")[[1]][1])
	Arrival_at_WG.m[i] <- as.numeric(strsplit(as.character(coe_trait$Arrival_of_WG[i]), "/")[[1]][2])
	Arrival_at_WG.y[i] <- as.numeric(strsplit(as.character(coe_trait$Arrival_of_WG[i]), "/")[[1]][3])
	# WG_departure
	WG_departure.d[i] <- as.numeric(strsplit(as.character(coe_trait$WG_departure[i]), "/")[[1]][1])
	WG_departure.m[i] <- as.numeric(strsplit(as.character(coe_trait$WG_departure[i]), "/")[[1]][2])
	WG_departure.y[i] <- as.numeric(strsplit(as.character(coe_trait$WG_departure[i]), "/")[[1]][3])
	# Colony_arrival
	Colony_arrival.d[i] <- as.numeric(strsplit(as.character(coe_trait$Colony_arrival[i]), "/")[[1]][1])
	Colony_arrival.m[i] <- as.numeric(strsplit(as.character(coe_trait$Colony_arrival[i]), "/")[[1]][2])
	Colony_arrival.y[i] <- as.numeric(strsplit(as.character(coe_trait$Colony_arrival[i]), "/")[[1]][3])
	# post_hatch
	post_hatch.d[i] <- as.numeric(strsplit(as.character(coe_trait$post_hatch[i]), "/")[[1]][1])
	post_hatch.m[i] <- as.numeric(strsplit(as.character(coe_trait$post_hatch[i]), "/")[[1]][2])
	post_hatch.y[i] <- as.numeric(strsplit(as.character(coe_trait$post_hatch[i]), "/")[[1]][3])
	# post_fledge
	post_fledge.d[i] <- as.numeric(strsplit(as.character(coe_trait$post_fledge[i]), "/")[[1]][1])
	post_fledge.m[i] <- as.numeric(strsplit(as.character(coe_trait$post_fledge[i]), "/")[[1]][2])
	post_fledge.y[i] <- as.numeric(strsplit(as.character(coe_trait$post_fledge[i]), "/")[[1]][3])
}
prior_lay.dd <- yday(dmy(coe_trait$prior_lay))
prior_hatch.dd <- yday(dmy(coe_trait$prior_hatch))
prior_fledge.dd <- yday(dmy(coe_trait$prior_fledge))
Colony_departure.dd <- yday(dmy(coe_trait$Colony_departure))
Arrival_at_WG.dd <- yday(dmy(coe_trait$Arrival_of_WG))
WG_departure.dd <- yday(dmy(coe_trait$WG_departure))
Colony_arrival.dd <- yday(dmy(coe_trait$Colony_arrival))
post_hatch.dd <- yday(dmy(coe_trait$post_hatch))
post_fledge.dd <- yday(dmy(coe_trait$post_fledge))

# cleanups for missing (pending) 2014 data
post_hatch.d[post_hatch.d == 2014] <- NA
post_fledge.d[post_fledge.d == 2014] <- NA


# convert dates to POSIX
coe_trait$prior_lay <- as.POSIXct(strptime(coe_trait$prior_lay, format = "%d/%m/%y"), tz = "GMT")
coe_trait$prior_hatch <- as.POSIXct(strptime(coe_trait$prior_hatch, format = "%d/%m/%y"), tz = "GMT")
coe_trait$prior_fledge <- as.POSIXct(strptime(coe_trait$prior_fledge, format = "%d/%m/%y"), tz = "GMT")
coe_trait$Colony_departure <- as.POSIXct(strptime(coe_trait$Colony_departure, format = "%d/%m/%y"), tz = "GMT")
coe_trait$Arrival_of_WG <- as.POSIXct(strptime(coe_trait$Arrival_of_WG, format = "%d/%m/%y"), tz = "GMT")
coe_trait$WG_departure <- as.POSIXct(strptime(coe_trait$WG_departure, format = "%d/%m/%y"), tz = "GMT")
coe_trait$Colony_arrival <- as.POSIXct(strptime(coe_trait$Colony_arrival, format = "%d/%m/%y"), tz = "GMT")
coe_trait$post_hatch <- as.POSIXct(strptime(coe_trait$post_hatch, format = "%d/%m/%y"), tz = "GMT")
coe_trait$post_fledge <- as.POSIXct(strptime(coe_trait$post_fledge, format = "%d/%m/%y"), tz = "GMT")

# change "SKIP?" and "2014"
#coe_trait$post_rs[which(coe_trait$post_rs == "SKIP?")] <- "SKIP"
#coe_trait$post_rs[which(coe_trait$post_rs == "2014")] <- NA

plot(coe_trait$prior_lay)
plot(coe_trait$prior_hatch)
plot(coe_trait$prior_fledge)
plot(coe_trait$Colony_departure)
plot(coe_trait$Arrival_of_WG)
plot(coe_trait$WG_departure)
plot(coe_trait$Colony_arrival)
plot(coe_trait$post_hatch)
plot(coe_trait$post_fledge)


rs_dates <- data.frame(yday(coe_trait$prior_lay),
						yday(coe_trait$prior_hatch),
						yday(coe_trait$prior_fledge),
						yday(coe_trait$Colony_departure),
						yday(coe_trait$Arrival_of_WG),
						yday(coe_trait$WG_departure),
						yday(coe_trait$Colony_arrival),
						yday(coe_trait$post_hatch),
						yday(coe_trait$post_fledge)
						)
						
rs_p_dates <- data.frame(yday(coe_trait$prior_lay),
						yday(coe_trait$prior_hatch),
						yday(coe_trait$prior_fledge),
						yday(coe_trait$Colony_departure),
						yday(coe_trait$Arrival_of_WG),
						yday(coe_trait$WG_departure),
						yday(coe_trait$Colony_arrival),
						yday(coe_trait$post_hatch),
						yday(coe_trait$post_fledge),
						as.double(post_rs_clean)
						)
						
						
#include RS and test correlation with features: H_0 there is different corr pattern: prior lay & prior htach not correl in same way with RS						
						
						
rs_dates <- rs_dates[complete.cases(rs_dates), ]					
#rs_p_dates <- rs_p_dates[complete.cases(rs_p_dates), ]					
colnames(rs_dates) <- c("prior lay","prior hatch","prior fledge","colony departure","arrival at WG","WG departure","colony arrival","post hatch","post fledge")
colnames(rs_p_dates) <- c("prior lay","prior hatch","prior fledge","colony departure","arrival at WG","WG departure","colony arrival","post hatch","post fledge","RS")

## http://personality-project.org/r/r.graphics.html
## the following code and figure is adapted from the help file for pairs 
##   put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
#first create a function (panel.cor)
     panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r = (cor(x, y))
         txt <- format(c(r, 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex * abs(r))
     }
     
# now use the function for the epi data. (see figure)
pdf("correl.pdf", width=12, height=8)
pairs(rs_dates, lower.panel=panel.smooth, upper.panel=panel.cor)
dev.off()
# run this one the entire rs_dates object, before running complete.cases()
pdf("correl_bis.pdf", width=12, height=8)
pairs(rs_dates, lower.panel=panel.smooth, na.action = na.omit, pch=20)
dev.off()

panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y, use="pairwise.complete.obs")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)

  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("P= ", txt2, sep = "")
  #if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

#library(corrplot)
#corrplot(rs_p_dates)
#library(corrgram)
#corrgram(rs_p_dates,lower.panel=panel.pts,text.panel=panel.cor)
#cor(rs_p_dates, use="pairwise.complete.obs")


pdf("correl_ter.pdf", width=12, height=8)
pairs(rs_p_dates, use="pairwise.complete.obs", lower.panel=panel.smooth, , upper.panel= panel.cor, pch=20)
dev.off()


length(coe_trait$post_rs[coe_trait$GLS_Year_end == 2010])
length(coe_trait$post_rs[coe_trait$GLS_Year_end == 2011])
length(coe_trait$post_rs[coe_trait$GLS_Year_end == 2012])
length(coe_trait$post_rs[coe_trait$GLS_Year_end == 2013])
length(coe_trait$post_rs[coe_trait$GLS_Year_end == 2014])



####################################################################
# adaptive boosting on dates recoded as number of days in the year #
####################################################################
post_rs_clean <- as.factor(as.character(coe_trait$post_rs))
boost.df <- data.frame(post_rs_clean, prior_lay.dd, prior_hatch.dd, prior_fledge.dd, Colony_departure.dd, 
											Arrival_at_WG.dd, WG_departure.dd, Colony_arrival.dd, post_hatch.dd, 
											post_fledge.dd)
colnames(boost.df) <- c("post_rs", "prior_lay_dd", "prior_hatch_dd", "prior_fledge_dd", "Colony_departure_dd",
											"Arrival_at_WG_dd", "WG_departure_dd", "Colony_arrival_dd", "post_hatch_dd",
											"post_fledge_dd")
# only keep rows for which RS is not NA
boost.df  <- boost.df[!is.na(boost.df$post_rs),]
RS.adaboost <- boosting(post_rs~., data= boost.df, mfinal = 5000, coeflearn = "Freund", control=rpart.control(minsplit=5))
#RS.adaboost
importance_sorted <- sort(RS.adaboost$importance, decreasing=T)
importance_sorted


# 10-fold cross-validation
RS.boostcv <- boosting.cv(post_rs~., data= boost.df, v=20, mfinal = 5000, coeflearn = "Freund", control=rpart.control(minsplit=5))
RS.boostcv[-1] #  0.3873874
#RS.boostcv <- boosting.cv(post_rs~., data= boost.df, mfinal = 10, coeflearn = "Breiman", control=rpart.control(minsplit=3))
#RS.boostcv <- boosting.cv(post_rs~., data= boost.df, mfinal = 100, coeflearn = "Zhu", control=rpart.control(minsplit=4))
#RS.boostcv[-1] # 0.407767

succRate <- 100 - 100*RS.boostcv[-1]$error


################################
# removing correlated features #
################################
boost.df <- data.frame(post_rs_clean, prior_hatch.dd, prior_fledge.dd, 
											Arrival_at_WG.dd, WG_departure.dd, post_hatch.dd, 
											post_fledge.dd)
colnames(boost.df) <- c("post_rs", "prior_hatch_dd", "prior_fledge_dd",
											"Arrival_at_WG_dd", "WG_departure_dd", "post_hatch_dd",
											"post_fledge_dd")
# only keep rows for which RS is not NA
boost.df  <- boost.df[!is.na(boost.df$post_rs),]
RS.adaboost <- boosting(post_rs~., data= boost.df, mfinal = 100, coeflearn = "Freund", control=rpart.control(minsplit=5))
#RS.adaboost
importance_sorted <- sort(RS.adaboost$importance, decreasing=T)
importance_sorted
# 10-fold cross-validation
RS.boostcv <- boosting.cv(post_rs~., data= boost.df, mfinal = 100, coeflearn = "Freund", control=rpart.control(minsplit=5))
RS.boostcv[-1] #  0.3873874




registerDoMC(cores=16)
NREPS <- 10
succRateVect <- c()
#for(i in 1:NREPS){
succRateVect <- foreach(i = 1:NREPS,.combine='rbind') %dopar% {
	print(i)
	RS.boostcv <- boosting.cv(post_rs~., data= boost.df, coeflearn = "Freund", control=rpart.control(minsplit=5))
	RS.boostcv[-1] #  0.3873874
	#succRateVect[i] <- 100 - 100*RS.boostcv[-1]$error
	return(100 - 100*RS.boostcv[-1]$error)
}

pdf("figure_S2.pdf", family="Times",width=6,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1)) 
#hist(succRateVect,20,xlab="Boosting success rate (10-fold cross-validation)",main="")
plot(locfit( ~ lp(succRateVect, nn=.8)),xlim=c(45,100),col="blue",xlab="Boosting success rate (10-fold cross-validation)",ylab="Density")
abline(v=50,col="red",lty=2)
dev.off()



pdf("figure_S1.pdf", family="Times",width=6,height=10)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,1)) 
plot((importance_sorted), xlab="Classifier features",type="l",col="gray",lwd=2,ylab="Relative boosting importance (%)")
text(seq(1:length(importance_sorted)),(importance_sorted), names(importance_sorted),cex=.7)
minor.tick(nx=2, ny=0, tick.ratio=1)
plot(cumsum(importance_sorted), xlab="Classifier features",type="l",col="gray",lwd=2,ylab="Cumulative boosting importance (%)")
text(seq(1:length(importance_sorted)), cumsum(importance_sorted), names(importance_sorted),cex=.7)
minor.tick(nx=2, ny=0, tick.ratio=1)
abline(h=50,col="red",lty=3); text(1,52,"50%",col="red",cex=.7)
abline(h= succRate,col="red",lty=2); text(1,(succRate+2),paste(succRate,"%"),col="red",cex=.7)
abline(h=90,col="red",lty=1); text(1,92,"90%",col="red",cex=.7)
dev.off()


pdf("figure_1.pdf", family="Times",width=6,height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1)) 
plot(cumsum(importance_sorted), xlab="Classifier features",type="l",col="gray",lwd=2,ylab="Cumulative boosting importance (%)")
text(seq(1:length(importance_sorted)), cumsum(importance_sorted), names(importance_sorted),cex=.7)
minor.tick(nx=2, ny=0, tick.ratio=1)
abline(h=50,col="red",lty=3); text(1,52,"50%",col="red",cex=.7)
abline(h= succRate,col="red",lty=2); text(1,(succRate+2),paste(succRate,"%"),col="red",cex=.7)
abline(h=90,col="red",lty=1); text(1,92,"90%",col="red",cex=.7)
dev.off()


##################
# visualizations #
##################
pdf("date_effects.pdf", family="Times", width=12, height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(3,3))
plot(boost.df$post_rs, boost.df$prior_lay_dd, main=paste("prior_lay_dd (I=",format(importance_sorted[[1]],digits=3),")",sep=""))
plot(boost.df$post_rs, boost.df$Colony_departure_dd, main=paste("Colony_departure_dd (I=",format(importance_sorted[[2]],digits=3),")",sep=""))
plot(boost.df$post_rs, boost.df$WG_departure_dd, main=paste("WG_departure_dd (I=",format(importance_sorted[[3]],digits=3),")",sep=""))
plot(boost.df$post_rs, boost.df$Arrival_at_WG_dd, main=paste("Arrival_at_WG_dd (I=",format(importance_sorted[[4]],digits=3),")",sep=""))
plot(boost.df$post_rs, boost.df$prior_fledge_dd, main=paste("prior_fledge_dd (I=",format(importance_sorted[[5]],digits=3),")",sep=""))
plot(boost.df$post_rs, boost.df$Colony_arrival_dd, main=paste("Colony_arrival_dd (I=",format(importance_sorted[[6]],digits=3),")",sep=""))
plot(boost.df$post_rs, boost.df$prior_hatch_dd, main=paste("prior_hatch_dd (I=",format(importance_sorted[[7]],digits=3),")",sep=""))
plot(boost.df$post_rs, boost.df$post_hatch_dd, main=paste("post_hatch_dd (I=",format(importance_sorted[[8]],digits=3),")",sep=""))
plot(boost.df$post_rs, boost.df$post_fledge_dd, main=paste("post_fledge_dd (I=",format(importance_sorted[[9]],digits=3),")",sep=""))
dev.off()

RS <- boost.df$post_rs
aov1 <- glm(boost.df$prior_lay_dd ~ RS); summary(aov1)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 135.13953    1.11436 121.271   <2e-16 ***
# RSEGG        -4.35693    1.88770  -2.308   0.0238 *  
# RSSKIP       -0.03953    2.56544  -0.015   0.9877    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 53.39693)
# 
#     Null deviance: 4201.4  on 75  degrees of freedom
# Residual deviance: 3898.0  on 73  degrees of freedom
#   (35 observations deleted due to missingness)
# AIC: 522.93

aov2 <- glm(boost.df$prior_hatch_dd ~ RS); summary(aov2)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  188.250      3.107  60.595   <2e-16 ***
# RSEGG         -7.490      5.009  -1.495    0.139    
# RSSKIP        -4.750      6.947  -0.684    0.496    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 386.0633)
# 
#     Null deviance: 28690  on 74  degrees of freedom
# Residual deviance: 27797  on 72  degrees of freedom
#   (36 observations deleted due to missingness)
# AIC: 664.48

aov3 <- glm(boost.df$prior_fledge_dd ~ RS); summary(aov3)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  246.625      4.963  49.689   <2e-16 ***
# RSEGG          5.217      7.467   0.699    0.488    
# RSSKIP        -8.125      9.152  -0.888    0.379    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 591.253)
# 
#     Null deviance: 30734  on 52  degrees of freedom
# Residual deviance: 29563  on 50  degrees of freedom
#   (58 observations deleted due to missingness)
# AIC: 493.58




# RS <- as.character(boost.df$post_rs)
# RS[RS == "SKIP"] <- 0
# RS[RS == "EGG"] <- 1
# RS[RS == "CHICK"] <- 2
# RS <- as.numeric(RS)
# aov0 <- aov(RS ~ boost.df$prior_lay_dd +  boost.df$prior_hatch_dd +  boost.df$prior_fledge_dd +  boost.df$Colony_departure_dd + 
# boost.df$Arrival_at_WG_dd +  boost.df$WG_departure_dd +  boost.df$Colony_arrival_dd +  boost.df$post_hatch_dd + 
# boost.df$post_fledge_dd, data=na.omit(boost.df))
# summary(aov0)
# mod.sel <- stepAIC(aov0)
# mod.sel$anova
# aov_best <- aov(RS ~ boost.df$prior_lay_dd + boost.df$prior_hatch_dd + boost.df$prior_fledge_dd + 
#     boost.df$Arrival_at_WG_dd + boost.df$WG_departure_dd + boost.df$Colony_arrival_dd + 
#     boost.df$post_hatch_dd, data=na.omit(boost.df))
# summary(aov_best)
# 
# plot.design(RS ~ boost.df$prior_lay_dd + boost.df$prior_hatch_dd + boost.df$prior_fledge_dd + 
#     boost.df$Arrival_at_WG_dd + boost.df$WG_departure_dd + boost.df$Colony_arrival_dd + 
#     boost.df$post_hatch_dd)


###########################
# 2010-14_summary_day_all #
###########################
#read files
sum_10 <- read.csv("2010_summary_day.csv", header=T)
sum_11 <- read.csv("2011_summary_day.csv", header=T)
sum_12 <- read.csv("2012_summary_day.csv", header=T)
sum_13 <- read.csv("2013_summary_day.csv", header=T)
sum_14 <- read.csv("2014_summary_day.csv", header=T)

#keep only first five columns
sum_10 <- sum_10[,1:5]
sum_11 <- sum_11[,1:5]
sum_12 <- sum_12[,1:5]
sum_13 <- sum_13[,1:5]
sum_14 <- sum_14[,1:5]

#rename columns
colnames(sum_10) <- c("ind","DayMonthYear","nb_0","nb_200","nb_1to199")
colnames(sum_11) <- c("ind","DayMonthYear","nb_0","nb_200","nb_1to199")
colnames(sum_12) <- c("ind","DayMonthYear","nb_0","nb_200","nb_1to199")
colnames(sum_13) <- c("ind","DayMonthYear","nb_0","nb_200","nb_1to199")
colnames(sum_14) <- c("ind","DayMonthYear","nb_0","nb_200","nb_1to199")

# check dates
plot(as.POSIXct(strptime(sum_10$DayMonthYear, format = "%d/%m/%Y"), tz = "GMT"))
plot(as.POSIXct(strptime(sum_11$DayMonthYear, format = "%d/%m/%Y"), tz = "GMT"))
plot(as.POSIXct(strptime(sum_12$DayMonthYear, format = "%d/%m/%Y"), tz = "GMT"))
plot(as.POSIXct(strptime(sum_13$DayMonthYear, format = "%d/%m/%Y"), tz = "GMT"))
plot(as.POSIXct(strptime(sum_14$DayMonthYear, format = "%d/%m/%Y"), tz = "GMT"))

#bind data
sum_all <- rbind(sum_10,sum_11,sum_12,sum_13,sum_14)
sum_all_dates <- as.POSIXct(strptime(sum_all$DayMonthYear, format = "%d/%m/%Y"), tz = "GMT")
plot(sum_all_dates) # one odd bird? [fb22822 -- keep it?]
sum_all <- cbind(sum_all, sum_all_dates)

levels(sum_all$ind)
nb_birds <- length(levels(sum_all$ind)) # 88 birds
nb_birds

# missing birds
levels(sum_all$ind)
coe_trait$Ring
setdiff(levels(sum_all$ind), coe_trait$Ring)
setdiff(coe_trait$Ring, levels(sum_all$ind))



#######################
# added on Sep 25, 2015
#######################

Isthmus <- c(0.59,
0.7,
0.69,
0.55,
0.6,
NA)

NorthHaven <- c(0.92,
0.86,
0.91,
1,
0.84,
0.63)

wilcox.test(Isthmus, NorthHaven, paired=T)

#########################
# added on July 15, 2015 
# sex effects
#########################

# 3 levels: f m UNKNOWN
# 1 is f; 2 is m; 3 is ?
sex <- character(length(sum_all$ind))
for(i in 1:length(sum_all$ind)){
	sex[i] <- coe_trait$sex[which(coe_trait$Ring == sum_all$ind[i])]
}

sex[is.na(sex)] <- "3"


# pdf("summary_day_all.pdf", family="Times",width=10,height=7)
# par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# for(i in levels(sum_all$ind)){
# 	plot(sum_all$sum_all_dates[sum_all$ind == i], sum_all$nb_0[sum_all$ind == i], main=paste("Bird ID: ",i," -- ",coe_trait$prior_rs[which(coe_trait$Ring == i & coe_trait$GLS_Year_start == format(min(sum_all$sum_all_dates[sum_all$ind == i]), "%Y"))],sep=""), type="l", xlab="Date", ylab="Total time on a given date",xaxt="n",col="red")
# 	axis.POSIXct(1, at= ceiling_date(sum_all$sum_all_dates[sum_all$ind == i], "month"), labels=format(sum_all$sum_all_dates[sum_all$ind == i], "%m/%d/%Y"),cex=.5)
# 	lines(sum_all$sum_all_dates[sum_all$ind == i], sum_all$nb_200[sum_all$ind == i],col="blue")
# 	lines(sum_all$sum_all_dates[sum_all$ind == i], sum_all$nb_1to199[sum_all$ind == i],col="black",lty=2)
# 	legend("topleft",c("Flight time","Rest time","Inbetween"),col=c("red","blue","black"),lty=c(1,1,2),cex=.7)
# 	abline(v=coe_trait$prior_lay[which(coe_trait$Ring == i)], lwd=2)
# 	abline(v=coe_trait$prior_hatch[which(coe_trait$Ring == i)], lwd=2)
# 	abline(v=coe_trait$prior_fledge[which(coe_trait$Ring == i)], lwd=2)
# 	abline(v=coe_trait$Colony_departure[which(coe_trait$Ring == i)], lwd=2)
# 	abline(v=coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)], lwd=2)
# 	abline(v=coe_trait$WG_departure[which(coe_trait$Ring == i)], lwd=2)
# 	abline(v=coe_trait$Colony_arrival[which(coe_trait$Ring == i)], lwd=2)
# 	abline(v=coe_trait$post_hatch[which(coe_trait$Ring == i)], lwd=2)
# 	abline(v=coe_trait$post_fledge[which(coe_trait$Ring == i)], lwd=2)
# 	text(coe_trait$prior_lay[which(coe_trait$Ring == i)],-25,"PL",offset = 0.0,pos=2,cex=.5)
# 	text(coe_trait$prior_hatch[which(coe_trait$Ring == i)],-25,"PH",offset = 0.0,pos=2,cex=.5)
# 	text(coe_trait$prior_fledge[which(coe_trait$Ring == i)],-25,"PF",offset = 0.0,pos=2,cex=.5)
# 	text(coe_trait$Colony_departure[which(coe_trait$Ring == i)],-25,"CD",offset = 0.0,pos=2,cex=.5)
# 	text(coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)],-25,"AWG",offset = 0.0,pos=2,cex=.5)
# 	text(coe_trait$WG_departure[which(coe_trait$Ring == i)],-25,"DWG",offset = 0.0,pos=2,cex=.5)
# 	text(coe_trait$Colony_arrival[which(coe_trait$Ring == i)],-25,"CA",offset = 0.0,pos=2,cex=.5)
# 	text(coe_trait$post_hatch[which(coe_trait$Ring == i)],-25,"PH",offset = 0.0,pos=2,cex=.5)
# 	text(coe_trait$post_fledge[which(coe_trait$Ring == i)],-25,"PL",offset = 0.0,pos=2,cex=.5)
# }
# dev.off()






##########################
# added August 24, 2015
# for JOS meeting 
##########################
library(locfit)

mean_0_ColDep2ArrWG <- tapply(sum_all$nb_0, sum_all$ind, mean)

#### aggregate data: ColDep2ArrWG  CHICK ####
vect_0_ColDep2ArrWG_CHICK <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_0[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]
				vect_0_ColDep2ArrWG_CHICK <- c(vect_0_ColDep2ArrWG_CHICK, y)
			}
		}
	}
}
mean_0_ColDep2ArrWG_CHICK <- mean(vect_0_ColDep2ArrWG_CHICK)
sd_0_ColDep2ArrWG_CHICK   <- sd(vect_0_ColDep2ArrWG_CHICK)

vect_200_ColDep2ArrWG_CHICK <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_200[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]
				vect_200_ColDep2ArrWG_CHICK <- c(vect_200_ColDep2ArrWG_CHICK, y)
			}
		}
	}
}
mean_200_ColDep2ArrWG_CHICK <- mean(vect_200_ColDep2ArrWG_CHICK)
sd_200_ColDep2ArrWG_CHICK   <- sd(vect_200_ColDep2ArrWG_CHICK)

vect_1to199_ColDep2ArrWG_CHICK <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_1to199[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]
				vect_1to199_ColDep2ArrWG_CHICK <- c(vect_1to199_ColDep2ArrWG_CHICK, y)
			}
		}
	}
}
mean_1to199_ColDep2ArrWG_CHICK <- mean(vect_1to199_ColDep2ArrWG_CHICK)
sd_1to199_ColDep2ArrWG_CHICK   <- sd(vect_1to199_ColDep2ArrWG_CHICK)


#### aggregate data: ColDep2ArrWG  EGG ####
vect_0_ColDep2ArrWG_EGG <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_0[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]
				vect_0_ColDep2ArrWG_EGG <- c(vect_0_ColDep2ArrWG_EGG, y)
			}
		}
	}
}
mean_0_ColDep2ArrWG_EGG <- mean(vect_0_ColDep2ArrWG_EGG)
sd_0_ColDep2ArrWG_EGG   <- sd(vect_0_ColDep2ArrWG_EGG)

vect_200_ColDep2ArrWG_EGG <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_200[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]
				vect_200_ColDep2ArrWG_EGG <- c(vect_200_ColDep2ArrWG_EGG, y)
			}
		}
	}
}
mean_200_ColDep2ArrWG_EGG <- mean(vect_200_ColDep2ArrWG_EGG)
sd_200_ColDep2ArrWG_EGG   <- sd(vect_200_ColDep2ArrWG_EGG)

vect_1to199_ColDep2ArrWG_EGG <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_1to199[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]
				vect_1to199_ColDep2ArrWG_EGG <- c(vect_1to199_ColDep2ArrWG_EGG, y)
			}
		}
	}
}
mean_1to199_ColDep2ArrWG_EGG <- mean(vect_1to199_ColDep2ArrWG_EGG)
sd_1to199_ColDep2ArrWG_EGG   <- sd(vect_1to199_ColDep2ArrWG_EGG)


#### aggregate data: ColDep2ArrWG  SKIP ####
vect_0_ColDep2ArrWG_SKIP <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_0[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "SKIP") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]
				vect_0_ColDep2ArrWG_SKIP <- c(vect_0_ColDep2ArrWG_SKIP, y)
			}
		}
	}
}
mean_0_ColDep2ArrWG_SKIP <- mean(vect_0_ColDep2ArrWG_SKIP)
sd_0_ColDep2ArrWG_SKIP   <- sd(vect_0_ColDep2ArrWG_SKIP)

vect_200_ColDep2ArrWG_SKIP <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_200[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "SKIP") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]
				vect_200_ColDep2ArrWG_SKIP <- c(vect_200_ColDep2ArrWG_SKIP, y)
			}
		}
	}
}
mean_200_ColDep2ArrWG_SKIP <- mean(vect_200_ColDep2ArrWG_SKIP)
sd_200_ColDep2ArrWG_SKIP   <- sd(vect_200_ColDep2ArrWG_SKIP)

vect_1to199_ColDep2ArrWG_SKIP <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_1to199[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "SKIP") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])]
				vect_1to199_ColDep2ArrWG_SKIP <- c(vect_1to199_ColDep2ArrWG_SKIP, y)
			}
		}
	}
}
mean_1to199_ColDep2ArrWG_SKIP <- mean(vect_1to199_ColDep2ArrWG_SKIP)
sd_1to199_ColDep2ArrWG_SKIP   <- sd(vect_1to199_ColDep2ArrWG_SKIP)



######

#### aggregate data: ArrWG2DepWG  CHICK ####
vect_0_ArrWG2DepWG_CHICK <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_0[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]
				vect_0_ArrWG2DepWG_CHICK <- c(vect_0_ArrWG2DepWG_CHICK, y)
			}
		}
	}
}
mean_0_ArrWG2DepWG_CHICK <- mean(vect_0_ArrWG2DepWG_CHICK)
sd_0_ArrWG2DepWG_CHICK   <- sd(vect_0_ArrWG2DepWG_CHICK)

vect_200_ArrWG2DepWG_CHICK <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_200[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]
				vect_200_ArrWG2DepWG_CHICK <- c(vect_200_ArrWG2DepWG_CHICK, y)
			}
		}
	}
}
mean_200_ArrWG2DepWG_CHICK <- mean(vect_200_ArrWG2DepWG_CHICK)
sd_200_ArrWG2DepWG_CHICK   <- sd(vect_200_ArrWG2DepWG_CHICK)

vect_1to199_ArrWG2DepWG_CHICK <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_1to199[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]
				vect_1to199_ArrWG2DepWG_CHICK <- c(vect_1to199_ArrWG2DepWG_CHICK, y)
			}
		}
	}
}
mean_1to199_ArrWG2DepWG_CHICK <- mean(vect_1to199_ArrWG2DepWG_CHICK)
sd_1to199_ArrWG2DepWG_CHICK   <- sd(vect_1to199_ArrWG2DepWG_CHICK)


#### aggregate data: ArrWG2DepWG  EGG ####
vect_0_ArrWG2DepWG_EGG <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_0[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]
				vect_0_ArrWG2DepWG_EGG <- c(vect_0_ArrWG2DepWG_EGG, y)
			}
		}
	}
}
mean_0_ArrWG2DepWG_EGG <- mean(vect_0_ArrWG2DepWG_EGG)
sd_0_ArrWG2DepWG_EGG   <- sd(vect_0_ArrWG2DepWG_EGG)

vect_200_ArrWG2DepWG_EGG <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_200[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]
				vect_200_ArrWG2DepWG_EGG <- c(vect_200_ArrWG2DepWG_EGG, y)
			}
		}
	}
}
mean_200_ArrWG2DepWG_EGG <- mean(vect_200_ArrWG2DepWG_EGG)
sd_200_ArrWG2DepWG_EGG   <- sd(vect_200_ArrWG2DepWG_EGG)

vect_1to199_ArrWG2DepWG_EGG <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_1to199[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]
				vect_1to199_ArrWG2DepWG_EGG <- c(vect_1to199_ArrWG2DepWG_EGG, y)
			}
		}
	}
}
mean_1to199_ArrWG2DepWG_EGG <- mean(vect_1to199_ArrWG2DepWG_EGG)
sd_1to199_ArrWG2DepWG_EGG   <- sd(vect_1to199_ArrWG2DepWG_EGG)


#### aggregate data: ArrWG2DepWG  SKIP ####
vect_0_ArrWG2DepWG_SKIP <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_0[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "SKIP") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]
				vect_0_ArrWG2DepWG_SKIP <- c(vect_0_ArrWG2DepWG_SKIP, y)
			}
		}
	}
}
mean_0_ArrWG2DepWG_SKIP <- mean(vect_0_ArrWG2DepWG_SKIP)
sd_0_ArrWG2DepWG_SKIP   <- sd(vect_0_ArrWG2DepWG_SKIP)

vect_200_ArrWG2DepWG_SKIP <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_200[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "SKIP") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]
				vect_200_ArrWG2DepWG_SKIP <- c(vect_200_ArrWG2DepWG_SKIP, y)
			}
		}
	}
}
mean_200_ArrWG2DepWG_SKIP <- mean(vect_200_ArrWG2DepWG_SKIP)
sd_200_ArrWG2DepWG_SKIP   <- sd(vect_200_ArrWG2DepWG_SKIP)

vect_1to199_ArrWG2DepWG_SKIP <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_1to199[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "SKIP") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]
				vect_1to199_ArrWG2DepWG_SKIP <- c(vect_1to199_ArrWG2DepWG_SKIP, y)
			}
		}
	}
}
mean_1to199_ArrWG2DepWG_SKIP <- mean(vect_1to199_ArrWG2DepWG_SKIP)
sd_1to199_ArrWG2DepWG_SKIP   <- sd(vect_1to199_ArrWG2DepWG_SKIP)



######

#### aggregate data: DepWG2ColArr  CHICK ####
vect_0_DepWG2ColArr_CHICK <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$WG_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_0[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				vect_0_DepWG2ColArr_CHICK <- c(vect_0_DepWG2ColArr_CHICK, y)
			}
		}
	}
}
mean_0_DepWG2ColArr_CHICK <- mean(vect_0_DepWG2ColArr_CHICK)
sd_0_DepWG2ColArr_CHICK   <- sd(vect_0_DepWG2ColArr_CHICK)

vect_200_DepWG2ColArr_CHICK <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$WG_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_200[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				vect_200_DepWG2ColArr_CHICK <- c(vect_200_DepWG2ColArr_CHICK, y)
			}
		}
	}
}
mean_200_DepWG2ColArr_CHICK <- mean(vect_200_DepWG2ColArr_CHICK)
sd_200_DepWG2ColArr_CHICK   <- sd(vect_200_DepWG2ColArr_CHICK)

vect_1to199_DepWG2ColArr_CHICK <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$WG_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_1to199[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				vect_1to199_DepWG2ColArr_CHICK <- c(vect_1to199_DepWG2ColArr_CHICK, y)
			}
		}
	}
}
mean_1to199_DepWG2ColArr_CHICK <- mean(vect_1to199_DepWG2ColArr_CHICK)
sd_1to199_DepWG2ColArr_CHICK   <- sd(vect_1to199_DepWG2ColArr_CHICK)


#### aggregate data: DepWG2ColArr  EGG ####
vect_0_DepWG2ColArr_EGG <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$WG_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_0[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				vect_0_DepWG2ColArr_EGG <- c(vect_0_DepWG2ColArr_EGG, y)
			}
		}
	}
}
mean_0_DepWG2ColArr_EGG <- mean(vect_0_DepWG2ColArr_EGG)
sd_0_DepWG2ColArr_EGG   <- sd(vect_0_DepWG2ColArr_EGG)

vect_200_DepWG2ColArr_EGG <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$WG_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_200[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				vect_200_DepWG2ColArr_EGG <- c(vect_200_DepWG2ColArr_EGG, y)
			}
		}
	}
}
mean_200_DepWG2ColArr_EGG <- mean(vect_200_DepWG2ColArr_EGG)
sd_200_DepWG2ColArr_EGG   <- sd(vect_200_DepWG2ColArr_EGG)

vect_1to199_DepWG2ColArr_EGG <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$WG_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_1to199[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				vect_1to199_DepWG2ColArr_EGG <- c(vect_1to199_DepWG2ColArr_EGG, y)
			}
		}
	}
}
mean_1to199_DepWG2ColArr_EGG <- mean(vect_1to199_DepWG2ColArr_EGG)
sd_1to199_DepWG2ColArr_EGG   <- sd(vect_1to199_DepWG2ColArr_EGG)


#### aggregate data: DepWG2ColArr  SKIP ####
vect_0_DepWG2ColArr_SKIP <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$WG_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_0[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "SKIP") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				vect_0_DepWG2ColArr_SKIP <- c(vect_0_DepWG2ColArr_SKIP, y)
			}
		}
	}
}
mean_0_DepWG2ColArr_SKIP <- mean(vect_0_DepWG2ColArr_SKIP)
sd_0_DepWG2ColArr_SKIP   <- sd(vect_0_DepWG2ColArr_SKIP)

vect_200_DepWG2ColArr_SKIP <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$WG_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_200[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "SKIP") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				vect_200_DepWG2ColArr_SKIP <- c(vect_200_DepWG2ColArr_SKIP, y)
			}
		}
	}
}
mean_200_DepWG2ColArr_SKIP <- mean(vect_200_DepWG2ColArr_SKIP)
sd_200_DepWG2ColArr_SKIP   <- sd(vect_200_DepWG2ColArr_SKIP)

vect_1to199_DepWG2ColArr_SKIP <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$WG_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_1to199[(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "SKIP") & (sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$WG_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				vect_1to199_DepWG2ColArr_SKIP <- c(vect_1to199_DepWG2ColArr_SKIP, y)
			}
		}
	}
}
mean_1to199_DepWG2ColArr_SKIP <- mean(vect_1to199_DepWG2ColArr_SKIP)
sd_1to199_DepWG2ColArr_SKIP   <- sd(vect_1to199_DepWG2ColArr_SKIP)







# barplot
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
	stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

y.mean <- c(
			mean_0_ColDep2ArrWG_CHICK, mean_1to199_ColDep2ArrWG_CHICK, mean_200_ColDep2ArrWG_CHICK,
			mean_0_ColDep2ArrWG_EGG, mean_1to199_ColDep2ArrWG_EGG, mean_200_ColDep2ArrWG_EGG,
			mean_0_ColDep2ArrWG_SKIP, mean_1to199_ColDep2ArrWG_SKIP, mean_200_ColDep2ArrWG_SKIP,

			mean_0_ArrWG2DepWG_CHICK, mean_1to199_ArrWG2DepWG_CHICK, mean_200_ArrWG2DepWG_CHICK,
			mean_0_ArrWG2DepWG_EGG, mean_1to199_ArrWG2DepWG_EGG, mean_200_ArrWG2DepWG_EGG,
			mean_0_ArrWG2DepWG_SKIP, mean_1to199_ArrWG2DepWG_SKIP, mean_200_ArrWG2DepWG_SKIP,

			mean_0_DepWG2ColArr_CHICK, mean_1to199_DepWG2ColArr_CHICK, mean_200_DepWG2ColArr_CHICK,
			mean_0_DepWG2ColArr_EGG, mean_1to199_DepWG2ColArr_EGG, mean_200_DepWG2ColArr_EGG,
			mean_0_DepWG2ColArr_SKIP, mean_1to199_DepWG2ColArr_SKIP, mean_200_DepWG2ColArr_SKIP
			
		)
y.sd <- c(
			sd_0_ColDep2ArrWG_CHICK, sd_1to199_ColDep2ArrWG_CHICK, sd_200_ColDep2ArrWG_CHICK,
			sd_0_ColDep2ArrWG_EGG, sd_1to199_ColDep2ArrWG_EGG, sd_200_ColDep2ArrWG_EGG,
			sd_0_ColDep2ArrWG_SKIP, sd_1to199_ColDep2ArrWG_SKIP, sd_200_ColDep2ArrWG_SKIP,

			sd_0_ArrWG2DepWG_CHICK, sd_1to199_ArrWG2DepWG_CHICK, sd_200_ArrWG2DepWG_CHICK,
			sd_0_ArrWG2DepWG_EGG, sd_1to199_ArrWG2DepWG_EGG, sd_200_ArrWG2DepWG_EGG,
			sd_0_ArrWG2DepWG_SKIP, sd_1to199_ArrWG2DepWG_SKIP, sd_200_ArrWG2DepWG_SKIP,

			sd_0_DepWG2ColArr_CHICK, sd_1to199_DepWG2ColArr_CHICK, sd_200_DepWG2ColArr_CHICK,
			sd_0_DepWG2ColArr_EGG, sd_1to199_DepWG2ColArr_EGG, sd_200_DepWG2ColArr_EGG,
			sd_0_DepWG2ColArr_SKIP, sd_1to199_DepWG2ColArr_SKIP, sd_200_DepWG2ColArr_SKIP
			
		)
yy <- matrix(y.mean,3,9,byrow=T)
y.sd <- matrix(y.sd,3,9,byrow=T)
barx <- barplot(yy, beside=T, names.arg=rep(c("C","E","S"), 9), axis.lty=1, xlab="RS / dates", ylab="Salt data (arbitrary units)")
#error.bar(barx,yy, y.sd)












# see http://colorbrewer2.org/ for colors
pdf("byDate_1.pdf", width=6, height=4)
plot(locfit( ~ lp(vect_200_ColDep2ArrWG)), col="#8856a7", lwd=2, xlab="Activity (salt data)", main="From colony departure to WG arrival")
lines(locfit( ~ lp(vect_1to199_ColDep2ArrWG)), col="#9ebcda", lwd=2)
lines(locfit( ~ lp(vect_0_ColDep2ArrWG)), col="#e0ecf4", lwd=2)
legend("topright", c("salt 0","salt 1-199","salt 200"), lty=1, lwd=2, col=c("#e0ecf4","#9ebcda","#8856a7"))
dev.off()

#### aggregate data: ArrWG2DepWG ####
vect_0_ArrWG2DepWG <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_0[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]
				vect_0_ArrWG2DepWG <- c(vect_0_ArrWG2DepWG, y)
			}
		}
	}
}
mean_0_ArrWG2DepWG <- mean(vect_0_ArrWG2DepWG)
sd_0_ArrWG2DepWG   <- sd(vect_0_ArrWG2DepWG)

vect_200_ArrWG2DepWG <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_200[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]
				vect_200_ArrWG2DepWG <- c(vect_200_ArrWG2DepWG, y)
			}
		}
	}
}
mean_200_ArrWG2DepWG <- mean(vect_200_ArrWG2DepWG)
sd_200_ArrWG2DepWG   <- sd(vect_200_ArrWG2DepWG)

vect_1to199_ArrWG2DepWG <- c()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]) > 0)){
				y <- sum_all$nb_1to199[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Arrival_of_WG[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$WG_departure[which(coe_trait$Ring == i)][year])]
				vect_1to199_ArrWG2DepWG <- c(vect_1to199_ArrWG2DepWG, y)
			}
		}
	}
}
mean_1to199_ArrWG2DepWG <- mean(vect_1to199_ArrWG2DepWG)
sd_1to199_ArrWG2DepWG   <- sd(vect_1to199_ArrWG2DepWG)

# see http://colorbrewer2.org/ for colors
pdf("byDate_2.pdf", width=6, height=4)
plot(locfit( ~ lp(vect_0_ArrWG2DepWG)), col="#e0ecf4", lwd=2, xlab="Activity (salt data)", main="From WG arrival to WG departure")
lines(locfit( ~ lp(vect_1to199_ArrWG2DepWG)), col="#9ebcda", lwd=2)
lines(locfit( ~ lp(vect_200_ArrWG2DepWG)), col="#8856a7", lwd=2)
legend("topright", c("salt 0","salt 1-199","salt 200"), lty=1, lwd=2, col=c("#e0ecf4","#9ebcda","#8856a7"))
dev.off()





################
# nb_0 results #
################

pdf("summary_day_all_by_year_nb0.pdf", family="Times",width=30,height=40)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(10,5))
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				#plot(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])], sum_all$nb_0[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])], main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""), type="l", xlab="Date", ylab="Total time on a given date",xaxt="n",col="red")
				#axis.POSIXct(1, at= ceiling_date(sum_all$sum_all_dates[sum_all$ind == i], "month"), labels=format(sum_all$sum_all_dates[sum_all$ind == i], "%m/%d/%Y"),cex=.5)
				dates <- sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				startdate <- dates[1]
				y <- sum_all$nb_0[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				y_ts <- xts(y, order.by  = dates)
				plot(y_ts, main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
			}
		}
	}
}
dev.off()

trend_chick <- list()
trend_egg <- list()
trend_skip <- list()
pdf("figure_Sn_nb0.pdf", family="Times",width=10,height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((!is.na(coe_trait$post_rs[which(coe_trait$Ring == i)])) & (length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				dates <- sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				startdate <- dates[1]
				y <- sum_all$nb_0[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				y_ts <- ts(y, frequency = 30, start = startdate)
				if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK"){
					trend_chick[[i]] <- decompose(y_ts,type = "additive")$trend
					plot(decompose(y_ts,type = "additive"))
					title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
					#plot(trend_chick[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
				}else{
					if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG"){
						trend_egg[[i]] <- decompose(y_ts,type = "additive")$trend
						plot(decompose(y_ts,type = "additive"))
						#plot(trend_egg[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
					}else{
						trend_skip[[i]] <- decompose(y_ts,type = "additive")$trend
						plot(decompose(y_ts,type = "additive"))
						title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						#plot(trend_skip[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
					}
				}
			}
		}
	}
}
dev.off()



pdf("trends_nb0.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(trend_chick[[1]], xlab="Time (arbitrary units)", ylab="Trend")
for(i in 2:length(trend_chick)){
	par(new=TRUE)
	plot(trend_chick[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="blue")
}
for(i in 2:length(trend_egg)){
	par(new=TRUE)
	plot(trend_egg[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="red")
}
for(i in 2:length(trend_skip)){
	par(new=TRUE)
	plot(trend_skip[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="black",lwd=2)
}
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()


trend_chick_mean <- sapply(trend_chick, mean, na.rm=T)
trend_egg_mean <- sapply(trend_egg, mean, na.rm=T)
trend_skip_mean <- sapply(trend_skip, mean, na.rm=T)
pdf("trends_nb0_mean.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
	plot(trend_chick_mean,type="l",ylim=c(0,2000),col="blue")
	lines(trend_egg_mean,col="red")
	lines(trend_skip_mean,lwd=2)
dev.off()





pdf("trends_nb0_means.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# chick
min_length_chick <- min(sapply(trend_chick, length))
mean_trend_chick <- numeric(min_length_chick)
ss_trend_chick <- numeric(min_length_chick)
sd_trend_chick <- numeric(min_length_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		mean_trend_chick[i] <- mean_trend_chick[i] + trend_chick[[j]][i]
		ss_trend_chick[i] <- ss_trend_chick[i] + (trend_chick[[j]][i])^2
	}
}
mean_trend_chick <- mean_trend_chick / length(trend_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		sd_trend_chick[i] <- sd_trend_chick[i] + (trend_chick[[j]][i] - mean_trend_chick[i])^2
	}
}
sd_trend_chick <- sqrt( sd_trend_chick / length(trend_chick))
# egg
min_length_egg <- min(sapply(trend_egg, length))
mean_trend_egg <- numeric(min_length_egg)
ss_trend_egg <- numeric(min_length_egg)
sd_trend_egg <- numeric(min_length_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		mean_trend_egg[i] <- mean_trend_egg[i] + trend_egg[[j]][i]
		ss_trend_egg[i] <- ss_trend_egg[i] + (trend_egg[[j]][i])^2
	}
}
mean_trend_egg <- mean_trend_egg / length(trend_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		sd_trend_egg[i] <- sd_trend_egg[i] + (trend_egg[[j]][i] - mean_trend_egg[i])^2
	}
}
sd_trend_egg <- sqrt( sd_trend_egg / length(trend_egg))
# skip
min_length_skip <- min(sapply(trend_skip, length))
mean_trend_skip <- numeric(min_length_skip)
ss_trend_skip <- numeric(min_length_skip)
sd_trend_skip <- numeric(min_length_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		mean_trend_skip[i] <- mean_trend_skip[i] + trend_skip[[j]][i]
		ss_trend_skip[i] <- ss_trend_skip[i] + (trend_skip[[j]][i])^2
	}
}
mean_trend_skip <- mean_trend_skip / length(trend_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		sd_trend_skip[i] <- sd_trend_skip[i] + (trend_skip[[j]][i] - mean_trend_skip[i])^2
	}
}
sd_trend_skip <- sqrt( sd_trend_skip / length(trend_skip))
#
plot(mean_trend_chick,col="blue",type="l",xlab="Relative time (days)",ylab="Mean trend",ylim=c(min(c(mean_trend_chick - sd_trend_chick, mean_trend_egg - sd_trend_egg, mean_trend_skip - sd_trend_skip),na.rm=T), max(c(mean_trend_chick + sd_trend_chick, mean_trend_egg + sd_trend_egg, mean_trend_skip + sd_trend_skip),na.rm=T)))
x_seq_chick <- seq(1:min_length_chick)
#segments(x_seq_chick, mean_trend_chick - sd_trend_chick, x_seq_chick, mean_trend_chick + sd_trend_chick,col="blue")
lines(mean_trend_chick - sd_trend_chick, col="blue", lty=2)
lines(mean_trend_chick + sd_trend_chick, col="blue", lty=2)
lines(mean_trend_egg,col="red")
x_seq_egg <- seq(1:min_length_egg)
#segments(x_seq_egg, mean_trend_egg - sd_trend_egg, x_seq_egg, mean_trend_egg + sd_trend_egg,col="red")
lines(mean_trend_egg - sd_trend_egg, col="red", lty=2)
lines(mean_trend_egg + sd_trend_egg, col="red", lty=2)
lines(mean_trend_skip,lwd=2)
x_seq_skip <- seq(1:min_length_skip)
#segments(x_seq_skip, mean_trend_skip - sd_trend_skip, x_seq_skip, mean_trend_skip + sd_trend_skip,col="black")
lines(mean_trend_skip - sd_trend_skip, col="black", lty=2)
lines(mean_trend_skip + sd_trend_skip, col="black", lty=2)
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()



pdf("ecdf_nb0.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
ecdf_mean_trend_chick <- ecdf(mean_trend_chick)
ecdf_mean_trend_egg <- ecdf(mean_trend_egg)
ecdf_mean_trend_skip <- ecdf(mean_trend_skip)
plot(ecdf_mean_trend_chick,col="blue",main="ECDFs",xlim=c(50,500))
lines(ecdf_mean_trend_egg,col="red")
lines(ecdf_mean_trend_skip)
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")

# manual computation of ECDFs...
sc <- sort(mean_trend_chick)
se <- sort(mean_trend_egg)
ss <- sort(mean_trend_skip)
nc <- length(sc); ne <- length(se); ns <- length(ss)
plot(sc, (1:nc)/nc, type = 's', ylim = c(0, 1), col="blue")    
lines(se, (1:ne)/ne, type = 's', col="red")    
lines(ss, (1:ns)/ns, type = 's', lwd=2)    
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")
dev.off()

# KS tests [Bonferonni required?]
ks.test(sc, se)
ks.test(sc, ss)

mean_trend_chick_0 <- mean_trend_chick
mean_trend_egg_0 <- mean_trend_egg
mean_trend_skip_0 <- mean_trend_skip


##################
# nb_200 results #
##################

pdf("summary_day_all_by_year_nb200.pdf", family="Times",width=30,height=40)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(10,5))
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				#plot(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])], sum_all$nb_200[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])], main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""), type="l", xlab="Date", ylab="Total time on a given date",xaxt="n",col="red")
				#axis.POSIXct(1, at= ceiling_date(sum_all$sum_all_dates[sum_all$ind == i], "month"), labels=format(sum_all$sum_all_dates[sum_all$ind == i], "%m/%d/%Y"),cex=.5)
				dates <- sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				startdate <- dates[1]
				y <- sum_all$nb_200[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				y_ts <- xts(y, order.by  = dates)
				plot(y_ts, main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
			}
		}
	}
}
dev.off()

trend_chick <- list()
trend_egg <- list()
trend_skip <- list()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((!is.na(coe_trait$post_rs[which(coe_trait$Ring == i)])) & (length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				dates <- sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				startdate <- dates[1]
				y <- sum_all$nb_200[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				y_ts <- ts(y, frequency = 30, start = startdate)
				if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK"){
					trend_chick[[i]] <- decompose(y_ts,type = "additive")$trend
					#plot(trend_chick[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
				}else{
					if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG"){
						trend_egg[[i]] <- decompose(y_ts,type = "additive")$trend
						#plot(trend_egg[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
					}else{
						trend_skip[[i]] <- decompose(y_ts,type = "additive")$trend
						#plot(trend_skip[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
					}
				}
			}
		}
	}
}

pdf("trends_nb200.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(trend_chick[[1]], xlab="Time (arbitrary units)", ylab="Trend")
for(i in 2:length(trend_chick)){
	par(new=TRUE)
	plot(trend_chick[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="blue")
}
for(i in 2:length(trend_egg)){
	par(new=TRUE)
	plot(trend_egg[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="red")
}
for(i in 2:length(trend_skip)){
	par(new=TRUE)
	plot(trend_skip[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="black",lwd=2)
}
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()


trend_chick_mean <- sapply(trend_chick, mean, na.rm=T)
trend_egg_mean <- sapply(trend_egg, mean, na.rm=T)
trend_skip_mean <- sapply(trend_skip, mean, na.rm=T)
pdf("trends_nb200_mean.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
	plot(trend_chick_mean,type="l",ylim=c(0,2000),col="blue")
	lines(trend_egg_mean,col="red")
	lines(trend_skip_mean,lwd=2)
dev.off()


pdf("trends_nb200_means.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# chick
min_length_chick <- min(sapply(trend_chick, length))
mean_trend_chick <- numeric(min_length_chick)
ss_trend_chick <- numeric(min_length_chick)
sd_trend_chick <- numeric(min_length_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		mean_trend_chick[i] <- mean_trend_chick[i] + trend_chick[[j]][i]
		ss_trend_chick[i] <- ss_trend_chick[i] + (trend_chick[[j]][i])^2
	}
}
mean_trend_chick <- mean_trend_chick / length(trend_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		sd_trend_chick[i] <- sd_trend_chick[i] + (trend_chick[[j]][i] - mean_trend_chick[i])^2
	}
}
sd_trend_chick <- sqrt( sd_trend_chick / length(trend_chick))
# egg
min_length_egg <- min(sapply(trend_egg, length))
mean_trend_egg <- numeric(min_length_egg)
ss_trend_egg <- numeric(min_length_egg)
sd_trend_egg <- numeric(min_length_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		mean_trend_egg[i] <- mean_trend_egg[i] + trend_egg[[j]][i]
		ss_trend_egg[i] <- ss_trend_egg[i] + (trend_egg[[j]][i])^2
	}
}
mean_trend_egg <- mean_trend_egg / length(trend_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		sd_trend_egg[i] <- sd_trend_egg[i] + (trend_egg[[j]][i] - mean_trend_egg[i])^2
	}
}
sd_trend_egg <- sqrt( sd_trend_egg / length(trend_egg))
# skip
min_length_skip <- min(sapply(trend_skip, length))
mean_trend_skip <- numeric(min_length_skip)
ss_trend_skip <- numeric(min_length_skip)
sd_trend_skip <- numeric(min_length_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		mean_trend_skip[i] <- mean_trend_skip[i] + trend_skip[[j]][i]
		ss_trend_skip[i] <- ss_trend_skip[i] + (trend_skip[[j]][i])^2
	}
}
mean_trend_skip <- mean_trend_skip / length(trend_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		sd_trend_skip[i] <- sd_trend_skip[i] + (trend_skip[[j]][i] - mean_trend_skip[i])^2
	}
}
sd_trend_skip <- sqrt( sd_trend_skip / length(trend_skip))
#
plot(mean_trend_chick,col="blue",type="l",xlab="Relative time (days)",ylab="Mean trend",ylim=c(min(c(mean_trend_chick - sd_trend_chick, mean_trend_egg - sd_trend_egg, mean_trend_skip - sd_trend_skip),na.rm=T), max(c(mean_trend_chick + sd_trend_chick, mean_trend_egg + sd_trend_egg, mean_trend_skip + sd_trend_skip),na.rm=T)))
x_seq_chick <- seq(1:min_length_chick)
#segments(x_seq_chick, mean_trend_chick - sd_trend_chick, x_seq_chick, mean_trend_chick + sd_trend_chick,col="blue")
lines(mean_trend_chick - sd_trend_chick, col="blue", lty=2)
lines(mean_trend_chick + sd_trend_chick, col="blue", lty=2)
lines(mean_trend_egg,col="red")
x_seq_egg <- seq(1:min_length_egg)
#segments(x_seq_egg, mean_trend_egg - sd_trend_egg, x_seq_egg, mean_trend_egg + sd_trend_egg,col="red")
lines(mean_trend_egg - sd_trend_egg, col="red", lty=2)
lines(mean_trend_egg + sd_trend_egg, col="red", lty=2)
lines(mean_trend_skip,lwd=2)
x_seq_skip <- seq(1:min_length_skip)
#segments(x_seq_skip, mean_trend_skip - sd_trend_skip, x_seq_skip, mean_trend_skip + sd_trend_skip,col="black")
lines(mean_trend_skip - sd_trend_skip, col="black", lty=2)
lines(mean_trend_skip + sd_trend_skip, col="black", lty=2)
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()

pdf("ecdf_nb200.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
ecdf_mean_trend_chick <- ecdf(mean_trend_chick)
ecdf_mean_trend_egg <- ecdf(mean_trend_egg)
ecdf_mean_trend_skip <- ecdf(mean_trend_skip)
plot(ecdf_mean_trend_chick,col="blue",main="ECDFs",xlim=c(100,750))
lines(ecdf_mean_trend_egg,col="red")
lines(ecdf_mean_trend_skip)
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")

# manual computation of ECDFs...
sc <- sort(mean_trend_chick)
se <- sort(mean_trend_egg)
ss <- sort(mean_trend_skip)
nc <- length(sc); ne <- length(se); ns <- length(ss)
plot(sc, (1:nc)/nc, type = 's', ylim = c(0, 1), col="blue")    
lines(se, (1:ne)/ne, type = 's', col="red")    
lines(ss, (1:ns)/ns, type = 's', lwd=2)    
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")
dev.off()

# KS tests [Bonferonni required?]
ks.test(sc, se)
ks.test(sc, ss)


mean_trend_chick_200 <- mean_trend_chick
mean_trend_egg_200 <- mean_trend_egg
mean_trend_skip_200 <- mean_trend_skip

#####################
# nb_1to199 results #
#####################

pdf("summary_day_all_by_year_nb1to199.pdf", family="Times",width=30,height=40)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(10,5))
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				#plot(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])], sum_all$nb_1to199[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])], main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""), type="l", xlab="Date", ylab="Total time on a given date",xaxt="n",col="red")
				#axis.POSIXct(1, at= ceiling_date(sum_all$sum_all_dates[sum_all$ind == i], "month"), labels=format(sum_all$sum_all_dates[sum_all$ind == i], "%m/%d/%Y"),cex=.5)
				dates <- sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				startdate <- dates[1]
				y <- sum_all$nb_1to199[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				y_ts <- xts(y, order.by  = dates)
				plot(y_ts, main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
			}
		}
	}
}
dev.off()

trend_chick <- list()
trend_egg <- list()
trend_skip <- list()
for(i in levels(sum_all$ind)){
	n_years <- length(which(coe_trait$Ring == i))
	for(year in 1:n_years){
		if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
			if((!is.na(coe_trait$post_rs[which(coe_trait$Ring == i)])) & (length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
				dates <- sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				startdate <- dates[1]
				y <- sum_all$nb_1to199[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
				y_ts <- ts(y, frequency = 30, start = startdate)
				if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK"){
					trend_chick[[i]] <- decompose(y_ts,type = "additive")$trend
					#plot(trend_chick[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
				}else{
					if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG"){
						trend_egg[[i]] <- decompose(y_ts,type = "additive")$trend
						#plot(trend_egg[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
					}else{
						trend_skip[[i]] <- decompose(y_ts,type = "additive")$trend
						#plot(trend_skip[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
					}
				}
			}
		}
	}
}

pdf("trends_nb1to199.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(trend_chick[[1]], xlab="Time (arbitrary units)", ylab="Trend")
for(i in 2:length(trend_chick)){
	par(new=TRUE)
	plot(trend_chick[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="blue")
}
for(i in 2:length(trend_egg)){
	par(new=TRUE)
	plot(trend_egg[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="red")
}
for(i in 2:length(trend_skip)){
	par(new=TRUE)
	plot(trend_skip[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="black",lwd=2)
}
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()


trend_chick_mean <- sapply(trend_chick, mean, na.rm=T)
trend_egg_mean <- sapply(trend_egg, mean, na.rm=T)
trend_skip_mean <- sapply(trend_skip, mean, na.rm=T)
pdf("trends_nb1to199_mean.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
	plot(trend_chick_mean,type="l",ylim=c(0,2000),col="blue")
	lines(trend_egg_mean,col="red")
	lines(trend_skip_mean,lwd=2)
dev.off()



pdf("trends_nb1to199_means.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# chick
min_length_chick <- min(sapply(trend_chick, length))
mean_trend_chick <- numeric(min_length_chick)
ss_trend_chick <- numeric(min_length_chick)
sd_trend_chick <- numeric(min_length_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		mean_trend_chick[i] <- mean_trend_chick[i] + trend_chick[[j]][i]
		ss_trend_chick[i] <- ss_trend_chick[i] + (trend_chick[[j]][i])^2
	}
}
mean_trend_chick <- mean_trend_chick / length(trend_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		sd_trend_chick[i] <- sd_trend_chick[i] + (trend_chick[[j]][i] - mean_trend_chick[i])^2
	}
}
sd_trend_chick <- sqrt( sd_trend_chick / length(trend_chick))
# egg
min_length_egg <- min(sapply(trend_egg, length))
mean_trend_egg <- numeric(min_length_egg)
ss_trend_egg <- numeric(min_length_egg)
sd_trend_egg <- numeric(min_length_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		mean_trend_egg[i] <- mean_trend_egg[i] + trend_egg[[j]][i]
		ss_trend_egg[i] <- ss_trend_egg[i] + (trend_egg[[j]][i])^2
	}
}
mean_trend_egg <- mean_trend_egg / length(trend_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		sd_trend_egg[i] <- sd_trend_egg[i] + (trend_egg[[j]][i] - mean_trend_egg[i])^2
	}
}
sd_trend_egg <- sqrt( sd_trend_egg / length(trend_egg))
# skip
min_length_skip <- min(sapply(trend_skip, length))
mean_trend_skip <- numeric(min_length_skip)
ss_trend_skip <- numeric(min_length_skip)
sd_trend_skip <- numeric(min_length_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		mean_trend_skip[i] <- mean_trend_skip[i] + trend_skip[[j]][i]
		ss_trend_skip[i] <- ss_trend_skip[i] + (trend_skip[[j]][i])^2
	}
}
mean_trend_skip <- mean_trend_skip / length(trend_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		sd_trend_skip[i] <- sd_trend_skip[i] + (trend_skip[[j]][i] - mean_trend_skip[i])^2
	}
}
sd_trend_skip <- sqrt( sd_trend_skip / length(trend_skip))
#
plot(mean_trend_chick,col="blue",type="l",xlab="Relative time (days)",ylab="Mean trend",ylim=c(min(c(mean_trend_chick - sd_trend_chick, mean_trend_egg - sd_trend_egg, mean_trend_skip - sd_trend_skip),na.rm=T), max(c(mean_trend_chick + sd_trend_chick, mean_trend_egg + sd_trend_egg, mean_trend_skip + sd_trend_skip),na.rm=T)))
x_seq_chick <- seq(1:min_length_chick)
#segments(x_seq_chick, mean_trend_chick - sd_trend_chick, x_seq_chick, mean_trend_chick + sd_trend_chick,col="blue")
lines(mean_trend_chick - sd_trend_chick, col="blue", lty=2)
lines(mean_trend_chick + sd_trend_chick, col="blue", lty=2)
lines(mean_trend_egg,col="red")
x_seq_egg <- seq(1:min_length_egg)
#segments(x_seq_egg, mean_trend_egg - sd_trend_egg, x_seq_egg, mean_trend_egg + sd_trend_egg,col="red")
lines(mean_trend_egg - sd_trend_egg, col="red", lty=2)
lines(mean_trend_egg + sd_trend_egg, col="red", lty=2)
lines(mean_trend_skip,lwd=2)
x_seq_skip <- seq(1:min_length_skip)
#segments(x_seq_skip, mean_trend_skip - sd_trend_skip, x_seq_skip, mean_trend_skip + sd_trend_skip,col="black")
lines(mean_trend_skip - sd_trend_skip, col="black", lty=2)
lines(mean_trend_skip + sd_trend_skip, col="black", lty=2)
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()

pdf("ecdf_nb1to199.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
ecdf_mean_trend_chick <- ecdf(mean_trend_chick)
ecdf_mean_trend_egg <- ecdf(mean_trend_egg)
ecdf_mean_trend_skip <- ecdf(mean_trend_skip)
plot(ecdf_mean_trend_chick,col="blue",main="ECDFs",xlim=c(650,875))
lines(ecdf_mean_trend_egg,col="red")
lines(ecdf_mean_trend_skip)
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")

# manual computation of ECDFs...
sc <- sort(mean_trend_chick)
se <- sort(mean_trend_egg)
ss <- sort(mean_trend_skip)
nc <- length(sc); ne <- length(se); ns <- length(ss)
plot(sc, (1:nc)/nc, type = 's', ylim = c(0, 1), col="blue")    
lines(se, (1:ne)/ne, type = 's', col="red")    
lines(ss, (1:ns)/ns, type = 's', lwd=2)    
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")
dev.off()

# KS tests [Bonferonni required?]
ks.test(sc, se)
ks.test(sc, ss)

mean_trend_chick_1_199 <- mean_trend_chick
mean_trend_egg_1_199 <- mean_trend_egg
mean_trend_skip_1_199 <- mean_trend_skip


######################################################

############################
# nb_0 results for females #
############################

trend_chick <- list()
trend_egg <- list()
trend_skip <- list()
pdf("figure_Sn_nb0_females.pdf", family="Times",width=10,height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
for(i in levels(sum_all$ind)){
	if(sex[which(sum_all$ind == i)[1]] == "1"){
		n_years <- length(which(coe_trait$Ring == i))
		for(year in 1:n_years){
			if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
				if((!is.na(coe_trait$post_rs[which(coe_trait$Ring == i)])) & (length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
					dates <- sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
					startdate <- dates[1]
					y <- sum_all$nb_0[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
					y_ts <- ts(y, frequency = 30, start = startdate)
					if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK"){
						trend_chick[[i]] <- decompose(y_ts,type = "additive")$trend
						plot(decompose(y_ts,type = "additive"))
						title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						#plot(trend_chick[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
					}else{
						if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG"){
							trend_egg[[i]] <- decompose(y_ts,type = "additive")$trend
							plot(decompose(y_ts,type = "additive"))
							#plot(trend_egg[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
							title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						}else{
							trend_skip[[i]] <- decompose(y_ts,type = "additive")$trend
							plot(decompose(y_ts,type = "additive"))
							title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
							#plot(trend_skip[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						}
					}
				}
			}
		}
	}
}
dev.off()



pdf("trends_nb0_females.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(trend_chick[[1]], xlab="Time (arbitrary units)", ylab="Trend")
for(i in 2:length(trend_chick)){
	par(new=TRUE)
	plot(trend_chick[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="blue")
}
for(i in 2:length(trend_egg)){
	par(new=TRUE)
	plot(trend_egg[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="red")
}
for(i in 2:length(trend_skip)){
	par(new=TRUE)
	plot(trend_skip[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="black",lwd=2)
}
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()


trend_chick_mean <- sapply(trend_chick, mean, na.rm=T)
trend_egg_mean <- sapply(trend_egg, mean, na.rm=T)
trend_skip_mean <- sapply(trend_skip, mean, na.rm=T)
pdf("trends_nb0_mean_females.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
	plot(trend_chick_mean,type="l",ylim=c(0,2000),col="blue")
	lines(trend_egg_mean,col="red")
	lines(trend_skip_mean,lwd=2)
dev.off()





pdf("trends_nb0_means_females.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# chick
min_length_chick <- min(sapply(trend_chick, length))
mean_trend_chick <- numeric(min_length_chick)
ss_trend_chick <- numeric(min_length_chick)
sd_trend_chick <- numeric(min_length_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		mean_trend_chick[i] <- mean_trend_chick[i] + trend_chick[[j]][i]
		ss_trend_chick[i] <- ss_trend_chick[i] + (trend_chick[[j]][i])^2
	}
}
mean_trend_chick <- mean_trend_chick / length(trend_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		sd_trend_chick[i] <- sd_trend_chick[i] + (trend_chick[[j]][i] - mean_trend_chick[i])^2
	}
}
sd_trend_chick <- sqrt( sd_trend_chick / length(trend_chick))
# egg
min_length_egg <- min(sapply(trend_egg, length))
mean_trend_egg <- numeric(min_length_egg)
ss_trend_egg <- numeric(min_length_egg)
sd_trend_egg <- numeric(min_length_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		mean_trend_egg[i] <- mean_trend_egg[i] + trend_egg[[j]][i]
		ss_trend_egg[i] <- ss_trend_egg[i] + (trend_egg[[j]][i])^2
	}
}
mean_trend_egg <- mean_trend_egg / length(trend_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		sd_trend_egg[i] <- sd_trend_egg[i] + (trend_egg[[j]][i] - mean_trend_egg[i])^2
	}
}
sd_trend_egg <- sqrt( sd_trend_egg / length(trend_egg))
# skip
min_length_skip <- min(sapply(trend_skip, length))
mean_trend_skip <- numeric(min_length_skip)
ss_trend_skip <- numeric(min_length_skip)
sd_trend_skip <- numeric(min_length_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		mean_trend_skip[i] <- mean_trend_skip[i] + trend_skip[[j]][i]
		ss_trend_skip[i] <- ss_trend_skip[i] + (trend_skip[[j]][i])^2
	}
}
mean_trend_skip <- mean_trend_skip / length(trend_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		sd_trend_skip[i] <- sd_trend_skip[i] + (trend_skip[[j]][i] - mean_trend_skip[i])^2
	}
}
sd_trend_skip <- sqrt( sd_trend_skip / length(trend_skip))
#
plot(mean_trend_chick,col="blue",type="l",xlab="Relative time (days)",ylab="Mean trend",ylim=c(min(c(mean_trend_chick - sd_trend_chick, mean_trend_egg - sd_trend_egg, mean_trend_skip - sd_trend_skip),na.rm=T), max(c(mean_trend_chick + sd_trend_chick, mean_trend_egg + sd_trend_egg, mean_trend_skip + sd_trend_skip),na.rm=T)))
x_seq_chick <- seq(1:min_length_chick)
#segments(x_seq_chick, mean_trend_chick - sd_trend_chick, x_seq_chick, mean_trend_chick + sd_trend_chick,col="blue")
lines(mean_trend_chick - sd_trend_chick, col="blue", lty=2)
lines(mean_trend_chick + sd_trend_chick, col="blue", lty=2)
lines(mean_trend_egg,col="red")
x_seq_egg <- seq(1:min_length_egg)
#segments(x_seq_egg, mean_trend_egg - sd_trend_egg, x_seq_egg, mean_trend_egg + sd_trend_egg,col="red")
lines(mean_trend_egg - sd_trend_egg, col="red", lty=2)
lines(mean_trend_egg + sd_trend_egg, col="red", lty=2)
lines(mean_trend_skip,lwd=2)
x_seq_skip <- seq(1:min_length_skip)
#segments(x_seq_skip, mean_trend_skip - sd_trend_skip, x_seq_skip, mean_trend_skip + sd_trend_skip,col="black")
lines(mean_trend_skip - sd_trend_skip, col="black", lty=2)
lines(mean_trend_skip + sd_trend_skip, col="black", lty=2)
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()



pdf("ecdf_nb0_females.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
ecdf_mean_trend_chick <- ecdf(mean_trend_chick)
ecdf_mean_trend_egg <- ecdf(mean_trend_egg)
ecdf_mean_trend_skip <- ecdf(mean_trend_skip)
plot(ecdf_mean_trend_chick,col="blue",main="ECDFs",xlim=c(50,500))
lines(ecdf_mean_trend_egg,col="red")
lines(ecdf_mean_trend_skip)
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")

# manual computation of ECDFs...
sc <- sort(mean_trend_chick)
se <- sort(mean_trend_egg)
ss <- sort(mean_trend_skip)
nc <- length(sc); ne <- length(se); ns <- length(ss)
plot(sc, (1:nc)/nc, type = 's', ylim = c(0, 1), col="blue")    
lines(se, (1:ne)/ne, type = 's', col="red")    
lines(ss, (1:ns)/ns, type = 's', lwd=2)    
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")
dev.off()

# KS tests [Bonferonni required?]
ks.test(sc, se)
ks.test(sc, ss)

sc_females <- sc
se_females <- se
ss_females <- ss

mean_trend_chick_0 <- mean_trend_chick
mean_trend_egg_0 <- mean_trend_egg
mean_trend_skip_0 <- mean_trend_skip


############################
# nb_0 results for males #
############################

trend_chick <- list()
trend_egg <- list()
trend_skip <- list()
pdf("figure_Sn_nb0_males.pdf", family="Times",width=10,height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
for(i in levels(sum_all$ind)){
	if(sex[which(sum_all$ind == i)[1]] == "2"){
		n_years <- length(which(coe_trait$Ring == i))
		for(year in 1:n_years){
			if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
				if((!is.na(coe_trait$post_rs[which(coe_trait$Ring == i)])) & (length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
					dates <- sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
					startdate <- dates[1]
					y <- sum_all$nb_0[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
					y_ts <- ts(y, frequency = 30, start = startdate)
					if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK"){
						trend_chick[[i]] <- decompose(y_ts,type = "additive")$trend
						plot(decompose(y_ts,type = "additive"))
						title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						#plot(trend_chick[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
					}else{
						if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG"){
							trend_egg[[i]] <- decompose(y_ts,type = "additive")$trend
							plot(decompose(y_ts,type = "additive"))
							#plot(trend_egg[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
							title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						}else{
							trend_skip[[i]] <- decompose(y_ts,type = "additive")$trend
							plot(decompose(y_ts,type = "additive"))
							title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
							#plot(trend_skip[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						}
					}
				}
			}
		}
	}
}
dev.off()



pdf("trends_nb0_males.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(trend_chick[[1]], xlab="Time (arbitrary units)", ylab="Trend")
for(i in 2:length(trend_chick)){
	par(new=TRUE)
	plot(trend_chick[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="blue")
}
for(i in 2:length(trend_egg)){
	par(new=TRUE)
	plot(trend_egg[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="red")
}
for(i in 2:length(trend_skip)){
	par(new=TRUE)
	plot(trend_skip[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="black",lwd=2)
}
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()


trend_chick_mean <- sapply(trend_chick, mean, na.rm=T)
trend_egg_mean <- sapply(trend_egg, mean, na.rm=T)
trend_skip_mean <- sapply(trend_skip, mean, na.rm=T)
pdf("trends_nb0_mean_males.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
	plot(trend_chick_mean,type="l",ylim=c(0,2000),col="blue")
	lines(trend_egg_mean,col="red")
	lines(trend_skip_mean,lwd=2)
dev.off()





pdf("trends_nb0_means_males.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# chick
min_length_chick <- min(sapply(trend_chick, length))
mean_trend_chick <- numeric(min_length_chick)
ss_trend_chick <- numeric(min_length_chick)
sd_trend_chick <- numeric(min_length_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		mean_trend_chick[i] <- mean_trend_chick[i] + trend_chick[[j]][i]
		ss_trend_chick[i] <- ss_trend_chick[i] + (trend_chick[[j]][i])^2
	}
}
mean_trend_chick <- mean_trend_chick / length(trend_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		sd_trend_chick[i] <- sd_trend_chick[i] + (trend_chick[[j]][i] - mean_trend_chick[i])^2
	}
}
sd_trend_chick <- sqrt( sd_trend_chick / length(trend_chick))
# egg
min_length_egg <- min(sapply(trend_egg, length))
mean_trend_egg <- numeric(min_length_egg)
ss_trend_egg <- numeric(min_length_egg)
sd_trend_egg <- numeric(min_length_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		mean_trend_egg[i] <- mean_trend_egg[i] + trend_egg[[j]][i]
		ss_trend_egg[i] <- ss_trend_egg[i] + (trend_egg[[j]][i])^2
	}
}
mean_trend_egg <- mean_trend_egg / length(trend_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		sd_trend_egg[i] <- sd_trend_egg[i] + (trend_egg[[j]][i] - mean_trend_egg[i])^2
	}
}
sd_trend_egg <- sqrt( sd_trend_egg / length(trend_egg))
# skip
min_length_skip <- min(sapply(trend_skip, length))
mean_trend_skip <- numeric(min_length_skip)
ss_trend_skip <- numeric(min_length_skip)
sd_trend_skip <- numeric(min_length_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		mean_trend_skip[i] <- mean_trend_skip[i] + trend_skip[[j]][i]
		ss_trend_skip[i] <- ss_trend_skip[i] + (trend_skip[[j]][i])^2
	}
}
mean_trend_skip <- mean_trend_skip / length(trend_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		sd_trend_skip[i] <- sd_trend_skip[i] + (trend_skip[[j]][i] - mean_trend_skip[i])^2
	}
}
sd_trend_skip <- sqrt( sd_trend_skip / length(trend_skip))
#
plot(mean_trend_chick,col="blue",type="l",xlab="Relative time (days)",ylab="Mean trend",ylim=c(min(c(mean_trend_chick - sd_trend_chick, mean_trend_egg - sd_trend_egg, mean_trend_skip - sd_trend_skip),na.rm=T), max(c(mean_trend_chick + sd_trend_chick, mean_trend_egg + sd_trend_egg, mean_trend_skip + sd_trend_skip),na.rm=T)))
x_seq_chick <- seq(1:min_length_chick)
#segments(x_seq_chick, mean_trend_chick - sd_trend_chick, x_seq_chick, mean_trend_chick + sd_trend_chick,col="blue")
lines(mean_trend_chick - sd_trend_chick, col="blue", lty=2)
lines(mean_trend_chick + sd_trend_chick, col="blue", lty=2)
lines(mean_trend_egg,col="red")
x_seq_egg <- seq(1:min_length_egg)
#segments(x_seq_egg, mean_trend_egg - sd_trend_egg, x_seq_egg, mean_trend_egg + sd_trend_egg,col="red")
lines(mean_trend_egg - sd_trend_egg, col="red", lty=2)
lines(mean_trend_egg + sd_trend_egg, col="red", lty=2)
lines(mean_trend_skip,lwd=2)
x_seq_skip <- seq(1:min_length_skip)
#segments(x_seq_skip, mean_trend_skip - sd_trend_skip, x_seq_skip, mean_trend_skip + sd_trend_skip,col="black")
lines(mean_trend_skip - sd_trend_skip, col="black", lty=2)
lines(mean_trend_skip + sd_trend_skip, col="black", lty=2)
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()



pdf("ecdf_nb0_males.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
ecdf_mean_trend_chick <- ecdf(mean_trend_chick)
ecdf_mean_trend_egg <- ecdf(mean_trend_egg)
ecdf_mean_trend_skip <- ecdf(mean_trend_skip)
plot(ecdf_mean_trend_chick,col="blue",main="ECDFs",xlim=c(50,500))
lines(ecdf_mean_trend_egg,col="red")
lines(ecdf_mean_trend_skip)
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")

# manual computation of ECDFs...
sc <- sort(mean_trend_chick)
se <- sort(mean_trend_egg)
ss <- sort(mean_trend_skip)
nc <- length(sc); ne <- length(se); ns <- length(ss)
plot(sc, (1:nc)/nc, type = 's', ylim = c(0, 1), col="blue")    
lines(se, (1:ne)/ne, type = 's', col="red")    
lines(ss, (1:ns)/ns, type = 's', lwd=2)    
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")
dev.off()

# KS tests [Bonferonni required?]
ks.test(sc, se)
ks.test(sc, ss)

sc_males <- sc
se_males <- se
ss_males <- ss

mean_trend_chick_0 <- mean_trend_chick
mean_trend_egg_0 <- mean_trend_egg
mean_trend_skip_0 <- mean_trend_skip


# test sex effects #
ks.test(sc_males, sc_females, exact=T)

ks.test(se_males, se_females, exact=T)

ks.test(ss_males, ss_females, exact=T)






############################
# nb_200 results for females #
############################

trend_chick <- list()
trend_egg <- list()
trend_skip <- list()
pdf("figure_Sn_nb200_females.pdf", family="Times",width=10,height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
for(i in levels(sum_all$ind)){
	if(sex[which(sum_all$ind == i)[1]] == "1"){
		n_years <- length(which(coe_trait$Ring == i))
		for(year in 1:n_years){
			if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
				if((!is.na(coe_trait$post_rs[which(coe_trait$Ring == i)])) & (length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
					dates <- sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
					startdate <- dates[1]
					y <- sum_all$nb_200[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
					y_ts <- ts(y, frequency = 30, start = startdate)
					if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK"){
						trend_chick[[i]] <- decompose(y_ts,type = "additive")$trend
						plot(decompose(y_ts,type = "additive"))
						title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						#plot(trend_chick[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
					}else{
						if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG"){
							trend_egg[[i]] <- decompose(y_ts,type = "additive")$trend
							plot(decompose(y_ts,type = "additive"))
							#plot(trend_egg[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
							title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						}else{
							trend_skip[[i]] <- decompose(y_ts,type = "additive")$trend
							plot(decompose(y_ts,type = "additive"))
							title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
							#plot(trend_skip[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						}
					}
				}
			}
		}
	}
}
dev.off()



pdf("trends_nb200_females.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(trend_chick[[1]], xlab="Time (arbitrary units)", ylab="Trend")
for(i in 2:length(trend_chick)){
	par(new=TRUE)
	plot(trend_chick[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="blue")
}
for(i in 2:length(trend_egg)){
	par(new=TRUE)
	plot(trend_egg[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="red")
}
for(i in 2:length(trend_skip)){
	par(new=TRUE)
	plot(trend_skip[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="black",lwd=2)
}
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()


trend_chick_mean <- sapply(trend_chick, mean, na.rm=T)
trend_egg_mean <- sapply(trend_egg, mean, na.rm=T)
trend_skip_mean <- sapply(trend_skip, mean, na.rm=T)
pdf("trends_nb200_mean_females.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
	plot(trend_chick_mean,type="l",ylim=c(0,2000),col="blue")
	lines(trend_egg_mean,col="red")
	lines(trend_skip_mean,lwd=2)
dev.off()





pdf("trends_nb200_means_females.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# chick
min_length_chick <- min(sapply(trend_chick, length))
mean_trend_chick <- numeric(min_length_chick)
ss_trend_chick <- numeric(min_length_chick)
sd_trend_chick <- numeric(min_length_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		mean_trend_chick[i] <- mean_trend_chick[i] + trend_chick[[j]][i]
		ss_trend_chick[i] <- ss_trend_chick[i] + (trend_chick[[j]][i])^2
	}
}
mean_trend_chick <- mean_trend_chick / length(trend_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		sd_trend_chick[i] <- sd_trend_chick[i] + (trend_chick[[j]][i] - mean_trend_chick[i])^2
	}
}
sd_trend_chick <- sqrt( sd_trend_chick / length(trend_chick))
# egg
min_length_egg <- min(sapply(trend_egg, length))
mean_trend_egg <- numeric(min_length_egg)
ss_trend_egg <- numeric(min_length_egg)
sd_trend_egg <- numeric(min_length_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		mean_trend_egg[i] <- mean_trend_egg[i] + trend_egg[[j]][i]
		ss_trend_egg[i] <- ss_trend_egg[i] + (trend_egg[[j]][i])^2
	}
}
mean_trend_egg <- mean_trend_egg / length(trend_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		sd_trend_egg[i] <- sd_trend_egg[i] + (trend_egg[[j]][i] - mean_trend_egg[i])^2
	}
}
sd_trend_egg <- sqrt( sd_trend_egg / length(trend_egg))
# skip
min_length_skip <- min(sapply(trend_skip, length))
mean_trend_skip <- numeric(min_length_skip)
ss_trend_skip <- numeric(min_length_skip)
sd_trend_skip <- numeric(min_length_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		mean_trend_skip[i] <- mean_trend_skip[i] + trend_skip[[j]][i]
		ss_trend_skip[i] <- ss_trend_skip[i] + (trend_skip[[j]][i])^2
	}
}
mean_trend_skip <- mean_trend_skip / length(trend_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		sd_trend_skip[i] <- sd_trend_skip[i] + (trend_skip[[j]][i] - mean_trend_skip[i])^2
	}
}
sd_trend_skip <- sqrt( sd_trend_skip / length(trend_skip))
#
plot(mean_trend_chick,col="blue",type="l",xlab="Relative time (days)",ylab="Mean trend",ylim=c(min(c(mean_trend_chick - sd_trend_chick, mean_trend_egg - sd_trend_egg, mean_trend_skip - sd_trend_skip),na.rm=T), max(c(mean_trend_chick + sd_trend_chick, mean_trend_egg + sd_trend_egg, mean_trend_skip + sd_trend_skip),na.rm=T)))
x_seq_chick <- seq(1:min_length_chick)
#segments(x_seq_chick, mean_trend_chick - sd_trend_chick, x_seq_chick, mean_trend_chick + sd_trend_chick,col="blue")
lines(mean_trend_chick - sd_trend_chick, col="blue", lty=2)
lines(mean_trend_chick + sd_trend_chick, col="blue", lty=2)
lines(mean_trend_egg,col="red")
x_seq_egg <- seq(1:min_length_egg)
#segments(x_seq_egg, mean_trend_egg - sd_trend_egg, x_seq_egg, mean_trend_egg + sd_trend_egg,col="red")
lines(mean_trend_egg - sd_trend_egg, col="red", lty=2)
lines(mean_trend_egg + sd_trend_egg, col="red", lty=2)
lines(mean_trend_skip,lwd=2)
x_seq_skip <- seq(1:min_length_skip)
#segments(x_seq_skip, mean_trend_skip - sd_trend_skip, x_seq_skip, mean_trend_skip + sd_trend_skip,col="black")
lines(mean_trend_skip - sd_trend_skip, col="black", lty=2)
lines(mean_trend_skip + sd_trend_skip, col="black", lty=2)
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()



pdf("ecdf_nb200_females.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
ecdf_mean_trend_chick <- ecdf(mean_trend_chick)
ecdf_mean_trend_egg <- ecdf(mean_trend_egg)
ecdf_mean_trend_skip <- ecdf(mean_trend_skip)
plot(ecdf_mean_trend_chick,col="blue",main="ECDFs",xlim=c(100,750))
lines(ecdf_mean_trend_egg,col="red")
lines(ecdf_mean_trend_skip)
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")

# manual computation of ECDFs...
sc <- sort(mean_trend_chick)
se <- sort(mean_trend_egg)
ss <- sort(mean_trend_skip)
nc <- length(sc); ne <- length(se); ns <- length(ss)
plot(sc, (1:nc)/nc, type = 's', ylim = c(0, 1), col="blue")    
lines(se, (1:ne)/ne, type = 's', col="red")    
lines(ss, (1:ns)/ns, type = 's', lwd=2)    
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")
dev.off()

# KS tests [Bonferonni required?]
ks.test(sc, se)
ks.test(sc, ss)

sc_females <- sc
se_females <- se
ss_females <- ss

mean_trend_chick_0 <- mean_trend_chick
mean_trend_egg_0 <- mean_trend_egg
mean_trend_skip_0 <- mean_trend_skip


############################
# nb_200 results for males #
############################

trend_chick <- list()
trend_egg <- list()
trend_skip <- list()
pdf("figure_Sn_nb200_males.pdf", family="Times",width=10,height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
for(i in levels(sum_all$ind)){
	if(sex[which(sum_all$ind == i)[1]] == "2"){
		n_years <- length(which(coe_trait$Ring == i))
		for(year in 1:n_years){
			if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
				if((!is.na(coe_trait$post_rs[which(coe_trait$Ring == i)])) & (length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
					dates <- sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
					startdate <- dates[1]
					y <- sum_all$nb_200[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
					y_ts <- ts(y, frequency = 30, start = startdate)
					if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK"){
						trend_chick[[i]] <- decompose(y_ts,type = "additive")$trend
						plot(decompose(y_ts,type = "additive"))
						title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						#plot(trend_chick[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
					}else{
						if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG"){
							trend_egg[[i]] <- decompose(y_ts,type = "additive")$trend
							plot(decompose(y_ts,type = "additive"))
							#plot(trend_egg[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
							title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						}else{
							trend_skip[[i]] <- decompose(y_ts,type = "additive")$trend
							plot(decompose(y_ts,type = "additive"))
							title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
							#plot(trend_skip[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						}
					}
				}
			}
		}
	}
}
dev.off()



pdf("trends_nb200_males.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(trend_chick[[1]], xlab="Time (arbitrary units)", ylab="Trend")
for(i in 2:length(trend_chick)){
	par(new=TRUE)
	plot(trend_chick[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="blue")
}
for(i in 2:length(trend_egg)){
	par(new=TRUE)
	plot(trend_egg[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="red")
}
for(i in 2:length(trend_skip)){
	par(new=TRUE)
	plot(trend_skip[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="black",lwd=2)
}
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()


trend_chick_mean <- sapply(trend_chick, mean, na.rm=T)
trend_egg_mean <- sapply(trend_egg, mean, na.rm=T)
trend_skip_mean <- sapply(trend_skip, mean, na.rm=T)
pdf("trends_nb200_mean_males.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
	plot(trend_chick_mean,type="l",ylim=c(0,2000),col="blue")
	lines(trend_egg_mean,col="red")
	lines(trend_skip_mean,lwd=2)
dev.off()





pdf("trends_nb200_means_males.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# chick
min_length_chick <- min(sapply(trend_chick, length))
mean_trend_chick <- numeric(min_length_chick)
ss_trend_chick <- numeric(min_length_chick)
sd_trend_chick <- numeric(min_length_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		mean_trend_chick[i] <- mean_trend_chick[i] + trend_chick[[j]][i]
		ss_trend_chick[i] <- ss_trend_chick[i] + (trend_chick[[j]][i])^2
	}
}
mean_trend_chick <- mean_trend_chick / length(trend_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		sd_trend_chick[i] <- sd_trend_chick[i] + (trend_chick[[j]][i] - mean_trend_chick[i])^2
	}
}
sd_trend_chick <- sqrt( sd_trend_chick / length(trend_chick))
# egg
min_length_egg <- min(sapply(trend_egg, length))
mean_trend_egg <- numeric(min_length_egg)
ss_trend_egg <- numeric(min_length_egg)
sd_trend_egg <- numeric(min_length_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		mean_trend_egg[i] <- mean_trend_egg[i] + trend_egg[[j]][i]
		ss_trend_egg[i] <- ss_trend_egg[i] + (trend_egg[[j]][i])^2
	}
}
mean_trend_egg <- mean_trend_egg / length(trend_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		sd_trend_egg[i] <- sd_trend_egg[i] + (trend_egg[[j]][i] - mean_trend_egg[i])^2
	}
}
sd_trend_egg <- sqrt( sd_trend_egg / length(trend_egg))
# skip
min_length_skip <- min(sapply(trend_skip, length))
mean_trend_skip <- numeric(min_length_skip)
ss_trend_skip <- numeric(min_length_skip)
sd_trend_skip <- numeric(min_length_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		mean_trend_skip[i] <- mean_trend_skip[i] + trend_skip[[j]][i]
		ss_trend_skip[i] <- ss_trend_skip[i] + (trend_skip[[j]][i])^2
	}
}
mean_trend_skip <- mean_trend_skip / length(trend_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		sd_trend_skip[i] <- sd_trend_skip[i] + (trend_skip[[j]][i] - mean_trend_skip[i])^2
	}
}
sd_trend_skip <- sqrt( sd_trend_skip / length(trend_skip))
#
plot(mean_trend_chick,col="blue",type="l",xlab="Relative time (days)",ylab="Mean trend",ylim=c(min(c(mean_trend_chick - sd_trend_chick, mean_trend_egg - sd_trend_egg, mean_trend_skip - sd_trend_skip),na.rm=T), max(c(mean_trend_chick + sd_trend_chick, mean_trend_egg + sd_trend_egg, mean_trend_skip + sd_trend_skip),na.rm=T)))
x_seq_chick <- seq(1:min_length_chick)
#segments(x_seq_chick, mean_trend_chick - sd_trend_chick, x_seq_chick, mean_trend_chick + sd_trend_chick,col="blue")
lines(mean_trend_chick - sd_trend_chick, col="blue", lty=2)
lines(mean_trend_chick + sd_trend_chick, col="blue", lty=2)
lines(mean_trend_egg,col="red")
x_seq_egg <- seq(1:min_length_egg)
#segments(x_seq_egg, mean_trend_egg - sd_trend_egg, x_seq_egg, mean_trend_egg + sd_trend_egg,col="red")
lines(mean_trend_egg - sd_trend_egg, col="red", lty=2)
lines(mean_trend_egg + sd_trend_egg, col="red", lty=2)
lines(mean_trend_skip,lwd=2)
x_seq_skip <- seq(1:min_length_skip)
#segments(x_seq_skip, mean_trend_skip - sd_trend_skip, x_seq_skip, mean_trend_skip + sd_trend_skip,col="black")
lines(mean_trend_skip - sd_trend_skip, col="black", lty=2)
lines(mean_trend_skip + sd_trend_skip, col="black", lty=2)
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()



pdf("ecdf_nb200_males.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
ecdf_mean_trend_chick <- ecdf(mean_trend_chick)
ecdf_mean_trend_egg <- ecdf(mean_trend_egg)
ecdf_mean_trend_skip <- ecdf(mean_trend_skip)
plot(ecdf_mean_trend_chick,col="blue",main="ECDFs",xlim=c(100,750))
lines(ecdf_mean_trend_egg,col="red")
lines(ecdf_mean_trend_skip)
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")

# manual computation of ECDFs...
sc <- sort(mean_trend_chick)
se <- sort(mean_trend_egg)
ss <- sort(mean_trend_skip)
nc <- length(sc); ne <- length(se); ns <- length(ss)
plot(sc, (1:nc)/nc, type = 's', ylim = c(0, 1), col="blue")    
lines(se, (1:ne)/ne, type = 's', col="red")    
lines(ss, (1:ns)/ns, type = 's', lwd=2)    
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")
dev.off()

# KS tests [Bonferonni required?]
ks.test(sc, se)
ks.test(sc, ss)

sc_males <- sc
se_males <- se
ss_males <- ss

mean_trend_chick_0 <- mean_trend_chick
mean_trend_egg_0 <- mean_trend_egg
mean_trend_skip_0 <- mean_trend_skip


# test sex effects #
ks.test(sc_males, sc_females, exact=T)
#	Two-sample Kolmogorov-Smirnov test
#
#data:  sc_males and sc_females
#D = 0.08768, p-value = 0.5866
#alternative hypothesis: two-sided

ks.test(se_males, se_females, exact=T)
#	Two-sample Kolmogorov-Smirnov test
#
#data:  se_males and se_females
#D = 0.56, p-value < 2.2e-16
#alternative hypothesis: two-sided

ks.test(ss_males, ss_females, exact=T)
#	Two-sample Kolmogorov-Smirnov test
#
#data:  ss_males and ss_females
#D = 0.12688, p-value = 0.1427
#alternative hypothesis: two-sided


############################
# nb_1to199 results for females #
############################

trend_chick <- list()
trend_egg <- list()
trend_skip <- list()
pdf("figure_Sn_nb1to199_females.pdf", family="Times",width=10,height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
for(i in levels(sum_all$ind)){
	if(sex[which(sum_all$ind == i)[1]] == "1"){
		n_years <- length(which(coe_trait$Ring == i))
		for(year in 1:n_years){
			if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
				if((!is.na(coe_trait$post_rs[which(coe_trait$Ring == i)])) & (length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
					dates <- sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
					startdate <- dates[1]
					y <- sum_all$nb_1to199[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
					y_ts <- ts(y, frequency = 30, start = startdate)
					if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK"){
						trend_chick[[i]] <- decompose(y_ts,type = "additive")$trend
						plot(decompose(y_ts,type = "additive"))
						title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						#plot(trend_chick[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
					}else{
						if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG"){
							trend_egg[[i]] <- decompose(y_ts,type = "additive")$trend
							plot(decompose(y_ts,type = "additive"))
							#plot(trend_egg[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
							title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						}else{
							trend_skip[[i]] <- decompose(y_ts,type = "additive")$trend
							plot(decompose(y_ts,type = "additive"))
							title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
							#plot(trend_skip[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						}
					}
				}
			}
		}
	}
}
dev.off()



pdf("trends_nb1to199_females.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(trend_chick[[1]], xlab="Time (arbitrary units)", ylab="Trend")
for(i in 2:length(trend_chick)){
	par(new=TRUE)
	plot(trend_chick[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="blue")
}
for(i in 2:length(trend_egg)){
	par(new=TRUE)
	plot(trend_egg[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="red")
}
for(i in 2:length(trend_skip)){
	par(new=TRUE)
	plot(trend_skip[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="black",lwd=2)
}
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()


trend_chick_mean <- sapply(trend_chick, mean, na.rm=T)
trend_egg_mean <- sapply(trend_egg, mean, na.rm=T)
trend_skip_mean <- sapply(trend_skip, mean, na.rm=T)
pdf("trends_nb1to199_mean_females.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
	plot(trend_chick_mean,type="l",ylim=c(0,2000),col="blue")
	lines(trend_egg_mean,col="red")
	lines(trend_skip_mean,lwd=2)
dev.off()





pdf("trends_nb1to199_means_females.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# chick
min_length_chick <- min(sapply(trend_chick, length))
mean_trend_chick <- numeric(min_length_chick)
ss_trend_chick <- numeric(min_length_chick)
sd_trend_chick <- numeric(min_length_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		mean_trend_chick[i] <- mean_trend_chick[i] + trend_chick[[j]][i]
		ss_trend_chick[i] <- ss_trend_chick[i] + (trend_chick[[j]][i])^2
	}
}
mean_trend_chick <- mean_trend_chick / length(trend_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		sd_trend_chick[i] <- sd_trend_chick[i] + (trend_chick[[j]][i] - mean_trend_chick[i])^2
	}
}
sd_trend_chick <- sqrt( sd_trend_chick / length(trend_chick))
# egg
min_length_egg <- min(sapply(trend_egg, length))
mean_trend_egg <- numeric(min_length_egg)
ss_trend_egg <- numeric(min_length_egg)
sd_trend_egg <- numeric(min_length_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		mean_trend_egg[i] <- mean_trend_egg[i] + trend_egg[[j]][i]
		ss_trend_egg[i] <- ss_trend_egg[i] + (trend_egg[[j]][i])^2
	}
}
mean_trend_egg <- mean_trend_egg / length(trend_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		sd_trend_egg[i] <- sd_trend_egg[i] + (trend_egg[[j]][i] - mean_trend_egg[i])^2
	}
}
sd_trend_egg <- sqrt( sd_trend_egg / length(trend_egg))
# skip
min_length_skip <- min(sapply(trend_skip, length))
mean_trend_skip <- numeric(min_length_skip)
ss_trend_skip <- numeric(min_length_skip)
sd_trend_skip <- numeric(min_length_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		mean_trend_skip[i] <- mean_trend_skip[i] + trend_skip[[j]][i]
		ss_trend_skip[i] <- ss_trend_skip[i] + (trend_skip[[j]][i])^2
	}
}
mean_trend_skip <- mean_trend_skip / length(trend_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		sd_trend_skip[i] <- sd_trend_skip[i] + (trend_skip[[j]][i] - mean_trend_skip[i])^2
	}
}
sd_trend_skip <- sqrt( sd_trend_skip / length(trend_skip))
#
plot(mean_trend_chick,col="blue",type="l",xlab="Relative time (days)",ylab="Mean trend",ylim=c(min(c(mean_trend_chick - sd_trend_chick, mean_trend_egg - sd_trend_egg, mean_trend_skip - sd_trend_skip),na.rm=T), max(c(mean_trend_chick + sd_trend_chick, mean_trend_egg + sd_trend_egg, mean_trend_skip + sd_trend_skip),na.rm=T)))
x_seq_chick <- seq(1:min_length_chick)
#segments(x_seq_chick, mean_trend_chick - sd_trend_chick, x_seq_chick, mean_trend_chick + sd_trend_chick,col="blue")
lines(mean_trend_chick - sd_trend_chick, col="blue", lty=2)
lines(mean_trend_chick + sd_trend_chick, col="blue", lty=2)
lines(mean_trend_egg,col="red")
x_seq_egg <- seq(1:min_length_egg)
#segments(x_seq_egg, mean_trend_egg - sd_trend_egg, x_seq_egg, mean_trend_egg + sd_trend_egg,col="red")
lines(mean_trend_egg - sd_trend_egg, col="red", lty=2)
lines(mean_trend_egg + sd_trend_egg, col="red", lty=2)
lines(mean_trend_skip,lwd=2)
x_seq_skip <- seq(1:min_length_skip)
#segments(x_seq_skip, mean_trend_skip - sd_trend_skip, x_seq_skip, mean_trend_skip + sd_trend_skip,col="black")
lines(mean_trend_skip - sd_trend_skip, col="black", lty=2)
lines(mean_trend_skip + sd_trend_skip, col="black", lty=2)
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()



pdf("ecdf_nb1to199_females.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
ecdf_mean_trend_chick <- ecdf(mean_trend_chick)
ecdf_mean_trend_egg <- ecdf(mean_trend_egg)
ecdf_mean_trend_skip <- ecdf(mean_trend_skip)
plot(ecdf_mean_trend_chick,col="blue",main="ECDFs",xlim=c(650,875))
lines(ecdf_mean_trend_egg,col="red")
lines(ecdf_mean_trend_skip)
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")

# manual computation of ECDFs...
sc <- sort(mean_trend_chick)
se <- sort(mean_trend_egg)
ss <- sort(mean_trend_skip)
nc <- length(sc); ne <- length(se); ns <- length(ss)
plot(sc, (1:nc)/nc, type = 's', ylim = c(0, 1), col="blue")    
lines(se, (1:ne)/ne, type = 's', col="red")    
lines(ss, (1:ns)/ns, type = 's', lwd=2)    
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")
dev.off()

# KS tests [Bonferonni required?]
ks.test(sc, se)
ks.test(sc, ss)

sc_females <- sc
se_females <- se
ss_females <- ss

mean_trend_chick_0 <- mean_trend_chick
mean_trend_egg_0 <- mean_trend_egg
mean_trend_skip_0 <- mean_trend_skip


############################
# nb_1to199 results for males #
############################

trend_chick <- list()
trend_egg <- list()
trend_skip <- list()
pdf("figure_Sn_nb1to199_males.pdf", family="Times",width=10,height=12)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
for(i in levels(sum_all$ind)){
	if(sex[which(sum_all$ind == i)[1]] == "2"){
		n_years <- length(which(coe_trait$Ring == i))
		for(year in 1:n_years){
			if(!is.na(coe_trait$Colony_departure[which(coe_trait$Ring == i)][year])){
				if((!is.na(coe_trait$post_rs[which(coe_trait$Ring == i)])) & (length(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]) > 0)){
					dates <- sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
					startdate <- dates[1]
					y <- sum_all$nb_1to199[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])]
					y_ts <- ts(y, frequency = 30, start = startdate)
					if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "CHICK"){
						trend_chick[[i]] <- decompose(y_ts,type = "additive")$trend
						plot(decompose(y_ts,type = "additive"))
						title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						#plot(trend_chick[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
					}else{
						if(coe_trait$post_rs[which(coe_trait$Ring == i)][year] == "EGG"){
							trend_egg[[i]] <- decompose(y_ts,type = "additive")$trend
							plot(decompose(y_ts,type = "additive"))
							#plot(trend_egg[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
							title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						}else{
							trend_skip[[i]] <- decompose(y_ts,type = "additive")$trend
							plot(decompose(y_ts,type = "additive"))
							title(paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
							#plot(trend_skip[[i]], col="red", main=paste("Bird ID: ",i," -- ",coe_trait$post_rs[which(coe_trait$Ring == i)][year],sep=""))
						}
					}
				}
			}
		}
	}
}
dev.off()



pdf("trends_nb1to199_males.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
plot(trend_chick[[1]], xlab="Time (arbitrary units)", ylab="Trend")
for(i in 2:length(trend_chick)){
	par(new=TRUE)
	plot(trend_chick[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="blue")
}
for(i in 2:length(trend_egg)){
	par(new=TRUE)
	plot(trend_egg[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="red")
}
for(i in 2:length(trend_skip)){
	par(new=TRUE)
	plot(trend_skip[[i]],xaxt="n",yaxt="n",xlab="",ylab="",col="black",lwd=2)
}
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()


trend_chick_mean <- sapply(trend_chick, mean, na.rm=T)
trend_egg_mean <- sapply(trend_egg, mean, na.rm=T)
trend_skip_mean <- sapply(trend_skip, mean, na.rm=T)
pdf("trends_nb1to199_mean_males.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
	plot(trend_chick_mean,type="l",ylim=c(0,2000),col="blue")
	lines(trend_egg_mean,col="red")
	lines(trend_skip_mean,lwd=2)
dev.off()





pdf("trends_nb1to199_means_males.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
# chick
min_length_chick <- min(sapply(trend_chick, length))
mean_trend_chick <- numeric(min_length_chick)
ss_trend_chick <- numeric(min_length_chick)
sd_trend_chick <- numeric(min_length_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		mean_trend_chick[i] <- mean_trend_chick[i] + trend_chick[[j]][i]
		ss_trend_chick[i] <- ss_trend_chick[i] + (trend_chick[[j]][i])^2
	}
}
mean_trend_chick <- mean_trend_chick / length(trend_chick)
for(i in 1:min_length_chick){
	for(j in 1:length(trend_chick)){
		sd_trend_chick[i] <- sd_trend_chick[i] + (trend_chick[[j]][i] - mean_trend_chick[i])^2
	}
}
sd_trend_chick <- sqrt( sd_trend_chick / length(trend_chick))
# egg
min_length_egg <- min(sapply(trend_egg, length))
mean_trend_egg <- numeric(min_length_egg)
ss_trend_egg <- numeric(min_length_egg)
sd_trend_egg <- numeric(min_length_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		mean_trend_egg[i] <- mean_trend_egg[i] + trend_egg[[j]][i]
		ss_trend_egg[i] <- ss_trend_egg[i] + (trend_egg[[j]][i])^2
	}
}
mean_trend_egg <- mean_trend_egg / length(trend_egg)
for(i in 1:min_length_egg){
	for(j in 1:length(trend_egg)){
		sd_trend_egg[i] <- sd_trend_egg[i] + (trend_egg[[j]][i] - mean_trend_egg[i])^2
	}
}
sd_trend_egg <- sqrt( sd_trend_egg / length(trend_egg))
# skip
min_length_skip <- min(sapply(trend_skip, length))
mean_trend_skip <- numeric(min_length_skip)
ss_trend_skip <- numeric(min_length_skip)
sd_trend_skip <- numeric(min_length_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		mean_trend_skip[i] <- mean_trend_skip[i] + trend_skip[[j]][i]
		ss_trend_skip[i] <- ss_trend_skip[i] + (trend_skip[[j]][i])^2
	}
}
mean_trend_skip <- mean_trend_skip / length(trend_skip)
for(i in 1:min_length_skip){
	for(j in 1:length(trend_skip)){
		sd_trend_skip[i] <- sd_trend_skip[i] + (trend_skip[[j]][i] - mean_trend_skip[i])^2
	}
}
sd_trend_skip <- sqrt( sd_trend_skip / length(trend_skip))
#
plot(mean_trend_chick,col="blue",type="l",xlab="Relative time (days)",ylab="Mean trend",ylim=c(min(c(mean_trend_chick - sd_trend_chick, mean_trend_egg - sd_trend_egg, mean_trend_skip - sd_trend_skip),na.rm=T), max(c(mean_trend_chick + sd_trend_chick, mean_trend_egg + sd_trend_egg, mean_trend_skip + sd_trend_skip),na.rm=T)))
x_seq_chick <- seq(1:min_length_chick)
#segments(x_seq_chick, mean_trend_chick - sd_trend_chick, x_seq_chick, mean_trend_chick + sd_trend_chick,col="blue")
lines(mean_trend_chick - sd_trend_chick, col="blue", lty=2)
lines(mean_trend_chick + sd_trend_chick, col="blue", lty=2)
lines(mean_trend_egg,col="red")
x_seq_egg <- seq(1:min_length_egg)
#segments(x_seq_egg, mean_trend_egg - sd_trend_egg, x_seq_egg, mean_trend_egg + sd_trend_egg,col="red")
lines(mean_trend_egg - sd_trend_egg, col="red", lty=2)
lines(mean_trend_egg + sd_trend_egg, col="red", lty=2)
lines(mean_trend_skip,lwd=2)
x_seq_skip <- seq(1:min_length_skip)
#segments(x_seq_skip, mean_trend_skip - sd_trend_skip, x_seq_skip, mean_trend_skip + sd_trend_skip,col="black")
lines(mean_trend_skip - sd_trend_skip, col="black", lty=2)
lines(mean_trend_skip + sd_trend_skip, col="black", lty=2)
legend("bottomleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(1,1,2),bty="n")
dev.off()



pdf("ecdf_nb1to199_males.pdf", family="Times",width=10,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
ecdf_mean_trend_chick <- ecdf(mean_trend_chick)
ecdf_mean_trend_egg <- ecdf(mean_trend_egg)
ecdf_mean_trend_skip <- ecdf(mean_trend_skip)
plot(ecdf_mean_trend_chick,col="blue",main="ECDFs",xlim=c(650,875))
lines(ecdf_mean_trend_egg,col="red")
lines(ecdf_mean_trend_skip)
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")

# manual computation of ECDFs...
sc <- sort(mean_trend_chick)
se <- sort(mean_trend_egg)
ss <- sort(mean_trend_skip)
nc <- length(sc); ne <- length(se); ns <- length(ss)
plot(sc, (1:nc)/nc, type = 's', ylim = c(0, 1), col="blue")    
lines(se, (1:ne)/ne, type = 's', col="red")    
lines(ss, (1:ns)/ns, type = 's', lwd=2)    
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")
dev.off()

# KS tests [Bonferonni required?]
ks.test(sc, se)
ks.test(sc, ss)

sc_males <- sc
se_males <- se
ss_males <- ss

mean_trend_chick_0 <- mean_trend_chick
mean_trend_egg_0 <- mean_trend_egg
mean_trend_skip_0 <- mean_trend_skip


# test sex effects #
ks.test(sc_males, sc_females, exact=T)
#	Two-sample Kolmogorov-Smirnov test
#
#data:  sc_males and sc_females
#D = 0.34468, p-value = 1.799e-08
#alternative hypothesis: two-sided

ks.test(se_males, se_females, exact=T)
#	Two-sample Kolmogorov-Smirnov test
#
#data:  se_males and se_females
#D = 0.63333, p-value < 2.2e-16
#alternative hypothesis: two-sided

ks.test(ss_males, ss_females, exact=T)
#	Two-sample Kolmogorov-Smirnov test
#
#data:  ss_males and ss_females
#D = 0.37011, p-value = 3.526e-10
#alternative hypothesis: two-sided


####################
# alternative idea #
####################

sum_all$nb_0
hist(sum_all$nb_0)
hist(sum_all$nb_0[sex=="1"])

t.test(sum_all$nb_0[sex=="1"],sum_all$nb_0[sex=="2"])
t.test(sum_all$nb_200[sex=="1"],sum_all$nb_200[sex=="2"])
t.test(sum_all$nb_1to199[sex=="1"],sum_all$nb_1to199[sex=="2"])

pdf("figure_Sn_Sex.pdf", family="Times", width=4,height=7)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(3,1))
plot(locfit( ~ lp(sum_all$nb_0[sex=="1"] )),col="red",xlab="Daily salt data (flying)")
lines(locfit( ~ lp(sum_all$nb_0[sex=="2"] )),col="blue")
legend("topright",c("Females","Males"),lwd=1,col=c("red","blue"),bty = "n")
x <- t.test(sum_all$nb_0[sex=="1"],sum_all$nb_0[sex=="2"])
text(0,.0001,paste0("t = ",x$statistic,", df = ",x$parameter,", P = ",x$p.value),adj = c(0,0),cex=.8)
plot(locfit( ~ lp(sum_all$nb_200[sex=="1"] )),col="red",xlab="Daily salt data (resting)")
lines(locfit( ~ lp(sum_all$nb_200[sex=="2"] )),col="blue")
x <- t.test(sum_all$nb_200[sex=="1"],sum_all$nb_200[sex=="2"])
text(0,.0001,paste0("t = ",x$statistic,", df = ",x$parameter,", P = ",x$p.value),adj = c(0,0),cex=.8)
plot(locfit( ~ lp(sum_all$nb_1to199[sex=="1"] )),col="red",xlab="Daily salt data (foraging)")
lines(locfit( ~ lp(sum_all$nb_1to199[sex=="2"] )),col="blue")
x <- t.test(sum_all$nb_1to199[sex=="1"],sum_all$nb_1to199[sex=="2"])
text(0,.0001,paste0("t = ",x$statistic,", df = ",x$parameter,", P = ",x$p.value),adj = c(0,0),cex=.8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1))
dev.off()

###########################
# figure 2 for MS


pdf("figure_2.pdf", family="Times", width=3,height=9)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(3,1))
# flight
sc <- sort(mean_trend_chick_0)
se <- sort(mean_trend_egg_0)
ss <- sort(mean_trend_skip_0)
nc <- length(sc); ne <- length(se); ns <- length(ss)
plot(sc, (1:nc)/nc, type = 's', ylim = c(0, 1), col="blue",xlim=c(min(sc,se,ss),max(sc,se,ss)),ylab="Empirical cumulative distribution: flight",xlab="Activity (arbitrary unit)")    
lines(se, (1:ne)/ne, type = 's', col="red")    
lines(ss, (1:ns)/ns, type = 's', lwd=2)    
legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")
# rest
sc <- sort(mean_trend_chick_200)
se <- sort(mean_trend_egg_200)
ss <- sort(mean_trend_skip_200)
nc <- length(sc); ne <- length(se); ns <- length(ss)
plot(sc, (1:nc)/nc, type = 's', ylim = c(0, 1), col="blue",xlim=c(min(sc,se,ss),max(sc,se,ss)),ylab="Empirical cumulative distribution: rest",xlab="Activity (arbitrary unit)")    
lines(se, (1:ne)/ne, type = 's', col="red")    
lines(ss, (1:ns)/ns, type = 's', lwd=2)    
#legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")
# other
sc <- sort(mean_trend_chick_1_199)
se <- sort(mean_trend_egg_1_199)
ss <- sort(mean_trend_skip_1_199)
nc <- length(sc); ne <- length(se); ns <- length(ss)
plot(sc, (1:nc)/nc, type = 's', ylim = c(0, 1), col="blue",xlim=c(min(sc,se,ss),max(sc,se,ss)),ylab="Empirical cumulative distribution: foraging",xlab="Activity (arbitrary unit)")    
lines(se, (1:ne)/ne, type = 's', col="red")    
lines(ss, (1:ns)/ns, type = 's', lwd=2)    
#legend("topleft",c("Chick","Egg","Skip"),col=c("blue","red","black"),lwd=c(2,2,2),bty="n")

dev.off()


###########################
# carryover effect tested #

nreps <- 100000
test_stat <- c()
chi2_stat <- c()
for(reps in 1: nreps){
	row1 <- t(rmultinom(1, size = rowsums[1], prob = c(1/3,1/3,1/3)))
	row2 <- t(rmultinom(1, size = rowsums[2], prob = c(1/3,1/3,1/3)))
	row3 <- t(rmultinom(1, size = rowsums[3], prob = c(1/3,1/3,1/3)))
	rand_mat <- matrix (c(row1, row2, row3), ncol = 3, byrow = T)
	# test stat
	test_stat[reps] <- sum(((mymat - rand_mat)^2)/(rand_mat))
	chi2_stat[reps] <- chisq.test(rand_mat)$statistic
}
hist(chi2_stat,100)
abline(v=chisq.test(mymat)$statistic)
hist(test_stat,100)


# NB. you actually have a transition matrix, not a contingency table...

mymat <- matrix (c(197, 82, 32, 49, 22, 8, 18, 3, 1), ncol = 3, byrow = TRUE)
tot_birds <- sum(mymat)
rowsums <- rowSums(mymat)



my_transition_mat <- mymat / rowsums
prob_skip2chick <- my_transition_mat[3,1]
nreps <- 100000
rand_prob_skip2chick <- c()
rand_prob_2chick <- c() # probs of success in following year
for(reps in 1: nreps){
	n1 <- runif(1,0,1)
	n2 <- runif(1,0,1)
	n3 <- runif(1,0,1)
	sum_n <- sum(c(n1,n2,n3))
	row1 <- c(n1/sum_n, n2/sum_n, n3/sum_n)
	n1 <- runif(1,0,1)
	n2 <- runif(1,0,1)
	n3 <- runif(1,0,1)
	sum_n <- sum(c(n1,n2,n3))
	row2 <- c(n1/sum_n, n2/sum_n, n3/sum_n)
	n1 <- runif(1,0,1)
	n2 <- runif(1,0,1)
	n3 <- runif(1,0,1)
	sum_n <- sum(c(n1,n2,n3))
	row3 <- c(n1/sum_n, n2/sum_n, n3/sum_n)
	rand_mat <- matrix (c(row1, row2, row3), ncol = 3, byrow = T)
	# test 1: on S2C transition
	if(rand_mat[3,1] > prob_skip2chick){
		rand_prob_skip2chick[reps] <- 1
	}else{
		rand_prob_skip2chick[reps] <- 0
	}
	# test 2: on probs of success
	if( (rand_mat[1,1] > my_transition_mat[1,1]) & (rand_mat[2,1] > my_transition_mat[2,1]) & (rand_mat[3,1] > my_transition_mat[3,1]) ){
		rand_prob_2chick[reps] <- 1
	}else{
		rand_prob_2chick[reps] <- 0
	}
}
sum(rand_prob_skip2chick) / nreps
sum(rand_prob_2chick) / nreps






library(DTMCPack)
library(locfit)
my_lambda <- eigen(my_transition_mat)
my_det <- det(my_transition_mat)
pi_dist <- statdistr(my_transition_mat)
nreps <- 100000
test_stat <- c(NA,NA,NA)
rand_lambda <- c(NA,NA,NA)
rand_det <- c()
for(reps in 1: nreps){
	n1 <- runif(1,0,1)
	n2 <- runif(1,0,1)
	n3 <- runif(1,0,1)
	sum_n <- sum(c(n1,n2,n3))
	row1 <- c(n1/sum_n, n2/sum_n, n3/sum_n)
	n1 <- runif(1,0,1)
	n2 <- runif(1,0,1)
	n3 <- runif(1,0,1)
	sum_n <- sum(c(n1,n2,n3))
	row2 <- c(n1/sum_n, n2/sum_n, n3/sum_n)
	n1 <- runif(1,0,1)
	n2 <- runif(1,0,1)
	n3 <- runif(1,0,1)
	sum_n <- sum(c(n1,n2,n3))
	row3 <- c(n1/sum_n, n2/sum_n, n3/sum_n)
	rand_mat <- matrix (c(row1, row2, row3), ncol = 3, byrow = T)
	rand_mat_pi_dist <- statdistr(rand_mat)
	test_stat <- rbind(test_stat, rand_mat_pi_dist)
	rand_det[reps] <- det(rand_mat)
	rand_lambda <- rbind(rand_lambda, Re(eigen(rand_mat)$values))
}

test_stat <- test_stat[-1,] # removes the first line, of NAs
rand_lambda <- rand_lambda[-1,] # removes the first line, of NAs


pdf("test_stationary_distr.pdf", family="Times",width=10,height=2)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,3))
plot(locfit( ~ lp(test_stat[,1], nn=.4)),col="blue",xlab="pi Chick")
abline(v= pi_dist[1],lty=2,col="blue")
legend("topright",paste("P = ",format(sum(test_stat[,1] > pi_dist[1])/nreps,digits=3),sep=""),bty="n")
plot(locfit( ~ lp(test_stat[,2], nn=.4)),col="red",xlab="pi Egg")
abline(v= pi_dist[2],lty=2,col="red")
legend("topright",paste("P = ",format(sum(test_stat[,2] < pi_dist[2])/nreps,digits=3),sep=""),bty="n")
plot(locfit( ~ lp(test_stat[,3], nn=.4)),col="black",xlab="pi Skip")
abline(v= pi_dist[3],lty=2,col="black")
legend("topright",paste("P = ",format(sum(test_stat[,3] < pi_dist[3])/nreps,digits=3),sep=""),bty="n")
dev.off()



library(markovchain)
statesNames <- c("Chick","Egg","Skip")
dimnames(my_transition_mat) <- list(statesNames,statesNames)
mcA <- new("markovchain", transitionMatrix = my_transition_mat)
absorbingStates(mcA)
transientStates(mcA)
is.irreducible(mcA)
period(mcA)
steadyStates(mcA) # same as pi_dist, of course


FP_Chick <- firstPassage  (mcA,"Chick",10) # longevity of 50 years?
FP_Egg <- firstPassage  (mcA,"Egg",10) # longevity of 50 years?
FP_Skip <- firstPassage  (mcA,"Skip",10) # longevity of 50 years?


pdf("FistPassageTimes.pdf", family="Times",width=10,height=2)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,3))
plot(FP_Chick[,1],type="l",col="blue",xlab="Time (in years)")
lines(FP_Chick[,2],type="l",col="red")
lines(FP_Chick[,3],type="l")
plot(FP_Egg[,1],type="l",col="blue",xlab="Time (in years)")
lines(FP_Egg[,2],type="l",col="red")
lines(FP_Egg[,3],type="l")
plot(FP_Skip[,1],type="l",col="blue",xlab="Time (in years)")
lines(FP_Skip[,2],type="l",col="red")
lines(FP_Skip[,3],type="l")
dev.off()


##########################################################################



save.image("coe_15.RData")
q(save = "no")



##########################################################################
# below is not used 
##########################################################################



# i <- "fb32249"
# 
# 
# 
# plot(y_ts, main=paste("Bird ID: ",i," -- ",coe_trait$prior_rs[which(coe_trait$Ring == i)][year],sep=""))
# 
# periodicity(y_ts)
# plot(decompose(y_ts,type = "additive"))
# axis(1, labels=time, at=seq(from= startdate, by=10, length.out=length(y)) )
# 
# par(new=TRUE)
# plot(sum_all$sum_all_dates[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])], sum_all$nb_0[(sum_all$ind == i) & (sum_all$sum_all_dates > coe_trait$Colony_departure[which(coe_trait$Ring == i)][year]) & (sum_all$sum_all_dates < coe_trait$Colony_arrival[which(coe_trait$Ring == i)][year])], main=paste("Bird ID: ",i," -- ",coe_trait$prior_rs[which(coe_trait$Ring == i)][year],sep=""), type="l", xlab="Date", ylab="Total time on a given date",xaxt="n",col="red")
# axis.POSIXct(1, at= ceiling_date(sum_all$sum_all_dates[sum_all$ind == i], "month"), labels=format(sum_all$sum_all_dates[sum_all$ind == i], "%m/%d/%Y"),cex=.5)
# 
# 
# ssp <- spectrum(y_ts)  
# per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
# reslm <- lm(y ~ sin(2*pi/per*t)+cos(2*pi/per*t))
# summary(reslm)
# plot(y ~ t,type="l")
# lines(fitted(reslm)~t,col="red",lty=2,lwd=2)
# 
# 
# 
# 
# ##########################################################################
# 
# 
# i = levels(sum_all$ind)[1]
# 
# # all data
# plot(sum_all$sum_all_dates[sum_all$ind == i], sum_all$nb_0[sum_all$ind == i], main=paste("Bird ID: ",i," -- ",coe_trait$prior_rs[which(coe_trait$Ring == i & coe_trait$GLS_Year_start == format(min(sum_all$sum_all_dates[sum_all$ind == i]), "%Y"))],sep=""), type="l", xlab="Date", ylab="Total time on a given date",xaxt="n",col="red")
# axis.POSIXct(1, at= ceiling_date(sum_all$sum_all_dates[sum_all$ind == i], "month"), labels=format(sum_all$sum_all_dates[sum_all$ind == i], "%m/%d/%Y"),cex=.5)
# 
# # data between date1 and date2 (picked randomly, but should come form your trait file)
# date1 <- sum_all$sum_all_dates[sum_all$ind == i][100]
# date2 <- sum_all$sum_all_dates[sum_all$ind == i][300]
# 
# pdf()
# plot(sum_all$sum_all_dates[sum_all$ind == i & sum_all$sum_all_dates > date1 & sum_all$sum_all_dates < date2], sum_all$nb_0[sum_all$ind == i & sum_all$sum_all_dates > date1 & sum_all$sum_all_dates < date2], main=paste("Bird ID: ",i," -- ",coe_trait$prior_rs[which(coe_trait$Ring == i & coe_trait$GLS_Year_start == format(min(sum_all$sum_all_dates[sum_all$ind == i]), "%Y"))],sep=""), type="l", xlab="Date", ylab="Total time on a given date",xaxt="n",col="red")
# axis.POSIXct(1, at= ceiling_date(sum_all$sum_all_dates[sum_all$ind == i], "month"), labels=format(sum_all$sum_all_dates[sum_all$ind == i], "%m/%d/%Y"),cex=.5)
# dev.off()
# 
# #
# y <- sum_all$nb_0[sum_all$ind == i & sum_all$sum_all_dates > date1 & sum_all$sum_all_dates < date2]
# t <- 1:length(sum_all$sum_all_dates[sum_all$ind == i & sum_all$sum_all_dates > date1 & sum_all$sum_all_dates < date2])
# ssp <- spectrum(y)  
# per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
# reslm <- lm(y ~ sin(2*pi/per*t)+cos(2*pi/per*t))
# summary(reslm)
# plot(y ~ t,type="l")
# lines(fitted(reslm)~t,col="red",lty=2,lwd=2)
# 
# #
# pdf()
# y_ts <- ts(y, , frequency = 12, start=date1)
# plot(decompose(y_ts,type = "additive"))
# dev.off()
# 
# 
# y_ts <- ts(y, , frequency = 24, start=date1)
# plot(decompose(y_ts,type = "additive"))
# plot(decompose(y_ts,type = "additive")$trend)


