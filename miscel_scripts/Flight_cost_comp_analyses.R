# Second script used in Guigueno M.F., Shoji A., Elliott K.H. and Aris-Brosou S. 2019. 
# Flight costs in volant vertebrates: a phylogenetically-controlled meta-analysis of birds and bats. 
# Comparative Biochemistry and Physiology - Part A: Molecular & Integrative Physiology. 235:193-201.
# https://doi.org/10.1016/j.cbpa.2019.06.003

library(ape)
library(robust)
library(caper)
library(ggplot2)
require(scales)

# load("180509_flight_cost_birds.RData")
# load("180509_flight_cost_bats.RData")

run_phyml <- 0

###############################
# estimate phylogenetic trees #
###############################

if(run_phyml){
	curdir <- getwd()
	setwd(paste0(curdir, "/split_analyses/2.birds"))
	system("phyml -i cox1_cytb_all.aln -b -4 --run_id GTR+I+G -m 012345 -f m -v e -c 4 -a e -s BEST --r_seed 1203 &")
	system("phyml -i cox1_cytb_all.aln -b 250 --run_id GTR+I+G_boot1000a -m 012345 -f m -v e -c 4 -a e -s BEST --r_seed 1235 &")
	system("phyml -i cox1_cytb_all.aln -b 250 --run_id GTR+I+G_boot1000b -m 012345 -f m -v e -c 4 -a e -s BEST --r_seed 1237 &")
	system("phyml -i cox1_cytb_all.aln -b 250 --run_id GTR+I+G_boot1000c -m 012345 -f m -v e -c 4 -a e -s BEST --r_seed 1239 &")
	system("phyml -i cox1_cytb_all.aln -b 250 --run_id GTR+I+G_boot1000d -m 012345 -f m -v e -c 4 -a e -s BEST --r_seed 1231 &")
	system("cat cox1_cytb_all.aln_phyml_boot_trees_GTR+I+G_boot1000a.txt cox1_cytb_all.aln_phyml_boot_trees_GTR+I+G_boot1000b.txt cox1_cytb_all.aln_phyml_boot_trees_GTR+I+G_boot1000c.txt cox1_cytb_all.aln_phyml_boot_trees_GTR+I+G_boot1000d.txt > boot.trees")


	
	setwd(paste0(curdir, "../../180122_flightCosts/split_analyses/1.bats"))
	system("phyml -i cox1_cytb_dups.phy -b -4 --run_id GTR+I+G -m 012345 -f m -v e -c 4 -a e -s BEST --r_seed 1203 &")
	system("phyml -i cox1_cytb_dups.phy -b 250 --run_id GTR+I+G_boot1000a -m 012345 -f m -v e -c 4 -a e -s BEST --r_seed 1235 &")
	system("phyml -i cox1_cytb_dups.phy -b 250 --run_id GTR+I+G_boot1000b -m 012345 -f m -v e -c 4 -a e -s BEST --r_seed 1237 &")
	system("phyml -i cox1_cytb_dups.phy -b 250 --run_id GTR+I+G_boot1000c -m 012345 -f m -v e -c 4 -a e -s BEST --r_seed 1239 &")
	system("phyml -i cox1_cytb_dups.phy -b 250 --run_id GTR+I+G_boot1000d -m 012345 -f m -v e -c 4 -a e -s BEST --r_seed 1231 &")
	system("cat cox1_cytb_dups.phy_phyml_boot_trees_GTR+I+G_boot1000*.txt > boot.trees")
	
	setwd(curdir)
}

###########################
# load("flight_cost.RData")
###########################

ml.tree.bats      	<- read.tree("../../180122_flightCosts/split_analyses/1.bats/cox1_cytb_dups.phy_phyml_tree_GTR+I+G.txt")
boot.trees.bats   	<- read.tree("../../180122_flightCosts/split_analyses/1.bats/boot.trees")
ml.tree.birds      	<- read.tree("cox1_cytb_dups.phy_phyml_tree_GTR+I+G.txt")
boot.trees.birds   	<- read.tree("2.birds/boot.trees")
traits_bats		    <- read.csv("../180122_flightCosts/data_bats.csv")
traits_birds	    <- read.csv("../Flight_model_birds.csv")



####################
# work on ML trees #
####################

# bats
ml.tree.bats.ori <- root(ml.tree.bats, c("Mus_terricolor","Myrmecobius_fasciatus","Tamias_sibiricus"),resolve.root=T)
ml.tree.bats <- drop.tip(ml.tree.bats.ori, c("Mus_terricolor","Myrmecobius_fasciatus","Tamias_sibiricus"))
ml.tree.bats$edge.length[ml.tree.bats$edge.length == 0] <- .0001
# BP colors
co <- c("black", "grey", "white")
p_bats <- character(length(ml.tree.bats$node.label))
p_bats[ml.tree.bats$node.label >= .90] <- co[1]
p_bats[ml.tree.bats$node.label < .90 & ml.tree.bats$node.label >= .70] <- co[2]
p_bats[ml.tree.bats$node.label < .70] <- co[3]
pdf("ml_tree_all_bats.pdf", width=4, height=6)
par (oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
plot(ladderize(ml.tree.bats), no.margin = T)
add.scale.bar()
nodelabels(pch = 21, bg = p_bats, cex = .8, node=1: ml.tree.bats$Nnode+Ntip(ml.tree.bats))
dev.off()


# birds
ml.tree.birds.ori <- root(ml.tree.birds, c("Rhea_americana","Dromaius_novaehollandiae","Struthio_camelus"),resolve.root=T)
ml.tree.birds <- drop.tip(ml.tree.birds.ori, c("Rhea_americana","Dromaius_novaehollandiae","Struthio_camelus"))
ml.tree.birds$edge.length[ml.tree.birds$edge.length == 0] <- .0001
# BP colors
co <- c("black", "grey", "white")
p_birds <- character(length(ml.tree.birds$node.label))
p_birds[ml.tree.birds$node.label >= .90] <- co[1]
p_birds[ml.tree.birds$node.label < .90 & ml.tree.birds$node.label >= .70] <- co[2]
p_birds[ml.tree.birds$node.label < .70] <- co[3]
pdf("ml_tree_all_birds.pdf", width=4, height=6)
par (oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
plot(ladderize(ml.tree.birds), no.margin = T, cex=.2)
add.scale.bar()
nodelabels(pch = 21, bg = p_birds, cex = .8, node=1: ml.tree.birds$Nnode+Ntip(ml.tree.birds))
dev.off()

# check species consitency between trees and traits: bats
mis_bat_tree   <- setdiff(ml.tree.bats$tip.label, unique(traits_bats$Species))
mis_bat_traits <- setdiff(unique(traits_bats$Species), ml.tree.bats$tip.label) 
# [1] "Pteropus_alecto"
traits_bats <- traits_bats[-which(traits_bats$Species %in% mis_bat_traits),]

# check species consitency between trees and traits: birds
mis_bird_tree   <- setdiff(ml.tree.birds$tip.label, unique(traits_birds$Species))
mis_bird_traits <- setdiff(unique(traits_birds$Species), ml.tree.birds$tip.label) 
# "Nymphicus_hollandicus"
traits_birds <- traits_birds[-which(traits_birds$Species %in% mis_bird_traits),]


#################################
save.image("180509_flight_cost.RData")
#################################


###############################################################################
# 3. Does efficiency depend on body mass, method or phylogeny (birds)?
###############################################################################


species2drop_tab2 <- setdiff(ml.tree.birds$tip.label, traits_birds$Species)
ml.tree_birds <- drop.tip(ml.tree.birds, species2drop_tab2)
# BP colors
co <- c("black", "grey", "white")
p_birds <- character(length(ml.tree_birds$node.label))
p_birds[ml.tree_birds$node.label >= .90] <- co[1]
p_birds[ml.tree_birds$node.label < .90 & ml.tree_birds$node.label >= .70] <- co[2]
p_birds[ml.tree_birds$node.label < .70] <- co[3]

pdf("ml_tree_birds.pdf", width=4, height=6)
par (oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
plot(ladderize(ml.tree_birds), no.margin = T, cex=.2)
add.scale.bar()
nodelabels(pch = 21, bg = p_birds, cex = .8, node=1: ml.tree_birds$Nnode+Ntip(ml.tree_birds))
dev.off()


mydata_birds <- data.frame(as.numeric(log10(traits_birds$BMR)), 
						as.numeric(log10(traits_birds$Mass)), 
						as.numeric(traits_birds$Mass), 
						as.numeric(traits_birds$Method), 
						as.numeric(log10(traits_birds$Efficiency)), 
						as.numeric(traits_birds$Efficiency), 
						as.numeric(traits_birds$Speed), 
						taxa=traits_birds$Species)
colnames(mydata_birds) <- c("LogMR","LogBM","BM","Meth", "LogEffic", "Effic", "Speed", "taxa")

ml.tree_birds$node.label <- NA
cdat_tab3 <- comparative.data(phy = ml.tree_birds, data = mydata_birds, names.col = "taxa", vcv = TRUE)

pgls_21 <- pgls(LogEffic ~ LogBM + Meth + Speed, cdat_tab3)
summary(pgls_21)

pgls_22 <- pgls(LogEffic ~ LogBM + Meth + Speed, cdat_tab3, lambda='ML')
summary(pgls_22)

pgls_21$aicc
pgls_22$aicc
pgls_22$aicc - pgls_21$aicc


# over the bootstrapped trees
NREPS <- length(boot.trees.bats)
b_aicc21 <- b_aicc22 <- b_lambda2 <- b_logBM21 <- b_logMeth21 <- b_logSpeed21 <- b_logBM22  <- b_logMeth22 <- b_logSpeed22 <- b_P_logMeth22 <- b_P_logBM22 <- b_P_logSpeed22 <- rep(NA, NREPS)

for(i in 1:NREPS){
		if(is.monophyletic(boot.trees.birds[[i]], c("Rhea_americana","Dromaius_novaehollandiae","Struthio_camelus"))){
			# boot trees
			b_ml.tree.birds <- root(boot.trees.birds[[i]], c("Rhea_americana","Dromaius_novaehollandiae","Struthio_camelus"),resolve.root=T)
			b_ml.tree.birds$edge.length[b_ml.tree.birds$edge.length == 0] <- .0001
			
			# drop species
			species2drop_tab2 <- setdiff(b_ml.tree.birds$tip.label, traits_birds$Species)
			ml.tree_birds_b <- drop.tip(b_ml.tree.birds, species2drop_tab2)
			
			# comparative data
			b_cdat_tab22 <- comparative.data(phy = ml.tree_birds_b, data = mydata_birds, names.col = "taxa", vcv = TRUE)
			
			# pgls
			b_pgls_21 <- pgls(LogEffic ~ LogBM + Meth + Speed, b_cdat_tab22)
			summary(b_pgls_21)
			b_logBM21[i] <- summary(b_pgls_21)[[5]][2,1]
			b_logMeth21[i] <- summary(b_pgls_21)[[5]][3,1]
			b_logSpeed21[i] <- summary(b_pgls_21)[[5]][4,1]
			b_aicc21[i] <- b_pgls_21$aicc
			
			try(b_pgls_22 <- pgls(LogEffic ~ LogBM + Meth + Speed, b_cdat_tab22, lambda='ML')) # optim() fails on some trees...
			summary(b_pgls_22)
			b_logBM22[i] <- summary(b_pgls_22)[[5]][2,1]
			b_logMeth22[i] <- summary(b_pgls_22)[[5]][3,1]
			b_logSpeed22[i] <- summary(b_pgls_22)[[5]][4,1]
			b_aicc22[i]  <- b_pgls_22$aicc
			b_P_logBM22[i] <- summary(b_pgls_22)[[5]][2,4]
			b_P_logSpeed22[i] <- summary(b_pgls_22)[[5]][4,4]
			b_P_logMeth22[i] <- summary(b_pgls_22)[[5]][3,4]
			# lambda
			b_lambda2[i] <- b_pgls_22$param[2]
		}
}

aicc_diff2 <- b_aicc22 - b_aicc21
sum(aicc_diff2 <= 0)/1000
phylo_p3 <- sum(aicc_diff2 < -2, na.rm=T)/1000
pgls_dat3 <- data.frame(b_aicc21, b_aicc22, aicc_diff2, b_lambda2, b_logBM21, b_logBM22, b_logSpeed21, b_logSpeed22, b_logMeth21, b_logMeth22)

ggplot() + scale_x_continuous(limits = c(min(aicc_diff2,na.rm=T), 0)) + labs(x = expression(paste(Delta,"AICc")), y = "Density") +
  geom_histogram(aes(x=(aicc_diff2)), binwidth = .15, fill="red", colour="red", data= pgls_dat3) + #geom_density(aes(x=(aicc_diff)), colour="red", data= pgls_dat2)
  geom_vline(xintercept = -2, colour="grey", linetype = "longdash") +
  geom_vline(xintercept = 2, colour="grey", linetype = "longdash") + 
  theme_bw() +
  annotate("text", x = min(aicc_diff2,na.rm=T)/2, y = 10, hjust = 0, size=4, label = paste("P(DAICc < -2) = ", phylo_p3))

ggsave("birds.pgls_delat_AICc_log10.pdf", width = 4, height = 4)


P_b_logBM21 <- sum(b_logBM22 >= 0, na.rm=T)/1000
P_b_logMeth21 <- sum(b_logMeth22 >= 0, na.rm=T)/1000
P_b_logSpeed21 <- sum(b_logSpeed22 >= 0, na.rm=T)/1000
P_b_P_logBM22 <- sum(b_P_logBM22 < .05, na.rm=T)/1000
P_b_P_logSpeed22 <- sum(b_P_logSpeed22 < .05, na.rm=T)/1000
P_b_P_logMeth22 <- sum(b_P_logMeth22 < .05, na.rm=T)/1000
ggplot() + theme_bw() + scale_x_continuous(limits = c(min(c(b_logBM22, b_logMeth22, b_logSpeed22),na.rm=T), max(c(b_logBM22, b_logMeth22, b_logSpeed22),na.rm=T))) + 
  labs(x = "Slopes from bootstrapped trees (LogEfficiency)", y = "Density") +
  geom_density(aes(x=(b_logBM22)), colour="red", data= pgls_dat3) +
  geom_density(aes(x=(b_logMeth22)), colour="orange", data= pgls_dat3) +
  geom_density(aes(x=(b_logSpeed22)), colour="blue", data= pgls_dat3) +
  annotate("text", x = -.1, y = c(500,450,400), hjust = 1, size=4, colour = c("blue","orange","red"), 
  		label = c(paste0("speed (PGLS): P(slope >= 0) = ", P_b_logSpeed21, "P(P<.05) = ", P_b_P_logSpeed22), 
  				  paste0("method (PGLS): P(slope >= 0) = ", P_b_logMeth21, "P(P<.05) = ", P_b_P_logSpeed22), 
  				  paste0("body mass (PGLS): P(slope >= 0) = ", P_b_logBM21, "P(P<.05) = ", P_b_P_logSpeed22))) 

ggsave("birds.pgls_slopes_log10.pdf", width = 4, height = 4)



#################################
save.image("180509_flight_cost_birds.RData")
#################################


###############################################################################
# 3. Does efficiency depend on body mass, method or phylogeny (bats)?
###############################################################################


#species2drop_tab2 <- setdiff(ml.tree.bats$tip.label, traits_bats$Species)
#ml.tree_bats <- drop.tip(ml.tree.bats, species2drop_tab2)
#species2drop_tab3 <- setdiff(traits_bats$Species, ml.tree.bats$tip.label)
#traits_bats <- traits_bats[-which(traits_bats$Species %in% species2drop_tab3),]
ml.tree_bats <- ml.tree.bats
# BP colors
co <- c("black", "grey", "white")
p_bats <- character(length(ml.tree_bats$node.label))
p_bats[ml.tree_bats$node.label >= .90] <- co[1]
p_bats[ml.tree_bats$node.label < .90 & ml.tree_bats$node.label >= .70] <- co[2]
p_bats[ml.tree_bats$node.label < .70] <- co[3]

pdf("ml_tree_bats.pdf", width=4, height=6)
par (oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,1))
plot(ladderize(ml.tree_bats), no.margin = T, cex=.4)
add.scale.bar()
nodelabels(pch = 21, bg = p_bats, cex = .8, node=1: ml.tree_bats$Nnode+Ntip(ml.tree_bats))
dev.off()


mydata_bats <- data.frame(as.numeric(log10(traits_bats$BMR)), 
						as.numeric(log10(traits_bats$Mass)), 
						as.numeric(((traits_bats$Mass))), 
						#as.numeric(Methods2_ave), 
						as.numeric(log10(traits_bats$Efficiency)), 
						as.numeric(traits_bats$Efficiency), 
						as.numeric(traits_bats$Speed), 
						taxa=traits_bats$Species)
colnames(mydata_bats) <- c("LogMR","LogBM","BM","LogEffic", "Effic", "Speed", "taxa")

ml.tree_bats$node.label <- NA
cdat_tab3 <- comparative.data(phy = ml.tree_bats, data = mydata_bats, names.col = "taxa", vcv = TRUE)

pgls_21 <- pgls(LogEffic ~ LogBM + Speed, cdat_tab3)
summary(pgls_21)

pgls_22 <- pgls(LogEffic ~ LogBM + Speed, cdat_tab3, lambda='ML')
summary(pgls_22)

pgls_21$aicc
pgls_22$aicc
pgls_22$aicc - pgls_21$aicc


# over the bootstrapped trees
NREPS <- length(boot.trees.bats)
b_aicc21 <- b_aicc22 <- b_lambda2 <- b_logBM21 <- b_logSpeed21 <- b_logBM22 <- b_logSpeed22 <- b_P_logBM22 <- b_P_logSpeed22 <- rep(NA, NREPS)

for(i in 1:NREPS){
		# boot trees
		if(is.monophyletic(boot.trees.bats[[i]], c("Mus_terricolor","Myrmecobius_fasciatus","Tamias_sibiricus"))){
			b_ml.tree.bats <- root(boot.trees.bats[[i]], c("Mus_terricolor","Myrmecobius_fasciatus","Tamias_sibiricus"),resolve.root=T)
			b_ml.tree.bats$edge.length[b_ml.tree.bats$edge.length == 0] <- .0001
			
			# drop species
			species2drop_tab2 <- setdiff(b_ml.tree.bats$tip.label, traits_bats$Species)
			ml.tree_bats_b <- drop.tip(b_ml.tree.bats, species2drop_tab2)
			
			# comparative data
			b_cdat_tab22 <- comparative.data(phy = ml.tree_bats_b, data = mydata_bats, names.col = "taxa", vcv = TRUE)
			
			# pgls
			b_pgls_21 <- pgls(LogEffic ~ LogBM + Speed, b_cdat_tab22)
			summary(b_pgls_21)
			b_logBM21[i] <- summary(b_pgls_21)[[5]][2,1]
			b_logSpeed21[i] <- summary(b_pgls_21)[[5]][3,1]
			b_aicc21[i] <- b_pgls_21$aicc
			
			try(b_pgls_22 <- pgls(LogEffic ~ LogBM + Speed, b_cdat_tab22, lambda='ML')) # optim() fails on some trees...
			summary(b_pgls_22)
			b_logBM22[i] <- summary(b_pgls_22)[[5]][2,1]
			b_logSpeed22[i] <- summary(b_pgls_22)[[5]][3,1]
			b_aicc22[i]  <- b_pgls_22$aicc
			b_P_logBM22[i] <- summary(b_pgls_22)[[5]][2,4]
			b_P_logSpeed22[i] <- summary(b_pgls_22)[[5]][3,4]
			# lambda
			b_lambda2[i] <- b_pgls_22$param[2]
		}
}

aicc_diff2 <- b_aicc22 - b_aicc21
sum(aicc_diff2 <= 0)/1000
phylo_p3 <- sum(aicc_diff2 < -2, na.rm=T)/1000
pgls_dat3 <- data.frame(b_aicc21, b_aicc22, aicc_diff2, b_lambda2, b_logBM21, b_logBM22, b_logSpeed21, b_logSpeed22, b_logMeth21, b_logMeth22)

ggplot() + scale_x_continuous(limits = c(min(aicc_diff2,na.rm=T), 0)) + 
  labs(x = expression(paste(Delta,"AICc")), y = "Density") +
  theme_bw() +
  geom_histogram(aes(x=(aicc_diff2)), binwidth = .15, fill="red", colour="red", data= pgls_dat3) + #geom_density(aes(x=(aicc_diff)), colour="red", data= pgls_dat2)
  geom_vline(xintercept = -2, colour="grey", linetype = "longdash") +
  geom_vline(xintercept = 2, colour="grey", linetype = "longdash") + 
  annotate("text", x = min(aicc_diff2,na.rm=T)/2, y = 10, hjust = 0, size=4, label = paste("P(DAICc < -2) = ", phylo_p3))

ggsave("bats.pgls_delat_AICc_log10.pdf", width = 4, height = 4)


P_b_logBM22 <- sum(b_logBM22 >= 0, na.rm=T)/1000
#P_b_logMeth22 <- sum(b_logMeth22 >= 0, na.rm=T)/1000
P_b_logSpeed22 <- sum(b_logSpeed22 >= 0, na.rm=T)/1000
P_b_P_logBM22 <- sum(b_P_logBM22 < .05, na.rm=T)/1000
P_b_P_logSpeed22 <- sum(b_P_logSpeed22 < .05, na.rm=T)/1000
ggplot() + theme_bw() + scale_x_continuous(limits = c(min(c(b_logBM22, b_logSpeed22),na.rm=T), max(c(b_logBM22, b_logSpeed22),na.rm=T))) + 
  labs(x = "Slopes from bootstrapped trees (LogEfficiency)", y = "Density") +
  geom_density(aes(x=(b_logBM22)), colour="red", data= pgls_dat3) +
  #geom_density(aes(x=(b_logMeth21)), colour="orange", data= pgls_dat3) +
  geom_density(aes(x=(b_logSpeed22)), colour="blue", data= pgls_dat3) +
  #scale_y_continuous(trans=log10_trans()) +
  #scale_y_continuous(limits=c(1,20)) +
  annotate("text", x = .5, y = c(120,
  								 #110,
  								 100), hjust = 1, size=4, colour = c("blue",
  																	 #"orange",
  																	 "red"), 
  		label = c(paste0("speed (PGLS): P(slope >= 0) = ", P_b_logSpeed22, "P(P<.05) = ", P_b_P_logSpeed22), 
  				  #paste0("method (PGLS): P(slope >= 0) = ", P_b_logMeth22), 
  				  paste0("body mass (PGLS): P(slope >= 0) = ", P_b_logBM22, "P(P<.05) = ", P_b_P_logBM22))) 

ggsave("bats.pgls_slopes_log10.pdf", width = 4, height = 4)


posGT50 <- which(aicc_diff2 > -100)
posLT50 <- which(aicc_diff2 < -100)
treesTG50 <- treesLT50 <- list()
#
for(i in 1:length(posGT50)){
	rtr <- root(boot.trees.bats[[posGT50[i]]], c("Mus_terricolor","Myrmecobius_fasciatus","Tamias_sibiricus"),resolve.root=T)
	treesTG50[[i]] <- drop.tip(rtr, species2drop_tab2)
}
treeGT50 <- consensus(treesTG50, p=.95)
plot(treeGT50)
#
for(i in 1:length(posLT50)){
	rtr <- root(boot.trees.bats[[posLT50[i]]], c("Mus_terricolor","Myrmecobius_fasciatus","Tamias_sibiricus"),resolve.root=T)
	treesLT50[[i]] <- drop.tip(rtr, species2drop_tab2)
}
treeLT50 <- consensus(treesLT50, p=.95)
plot(treeLT50)

pdf("ml_tree_bats_1_2.pdf", width=18, height=6)
par (oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(1,2))
plot(treeLT50, main="< -100")
plot(treeGT50, main="> -100")
add.scale.bar()
dev.off()


#################################
save.image("180509_flight_cost_bats.RData")
#################################




# added 180706 #
LogE_birds <- pic(mydata_birds$LogEffic, ml.tree_birds)
LogM_birds <- pic(mydata_birds$LogBM, ml.tree_birds)
LogS_birds <- pic(mydata_birds$Speed, ml.tree_birds)
m_birds <- nls(LogE_birds ~ b * LogM_birds + (a + g * LogS_birds), start=list(a=5.44, b=0.64, g=-0.12))
summary(m_birds)
cor(LogE_birds, predict(m_birds))

# filtering out S < 2m/s?
which(mydata_birds$Speed < log10(2))
# integer(0)


LogE_bats <- pic(mydata_bats$LogEffic, ml.tree_bats)
LogM_bats <- pic(mydata_bats$LogBM, ml.tree_bats)
LogS_bats <- pic(mydata_bats$Speed, ml.tree_bats)
m_bats <- nls(LogE_bats ~ b * LogM_bats + (a + g * LogS_bats), start=list(a=2.94, b=-.32, g=-.18))
summary(m_bats)
cor(LogE_bats, predict(m_bats))

# filtering out S < 2m/s?
which(mydata_bats$Speed < log10(2))
# integer(0)
									###########
									# the end #
									###########
								
