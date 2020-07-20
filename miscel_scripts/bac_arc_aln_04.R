# This script was used in Colby et al. (revisions) to compare the 
# bacterial and archeal trees, with the full origianl alignments
# that contained missing data, and trimmed alignments where sites
# with > 50% indels where eliminated. The metric used to quantify
# the distance btw the full and half-gapped trees is the RF distance
# Because of the size of these alignments, bootstrapped trees could not 
# be used (too expensive), so we resorted to estimating a null
# distribution by finding SPR neighbors of the half-gapped trees.


library(seqinr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ape)
library(phangorn)
library(foreach)
library(doMC)
registerDoMC(cores=10)


#############
# WHOLE ALN #
#############

					###############
					# acrcheo aln #
					###############

arc <- read.fasta("arc_concat_edited_trimal_gappy.fa", seqtype = "AA")

# a simple count
word1count_arc <- NULL
for(i in 1:length(arc)){
	word1count_arc <- rbind(word1count_arc, count(arc[[i]], word = 1, alphabet = s2c(paste(a(aaa())[-1], "-", collapse=""))) / length(arc[[i]]) )
}
word1count_arc <- word1count_arc[,-1]

m_word1count_arc <- melt(word1count_arc)
names(m_word1count_arc) <- c("row", "AA", "Proportions")


box_arc <- ggplot(m_word1count_arc, aes(AA, Proportions)) + geom_boxplot() + theme_bw()


					################
					# bacterio aln #
					################

bac <- read.fasta("bac_concat_aln_50_trimal_gappy.fa", seqtype = "AA")

# a simple count
word1count_bac <- NULL
for(i in 1:length(bac)){
	word1count_bac <- rbind(word1count_bac, count(bac[[i]], word = 1, alphabet = s2c(paste(a(aaa())[-1], "-", collapse=""))) / length(bac[[i]]) )
}
word1count_bac <- word1count_bac[,-1]

m_word1count_bac <- melt(word1count_bac)
names(m_word1count_bac) <- c("row", "AA", "Proportions")


box_bac <- ggplot(m_word1count_bac, aes(AA, Proportions)) + geom_boxplot() + theme_bw()

# corr
word1count_arc_med <- tapply(m_word1count_arc$Proportions, m_word1count_arc$AA, median)
word1count_bac_med <- tapply(m_word1count_bac$Proportions, m_word1count_bac$AA, median)
word1count_arc_sd <- tapply(m_word1count_arc$Proportions, m_word1count_arc$AA, sd)
word1count_bac_sd <- tapply(m_word1count_bac$Proportions, m_word1count_bac$AA, sd)
word1count_arc_Q <- tapply(m_word1count_arc$Proportions, m_word1count_arc$AA, quantile,  probs = c(.1, .9))
word1count_bac_Q <- tapply(m_word1count_bac$Proportions, m_word1count_bac$AA, quantile,  probs = c(.1, .9))

word1count_arc_Q10 <- word1count_arc_Q90 <- word1count_bac_Q10 <- word1count_bac_Q90 <- NULL
for(i in 1:length(word1count_arc_Q)){
	word1count_arc_Q10[i] <- word1count_arc_Q[[i]][1]
	word1count_arc_Q90[i] <- word1count_arc_Q[[i]][2]
	word1count_bac_Q10[i] <- word1count_bac_Q[[i]][1]
	word1count_bac_Q90[i] <- word1count_bac_Q[[i]][2]
}

word1counts <- cbind(word1count_arc_med, word1count_bac_med, 
					word1count_arc_Q10, word1count_arc_Q90, 
					word1count_bac_Q10, word1count_bac_Q90 
					)
word1counts <- data.frame(word1counts)
colnames(word1counts) <- c("Median_Proportions_Archeae", "Median_Proportions_Bacteria", 
					"Q10_Proportions_Archeae", 
					"Q90_Proportions_Archeae", 
					"Q10_Proportions_Bacteria", 
					"Q90_Proportions_Bacteria")

corr_bac_arc <- ggplot(word1counts, aes(Median_Proportions_Bacteria, Median_Proportions_Archeae)) + geom_point() + theme_bw() +
			geom_segment(word1counts, mapping=aes(x= Median_Proportions_Bacteria , y= Q10_Proportions_Archeae, 
						xend= Median_Proportions_Bacteria, yend= Q90_Proportions_Archeae)) +
			geom_segment(word1counts, mapping=aes(x= Q10_Proportions_Bacteria, y= Median_Proportions_Archeae, 
						xend= Q90_Proportions_Bacteria, yend= Median_Proportions_Archeae)) +
			geom_abline(slope=1, intercept=0, color="gray", linetype="longdash") +
			geom_text(aes(label = rownames(word1counts)), hjust = rnorm(1, .05, .01), nudge_y = rnorm(1, .05, .01), size = 2) 


pdf("boxplots_word1_full.pdf", width=6, height=12)
grid.arrange( box_bac, box_arc, corr_bac_arc,
	nrow = 3
)
dev.off()




###############
# GBLOCKS ALN #
###############

					###############
					# acrcheo aln #
					###############

system("gblocks arc_concat_edited_trimal_gappy.fa -t=p -e=-half -b5=h")
arc50 <- read.fasta("arc_concat_edited_trimal_gappy.fa-half", seqtype = "AA")

# a simple count
word1count_arc50 <- NULL
for(i in 1:length(arc50)){
	word1count_arc50 <- rbind(word1count_arc50, count(arc50[[i]], word = 1, alphabet = s2c(paste(a(aaa())[-1], "-", collapse=""))) / length(arc[[i]]) )
}
word1count_arc50 <- word1count_arc50[,-1]

m_word1count_arc50 <- melt(word1count_arc50)
names(m_word1count_arc50) <- c("row", "AA", "Proportions")


box_arc50 <- ggplot(m_word1count_arc50, aes(AA, Proportions)) + geom_boxplot() + theme_bw()


					################
					# bacterio aln #
					################

system("gblocks bac_concat_aln_50_trimal_gappy.fa -t=p -e=-half -b5=h")
bac50 <- read.fasta("bac_concat_aln_50_trimal_gappy.fa-half", seqtype = "AA")

# a simple count
word1count_bac50 <- NULL
for(i in 1:length(bac50)){
	word1count_bac50 <- rbind(word1count_bac50, count(bac50[[i]], word = 1, alphabet = s2c(paste(a(aaa())[-1], "-", collapse=""))) / length(bac[[i]]) )
}
word1count_bac50 <- word1count_bac50[,-1]

m_word1count_bac50 <- melt(word1count_bac50)
names(m_word1count_bac50) <- c("row", "AA", "Proportions")


box_bac50 <- ggplot(m_word1count_bac50, aes(AA, Proportions)) + geom_boxplot() + theme_bw()

# corr
word1count_arc50_med <- tapply(m_word1count_arc50$Proportions, m_word1count_arc50$AA, median)
word1count_bac50_med <- tapply(m_word1count_bac50$Proportions, m_word1count_bac50$AA, median)
word1count_arc50_sd <- tapply(m_word1count_arc50$Proportions, m_word1count_arc50$AA, sd)
word1count_bac50_sd <- tapply(m_word1count_bac50$Proportions, m_word1count_bac50$AA, sd)
word1count_arc50_Q <- tapply(m_word1count_arc50$Proportions, m_word1count_arc50$AA, quantile,  probs = c(.1, .9))
word1count_bac50_Q <- tapply(m_word1count_bac50$Proportions, m_word1count_bac50$AA, quantile,  probs = c(.1, .9))

word1count_arc50_Q10 <- word1count_arc50_Q90 <- word1count_bac50_Q10 <- word1count_bac50_Q90 <- NULL
for(i in 1:length(word1count_arc50_Q)){
	word1count_arc50_Q10[i] <- word1count_arc50_Q[[i]][1]
	word1count_arc50_Q90[i] <- word1count_arc50_Q[[i]][2]
	word1count_bac50_Q10[i] <- word1count_bac50_Q[[i]][1]
	word1count_bac50_Q90[i] <- word1count_bac50_Q[[i]][2]
}

word1counts <- cbind(word1count_arc50_med, word1count_bac50_med, 
					word1count_arc50_Q10, word1count_arc50_Q90, 
					word1count_bac50_Q10, word1count_bac50_Q90 
					)
word1counts <- data.frame(word1counts)
colnames(word1counts) <- c("Median_Proportions_Archeae", "Median_Proportions_Bacteria", 
					"Q10_Proportions_Archeae", 
					"Q90_Proportions_Archeae", 
					"Q10_Proportions_Bacteria", 
					"Q90_Proportions_Bacteria")

corr_bac50_arc50 <- ggplot(word1counts, aes(Median_Proportions_Bacteria, Median_Proportions_Archeae)) + geom_point() + theme_bw() +
			geom_segment(word1counts, mapping=aes(x= Median_Proportions_Bacteria , y= Q10_Proportions_Archeae, 
						xend= Median_Proportions_Bacteria, yend= Q90_Proportions_Archeae)) +
			geom_segment(word1counts, mapping=aes(x= Q10_Proportions_Bacteria, y= Median_Proportions_Archeae, 
						xend= Q90_Proportions_Bacteria, yend= Median_Proportions_Archeae)) +
			geom_abline(slope=1, intercept=0, color="gray", linetype="longdash") +
			geom_text(aes(label = rownames(word1counts)), hjust = rnorm(1, .005, .001), nudge_y = rnorm(1, .005, .001), size = 2) 


pdf("boxplots_word1_50.pdf", width=6, height=12)
grid.arrange( box_bac50, box_arc50, corr_bac50_arc50,
	nrow = 3
)
dev.off()






############
# FASTTREE #
############

# full
system("fasttree -wag arc_concat_edited_trimal_gappy.fa > arc_concat_edited_trimal_gappy.tre&")
system("fasttree -wag bac_concat_aln_50_trimal_gappy.fa > bac_concat_aln_50_trimal_gappy.tre&")

# half
system("fasttree -wag arc_concat_edited_trimal_gappy.fa-half > arc_concat_edited_trimal_gappy-half.tre&")
system("fasttree -wag bac_concat_aln_50_trimal_gappy.fa-half > bac_concat_aln_50_trimal_gappy-half.tre&")

# arc
arc_tr <- read.tree("arc_concat_edited_trimal_gappy.tre")
bin_pos_arc <- grep("^Bin", arc_tr$tip.label)
tmp_root_arc <- arc_tr$tip.label[bin_pos_arc][1]
arc_tr <- root(arc_tr, tmp_root_arc, resolve.root=T)

arc_tr50 <- read.tree("arc_concat_edited_trimal_gappy-half.tre")
bin_pos_arc50 <- grep("^Bin", arc_tr50$tip.label)
arc_tr50 <- root(arc_tr50, tmp_root_arc, resolve.root=T)

pdf("trees_arc.pdf", width=10, height=10)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,2))

#plot(arc_tr, type="u", show.tip.label=F, no.margin=T)
plot(arc_tr, direction="r", show.tip.label=F, no.margin=T)
text(1.3, .1, "Archeae -- full aln", font = 2)
tiplabels("o", bin_pos_arc, frame = "n", adj = -2.5, col="red")
add.scale.bar()

#plot(arc_tr50, type="u", show.tip.label=F, no.margin=T)
plot(arc_tr50, direction="l", show.tip.label=F, no.margin=T)
text(.9, .1, "Archeae -- 1/2 gap aln", font = 2)
tiplabels("o", bin_pos_arc50, frame = "n", adj = 2.5, col="red")
add.scale.bar()

dev.off()





# bac
bac_tr <- read.tree("bac_concat_edited_trimal_gappy.tre")
bin_pos_bac <- grep("^Bin", bac_tr$tip.label)
tmp_root_bac <- bac_tr$tip.label[bin_pos_bac][1]
bac_tr <- root(bac_tr, tmp_root_bac, resolve.root=T)

bac_tr50 <- read.tree("bac_concat_edited_trimal_gappy-half.tre")
bin_pos_bac50 <- grep("^Bin", bac_tr50$tip.label)
bac_tr50 <- root(bac_tr50, tmp_root_bac, resolve.root=T)

pdf("trees_bac.pdf", width=10, height=10)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,2))

#plot(bac_tr, type="u", show.tip.label=F, no.margin=T)
plot(bac_tr, direction="r", show.tip.label=F, no.margin=T)
text(1.3, .1, "Bacteria -- full aln", font = 2)
tiplabels("o", bin_pos_bac, frame = "n", adj = -2.5, col="red")
add.scale.bar()

#plot(bac_tr50, type="u", show.tip.label=F, no.margin=T)
plot(bac_tr50, direction="l", show.tip.label=F, no.margin=T)
text(.9, .1, "Bacteria -- 1/2 gap aln", font = 2)
tiplabels("o", bin_pos_bac50, frame = "n", adj = 2.5, col="red")
add.scale.bar()

dev.off()



#################################
# NNI RF distance distributions #
#################################

# arc 
arc_tr <- read.tree("arc_concat_edited_trimal_gappy.tre")
arc_tr50 <- read.tree("arc_concat_edited_trimal_gappy-half.tre")

RF_dist_obs <- as.numeric(dist.topo(multi2di(arc_tr), multi2di(arc_tr50)))

Nreps <- 10000
RF_dist_rSPR_arc <- NULL
counter_arc <- 0

RF_dist_rSPR_arc <- foreach(i = 1: Nreps, .combine='rbind')%dopar%{
	if(!(i %% 25)){
		print(paste0("Now doing ", i, "/", Nreps))
	}
	#tmp <- as.numeric(dist.topo(multi2di(arc_tr), rSPR(multi2di(arc_tr50), 2)))
	tmp <- as.numeric(dist.topo(multi2di(arc_tr), rSPR(multi2di(arc_tr50), 1)))
	if(tmp < RF_dist_obs){
		counter_arc <- counter_arc + 1
	}
	return(tmp)
}

RF_dist_rSPR_arc_df <- data.frame(RF_dist_rSPR_arc)

h1 <- ggplot(RF_dist_rSPR_arc_df, aes(x=RF_dist_rSPR_arc)) + geom_histogram() + theme_bw() +
		labs(x="RF distance after SPR", y = "Distribution") +
		scale_x_continuous(n.breaks = 10) + 
		geom_vline(xintercept = RF_dist_obs, col="red", linetype="longdash") +
		geom_text(x=675, y=1200, label=paste0("P < ", format(round(counter_arc / Nreps, 2), nsmall = 2)), col="red")

save(RF_dist_rSPR_arc_df, counter_arc, file="RF_dist_rSPR_arc_df.RData")
save.image("bac_arc_aln_04.RData")



# bac 
bac_tr <- read.tree("bac_concat_edited_trimal_gappy.tre")
bac_tr50 <- read.tree("bac_concat_edited_trimal_gappy-half.tre")

RF_dist_obs <- as.numeric(dist.topo(multi2di(bac_tr), multi2di(bac_tr50)))

Nreps <- 10000
RF_dist_rSPR_bac <- NULL
counter_bac <- 0

RF_dist_rSPR_bac <- foreach(i = 1: Nreps, .combine='rbind')%dopar%{
	if(!(i %% 25)){
		print(paste0("Now doing ", i, "/", Nreps))
	}
	#tmp <- as.numeric(dist.topo(multi2di(bac_tr), rSPR(multi2di(bac_tr50), 2)))
	tmp <- as.numeric(dist.topo(multi2di(bac_tr), rSPR(multi2di(bac_tr50), 1)))
	if(tmp < RF_dist_obs){
		counter_bac <- counter_bac + 1
	}
	return(tmp)
}

RF_dist_rSPR_bac_df <- data.frame(RF_dist_rSPR_bac)

h2 <- ggplot(RF_dist_rSPR_bac_df, aes(x=RF_dist_rSPR_bac)) + geom_histogram() + theme_bw() +
		labs(x="RF distance after SPR", y = "Distribution") + 
		scale_x_continuous(n.breaks = 10) + 
		geom_vline(xintercept = RF_dist_obs, col="red", linetype="longdash") +
		geom_text(x=690, y=100, label=paste0("P < ", format(round(counter_bac / Nreps, 2), nsmall = 2)), col="red")

save(RF_dist_rSPR_bac_df, counter_bac, file="RF_dist_rSPR_bac_df.RData")



pdf("RF_dist_distrib.pdf", width=6, height=4)
grid.arrange( h1, #h2,
	nrow = 1
)
dev.off()


save.image("bac_arc_aln_04.RData")
q("no")







