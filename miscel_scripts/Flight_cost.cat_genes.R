# First script used in Guigueno M.F., Shoji A., Elliott K.H. and Aris-Brosou S. 2019. 
# Flight costs in volant vertebrates: a phylogenetically-controlled meta-analysis of birds and bats. 
# Comparative Biochemistry and Physiology - Part A: Molecular & Integrative Physiology. 235:193-201.
# https://doi.org/10.1016/j.cbpa.2019.06.003

library(ape)



# read aln files
cox1 <- read.dna("cox1/cox1.aln",format="fasta")
cytb <- read.dna("cytB/cytB.aln",format="fasta")

# get seq names
n_cox1 <- attributes(cox1)$dimnames[[1]]
n_cytb <- attributes(cytb)$dimnames[[1]]

# cleanup seq names
# cox1
sp_name <- NULL
for(i in 1:length(n_cox1)){
	#if(!grepl("_", n_cox1[i])){
	if(nchar(n_cox1[i]) > 30){
		tmp_name <- strsplit(n_cox1[i], " ")
		sp_name[i] <- paste0(tmp_name[[1]][2], "_", tmp_name[[1]][3])
	}else{
		sp_name[i] <- n_cox1[i]
	}
}
attributes(cox1)$dimnames[[1]] <- sp_name
# cytb
sp_name <- NULL
for(i in 1:length(n_cytb)){
	#if(!grepl("_", n_cytb[i])){
	if(nchar(n_cytb[i]) > 30){
		tmp_name <- strsplit(n_cytb[i], " ")
		sp_name[i] <- paste0(tmp_name[[1]][2], "_", tmp_name[[1]][3])
	}else{
		sp_name[i] <- n_cytb[i]
	}
}
attributes(cytb)$dimnames[[1]] <- sp_name

# update seq names
n_cox1 <- attributes(cox1)$dimnames[[1]]
n_cytb <- attributes(cytb)$dimnames[[1]]

# ID missing sequences in each alignment
intersect(n_cytb, n_cox1)
missing_cox1 <- setdiff(n_cytb, n_cox1)
missing_cytb <- setdiff(n_cox1, n_cytb)

				#############################
				# fill in missing sequences #
				#############################
# cox1
write.dna(cox1, "cox1_tmp.aln", format = "fasta")
for(i in 1:length(missing_cox1)){
	cat(paste0("
	> ", missing_cox1[i], "
	", paste(rep("N", dim(cox1)[2]), collapse="")) , file = "cox1_tmp.aln", sep = "", append = T)
}
# cytb
write.dna(cytb, "cytb_tmp.aln", format = "fasta")
for(i in 1:length(missing_cytb)){
	cat(paste0("
	> ", missing_cytb[i], "
	", paste(rep("N", dim(cytb)[2]), collapse="")) , file = "cytb_tmp.aln", sep = "", append = T)
}

# reload aln files
cox1 <- read.dna("cox1_tmp.aln",format="fasta")
cytb <- read.dna("cytb_tmp.aln",format="fasta")

# write concatenated alignment
cat_aln <- cbind(cox1, cytb)
write.dna(cat_aln, "cox1_cytb.aln")
write.dna(cat_aln, "cox1_cytb.fas", format="fasta")

				#########################################################
				# duplicate sequences as per # of samples in trait file #
				#########################################################
cat_aln <- read.dna("cox1_cytb.aln")
n_cat_aln <- attributes(cat_aln)$dimnames[[1]]

traits <- read.csv("../Flight_model_birds.csv", header=T)
sp_list <- traits$Species

# "cleanup" sp names
sp_list_clean <- NULL
for(i in 1:length(sp_list)){
	sp_list_clean[i] <- sub(" ", "_", as.character(sp_list[i]))
}

sp_list_clean_unique <- unique(sp_list_clean)

# sanity check...
setdiff(n_cat_aln, sp_list_clean_unique)
setdiff(sp_list_clean_unique, n_cat_aln)


tab_sp_list_clean <- table(sp_list_clean)
dup_list <- tab_sp_list_clean[tab_sp_list_clean > 1]
system("mkdir tmp_aln")

for(i in 1:length(dup_list)){
	nrep <- spname <- NULL
	nrep <- as.numeric(dup_list[i]) - 1
	spname <- names(dup_list[i])
	pos <- which(n_cat_aln == spname)
	for(j in 1:nrep){
		attributes(cat_aln)$dimnames[[1]][pos] <- paste0(spname, j)
		write.dna(cat_aln[pos,], paste0("tmp_aln/", spname, j,".aln"), format="fasta")
	}
}

system("cat cox1_cytb.fas tmp_aln/* > cox1_cytb_all.fas")
cat_aln_all <- read.dna("cox1_cytb_all.fas",format="fasta")
write.dna(cat_aln_all, "cox1_cytb_all.aln")



