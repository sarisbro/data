# reads blast output format 6 and assigns taxonomy
# contigs assembled by Trinity in a metagenomics project are in subfolder trinity/


# download and uncompress in current folder:
# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.Z

library(foreach)
library(doMC)
registerDoMC(cores=16)

#################################################################################################
#################################################################################################
#################################################################################################


						######################
						# reads BLAST output #
						######################
assembly <- "trinity"
out_files <- list.files(paste0(assembly, "/"), pattern = ".out")

seq_count <- 0
seq_name <- NULL
blast_tab <- NULL 

#for(i in rev(out_files)[1:1]){
for(i in rev(out_files)){
	print(paste0("Now doing ", i))
	ln <- NULL
	ln <- readLines(paste0(assembly, "/", i))
	nbln <- length(ln)
		
	# we read line by line
	start_reading <- query_found <- 0
	skip_line <- 1
	for(j in 1:nbln){
		
		# new query?
		if(grepl("^Query=", ln[j])){
			seq_count <- seq_count + 1
			seq_name[seq_count] <- unlist(strsplit(ln[j], " "))[2]
			blast_tmp <- NULL
		}
		
		# start reading
		if(grepl("^Sequences producing significant alignments", ln[j])){
			start_reading <- 1
		}else{
			# only lines starting w/ a letter
			if(start_reading & skip_line){
				skip_line <- 0
			}else{
				if(start_reading & grepl("^[A-Z]", ln[j]) & !skip_line){
					blast_tmp <- rbind(blast_tmp, ln[j])
				}
				
				if(start_reading & !grepl("^[A-Z]", ln[j]) & !skip_line){
					start_reading <- 0
					skip_line <- 1
					blast_tab[[seq_count]] <- blast_tmp
				}
			}
			
		}
		
		
	}

}

save(blast_tab, file="blast_tab.RData")
save.image("parse_blastn_res.RData")

#################################################################################################
#################################################################################################
#################################################################################################



						#########################
						# parse in tabular form #
						#########################


nseq <- length(seq_name)
blast_tabular <- NULL

for(i in 1:nseq){
	if(!(i %% 1000)){
		print(paste0("Now doing ", i, " -- out of ", nseq))
	}
	tmp <- blast_tab[[i]]
	newline <- NULL
	if(length(tmp) > 0){
		for(j in 1:length(tmp)){
			curline <- GI <- Species <- Evalue <- NULL
			curline <- tmp[j]
			curline_split <- unlist(strsplit(curline, "\\s+"))
			GI <- curline_split[1]
			Genus <- curline_split[2]
			Species <- paste(curline_split[2], curline_split[3])
			Evalue <- curline_split[length(curline_split)]
			newline <- rbind(newline, cbind(GI = GI, Genus = Genus, Species = Species, Evalue = Evalue))
		}
		
	}
	blast_tabular[[i]] <- newline
}

save(blast_tabular, file="blast_tabular.RData")
save.image("parse_blastn_res.RData")


#################################################################################################
#################################################################################################
#################################################################################################


						#####################
						# ID super kingdoms #
						#####################

# WARNING -- this chunk is *very* slow, with a large memory footprint

tax_names <- readLines("taxdump/names.dmp")
tax_nodes <- read.table("taxdump/nodes.dmp", sep="|")

admiss_tax_ranks <- c("clade", "phylum", "class", "order", "family")

# recursively finds taxonomy -- all the way to superkingdom
find_superKing <- function(cur_tax_ID, whole_tax){
	pos_names2 <- NULL
	pos_names2 <- which(tax_nodes[,1] == cur_tax_ID)
	parent_1 <- tax_nodes[pos_names2,2][1]
	tax_1 <- gsub("\t", "", as.character(tax_nodes[pos_names2,3][1]))
	if(tax_1 == "superkingdom"){
		pos <- tmp <- NULL
		pos <- grep(paste0("^", cur_tax_ID, "\t"), tax_names)[1]
		tmp <- tax_names[pos]
		tmp <- gsub("\\\"", "", tmp)
		tmp <- gsub("\\t", "", tmp)
		tmp <- unlist(strsplit(tmp, "\\|"))[2]
		return(c(tmp, paste(tmp, whole_tax, sep=";")))
	}else{
		if(tax_1 %in% admiss_tax_ranks){
			parent_1_name <- grep(paste0("^", parent_1, "\t"), tax_names)[1]
			parent_1_name <- tax_names[parent_1_name]
			parent_1_name <- gsub("\\\"", "", parent_1_name)
			parent_1_name <- gsub("\\t", "", parent_1_name)
			parent_1_name <- unlist(strsplit(parent_1_name, "\\|"))[2]
		}else{
			parent_1_name <- ""
		}
		find_superKing(as.character(parent_1), paste(parent_1_name, whole_tax, sep=";"))
	}
}



SupKing_lst <- list()

for(i in 1:nseq){
#for(i in 3:nseq){
	tmp_tab <- blast_tabular[[i]]
	print(paste0("Now doing ", i, " -- out of ", nseq))

	if(length(tmp_tab) > 0){
		sp_list <- tmp_tab[,3]
		sp_list <- gsub(",", "", sp_list)
		sp_list <- gsub("Uncultured", "", sp_list)
		sp_list <- gsub(" sp.", "", sp_list)
		sp_list <- gsub("^ ", "", sp_list)
		supKing <- NULL
		supKing <- foreach(j = 1:length(sp_list), .combine='rbind')%dopar%{
			# find taxon ID
			pos_names <- NULL
			pos_names <- grep(sp_list[j], tax_names)
			if(length(pos_names) > 0){
				tax_ID <- unlist(strsplit(tax_names[pos_names[1]], "\\|"))[1]
				tax_ID <- gsub("\t", "", tax_ID)
				#print(find_superKing(tax_ID, sp_list[j]))
				return(find_superKing(tax_ID, sp_list[j]))
			}else{
				return("Unidentified")
			}
		}
	}
	supKing[,2] <- gsub(";+", ";", supKing[,2])
	SupKing_lst[[i]] <- cbind(superkingdom = as.character(supKing[,1]), taxonomy = as.character(supKing[,2]), Evalue = tmp_tab[,4])
	save(SupKing_lst, file="SupKing_lst.RData")
}



#################################################################################################
#################################################################################################
#################################################################################################


save.image("parse_blastn_res.RData")
q("no")


#################################################################################################
#################################################################################################
#################################################################################################










































