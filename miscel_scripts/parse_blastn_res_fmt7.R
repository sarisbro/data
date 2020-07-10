# reads blast output format 7
# contigs assembled by Trinity in a metagenomics project are in subfolder trinity/


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
blast_tab <- blast_tmp <- NULL 

for(i in rev(out_files)){
	print(paste0("Now doing ", i))
	ln <- NULL
	ln <- readLines(paste0(assembly, "/", i))
	nbln <- length(ln)
		
	# we read line by line
	start_reading <- query_found <- 0
	for(j in 1:nbln){
		
		# new query?
		if(grepl("^\\# Query:", ln[j])){
			seq_count <- seq_count + 1
			seq_name[seq_count] <- unlist(strsplit(ln[j], " "))[3]
			blast_tmp <- NULL
		}
		
		# start reading
		if(grepl(" hits found", ln[j])){
			start_reading <- 1
		}else{
			# lines w/ blasts res
			if(start_reading & grepl("^[A-Z]", ln[j])){
			blast_tmp <- rbind(blast_tmp, ln[j])
			}
			# first line of next result: stop reading
			if(start_reading & grepl("^\\#", ln[j]) ){
			start_reading <- 0
			blast_tab[[seq_count]] <- blast_tmp
			}
			
		}
		
	}

}

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
			curline <- SK <- tax_ID <- PC_ID <- BS <- Species <- Evalue <- NULL
			curline <- tmp[j]
			curline_split <- unlist(strsplit(curline, "\\t"))
			SK 		<- curline_split[1]
			tax_ID	<- curline_split[2]
			Species	<- curline_split[3]
			PC_ID	<- curline_split[6]
			BS		<- curline_split[7]
			Evalue	<- curline_split[8]
			newline <- rbind(newline, cbind(superkingdom = SK, Species = Species,
						tax_ID = tax_ID, PC_ID = PC_ID, BitScore = BS, Evalue = Evalue))
		}
		
	}
	blast_tabular[[i]] <- newline
}





#################################################################################################
#################################################################################################
#################################################################################################


save.image("parse_blastn_res.RData")
q("no")


#################################################################################################
#################################################################################################
#################################################################################################










































