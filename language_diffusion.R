migrate <- function(ch,topfile) {
	# find where to start
	areafiles <- c("AFRICA.txt", "AMERICAS.txt", "AUSTRALIA.txt", "EURASIA.txt", "OCEANIA.txt")
	# specify the proportion of populated places in each area out of all ~2.9 million
	prop.p <- c(0.1385, 0.1771, 0.0034, 0.6001, 0.0809) 
	area <- sample(1:5, 1, replace=TRUE, prop.p)
	geodata <- read.table(file=areafiles[area])
	line <- sample(length(geodata[,1]), 1)
	startdata <- geodata[line,]
	code <- as.vector(startdata[1,3])
	codes <- read.table("CountryCodes.txt", quote="", stringsAsFactors=FALSE, header=FALSE, sep="\t")
	country <- codes[match(code, codes[,1]),2]
	orilat <- startdata[,1]
	orilong <- startdata[,2]
	geodata <- geodata[,-3]
	cat("The homeland is in ", country, " at the coordinates ", orilat, ", ", orilong, "\n", sep="")

	# migrate
	# make a table where rows are steps, columns are taxa, and 1's is presence of a lineage at a step 
	first <- gsub(".top","",topfile)
	top <- read.table(topfile)
	first.step <- min(top[,3])
	last.step <- as.numeric(strsplit(readLines(paste(first,".par",sep=""))[8], " ")[[1]][3])
	no.taxa <- as.numeric(strsplit(first, "_")[[1]][2])
	taxa <- sort(unique(c(top[,1],top[,2])))
	lats <- matrix(NA, nrow=last.step-first.step+1, ncol=no.taxa, dimnames=list(first.step:last.step,taxa))
	# find rownumber for first occurrence of a taxon and put 1's in the matrix starting from the 
	# step where a taxon is born
	for (i in 1:length(taxa)) {
		if (i < 3) {
			lats[1:length(lats[,1]),i] <- 1
		}
		if (i >= 3) {
			topline.born <- which(top[,2]==taxa[i], arr.ind = TRUE)
			when.born <- top[topline.born,3]
			lats[((when.born-first.step+1):length(lats[,1])),i] <- 1
		}
	}
	longs <- lats

	walk <- function(orilat,orilong,geodata,steps,ch) {
		m <- list()
		m[[1]] <- c(orilat,orilong)
		cat("\ndoing step\n")
		for (i in 1:steps) {
			cat(i," ")
			if ( round(i/10, 0)==i/10 ) {
				cat("\n")
			}
			m[[i+1]] <- new.coor(geodata,m[[i]],ch)
		}
		mt <- matrix(ncol=3,nrow=length(m))
		for (j in 1:length(m)) {
			mt[j,1] <- m[[j]][1]
			mt[j,2] <- m[[j]][2]
		}
		return(mt)
	}

	new.coor <- function(geodata,coor,ch) {
		add <- .01
		rg <- c()
		while ( length(rg) < 10000 ) {
			i1 <- which(geodata[,1] > coor[1] - add & geodata[,1] < coor[1] + add)
			i2 <- which(geodata[,2] > coor[2] - add & geodata[,2] < coor[2] + add)
			rg <- intersect(i1,i2)
			add <- add + .1
			if ( length(rg) > ch ) {
				break
			}
		}
		choice <- sample(rg, 1)
		nc <- c(geodata[choice,c(1:2)][1,1],geodata[choice,c(1:2)][1,2])
		return(nc)
	}

	w1 <- walk(orilat,orilong,geodata,length(lats[,1])-1,ch)
	lats[,1] <- w1[,1]
	longs[,1] <- w1[,2]
	w2 <- walk(lats[1,1],longs[1,1],geodata,length(lats[,1])-1,ch)
	lats[,2] <- w2[,1]
	longs[,2] <- w2[,2]
	rm(i,w1,w2)
	if ( length(lats[1,]) > 2 ) {
		for (i in 3:length(lats[1,])) {
			where_start <- min(as.vector(which(lats[,i]==1)))
			how_many <- length(lats[,1])-where_start
			mother <- as.character(top[which(top[,2]==as.numeric(colnames(lats)[i])),1])
			birthdate <- top[which(top[,2]==as.character(colnames(lats)[i])),3]
			newlat <- lats[as.character(birthdate),as.character(mother)]
			newlong <- longs[as.character(birthdate),as.character(mother)]
			w3 <- walk(newlat,newlong,geodata,how_many,ch)
			lats[where_start:length(lats[,1]),i] <- w3[,1]
			longs[where_start:length(longs[,1]),i] <- w3[,2]
		}
	}
	write.table(lats, file=paste(first,".lats",sep=""), sep="\t")
	write.table(longs, file=paste(first,".longs",sep=""), sep="\t")
}

