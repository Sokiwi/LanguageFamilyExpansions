library(OpenStreetMap)
library(rgeos)
library(maptools)
library(rJava)
library(argosfilter)

map <- function(first) {
	md <- mapfam(first)
	map.coordinates(md)
}

mapfam <- function(first) {
	lafile <- paste(first, ".lats", sep="")
	lofile <- paste(first, ".longs", sep="")
	alllat <- read.table(file=lafile)
	alllong <- read.table(file=lofile)
	orilat <- alllat[1,1]
	orilong <- alllong[1,1]
	termlat <- as.vector(as.matrix(alllat[length(alllat[,1]),]))
	termlong <- as.vector(as.matrix(alllong[length(alllong[,1]),]))
	alllat <- as.vector(as.matrix(alllat))
	alllong <- as.vector(as.matrix(alllong))
	if ( length(which(is.na(alllat))) > 0 ) {
		alllat <- alllat[-which(is.na(alllat))]
		alllong <- alllong[-which(is.na(alllong))]
	}
	md <- list(orilat, orilong, alllat, alllong, termlat, termlong, first) # md stands for migration data
	return(md)
}

map.coordinates <- function(md) {
	# get max min long and lat
	lat_range <- range(md[[3]])
	lng_range <- range(md[[4]])

	# add a frame of 10% around the points
	lat_extend <- 0.1 * diff(lat_range)
	lng_extend <- 0.1 * diff(lng_range)
	lat_range <- c(lat_range - lat_extend, lat_range + lat_extend)
	lng_range <- c(lng_range - lng_extend, lng_range + lng_extend)

	# get openstreetmap data - see ?openmap for "type" to get the list of map styles
	omap <- openmap(c(max(lat_range), min(lng_range)), c(min(lat_range), max(lng_range)), type="bing")

	# get map projection
	omap_proj <- omap$tiles[[1]]$projection

	# create a spatial object of read data
	data_sp <- SpatialPoints(cbind(md[[4]],md[[3]]))
	data_hl <- SpatialPoints(cbind(md[[4]][1],md[[3]][1]))
	data_term <- SpatialPoints(cbind(md[[6]],md[[5]]))
	proj4string(data_sp) <- CRS("+proj=longlat")
	proj4string(data_hl) <- CRS("+proj=longlat")
	proj4string(data_term) <- CRS("+proj=longlat")

	# transform the projection of data read to the one used by openstreet
	data_sp_proj <- spTransform(data_sp, omap_proj)
	data_hl_proj <- spTransform(data_hl, omap_proj)
	data_term_proj <- spTransform(data_term, omap_proj)

	# plot
	png(filename=paste(md[[7]],".png", sep=""))
	plot(omap)
	plot(data_sp_proj, pch=21, bg="white", add=T)
	plot(data_term_proj, pch=21, cex=1.8, bg="red", add=T)
	plot(data_hl_proj, pch=22, cex=1.8, bg="green", add=T)
	dev.off()
}
