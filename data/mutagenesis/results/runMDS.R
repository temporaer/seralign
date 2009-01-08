library(vegan)
path = "./"
files=dir( path=path, pattern=".*csv" )
act<<- read.delim("../188/is_active.txt", header=FALSE, sep=" ")
for (f in files){
	csvfile = paste(path, f, sep="/")
	simmat <<- read.csv (file=csvfile,  na.strings = "NA", nrows = -1, skip = 0, check.names = TRUE, strip.white = FALSE, blank.lines.skip = TRUE)
	nam = names(simmat)[-1] 

	z <- c(); for ( x in 1:length(nam) ) z<-c(z, if(act[["V2"]][which(act[["V1"]]==nam[x])]) "X" else "O")

	m <- as.matrix(simmat[nam])
	m<-matrix(as.numeric(m), ncol=length(nam))
#	m<-exp(-m)
	m[is.na(m)] <- 0
	w <- metaMDS(m,zerodist="add")
	png(paste(csvfile, "png", sep="."))
	plot(w$points, type="n", main=paste("2D embedding of distances:",csvfile), xlab="axis 1", ylab="axis 2")
	points(w$points[which(z=="X"),],col="blue")
	points(w$points[which(z=="O"),],col="red")
	dev.off()
}
