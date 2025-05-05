library(GEOquery)
if (file.exists("data/GSE26117-columns.csv")) {
	coldata <- read.csv("data/GSE26117-columns.csv")
}
else {
	gds <- getGEO("GSE26117", destdir="data/GEO")

	gse <- gds[[1L]]
	coldata <- pData(gse)
}





