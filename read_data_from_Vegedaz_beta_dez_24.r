## Get data from vegedaz23 (Beta-Version vom Dez 2023)
# path should lead to the DB folder in the VegedazQtrelease directory

read.vegedaz23.data <- function(path){

	artlist <- read.delim(paste(path,"SpeciesListCH.txt", sep=""), sep = "@", header = F)
	## variable names are specified in the first rows
	n <- do.call(rbind,strsplit(artlist[2:10,1], "\t"))[,2]
	artlist <- artlist[-1:-10,]
	## recreate data.frame
	artlist <- strsplit(artlist, "\t")
	artlist <- as.data.frame(do.call(rbind,artlist))

	names(artlist) <- n



	# Prepare Flora Indicativa
	zeilan <- read.table(paste(path,"IndicatorValuesCH.txt", sep=""), sep="\t", header = T)

	# prepare local names
	localnames <- read.delim(paste(path,"LocalNamesCH.txt", sep=""), sep = "\t", header = F)
	names(localnames) <- t(localnames[1,])
	localnames <- localnames[-1,-5]

		# Prepare RedList information
	redlist <- read.delim(paste(path,"RedListCH.txt", sep=""), sep = "\t", header = F)
	n <- t(redlist[2,])
	redlist <- as.data.frame(redlist[-1:-2,])
	names(redlist) <- n
	redlist[,ncol(redlist)] <- NULL


	# prepare names Euromed
	euromed <- read.delim(paste(path,"SpeciesListEuromed.txt", sep=""), sep = "@", header = F)
	n <- do.call(rbind,strsplit(euromed[2:10,1], "\t"))[,2]
	euromed <- euromed[-1:-10,]
	## recreate data.frame
	euromed <- strsplit(euromed, "\t")
	euromed <- as.data.frame(do.call(rbind,euromed))
	names(euromed) <- n

	# Prepare Eunis
	zeieu <- read.table(paste(path,"IndicatorValuesEunis.txt", sep=""), sep="\t", header = T)




	# merge Landolt
	## artlist and zeilan
	zeilan <-    merge(artlist, zeilan, by.x = "Taxon_plain", by.y = "Name_plain", all=T)
	## plus redlist
	zeilan <-    merge(zeilan, redlist, by.x = "Taxon_plain", by.y = "Name", all=T)

	## remove Taxa that have neither redlist information nor indicator values
	test <- matrix(1, nrow=nrow(zeilan), ncol=ncol(zeilan))
	test[is.na(zeilan)] <- 0

	## add plus local names
	zeilan <- merge(zeilan, localnames, by.x = "Taxon_plain", by.y = "Name", all.x = TRUE, all.y = FALSE)
	zeilan$Latin <- zeilan$Taxon_plain


	# merge Eunis
	zeieu <- merge(euromed, zeieu, by.x="Taxon_plain", by.y = "TaxonConcept", all=T)
	zeieu$Latin <- zeieu$Taxon_plain

	return(list(landolt = zeilan, eunis = zeieu))
	# rm(artlist, ind, con, i, l, n)
}


# d <- read.vegedaz23.data(path)
