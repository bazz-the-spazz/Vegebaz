

# Function to calculate weighted or mean indicator values
## data (abundance or presence/absence (0, 1)) has to be in the form: columns are species, rows are plots. Optionally: rownames are the names of the plots




get.indicator.value <- function(data, corrected.names, value, weighted="sqrt", source, na.rm=TRUE, method="mean", propose.alternatives=T, propose.alternatives.full=F, diversities=F, sozio=NULL, choose.names.sozz=TRUE){

	# Check if Latin in source
	if(!("Latin" %in% names(source))){
		if(("TaxonConcept" %in% names(source))){ source$Latin <- source$TaxonConcept} else{ # The dataframe of dengler has the variable TaxonConcep which is the Latin name of the species
			warning( paste('Latin is not in source!', sep=""), immediate. = F, call. = TRUE)
		}
	}


	# Check if value is in the source
	if(!(value %in% names(source))) warning( paste('"', value, '" is not in source!', sep=""), immediate. = F, call. = TRUE)

	# Check if value is a numeric
	if(!(is.numeric(source[,value]))) {
		cat(paste('"', value, '" is not numeric!\n', sep="") )
		do.calculations <- F} else do.calculations <- T

		# Manage Data
		# Use correct names
		if(!missing(corrected.names)){
			if(length(corrected.names) != ncol(data)){
				if(length(corrected.names) > ncol(data)){
					warning("To many species in corrected.names", call. = F)
				} else {
					warning("To few species in corrected.names", call. = F)
				}
			}
			names(data) <- corrected.names
		}

		# change NAs to zero and backup
		data[is.na(data)] <- 0
		data.bak <- data

		# subset the source data
		source.bak  <- X <- source
		rownames(source) <- source$Latin
		source <- source[names(data),value]  # values in correct order


		if(diversities) d.diversity.backup <- data # special backup for diversity

		# for propose.alternatives in case of missing values
		N <- names(data)
		if((propose.alternatives | propose.alternatives.full) & TRUE %in% is.na(source) ){

			nams <- namso <-  names(data)[is.na(source)]
			if(propose.alternatives.full==FALSE)  nams <- namso <- nams[nams %in% source.bak$Latin] # propose.alternatives will only check for species that previously have been corrected. propose.alternatives.full will search alternatives for all species without a value (also for 'cf's etc)

			if(length(nams)>0){

				nams <- paste(nams , " ") # add blank at the end
				nams <- gsub(" ", "  ", nams) # double all the blanks in the names
				nams <- gsub("cf\\.", "cf ", nams) # double all the blanks in the names
				nams <- gsub("sp\\.", "sp ", nams) # double all the blanks in the names
				for(i in c( "_", " cf ", " sp ", "   ", "  ", "  ")) nams <- gsub(i, " ", nams) # remove all the points, underscores, "cf", "sp", tripple and double blanks from the names
				nams <- substr(nams, start = 1, stop = nchar(nams)-1) # remove the blank from the end

				X$g <- sub(" .*", "", X$Latin) # get Genus of all Species
				X$s <- sub(" .*", "", substr(stop   = 100, start = nchar(X$g)+2, x = X$Latin))
				X$gs <- paste(X$g, X$s)



				for(i in 1:length(nams)){ #names(data)[is.na(source)]){ # loop for each species with with missing values


					#first get all entrys of the same species which have the data
					x <- X$gs[X$Latin==nams[i]]
					if(length(x)==0) {
						x <- X$gs[X$gs==nams[i]]
						if(length(x)==0) {
							x <- unique(X$gs[agrep(nams[i], X$gs)])
						}
					}

					alternative <- character()
					if(length(x)==1) alternative <- X[X$gs==x & !is.na(X[,value]), "Latin" ]
					if(length(x)>1) alternative <- X[X$g %in% sub(" .*", "", x)  & !is.na(X[,value]), "Latin" ]



					if(length(alternative)>0){
						I <-  0
						answer <- -1
						while(!(answer %in% 1:(length(alternative)+1) ) | answer==-1 ){
							if(I>0) cat(paste("\n Only NUMBERS between 1 and",length(alternative)+1,"allowed.\n"))
							cat(paste("\n",namso[i], ": No value '", value, "'. Choose other species:\n",
												paste(
													paste(1:length(alternative), ". ", alternative, " (", X[X$Latin %in% alternative, value], ")\n", sep="")
													, collapse = ""),
												paste(length(alternative)+1, ". keep '", namso[i],"'.\n", sep = "" )
												, sep=""))
							answer <- readline(prompt = paste("Choose number between 1 and ",length(alternative)+1," (or zero):", sep=""))
							I <- 1
							if(answer %in% c("zero", 0) ) answer <- length(alternative)+1
						}
						if(as.numeric(answer) <= length(alternative)) {
							source[names(data)==namso[i]] <-  X[X$Latin==alternative[as.numeric(answer)], value]
							N[N==namso[i]] <- alternative[as.numeric(answer)]
						}
					}

				}
			}
		}



		r2 <- data.frame(species=names(data), value=source)
		if(!identical(names(data), N)) r2$used.species <- N # if you used alternatives, report them

		## Function to calculate diversities
		diversitee <- function(x, q, margin=1, effective=TRUE, na.rm=TRUE, se=F, allow.incorrect.beta=FALSE){
			# based on Jost 2006 & 2007 calculate the (effective) diversity of a community. Loosely adapted from the diversity function in vegan.

			# Organize data
			x <- drop(as.matrix(x))
			if (!is.numeric(x)) stop("input data must be numeric")
			if (any(x < 0, na.rm = TRUE)) stop("input data must be non-negative")

			# Function to Calculate Effective Diversity
			f <- function(x, Q, na=na.rm, raw=FALSE){

				# Get Proportional Data
				x <- x[x>0] # Zeros are ignored
				if(na) x <- x[!is.na(x)]
				p <- x/sum(x)

				if(NA %in% p){
					data <- NA
				} else {
					if(Q==1){
						data <- exp(-sum(p*log(p))) #q=1 is not definded, calculate Shannon
						if(raw) data <- -sum(p*log(p)) # needed for Alpha
					} else {
						data <- sum(p^Q)^(1/(1-Q)) # Else use Jost's universal formula
					}
				}
				return(data)
			}


			# Apply the before defined function f to the data(frame)
			if(length(dim(x))>1){
				D <- apply(x, MARGIN = margin, FUN = function(x)f(x, Q=q) )
			} else {
				D <- f(x,Q=q)
			}


			# From Jost 2007 :PARTITIONING DIVERSITY INTO INDEPENDENT ALPHA AND BETA COMPONENTS
			# True Alpha diversity
			if(length(dim(x))>1){
				## First Calculate weights
				w <- apply(x, MARGIN = margin, FUN = function(x)f(x, Q = 0))
				w <- w/sum(w) # w is the statistical weight of Community j (usually the number of individuals in Community j divided by the total number of individuals in the region)

				if(q!=1){
					if(length(unique(w))==1 | allow.incorrect.beta){ # if weights w are not equal only Shannon might be used!
						sums <- D^(1/(1/(1-q)))
						A <- ( sum((w^q * sums))/ sum(w^q) )^(1/(1-q)) # eq 11a
					} else {
						A <- NA
					}
				}
				if(q==1){ # use formula 11b for q=1 (Shannon)
					lamdas <- apply(x, MARGIN = margin, FUN = function(x)f(x, Q=q, raw=TRUE) )
					A <- exp(sum(w*lamdas))
				}

				# Gamma? the numbers equivalent of the diversity index of the pooled samples
				G <- f( t(colSums(x/rowSums(d, na.rm = T), na.rm = T))/nrow(x), Q=q)



				# Beta (Following from eq9)
				B <- G/A



				# Standard Error

				if(se & (q==1 | allow.incorrect.beta)){
					se.a <- se.b <- se.g <- numeric()
					for(i in 1:(nrow(x)-1)){
						for(j in (i+1):nrow(x)){
							xx <- x[c(i,j),]
							w <- apply(xx, MARGIN = margin, FUN = function(x)f(x, Q = 0))
							w <- w/sum(w)
							se.a <-  c(se.a, exp(sum(w*apply(xx, MARGIN = margin, FUN = function(x)f(x, Q=q, raw=TRUE) ))))
							se.g <- c(se.g, f( t(apply(xx, MARGIN = ifelse(margin==1,2,1) , sum, na.rm=T)) ,Q=q))
							se.b <- c(se.b, se.g[length(se.g)]/se.a[length(se.a)])
						}
					}
					se.a <- sd(se.a)/sqrt(ifelse(margin==1, nrow(x), ncol(x)))
					se.b <- sd(se.b)/sqrt(ifelse(margin==1, nrow(x), ncol(x)))
					se.g <- sd(se.g)/sqrt(ifelse(margin==1, nrow(x), ncol(x)))
				}



			}


			# Transform to non-effective Indices (if needed)
			if(effective==FALSE){
				if(q==0){
					D <- D
					I <- "Species Richness"
				}
				if(q==1){
					D <- log(D)
					I <- "Shannon entropy"
				}
				if(q==2){
					D <- 1-(1/D)
					I <- "Gini-Simpson index"
				}
				if(!(q %in% 0:2)) {
					Value <- NA
					I <- NA
				}
				return(list(Value=D, index=I, gamma=G, alpha=A, beta=B))
			} else {
				if(se){
					return(list(D=D, gamma=G, alpha=A, beta=B, se.g=se.g, se.a=se.a, se.b=se.b))
				} else {
					return(list(D=D, gamma=G, alpha=A, beta=B))
				}
			}
		}

	if(is.data.frame(sozio)){
		sozfun <- function(d, sozz, newname=choose.names.sozz){ # The sozfun-function adds up the values of but doubles them when the species is more than 10% in relative abundance
			sp <- names(d)[names(d) %in% rownames(sozio)]
			cat(paste(c("\n These species are in the Sozio data\n\n", sp, "\n\n"), collapse = " : "))
			
			if(newname){
				nsp <- names(d)[!(names(d) %in% rownames(sozio))]
				tmp <- sozz
				tmp$Latin <- rownames(sozz)
				x  <- choose.name(names = nsp, data = tmp)
				translation <- data.frame(original_name=nsp, new_name=x)
				names(d)[!(names(d) %in% rownames(sozio))] <- x
				sp <- names(d)[names(d) %in% rownames(sozio)]
				cat(paste(c("\n These species are now in the Sozio data\n\n", sp, "\n\n"), collapse = " : "))
			} else translation <- "no additional Translation"
			
			
			d <- d/rowSums(d, na.rm = T) # get relative abundance
			l <- list()
			for(i in rownames(d)){ # loop for all plots
				l[[i]] <- sozio[0,]
				if(sum(d[i,sp], na.rm=T)>0){
					for(j in sp){ # loop for all species in plot
						if(d[i,j]>0.1) {   # is species more than 10%
							l[[i]] <- rbind(l[[i]], sozio[j,]*2) 
						} else {
							if(d[i,j]>0) l[[i]] <- rbind(l[[i]], sozio[j,]) 
						}
					}
					x <- data.frame(TypoCH = names(colSums(l[[i]], na.rm = T)), score=(colSums(l[[i]], na.rm = T)))
					x <- x[order(x$score, decreasing = T),]
					x <- x[x$score>0,]
					rownames(x) <- NULL
					l[[i]] <- x } else 					l[[i]] <- "no species in indicator list"
			}
			return(list(sozz=l, translation=translation))
		}
		sozR <- sozfun(d = d, sozz = sozio)
		# sozR
	}
	


		if(do.calculations){
			if(weighted=="normal" | weighted=="sqrt") {   # when the weighted values are needed, make that plots add up to 1. Else make the data presence absence.
				data.w <- data/rowSums(data, na.rm = T) # wheighted
				data.w2 <- data/rowSums(data[,which(!is.na(source))], na.rm = T)  # wheighted but ignore plants without Indicator value
				data.w2[, is.na(source)] <- NA
				data <- data.w2
				if(weighted=="sqrt") data <- sqrt(data.w2)
				weighted <- TRUE
				sqrtttt <- TRUE
			}  else {
				data[data>0 & !is.na(data)] <- 1
				data[data==0 | is.na(data)] <- NA
				weighted <- FALSE
			}
			R <- as.numeric()
			D <- as.character()
			for(i in 1:nrow(data)){ # loop for each plot
				if(method=="mean"){
					if(weighted) R[i] <- sum(t(data[i,])*source,na.rm=na.rm) else R[i] <- mean(t(data[i,])*source,na.rm=na.rm) # when weighted sum up the weighted components Else take the average of the occuring species
				}
				if(method=="sd"){ # if standard deviation is chosen
					if(weighted) R[i] <- sd(t(data[i,])*source,na.rm=na.rm) else R[i] <- sd(t(data[i,])*source,na.rm=na.rm)
				}
			}

			names(R) <-  rownames(data)


			Return <- (list(value=paste(ifelse(weighted, ifelse(sqrtttt, "sqrt weighted", "weighted"), ""), method, value), plots=R, species=r2))


		} else   Return <- (list(value=value, species=r2))

		if(diversities){
			sr <- diversitee(x = d.diversity.backup, q = 0)$D
			Return <- append(Return, list("Species richness"=sr))
			if(weighted){
				esh <- diversitee(x = d.diversity.backup, q = 1)
				Return <- append(Return, list("effective Shannon diversity"=esh$D, "Gamma diversity"=esh$gamma, "Beta diversity"=esh$beta))
			}
		}

		if(is.data.frame(sozio)){ Return <- append(Return, sozR) }

		return(Return)
}




















