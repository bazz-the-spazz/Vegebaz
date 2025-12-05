
# # Example

# 1. Read data from either Vegedaz, Vegedaz 2024, or Veg.
# Run the scripts to get the functions and the source data:
source("Calculate_indicator_values.r")
source("Correct_species_names.r")
source("read_data_from_Vegedaz.r")
source("read_data_from_Vegedaz_2024.r")
source("read_data_from_Veg_2015.r")
source("transpose_data.r")


# 1.1
# read data from Vegedaz 2021 (Flora Indicativa, Landolt)
# path=("~/.wine/drive_c/Program Files/Vegedaz/Daten/") # on Linux
# path=("C:/Program Files/Vegedaz/Daten/") # on Windows
path=("./Daten/") # or wherever Vegedaz is located on Windows
vegedaz21 <- read.vegedaz.data(path)


# 1.2 get data from Vegedaz Beta Dez 24
# read data from Vegedaz 2024 (Landolt and Eunis(Dengler et al 2023))
path=("./DB/") # The path should lead to the 'DB' folder in the program folder
vegedaz24 <- read.vegedaz24(path = path)

# 1.3
# get the data from Dengler etal 2023 (Ecological Indicator Values for Europe (EIVE) 1.0) https://vcs.pensoft.net/article/98324/
# 	supplementary material 8
require(openxlsx)
# dengler <- read.xlsx("vegetation_classification_and_survey-004-007-g008.xlsx", sheet=2) # From locally saved copy1
dengler <- read.xlsx("https://vcs.pensoft.net/article/98324/download/suppl/38/", sheet = 2) # Directly from the Net
dengler$Latin <- dengler$TaxonConcept

# .4
# for Veg (Flora Helvetica 2014), extract the Zip file and put the file "Zeigerliste.txt" into your work directory.
floraH <- getVegData(path="")

# 1.5
## choose your data source to proceed
source <- dengler   # the options are 'vegedaz21$landolt', 'vegedaz21$indicativa' from Vegedaz21, 'vegedaz24$landolt' or 'vegedaz24$eunis, dengler, or 'floraH' from Veg


# 1.6 Get Soziodata from https://www.envidat.ch/#/metadata/modified-typoch
sozio <- read.csv("https://www.envidat.ch/dataset/bcf50c9d-3b31-4aa7-a774-57a8f6f60ed7/resource/47adb4a8-f089-4622-b525-6dff0da104d9/download/typoch_modified-list.csv", header = F, row.names = NULL, sep = ";")

for(i in 1:ncol(sozio)) names(sozio)[i] <- paste(sozio[1,i], sozio[2,i]) # name the data
for(i in 2:ncol(sozio)) sozio[,i] <- as.numeric(sozio[,i]) # make numeric
sozio <- sozio[!(sozio[,1] %in% c("", "TypoCode", "Vegetation type ", "Structural variables", "Bryophytes", "Vascular plants")),]

sozio[,1][duplicated(sozio[,1])] # there are duplicated species with different scores
sozio[sozio[,1] %in% sozio[,1][duplicated(sozio[,1])],]
sozio <- sozio[duplicated(sozio[,1])==FALSE,] #delete them
rownames(sozio) <- sozio[,1]
sozio[,1] <- NULL



# 2.
# for the exercise create a random dataframe (with species as rows and plots as columns)

# 2.1
species <- c("Daucus carota",  "Lathyrus pratensis", "Scorzoneroides autumnalis", "Aegopodium podagraria cf", "Heracleum sphondyllum sp", "Erigeron annuus")
d <- data.frame(species=species, plotA= runif(length(species)), plotB= runif(length(species)), plotC= runif(length(species)))
d

# 2.2 transpose the data.frame for the analyses (species as columns and plots as rows)
d <- transpose.data(data = d)
rownames(d) <- d[,1]  # plots as rownames
d <- d[,-1]						# Remove plot column
d
# introduce zeros
for( i in 1:3) d[sample(nrow(d), 1), sample(ncol(d), 1)] <- 0

# 3.
# Use the choose.name()-function to correct the species names according to the chosen source
corrected.names <- choose.name(names = names(d), data = source, write.tmp.file = T)
corrected.names

## if you're tired of choosing species names you can type 'pause' and later resume the task with:
# choose.name(names = names(d), data = source, continue.after.pause = corrected.names)


# 4.
# Use the get.indicator.value()-function to calculate e.g. mean Temperaturzahl for the plots
indi <- get.indicator.value(
	data=d,
	corrected.names = corrected.names,
	value = "EIVEres-T",
	weighted = "sqrt",
	source = source,
	na.rm = T,
	propose.alternatives = T,
	diversities = T,
	sozio = sozio

)



indi

# 'data=' is your data.frame, column names are the corrected species names and the columns only contain numbers.
# 'corrected.names=' the string of names received from the choose.name-function in the same order as the columns of data. (optional)
# 'value=' is the indicator value you are interested in. Need not to be numeric. There are many options. Check out names(vegedaz$indicativa).
# 'weighted=' should the indicator value be weighted by the abundance of the species? Options: "sqrt", "normal", "none"
# 'source='  is the data that contains the indicator value for each species.
# 'na.rm=' should missing indicator values be ignored?
# 'method=' do you want the average value ("mean") or the standard deviation ("sd")?
# 'propose.alternatives=TRUE' sometimes a indicator value is missing but another subspecies might have it.
# 'propose.alternatives.full=TRUE' similar to the choose.name-function, the algorhythm is looking for similar sounding species.
# 'sozio=' soziodata to get likely TypoCH community.
# 'diversities = TRUE' Species richness, effective Shannon diversity, and Beta-Diversity are calculated.

