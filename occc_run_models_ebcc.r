
# read settings from shell
args <- commandArgs()
n <- args[length(args)]
index <- as.numeric(n)

set.seed(index+as.numeric(proc.time())[3])
seed <- .Random.seed
# run the models

#set to show whether atlas data or pseudoabsence data are used
atlas <- TRUE

setwd("/home/pearman/lud11_docs/wsl_research/projects/OcCC/r_code")

source("occc_run4models_function.r")

# set working directory to directory with species data structures

setwd("/home/pearman/lud11_docs/wsl_research/projects/OcCC/r_code/birds/ebcc")

# load the species and climate data

load("EBCC_25OBS")
data.a <- ebcc.25obs[-1]
load("EBCC_CLIMATE")
ebcc.climate <- ebcc.climate[-1]

#for (i in 12:19){
#  eval(parse(text=paste("ebcc.climate$bio",i," <- ebcc.climate$bio",i,"/10",sep="")))
#}

#load weights for climate polygons, must be vector on 1's if not used
#and the length of this vector must be as long as dim(data.a)[1]


###################  specify parameter sets here#################################
vec1 <- c("bio3","bio5","bio15","bio18")   #ok
vec2 <- c("bio1","bio2","bio8","bio17")    #ok 
vec3 <- c("bio2","bio9","bio18","bio19")   #ok 
vec4 <- c("bio5","bio6","bio14","bio15")   #ok
vec5 <- c("bio2","bio8","bio14","bio19")   #ok

# here you specify the species names as a single value or a vector.
# resdir, below, could potentially be a vector, for example if models were to built with two different datasets
# otherwise, each of the species directories will be set up in the results directory
species <- names(data.a)[index]
resdir <- species

# assign climate to variables

data1 <- cbind(data.a[species],ebcc.climate)


############################ SET WEIGHTS HERE ################################


#run <- "models_no_weights"
#climate.weights <- rep(1,dim(data1)[1])

load("/home/pearman/lud11_docs/wsl_research/projects/OcCC/data/weights_eoa/CLIMATE_WEIGHTS")

#climate.weights <- climate.weights$linear_wei
#run <- "models_lin_weights"

climate.weights <- climate.weights$log_weight
run <- "models_log_weights"


##  if there are more than 2000 zeros, take a sample of 2000 of the absence data and set atlas to FALSE
if ((length(which(data1[species]==0))>2000)==TRUE){
  atlas=FALSE
  rows.0 <- sample(which(data1[species]==0),2000,replace=FALSE)
  rows.1 <- which(data1[species]==1)
  data1 <- data1[c(rows.1,rows.0),]
  climate.weights <- climate.weights[c(rows.1,rows.0)]
}


# calculation of number of occurrences and the prevalence of each lineage or species
occurrences <- sum(data1[species])


# make sure that these next two items have the same structure
vars <- list(vec1,vec2,vec3,vec4,vec5)
para.names <- c("vec1","vec2","vec3","vec4","vec5")

# here is where the data set is specified
data <- data1

# specify the directory into which all the species directories should go
setwd(paste("/home/pearman/lud11_docs/wsl_research/projects/OcCC/r_code/birds/ebcc/",run,sep=""))


# these directories are necessary for smooth running of maxent
# full path names are needed for ME to work on Mac.  Maybe on windoze too.
me.dir1 <- paste(getwd(),"/",resdir,"/me_main",sep="")
me.dir2 <- paste(getwd(),"/",resdir,"/me_cv",sep="")

                       
                       
                       
# if there are a bunch of species, here you should set up a loop to run through the vector of species (or lineage) names
run4models(data=data,
                     species=species,
                     paraset=vars,
                     para.names=para.names,
                     models=c('me','glm','gam','gbm'),
                     atlas=atlas,
                     climate.weights=climate.weights,
                     k=10,
                     resdir=resdir,
                     me.dir1=me.dir1,
                     me.dir2=me.dir2,
                     gbm.args=list(n.trees=10000,shrinkage.bound=c(0.01,0.0005),interaction.depth=3,maxiter=15),
                     maxent.args=c('-J','-P','threshold=false','product=false','quadratic=false','linear=false'),
                     seed=seed,
                     plot=FALSE,
                     help=FALSE)


cat(paste(species," ran successfully\n",sep=""),file="calibration_logfile.log",append=TRUE)

