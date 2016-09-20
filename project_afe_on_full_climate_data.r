# read settings from shell
args <- commandArgs()
n <- args[length(args)]
index <- as.numeric(n)


#.libPaths(c('/home/pearman/R_library2','/opt/cluster/WSLSoftware/R-2.12.2/lib64/R/library'))
library(PresenceAbsence)
library(shapefiles)
library(maptools)
library(foreign)
#gpclibPermit()
library(gbm)
library(gam)
library(raster)
library(dismo)

run <- "models_lin_weights"
#  NOTE:  You have to know how many species there are, so that you can tell the cluster. Use the check
#         code immediately below

# function to make a brick out of layers in a netcdf file
ncdf2brick<-function(ncdf,varnames){
   tmp<-brick()
   for(i in varnames){
      tmp<-addLayer(tmp,raster(ncdf,varname=i))
   }
   tmp
}




species.list <- read.csv("/home/pearman/lud11_docs/wsl_research/projects/OcCC/r_code/plants/afe/models_lin_weights/auc_summary_across_scales.csv")[-1]

species.list <- species.list[which(species.list$auc.1km >= 0.70),]

species <- as.character(species.list$afe.code[index])
k <- species

models <- c("gam","glm","gbm","me")
para.sets <- c("vec1","vec2","vec3","vec4","vec5")  # adjust this


# collect auc values and weights

auc.weights <-as.data.frame(matrix(data=NA,nrow=20,ncol=5))
names(auc.weights) <- c("model","para.set","auc","weight","max.tss")
class(auc.weights$auc) <- "numeric"


setwd(paste("/home/pearman/lud11_docs/wsl_research/projects/OcCC/r_code/plants/afe/",run,"/",species,"/",species,sep=""))
mat.row <- 0
for(i in 1:length(models)){
  stat.mat <- read.csv(paste("/home/pearman/lud11_docs/wsl_research/projects/OcCC/r_code/plants/afe/",run,"/",species,"/",species,"/summary_stat_models/sum.stat.",models[i],"_",species,".csv",sep=""))[-1]
  for(j in 1:length(para.sets)){
    mat.row <- mat.row+1
    auc.weights[mat.row,1] <- models[i]
    auc.weights[mat.row,2] <- para.sets[j]
    auc.weights[mat.row,3] <- stat.mat[j,1]
    auc.weights[mat.row,4] <- NA
    y=load(paste(species,"_",models[i],"_",para.sets[j],sep=""))
    fname <- ifelse(models[i]=="me",paste(models[i],".max.tss",sep=""),paste(models[i],".tss.max",sep=""))
    eval(parse(text=paste("auc.weights[mat.row,5] <- ",fname,sep="")))
    rm(list=y)

  }
}
auc.mod <- auc.weights$auc - rep(0.5,dim(auc.weights)[1])
auc.weights$weight <- auc.mod/sum(auc.mod)
    
data.set.a <- c('world_bio_ch_1km_aea','dmi_1130','dmi_4160','dmi_8100','en_ccsm3_1130','en_ccsm3_4160','en_ccsm3_8100','en_echam5_1130','en_echam5_4160','en_echam5_8100',
                                       'knmi_1130','knmi_4160','knmi_8100','meto_1130','meto_4160','meto_8100','mpi_1130','mpi_4160','mpi_8100')

data.set <- paste(data.set.a,".nc",sep="")
###################


#auc.weights <-as.data.frame(matrix(data=NA,nrow=20,ncol=5))
#names(auc.weights) <- c("model","para.set","auc","weight","max.tss")
#class(auc.weights$auc) <- "numeric"



wkdir <- paste("/home/pearman/lud11_docs/wsl_research/projects/OcCC/r_code/plants/afe/",run,"/",species,sep="")
setwd(wkdir)
system("mkdir -p ch_projections")

z <- 0
for (m in data.set) {
  z <- z+1

  #prepare predictor variables and species vector, specify a vector of datasets that have current and future climate values
  setwd("/home/pearman/lud11_docs/wsl_research/projects/OcCC/data/climate/ch2014")


  varnames <- paste("bio_",1:19,sep="")
  varnames2 <- paste("bio",1:19,sep="")

  # load the climate data set here and make a raster brick
  data.brick <- ncdf2brick(m,varnames)

  projection(data.brick) <- CRS("+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs")
  names(data.brick) <- varnames2
  predict.data <- rasterToPoints(data.brick,spatial=TRUE)

 

  setwd(paste("/home/pearman/lud11_docs/wsl_research/projects/OcCC/r_code/plants/afe/",run,"/",k,sep=""))

  for(i in models){
    
    for(j in para.sets){
      
      print(paste("now i = ",i," and j = ",j," and k = ",k,sep=""))
      y <- eval(parse(text=paste("load('./",k,"/",k,"_",i,"_",j,"',.GlobalEnv)",sep="")))
      if((i=="gam")==TRUE) mod.obj <- gam.step
      if((i=="glm")==TRUE) mod.obj <- glm.step
      if((i=="gbm")==TRUE) mod.obj <- gbm.model
      if((i=="me")==TRUE)  mod.obj <- me

      if((i=="gbm")==TRUE){
        pred <- predict.gbm(mod.obj,predict.data,n.trees=best.itr1,type="response")
        eval(parse(text=paste("predict.data@data$",i,j," <- pred",sep="")))
      }

      if((i=="gam")==TRUE){
        pred <- predict.gam(mod.obj,predict.data,type="response")
        eval(parse(text=paste("predict.data@data$",i,j," <- pred",sep="")))
      }

      if((i=="glm")==TRUE){
        pred <- predict(mod.obj,predict.data,type="response")
        eval(parse(text=paste("predict.data@data$",i,j," <- pred",sep="")))
      }

      if((i=="me")==TRUE){
        system(paste("mkdir -p ",species,"/temp.path",sep=""))
        mod.obj@path <- paste(species,"/temp.path",sep="")
        pred <- dismo::predict(mod.obj,predict.data,type="response",na.rm=FALSE)
        eval(parse(text=paste("predict.data@data$",i,j," <- pred",sep="")))
      }

      rm(list=y) # remove the model that was used
      
    } # end loop for para.sets
  }  #end loop for models
  
  print(paste("now making weighted prediction for k = ",k,sep=""))
  predict.data@data$mean <- as.matrix(predict.data@data[,20:39])%*%auc.weights$weight
  
  eval(parse(text=paste("tss.data <- read.csv('/home/pearman/lud11_docs/wsl_research/projects/OcCC/r_code/plants/afe/",run,"/",species,"/validation/",species,"_validation.csv')[-1]",sep="")))
  tssmax <- tss.data[4,3]
  predict.data@data$bi.pred <- as.numeric(predict.data@data$mean >= tssmax)

  wkdir <- paste("/home/pearman/lud11_docs/wsl_research/projects/OcCC/r_code/plants/afe/",run,"/",species,"/ch_projections",sep="")
  setwd(wkdir)

  # name the saved object which has the projection for rcm and time slice
   file.name <- paste("en_results_",species,"_",data.set.a[z],sep="")  
  save(predict.data,file=file.name)

  # save a tiff of the projection as well, for the rcm and time slice
  file.name <- paste("pred_",species,"_",data.set.a[z],".tif",sep="")
  dummy <- data.brick[[1]]
  prediction <- rasterize(predict.data,y=dummy,field='bi.pred')
  writeRaster(prediction,filename=file.name,format="GTiff",overwrite=TRUE)
  
  
}  # end loop for datasets

wkdir <- paste("/home/pearman/lud11_docs/wsl_research/projects/OcCC/r_code/plants/afe/",run,sep="")
setwd(wkdir)

cat(paste(species," ran successfully\n",sep=""),file="logfile.log",append=TRUE)

########################################################




