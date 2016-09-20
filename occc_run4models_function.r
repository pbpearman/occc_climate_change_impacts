###############################################################################
#                                                                             #
#=============================================================================#
#           FUNCTION TO RUN FOUR MODELS (GBM, GAM, GLM, MAXENT)               #
#               ON A NUMBER OF SPECIES AND PARAMETER SETS                     #
#=============================================================================#
#                                                                             #
#                                                                             #
#   WRITTEN BY RAFAEL W†EST, DEC. 2011, WSL BIRMENSDORF, SWITZERLAND          #
#   rafael.wueest@wsl.ch                                                      #
#                                                                             #
#   original code was provided by Peter B. Pearman, WSL Birmensdorf and       #
#   Manuela D'Amen, Rom, Italy                                                #
#   and transformed into a function by Rafael WŸest                           #
#                                                                             #
###############################################################################

run4models<-function(data=data.frame(),
                     species=character(),
                     paraset=list(),
                     para.names=NULL,
                     models=c('gbm','me','gam','glm'),
                     atlas=TRUE,
                     climate.weights=NULL,
                     k=5,
                     resdir='results',
                     me.dir1 = NULL,
                     me.dir2 = NULL,
                     gbm.args=list(n.trees=10000,shrinkage.bound=c(0.0005,0.003),interaction.depth=3,maxiter=5),
                     maxent.args=c('-J','-P','threshold=false','product=false','quadratic=false','linear=false'),
                     seed=NULL,
                     plot=FALSE,
                     help=FALSE){
   
   ### helptext
   helptext='
   ----------------------------------------------------------------------------------
   This function runs a set of models using a set of parameterizations and requires
   the following input data...
   
   data: data.frame with both the environmental and species occurrence information
   
   species: Character vector with names of species to be modeled.
   The species names must match a column name in data!
   
   paraset: a vector or list, each element of which gives the names of variables
   to be used. The variable names must match a column name in data!
   
   models: a character vector specifying the models to be run. Model names must
   match any (or several) of \'gbm\', \'me\', \'gam\', \'glm\'. Defaults to all.
   
   k: indicating how many partitions will be used in k-fold cross validation,
   defaults to 5.
   
   resdir: directory of where to store the results (character vector),
   default is \'results\'. Within resdir, folders will be set up for each species
   and one directory for the summary tables.
   
   The function also requires parameter settings for all four models!
   They have to be provided in list format, the defaults are as described below.
   
   gbm.args: list of gbm parameters of the form
   list(n.trees=10000,shrinkage.bound=c(0.0005,0.003),interaction.depth=3,maxiter=5)
   
   maxent.args: character vector of parameters settings for maxent (exluding
   prevalence). Defaults to
   c(\'-J\',\'-P\',\'threshold=false\',\'product=false\',\'quadratic=false\',\'linear=false\')
   
   seed: this will set the seed to make results reproduceable. Default is NULL.
   
   plot: default is FALSE, will plot gbm results if TRUE.
   ----------------------------------------------------------------------------------\n'
   #if(help) {cat(helptext);stop('This is not an error, but just a description of the function\'s parameters!')}
   
   
   ### check for results directory
 #  if(resdir %in% list.files()){
 #     cat(paste('The directrory',resdir,'for your results already exists. Do you want to continue (yes|no)?'))
 #     answer<-scan(what='character',n=1)
 #     if(!answer %in% c('y','Yes','YES','yes','Y')) stop('Please rename your results directory!')
 #  } else{
      system(paste('mkdir -p',resdir))
 #  }
   system(paste('mkdir -p',me.dir1,sep=" "))
   system(paste('mkdir -p',me.dir2,sep=" "))
   ### libraries needed
   library(foreign)
   library(PresenceAbsence)
   library(dismo)
   library(gam)
   library(gbm)
   .Random.seed <- seed
   
   ### tss function to be compiled
   tss <- function(observed,prediction){
   levs <- seq(from=0.001,to=0.999,by=0.001)
   obsnum <- paste("obs_",1:length(observed),sep="")
   data1 <- data.frame(obsnum,observed,prediction)
   val <- numeric()
   sens <- numeric()
   spec <- numeric()
   thresh <- NA
   for(i in 1:length(levs)){
      confmat <- cmx(data1,threshold=levs[i])
      sens <- c(sens,as.numeric(sensitivity(confmat)[[1]]))
      spec <- c(spec,as.numeric(specificity(confmat)[[1]]))
   }
   val <- sens + spec -rep(1,length(levs))
   to.average <- which(val==max(val))
   thresh <- mean(levs[to.average])
   return(thresh)
   }
   
   
   ### define elements of the function
   ## seed if needed
   if(!is.null(seed)) set.seed(seed)
   ## create list if paraset is vector
#   if(!is.list(paraset)) {
#      paraset<-list(paraset);cat('The function assumes you only speciefed one set of parameters!\nPlease confirm with \'yes\' if this is correct.\n')
#      answer<-scan(what='character',n=1); if(!answer %in% c('yes','Yes','YES','y','Y')) stop('Please supply a valid parameter set as list or vector.')
#   }
   ## setting some parameters
   # order of the data
   order<-row.names(data)
   # occurrences/prevalence
   occurrences<-apply(data.frame(data[,species]),2,sum)
   names(occurrences)<-species
   prevalences<-occurrences/nrow(data)
   absences<-nrow(data)-occurrences
   # naming vectors
   names.cv.tab<-c(paste('rep',1:k,sep=''))
   names.eval.stat <- c('auc', 'pcc', 'sensitivity', 'specificity')
   
   
##############
# species loop
##############
   
   for(spec in species){
      cat(paste('\n\n\n#####################\n### ',spec,'\n#####################',sep=''))
      
      ### create directory for results
      direc <- paste(resdir,"/",spec,sep="")
      system(paste('mkdir -p ',direc,sep=''))
      
      ### create summary tables
      ## gbm
      summary.sp.stat.gbm<-data.frame(matrix(NA,nrow=length(paraset),ncol=4))
      dimnames(summary.sp.stat.gbm)<-list(paste('para',1:length(paraset),sep=''),names.eval.stat)
      ## gam
      summary.sp.stat.gam<-data.frame(matrix(NA,nrow=length(paraset),ncol=4))
      dimnames(summary.sp.stat.gam)<-list(paste('para',1:length(paraset),sep=''),names.eval.stat)
      ## glm
      summary.sp.stat.glm<-data.frame(matrix(NA,nrow=length(paraset),ncol=4))
      dimnames(summary.sp.stat.glm)<-list(paste('para',1:length(paraset),sep=''),names.eval.stat)
      ## me
      summary.sp.stat.me<-data.frame(matrix(NA,nrow=length(paraset),ncol=4))
      dimnames(summary.sp.stat.me)<-list(paste('para',1:length(paraset),sep=''),names.eval.stat)
      
      ### settings
      ## wanted prevalence
      prev.want<-ifelse(atlas,prevalences[spec],0.5)
      ## weights
      weights<-data[,spec]
      weights[weights==0]<-(occurrences[spec]*(1-prev.want))/(prev.want*absences[spec])
      if((is.null(climate.weights)==FALSE)) weights <- weights*climate.weights
      names(weights)<-row.names(data)
      
      
####################
# parameter-set loop
####################
      
      for(i in 1:length(paraset)){
         para<-paraset[[i]]
         cat(paste('\n\nParameter set ',i,'\n###############\n',sep=''))
         
         ### data set for cross validation (will be used for every model)
         data.ana<-data[,c(para,spec)]
         data.ana$folds.id<-kfold(data,k=k,by=data[,spec])
         
         
#---------------
# fit GBM models
#---------------
         if('gbm' %in% models){
            cat(paste('\n-------------------\n||| starting GBM    - parameter set ',i,' of ',spec,'\n-------------------',sep=''))
            
            ## formula creation
            form.gbm<-paste(spec,'~',sep='')
            for(s in 1:length(para)) {
               if(s==1) form.gbm<-paste(form.gbm,para[s],sep='') else form.gbm<-paste(form.gbm,para[s],sep='+')
            }
            
            ## results vector creation
            auc.vec.gbm<-numeric(k)
            pcc.vec.gbm<-numeric(k)
            sensitivity.gbm<-numeric(k)
            specificity.gbm<-numeric(k)
            
            ## setup of shrinkage optimization
            shrink<-mean(gbm.args$shrinkage.bound)
            shrinknew<-shrink
            test<-TRUE;testmin<-FALSE;testmax<-FALSE
            best.ntrees<-integer(k+1)
            counter<-0
            
            ## optimization iteration
            while(test&counter<=gbm.args$maxiter){
               counter<-counter+1
               # set up the new shrinkage if needed
               if(testmin){
                  ifelse(counter==gbm.args$maxiter,shrinknew<-min(gbm.args$shrinkage.bound),shrinknew<-mean(c(shrink,min(gbm.args$shrinkage.bound))))
               }
               if(testmax){
                  ifelse(counter==gbm.args$maxiter,shrinknew<-max(gbm.args$shrinkage.bound),shrinknew<-mean(c(shrink,max(gbm.args$shrinkage.bound))))
               }
               if(counter>gbm.args$maxiter) stop(paste('The process to find a good shrinkage value did not succeed within ',gbm.args$maxiter,' iterations.\nPlease adjust the bounds or the number of in gbm.args.',sep=''))
               shrink<-shrinknew
               # report on progress of optimization
               cat(paste('\n\nIteration ',counter,' of finding good shrinkage for GBM started. shrinkage=',shrink,'\n',sep=''))
               # create appropriate training data set for gbm
               training<-data[order(kfold(data,k=2,by=data[,spec])),]
               # fit the model
               gbm.model<-gbm(as.formula(form.gbm),data=training,weights=weights[row.names(training)],distribution='bernoulli',
                                 n.trees=gbm.args$n.trees,shrinkage=shrink,interaction.depth=gbm.args$interaction.depth,
                                 bag.fraction=0.5,train.fraction=1,cv.folds=0,verbose=FALSE)
               best.itr1<-suppressWarnings(gbm.perf(gbm.model,method='OOB',oobag.curve=plot,overlay=TRUE,plot.it=plot))# estimate number of trees needed for prediction
               # check conditions
               testmin<-best.itr1<2000;testmax<-best.itr1>(0.8*gbm.args$n.trees)
               test<-testmin|testmax
               if(test) next
               # predict the model 
               gbm.pred<-predict.gbm(gbm.model,newdata=training,n.trees=best.itr1,type='response')
               gbm.tss.max<-tss(data[spec],gbm.pred)
               gbm.binom<-as.integer(gbm.pred>gbm.tss.max)
               # k-fold cross validation
               cat('=> CrossValidation ')
               for(v in 1:k) {
                  cat(paste(' > ',v,sep=''))
                  train.data<-data.ana[which(data.ana$folds.id!=v),]
                  test.data<-data.ana[which(data.ana$folds.id==v),]
                  pres<-sum(train.data[,spec]) 
                  absence<-nrow(train.data)-pres
                  train.data$weights1<-train.data[,spec]
                  train.data$weights1[train.data$weights1==0]<-((pres*(1-prev.want))/(prev.want*absence))
                  training.cv<-train.data[order(kfold(train.data,k=2,by=train.data[,spec])),]
                  train.gbm1<-gbm(as.formula(form.gbm),data=training.cv,weights=weights1,distribution='bernoulli',
                                  n.trees=gbm.args$n.trees,shrinkage=shrink,interaction.depth=gbm.args$interaction.depth,
                                  bag.fraction=0.5,train.fraction=1,cv.folds=0,verbose=FALSE)     
                  best.itr1.cv<-suppressWarnings(gbm.perf(train.gbm1,method='OOB',oobag.curve=plot,overlay=TRUE,plot.it=plot))
                  # check conditions
                  testmin<-best.itr1.cv<2000;testmax<-best.itr1.cv>(0.8*gbm.args$n.trees)
                  test<-testmin|testmax
                  if(test) break
                  gbm.pred1<-predict.gbm(train.gbm1,test.data,best.itr1.cv,type='response')
                  gbm.id<-c(1:length(dim(test.data)[1]))
                  pred.occ1<-data.frame(gbm.id,test.data[,spec],gbm.pred1)
                  auc.vec.gbm[v]<-auc(pred.occ1)$AUC 
                  gbm.max.tss <-tss(pred.occ1[,2],pred.occ1[,3])
                  gbm.cmx<-cmx(pred.occ1,threshold=gbm.max.tss)
                  pcc.vec.gbm[v]<-pcc(gbm.cmx)$PCC
                  sensitivity.gbm[v] <-sensitivity(gbm.cmx)$sensitivity
                  specificity.gbm[v] <-specificity(gbm.cmx)$specificity     
               }
               
            }
            
            ## output preparation
            eval.table.cv<-data.frame(rbind(auc.vec.gbm,pcc.vec.gbm,sensitivity.gbm,specificity.gbm),row.names=names.eval.stat)
            names(eval.table.cv)<-names.cv.tab
            xval.auc.gbm<-mean(auc.vec.gbm)
            xval.pcc.gbm<-mean(pcc.vec.gbm)
            xval.sensitivity.gbm<-mean(sensitivity.gbm)
            xval.specificity.gbm<-mean(specificity.gbm)
            mean.eval<-data.frame(xval.auc.gbm,xval.pcc.gbm,xval.sensitivity.gbm,xval.specificity.gbm,row.names=para.names[i])
            names(mean.eval)<-names.eval.stat
            summary.sp.stat.gbm[i,]<-mean.eval
            gbm.out.data<-data.frame(data,gbm.pred=gbm.pred,gbm.binom=gbm.binom)
            
            ## writing files into results directory
            foreign::write.dbf(gbm.out.data,paste(direc,'/',spec,'_curr.pred_gbm_',para.names[i],'.dbf',sep=''))
            write.csv(eval.table.cv,paste(direc,'/',spec,'_rep.eval_gbm_',para.names[i],'.csv',sep=''))
            write.csv(mean.eval,paste(direc,'/',spec,'_mean.cv.eval_gbm_',para.names[i],'.csv',sep='')) 
            oblist.gbm<-c('gbm.model','best.itr1','gbm.tss.max','gbm.out.data','xval.auc.gbm','xval.pcc.gbm','xval.sensitivity.gbm','xval.specificity.gbm')
            file.name<-paste(direc,'/',spec,'_gbm_',para.names[i],sep='')
            save(list=oblist.gbm,file=file.name)
            rm(gbm.model,best.itr1,gbm.tss.max,gbm.out.data,xval.auc.gbm,xval.pcc.gbm,xval.sensitivity.gbm,xval.specificity.gbm,eval.table.cv,mean.eval)
            
         }
#------------------
#||| fit GAM models
#------------------
         if('gam' %in% models){
            cat(paste('\n-------------------\n||| starting GAM    - parameter set ',i,' of ',spec,'\n-------------------',sep=''))
            
            ## formula creation
            form.gam<-paste(spec,'~',sep='')
            scope.gam<-list()
            for(s in 1:length(para)) {
               if(s==1) form.gam<-paste(form.gam,'s(',para[s],',2)',sep='') else form.gam<-paste(form.gam,'+','s(',para[s],',2)',sep='')
               scope.gam[[s]]<-as.formula(paste('~1+s(',para[s],',2)',sep=''))
            }
            
            ## results vector creation
            auc.vec.gam<-numeric(k)
            pcc.vec.gam<-numeric(k)
            sensitivity.gam<-numeric(k)
            specificity.gam<-numeric(k)
            
            ## fit the model
            train.gam<-gam(as.formula(form.gam),family='binomial',data=data,weights=weights)
            environment(step.gam) <- environment()
            gam.step<-step.gam(train.gam,scope=scope.gam,direction='backward',trace=FALSE)
            # predict the model
            gam.pred<-predict.gam(gam.step,newdata=data,type='response')
            gam.tss.max<-tss(data[spec],gam.pred)
            gam.binom<-as.integer(gam.pred>gam.tss.max)
            # store formula for further use
            form.gam2<-gam.step$formula
            
            ## k-fold cross validation
            cat('\n=> CrossValidation ')
            for(v in 1:k){
               cat(paste(' > ',v,sep=''))
               train.data<-data.ana[which(data.ana$folds.id!=v),]            
               test.data<-data.ana[which(data.ana$folds.id==v),]
               pres<-sum(train.data[,spec])
               absence<-dim(train.data)[1]-pres
               weights.gam<-train.data[,spec]
               weights.gam[weights.gam==0]<-((pres*(1-prev.want))/(prev.want*absence))
               train.data$weights.gam<-weights.gam
               train.gam1<-gam(form.gam2,family='binomial',data=train.data,weights=weights.gam)
               gam.pred1<-predict.gam(train.gam1,test.data,type='response')
               gam.id<-1:nrow(test.data)
               pred.occ.gam1<-data.frame(gam.id,test.data[,spec],gam.pred1)
               auc.vec.gam[v]<-auc(pred.occ.gam1)$AUC
               gam.max.tss<-tss(pred.occ.gam1[,2],pred.occ.gam1[,3])
               gam.cmx<-cmx(pred.occ.gam1,threshold=gam.max.tss)
               pcc.vec.gam[v]<-pcc(gam.cmx)$PCC
               sensitivity.gam[v]<-sensitivity(gam.cmx)$sensitivity
               specificity.gam[v]<-specificity(gam.cmx)$specificity
            }
            
            ## output preparation
            eval.table.cv.gam=data.frame(rbind(auc.vec.gam,pcc.vec.gam,sensitivity.gam,specificity.gam),row.names=names.eval.stat)
            names(eval.table.cv.gam)<-names.cv.tab
            xval.auc.gam<-mean(auc.vec.gam)
            xval.pcc.gam<-mean(pcc.vec.gam)
            xval.sensitivity.gam<-mean(sensitivity.gam)
            xval.specificity.gam<-mean(specificity.gam)
            mean.eval.gam<-data.frame(xval.auc.gam,xval.pcc.gam,xval.sensitivity.gam,xval.specificity.gam,row.names=paste('para',i,sep=''))
            names(mean.eval.gam)<-names.eval.stat
            summary.sp.stat.gam[i,]<-mean.eval.gam
            gam.out.data<-data.frame(data,gam.pred=gam.pred,gam.binom=gam.binom)
            
            ## writing files into results directory
            foreign::write.dbf(gam.out.data,paste(direc,'/',spec,'_curr.pred_gam_',para.names[i],'.dbf',sep=''))
            write.csv(eval.table.cv.gam,paste(direc,'/',spec,'_rep.eval_gam_',para.names[i],'.csv',sep=''))
            write.csv(mean.eval.gam,paste(direc,'/',spec,'_mean.cv.eval_gam_',para.names[i],'.csv',sep=''))
            oblist.gam<-c('gam.step','gam.tss.max','gam.out.data','xval.auc.gam','xval.pcc.gam','xval.sensitivity.gam','xval.specificity.gam')
            file.name<-paste(direc,'/',spec,'_gam_',para.names[i],sep='')
            save(list=oblist.gam,file=file.name)
            rm(gam.step,gam.tss.max,gam.out.data,xval.auc.gam,xval.pcc.gam,xval.sensitivity.gam,xval.specificity.gam,eval.table.cv.gam,mean.eval.gam)
         }
         
#------------------
#||| fit GLM models
#------------------
         if('glm' %in% models){
            cat(paste('\n-------------------\n||| starting GLM    - parameter set ',i,' of ',spec,'\n-------------------',sep=''))
            
            ## formula creation
            form.glm<-paste(spec,'~',sep='')
            for(s in 1:length(para)) {
               if(s==1) form.glm<-paste(form.glm,para[s],sep='') else form.glm<-paste(form.glm,para[s],sep='+')
               form.glm<-paste(form.glm,'+I(',para[s],'^2)',sep='')
            }
            
            ## results vector creation
            auc.vec.glm<-numeric(k)
            pcc.vec.glm<-numeric(k)
            sensitivity.glm<-numeric(k)
            specificity.glm<-numeric(k)
            
            ## fit the model
            train.glm<-suppressWarnings(glm(as.formula(form.glm),family='binomial',data=data,weights=weights))
            glm.step<-suppressWarnings(step(train.glm,trace=0,direction='backward'))
            # predict the model
            glm.pred<-predict.glm(glm.step,newdata=data,type='response')
            glm.tss.max<-tss(data[,spec],glm.pred)
            glm.binom<-as.integer(glm.pred>glm.tss.max)
            # store the formula for further use
            form.glm2<-glm.step$formula
            
            ## k-fold cross validation
            cat('\n=> CrossValidation ')
            for(v in 1:k){
               cat(paste(' > ',v,sep=''))
               train.data<-data.ana[which(data.ana$folds.id!=v),]           
               test.data<-data.ana[which(data.ana$folds.id==v),]
               pres<-sum(train.data[,spec]) 
               absence<-dim(train.data)[1]-pres
               weights.glm<-train.data[,spec]
               weights.glm[weights.glm==0]<-(pres*(1-prev.want))/(prev.want*absence)
               train.data$weights.glm<-weights.glm
               train.glm1<-suppressWarnings(glm(form.glm2,family='binomial',data=train.data,weights=weights.glm))
               glm.pred1<-predict.glm(train.glm1,newdata=test.data,type='response')
               glm.id<-1:nrow(test.data)
               pred.occ1<-data.frame(glm.id,test.data[,spec],glm.pred1)
               auc.vec.glm[v]<-auc(pred.occ1)$AUC 
               glm.max.tss <-tss(pred.occ1[,2],pred.occ1[,3])
               glm.cmx<-cmx(pred.occ1,threshold=glm.max.tss)
               pcc.vec.glm[v]<-pcc(glm.cmx)$PCC
               sensitivity.glm[v] <-sensitivity(glm.cmx)$sensitivity
               specificity.glm[v] <-specificity(glm.cmx)$specificity
            }
            
            ## output preparation
            eval.table.cv.glm<-data.frame(rbind (auc.vec.glm,pcc.vec.glm,sensitivity.glm,specificity.glm),row.names=names.eval.stat)
            names(eval.table.cv.glm)<-names.cv.tab
            xval.auc.glm<-mean(auc.vec.glm)
            xval.pcc.glm<-mean(pcc.vec.glm)
            xval.sensitivity.glm<-mean(sensitivity.glm)
            xval.specificity.glm<-mean(specificity.glm)
            mean.eval.glm<-data.frame(xval.auc.glm,xval.pcc.glm,xval.sensitivity.glm,xval.specificity.glm,row.names=paste('para',i,sep=''))
            names(mean.eval.glm)<-names.eval.stat
            summary.sp.stat.glm[i,]<-mean.eval.glm
            glm.out.data<-data.frame(data,glm.pred=glm.pred,glm.binom=glm.binom)
            
            ## writing files into results directory
            foreign::write.dbf(glm.out.data,paste(direc,'/',spec,'_curr.pred_glm_',para.names[i],'.dbf',sep=''))
            write.csv(eval.table.cv.glm,paste(direc,'/',spec,'_rep.eval_glm_',para.names[i],'.csv',sep=''))
            write.csv(mean.eval.glm,paste(direc,'/',spec,'_mean.cv.eval_glm_',para.names[i],'.csv',sep=''))
            oblist<-c('glm.step','glm.tss.max','glm.out.data','xval.auc.glm','xval.pcc.glm','xval.sensitivity.glm','xval.specificity.glm')
            file.name<-paste(direc,'/',spec,'_glm_',para.names[i],sep='')
            save(list=oblist,file=file.name)
            rm(glm.step,glm.tss.max,glm.out.data,xval.auc.glm,xval.pcc.glm,xval.sensitivity.glm,xval.specificity.glm,eval.table.cv.glm,mean.eval.glm)
         }
         
#---------------------
#||| fit maxent models
#---------------------
         if('me' %in% models){
            cat(paste('\n-------------------\n||| starting maxent - parameter set ',i,' of ',spec,'\n-------------------',sep=''))
            
            # check for whether data is atlas data and set maxent arguments accordingly
            if(atlas){
               me.args<-c(maxent.args,paste('defaultprevalence',prevalences[spec],sep='='))
            } else {
               me.args<-maxent.args
            }
            
            ## results vector creation
            auc.vec<-numeric(k)
            pcc.vec.me<-numeric(k)
            sensitivity.me<-numeric(k)
            specificity.me<-numeric(k)
            
            ## fit the model
            background<-data[,para]
            occ<-data[,spec]
            me<-maxent(background,occ,args=me.args,path=me.dir1)  # ,path= to-main-model-info
            # predict the model
            prediction.me<-predict(me,background)
            id<-c(1:nrow(background))
            me.out<-data.frame(id,occ,prediction.me)
            me.max.tss<-tss(me.out[,2],me.out[,3])
            me.binom<-logical(length(id))
            me.binom<-as.integer(prediction.me>me.max.tss)
            
            ## k-fold cross validation
            cat('\n=> CrossValidation ')
            for(v in 1:k){
               cat(paste(' > ',v,sep=''))
               train.data<-data.ana[which(data.ana$folds.id!=v),]
               train.data.back<-train.data[,para]
               train.data.occ<-train.data[,spec]
               test.data<-data.ana[which(data.ana$folds.id==v),]
               test.data.back<-test.data[,para]
               test.data.occ<-test.data[,spec]
               me.auc<-maxent(train.data.back,train.data.occ,args=me.args,path=me.dir2) # ,path= to-cross-validation-models-which-can-be-overwritten
               pred<-predict(me.auc,test.data.back)
               me.id<-c(1:length(test.data.occ))
               pred.occ<-data.frame(me.id,test.data.occ,pred)
               auc.vec[v]<-auc(pred.occ)$AUC
               me.max.tss <-tss(pred.occ[,2],pred.occ[,3])
               me.cmx<-cmx(pred.occ,threshold=me.max.tss)
               pcc.vec.me[v]<-pcc(me.cmx)$PCC
               sensitivity.me[v] <-sensitivity(me.cmx)$sensitivity
               specificity.me[v] <-specificity(me.cmx)$specificity		
            }
            
            ## output preparation
            eval.table.cv.me = data.frame(rbind (auc.vec,pcc.vec.me,sensitivity.me,specificity.me),row.names = names.eval.stat)
            names(eval.table.cv.me)<-names.cv.tab
            xval.auc.me<-mean(auc.vec)
            xval.pcc.me<-mean(pcc.vec.me)
            xval.sensitivity.me<-mean(sensitivity.me) 
            xval.specificity.me<-mean(specificity.me)	
            mean.eval.me<-data.frame(xval.auc.me,xval.pcc.me,xval.sensitivity.me,xval.specificity.me,row.names = paste('para',i,sep=''))                       
            names(mean.eval.me)<-names.eval.stat
            summary.sp.stat.me[i,]<-mean.eval.me
            me.out.data<-data.frame(data,me.pred=prediction.me,me.binom=me.binom)
            
            ## writing files into results directory
            foreign::write.dbf(me.out.data,paste(direc,'/',spec,'_curr.pred_me_',para.names[i],'.dbf',sep=''))
            write.csv(eval.table.cv.me,paste(direc,'/',spec,'_rep.eval_me_',para.names[i],'.csv',sep=''))
            write.csv(mean.eval.me,paste(direc,'/',spec,'_mean.cv.eval_me_',para.names[i],'.csv',sep=''))
            oblist<-c('me','me.max.tss','me.out.data','xval.pcc.me','xval.sensitivity.me','xval.specificity.me')
            file.name<-paste(direc,'/',spec,'_me_',para.names[i],sep='')
            save(list=oblist,file=file.name)
            rm(me,me.max.tss,me.out.data,xval.pcc.me,xval.sensitivity.me,xval.specificity.me,eval.table.cv.me,mean.eval.me)
         }
      } # END parameter-set loop
   
   ### writing summary files into results directory
   
   ## create directory
     direc2 <- paste(direc,'/summary_stat_models/',sep='')
     system(paste('mkdir -p ',direc2,sep=""))
     for(u in models){
        write.csv(get(paste('summary.sp.stat.',u,sep='')),file=paste(direc2,'/sum.stat.',u,'_',spec,'.csv',sep=""))
     }
   } #END species loop
}

