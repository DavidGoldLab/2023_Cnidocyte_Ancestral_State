#!/usr/bin/Rscript

#SBATCH -D output/ # Project directory
#SBATCH --mail-type=ALL # types of emails to send
#SBATCH -o # path/to/output/logs/stochmap-example-o%j.txt
#SBATCH -e # path/to/error/logs/stochmap-example-e%j.txt
#SBATCH --mail-user= example@university.edu # Email to which notifications will be sent
#SBATCH -J sm-example # job name
#SBATCH -t 70:00:00 # wall time limit

# Run the following in the primary environment before starting script
#export R_LIBS=~/tools/Rlibs

type="acrophore" 
pmodel="ER"
cat=1

nSim = 10000 


######################################
## 0) Load environment and packages ##
######################################

# Load packages
library(geiger)
library(phytools)
library(corHMM)

# Set up directories
root="projects/cnidome/" 
date=gsub(paste0(format(Sys.Date(), "%Y"),"-"),"",Sys.Date())

# Initialize useful functions
set_class <- function(x, class, add=c("overwrite", "prepend", "append")){
  add = match.arg(add)
  class <- switch(add,
                  overwrite = class,
                  prepend = unique(c(class, class(x))),
                  append = unique(c(class(x), class)))
  if(length(class) == 0) # for class == ""
    unclass(x)
  else
    `class<-`(x, class)
} # credit: https://rdrr.io/github/thackl/thacklr/src/R/set_class.R

as.multiPhylo.list <- function(x, ...){
  if(!all(sapply(x, function(y) inherits(y, "phylo"))))
    stop("Need a list of phylo objects")
  set_class(x, "multiPhylo")
} # credit: https://rdrr.io/github/thackl/thacklr/src/R/phylo.R



##########################################
## 1) Load data and corHMM model object ##
##########################################

# Load tree
pruned.tree<-read.nexus(paste0("data/cnidaria-final477.tre"))

# Format data for corHMM
fitHRM.data<-read.table(paste0("data/character-data_477/",type,"_v3.txt"),header=TRUE,row.names=1,stringsAsFactors=TRUE)
if(name.check(pruned.tree,fitHRM.data)=="OK"){state.p<-setNames(fitHRM.data$State,rownames(fitHRM.data))} else{print("tip disagreement")}
corhmm.data<-data.frame(Genus_species=names(state.p),state=as.numeric(state.p)); rm(state.p,fitHRM.data)

## Reading back in corHMM model object
if(cat==2){
corHMM.model<-readRDS(file=paste0("output/",mt_date,"/model-object/cor.",pmodel,cat,"cat.",type,".RData"))
} else {corHMM.model<-readRDS(file=paste0("output/",mt_date,"/model-object/cor.",pmodel,".",type,".RData"))}



###############################
## 2) Run stochastic mapping ##
###############################

### Running with single best model using corHMM makeSimmap

# Set up model
phy = corHMM.model$phy
model = corHMM.model$solution
model[is.na(model)] <- 0
diag(model) <- -rowSums(model)

if(cat==2){root.p=c(1,0,0,0)
} else {root.p=c(1,0)}

# Run makeSimmap
sim.obj<-makeSimmap(tree = phy, data = corhmm.data, model = model, rate.cat=cat, nSim = nSim, root.p=root.p)
  
# Describe:
sim.obj.phylo<-as.multiPhylo.list(sim.obj) 
sim.obj.summ<-describe.simmap(sim.obj.phylo)
# change ace colnames to match presence/absence representation
  if(ncol(sim.obj.summ[["ace"]])==4){colnames(sim.obj.summ[["ace"]])<-c("0*","1*","0","1")} else {colnames(sim.obj.summ[["ace"]])<-c("0","1")}
saveRDS(sim.obj.summ,file=paste0("output/",date,"/simmap-object/sim",nSim,".",pmodel,cat,"cat.summ.",type,".tre"))

#######################
## 3) Plot summaries ##
#######################

# Set up color labels
if(cat==2) {
cols<-setNames(c("red","blue","lightcoral","skyblue"),c("2","1","2*","1*"))
colsL<-setNames(c("red","blue"),c("Presence","Absence"))
} else {
cols<-setNames(c("red","blue"),c("2","1"))
colsL<-setNames(c("red","blue"),c("Presence","Absence"))
}

# Plot node values
    pdf(paste0("output/",date,"/sim",nSim,".",pmodel,cat,"cat.",type,"-nodeSummary.pdf"), width=6, height=60)
    plot(sim.obj.summ,fsize = 0.3,colors=cols, lwd=0.5, offset=0.5,cex=0.25,setEnv=T)
    tiplabels(pie=corhmm.data$state,piecol=cols,cex=0.25,offset=5) 
    nodelabels(node=1:pruned.tree$Nnode+Ntip(pruned.tree),pie=sim.obj.summ$ace,piecol=rev(cols),cex=0.3)
    add.simmap.legend(colors=colsL,prompt=FALSE,x=0.9*par()$usr[1],y=0.1*par()$usr[3],fsize=0.5)
    dev.off()


# Plot density map (summarize on branches)
    pdf(paste0("output/",date,"/sim",nSim,".",pmodel,cat,"cat.",type,"-densityMap.pdf"), width=6, height=60)
    densityMap(sim.obj.phylo,fsize = 0.3,colors=cols, lwd=1, offset=0.5,cex=0.25,legend=F,setEnv=T)
    tiplabels(pie=corhmm.data$state,piecol=cols,cex=0.2,offset=9) 
    add.simmap.legend(colors=colsL,prompt=FALSE,x=0.9*par()$usr[1],y=0.1*par()$usr[3],fsize=0.5)
    dev.off()
    
