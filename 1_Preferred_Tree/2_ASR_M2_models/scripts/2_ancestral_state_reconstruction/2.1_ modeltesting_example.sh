#!/usr/bin/Rscript

#SBATCH -D output/ # Project directory
#SBATCH --mail-type=ALL # types of emails to send
#SBATCH -o # path/to/output/logs/modeltest-example-o%j.txt
#SBATCH -e # path/to/error/logs/modeltest-example-e%j.txt
#SBATCH --mail-user= example@university.edu # Email to which notifications will be sent
#SBATCH -J mt-example # job name
#SBATCH -t 70:00:00 # wall time limit


###############################
## 0) Load packages and data ##
###############################

# Run the following in the bash environment before starting script to allow slurm to locate R packages:
#export R_LIBS=~/tools/Rlibs #point path to R libraries in your account

# Recommend creating a script for each cnidocyte type and running in parallel on an hpc; 
	# alternatively the following can be copied into an R script and run locally:

type="acrophore"

dir.create(file.path('output',date), recursive = TRUE)
dir.create(file.path('output',date,'model-object'), recursive = TRUE)

root="output/" 
date=gsub(paste0(format(Sys.Date(), "%Y"),"-"),"",Sys.Date())

# load packages
library(geiger)
library(phytools)
library(corHMM)

pruned.tree<-read.nexus("data/2_asr/cnidaria-final_477.tre")

# format data for fitMk
fitHRM.data<-read.table(paste0("data/2_asr/character-data_477/",type,"_v3.txt"),header=TRUE,row.names=1,stringsAsFactors=TRUE)
# Make sure tree and data match
name.check(pruned.tree,fitHRM.data)
# Extract discrete character
state<-setNames(fitHRM.data$State,rownames(fitHRM.data))
# Create "states" vector
pi.states<-sort(unique(state))
# Create "states" factor for HRM 
cat2states <- c("0", "0*", "1", "1*")
pi.cat2.states <- factor(cat2states)

# Format data for corHMM
if(name.check(pruned.tree,fitHRM.data)=="OK"){state.p<-setNames(fitHRM.data$State,rownames(fitHRM.data))} else{print("tip disagreement")}
corhmm.data<-data.frame(Genus_species=names(state.p),state=as.numeric(state.p)); rm(state.p)



########################################
## 1) Run models on fitHRM and corHMM ##
########################################

# Run corHMM with 1 and 2 rate categories to generate Q matrix with user-input root probability
# Comparison of BF, AIC, AICc and -logLik to compare best fit model for each cnidocyte

## NCAT = 1, root fixed

cor.ER<-corHMM(phy=pruned.tree,data=corhmm.data,rate.cat=1,nstarts=20,model="ER",node.states="marginal",root.p=c(1,0),n.cores=2)
cor.ARD<-corHMM(phy=pruned.tree,data=corhmm.data,rate.cat=1,nstarts=20,model="ARD",node.states="marginal",root.p=c(1,0),n.cores=2)


## NCAT = 2, root fixed

cor.ER2cat<-corHMM(phy=pruned.tree,data=corhmm.data,rate.cat=2,nstarts=20,model="ER",node.states="marginal",root.p=c(1,0,0,0),n.cores=2)
cor.ARD2cat<-corHMM(phy=pruned.tree,data=corhmm.data,rate.cat=2,nstarts=20,model="ARD",node.states="marginal",root.p=c(1,0,0,0),n.cores=2)


# Save model objects
saveRDS(cor.ER, file=paste0("output/",date,"/model-object/cor.ER.",type,".RData"))
saveRDS(cor.ARD, file=paste0("output/",date,"/model-object/cor.ARD.",type,".RData"))
saveRDS(cor.ER2cat, file=paste0("output/",date,"/model-object/cor.ER2cat.",type,".RData"))
saveRDS(cor.ARD2cat, file=paste0("output/",date,"/model-object/cor.ARD2cat.",type,".RData"))

