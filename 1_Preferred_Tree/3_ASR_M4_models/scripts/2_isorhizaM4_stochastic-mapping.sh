#!/usr/bin/Rscript

#SBATCH -D /home/ncsierra/projects/cnidome/asr
#SBATCH --mail-type=ALL 
#SBATCH -o /home/ncsierra/projects/cnidome/slurm-logs/M4isorhiza-o%j.txt
#SBATCH -e /home/ncsierra/projects/cnidome/slurm-logs/M4isorhiza-e%j.txt
#SBATCH --mail-user=ncsierra@ucdavis.edu  # Email to which notifications will be sent
#SBATCH -J M4isorhiza
#SBATCH -t 70:00:00


# Run the following in the primary environment before starting script
#export R_LIBS=~/tools/Rlibs


################
## 3.0) Setup ##
################

### **The following is best run on an hpc cluster
	## upload 451 pruned tree and isorhizaM4 dataset to run

library(corHMM)
library(phytools)
require(ape)

mkdir("Supplemental-materials/isorhizaM4/simmap-object")

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

# Read in data
isorhiza.tree<-read.nexus("Supplemental-materials/isorhizaM4/isorhizaM4-final_451.tre")
isorhizaM4<-read.table("Supplemental-materials/isorhizaM4/isorhizaM4.txt",sep="\t",header=TRUE)

# Selected ER2cat as best model
cor.ER.2cat.fixed<-readRDS(file=paste0("Supplemental-materials/isorhizaM4/model-object/cor.ER.2cat.fixed.RData"))

# Set up model
model = cor.ER.2cat.fixed$solution
model[is.na(model)] <- 0
diag(model) <- -rowSums(model)

#############################
## 3.1) Stochastic Mapping ##
#############################

# Run stochmap
sim.obj<-makeSimmap(tree = pruned.tree, data = isorhizaM4, model = model, rate.cat=2, nSim = 10000, root.p=c(1,0,0,0,0,0,0,0))

#Describe
sim.obj.phylo<-as.multiPhylo.list(sim.obj) 
sim.obj.summ<-describe.simmap(sim.obj.phylo)

# Rewrite ace value to legible name and save
colnames(sim.obj.summ[["ace"]])<-c("0","1","2","3","0*","1*","2*","3*")

# Save asr output
saveRDS(sim.obj, file=paste0("Supplemental-materials/isorhizaM4/simmap-object/sim10000.ER2cat.isorhizaM4.tre"))
saveRDS(sim.obj.summ, file=paste0("Supplemental-materials/isorhizaM4/simmap-object/sim10000.ER2cat.summ.isorhizaM4.tre"))

###############################
## 3.2) Plot Stochmap Output ##
###############################

# Set up color labels - depends on model setup
# M4 (8 colors)
cols<-setNames(c("blue","tan","red","purple","skyblue","lemonchiffon","lightcoral","plum1"),c("0","1","2","3","0*","1*","2*","3*"))
cols2<-setNames(c("blue","tan","red","purple","skyblue","lemonchiffon","lightcoral","plum1"),c("1","2","3","4","1*","2*","3*","4*"))
colsL<-setNames(c("blue","tan","red","purple"),c("Absence","Atrichous","Trichous","Both"))


# Plot M4 - following method for plotting from: http://blog.phytools.org/2020/06/plotting-states-of-discrete-character.html
val<-as.factor(setNames(isorhizaM4$State,isorhizaM4$Species))

pdf(paste0("Supplemental-materials/isorhizaM4/sim10000.ER2cat.isorhizaM4-nodeSummary.pdf"), width=6, height=60)
plotTree(pruned.tree, fsize = 0.4, lwd=0.5, cex=0.25, offset=0.5,setEnv=T)
VAL<-to.matrix(val,levels(val))
tiplabels(pie=VAL,piecol=cols,cex=0.25,offset=5)
nodelabels(node=1:pruned.tree$Nnode+Ntip(pruned.tree),pie=sim.obj.summ$ace,piecol=cols2,cex=0.3)
dev.off()

