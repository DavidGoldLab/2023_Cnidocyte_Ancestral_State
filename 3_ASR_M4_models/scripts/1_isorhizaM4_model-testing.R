
######################
## 2) Model Testing ##
######################

mkdir("Supplemental-materials/isorhizaM4/model-object")

# Read in data
isorhiza.tree<-read.nexus("Supplemental-materials/isorhizaM4/isorhizaM4-final_451.tre")
isorhizaM4<-read.table("Supplemental-materials/isorhizaM4/isorhizaM4.txt",sep="\t",header=TRUE)

# BUILD MODELS

# NCAT = 1
LegendAndRate<-getStateMat4Dat(isorhizaM4)

ARD_R1 <- dropStateMatPars(LegendAndRate$rate.mat, c(3, 10)) # drop 1->4 and 4->1
pars2equal <- list(c(1, 3), c(2, 6), c(4,7), c(5,9), c(8,10))
ER_R1 <- equateStateMatPars(ARD_R1, pars2equal)

# NCAT = 2 
ERR2<-matrix(c(0,6,7,0,6,0,8,9,7,8,0,10,0,9,10,0),nrow=4) # create second with same symmetrical structure, just different parameters
# Set transition between rate cats
RateClassMat <- getRateCatMat(2) 
# Now group all into a model; the first element of the list corresponds to R1, the second to R2, etc.
StateMatsER <- list(ER_R1, ERR2)
ER_R2 <- getFullMat(StateMatsER, RateClassMat)

# Set up R2 for ARD
ARDR2<-matrix(c(0,11,12,0,16,0,13,14,17,18,0,15,0,19,20,0),nrow=4) # create second with same asymmetrical structure, just different parameters
# Set transition between rate cats
RateClassMat <- getRateCatMat(2) # same as other one
# Now group all into a model; the first element of the list corresponds to R1, the second to R2, etc.
StateMatsARD <- list(ARD_R1, ARDR2)
ARD_R2 <- getFullMat(StateMatsARD, RateClassMat)

# RUN MODELS

cor.ER.1cat.fixed<-corHMM(phy=isorhiza.tree,data=isorhizaM4,rate.cat=1,root.p=c(1,0,0,0),nstarts=20,rate.mat=ER_R1,node.states="marginal",n.cores=2)
cor.ARD.1cat.fixed<-corHMM(phy=isorhiza.tree,data=isorhizaM4,rate.cat=1,root.p=c(1,0,0,0),nstarts=20,rate.mat=ARD_R1,node.states="marginal",n.cores=2)
cor.ER.2cat.fixed<-corHMM(phy=isorhiza.tree,data=isorhizaM4,rate.cat=2,root.p=c(1,0,0,0,0,0,0,0),nstarts=20,rate.mat=ER_R2,node.states="marginal",n.cores=2)
cor.ARD.2cat.fixed<-corHMM(phy=isorhiza.tree,data=isorhizaM4,rate.cat=2,root.p=c(1,0,0,0,0,0,0,0),nstarts=20,rate.mat=ARD_R2,node.states="marginal",n.cores=2)
beep(3)

# Save output
saveRDS(cor.ER.1cat.fixed, file=paste0("Supplemental-materials/isorhizaM4/model-object/cor.ER.1cat.fixed.RData"))
saveRDS(cor.ARD.1cat.fixed, file=paste0("Supplemental-materials/isorhizaM4/model-object/cor.ARD.1cat.fixed.RData"))
saveRDS(cor.ER.2cat.fixed, file=paste0("Supplemental-materials/isorhizaM4/model-object/cor.ER.2cat.fixed.RData"))
saveRDS(cor.ARD.2cat.fixed, file=paste0("Supplemental-materials/isorhizaM4/model-object/cor.ARD.2cat.fixed.RData"))

# Model choice was performed using AIC
cor.ER.1cat.fixed$AIC
cor.ARD.1cat.fixed$AIC
cor.ER.2cat.fixed$AIC
cor.ARD.2cat.fixed$AIC

cor.ER.1cat.fixed$loglik
cor.ARD.1cat.fixed$loglik
cor.ER.2cat.fixed$loglik
cor.ARD.2cat.fixed$loglik

