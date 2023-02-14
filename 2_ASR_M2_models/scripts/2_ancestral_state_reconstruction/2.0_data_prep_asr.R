
#############################
# PRUNE TREE TO SPP W/ DATA #
#############################
## The asr techniques we used do not code for missing states, so species for whom no data can be found 
## must be removed from the dataset for the ancestral state reconstruction 


## Prune new consensus molecular clock 

tree<-read.nexus("data/2_asr/cnidaria-final_908.con.tre")
keep<-scan(file="data/2_asr/species477.txt",what = character(), sep = "\t")
pruned_tree<-keep.tip(tree,keep)

# Save new tree file
writeNexus(pruned_tree, "data/2_asr/cnidaria-final477.tre")


#################################
# CREATE AN M2 CHARACTER MATRIX #
#################################

# Read in raw data to data frame
data<-read.table(file="data/0_data_cleanup/cnidome_v3.txt",header=T, sep = "\t")
# Read in tree
tree<-read.nexus("data/cnidaria-final_477.tre")

# A k=2 model where 
# 0 = ABSENCE in species 
# 1 = PRESENCE in species 
# ? = NO DATA in species

# Get unique cnidocytes to generate M2 file for each
cnidocytes=as.vector(unique(data[data$Cnidocyte!="","Cnidocyte"]))

for (type in cnidocytes){
  # Make an empty df, populate states$Species with all species on tree
  states <- setNames(data.frame(matrix(ncol=2, nrow=length(tree$tip.label))), c("Species", "State"))
  states[,1]<-sort(tree$tip.label)
  
  # Populate states$State with cnidocyte present(1), absent(0), no data (?)
  spp_with <- as.vector(unique(data[which(data$Cnidocyte == type),7]))
  states$State<-ifelse(states$Species %in% spp_with,1,0)
  
  # Save as a txt file
  write.table(states, file=paste0("data/2_asr/character-data_477/",type,"_v3.txt"), sep="\t",quote=F, row.names=F)
  
 
}; rm(spp_with,type,states,)
