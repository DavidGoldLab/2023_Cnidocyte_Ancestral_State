

###############################
## 0) Load packages and data ##
###############################
  
  library(phytools)
  library(corHMM)
  library(beepr) # optional

  # Open data
  table<-read.table(file="Supplemental-materials/data/character_04-19.txt", sep="\t",header=TRUE)
  # Open species list of 477
  spp477<-scan("Supplemental-materials/data/species477.txt", character(), quote = "")
  
#################################
## 1) Make isorhiza M4 dataset ##
#################################
  
  # SUBSET DATA BY TRICH TYPE
  
  # Subset isorhiza
  isorhiza.all<- subset(table, Cnidocyte == "isorhiza")
  # Subset only those with trich definitions
  isorhiza <- isorhiza.all[grep("trichous", isorhiza.all$Descriptor), ]
    # Check - should be true
    sum(isorhiza$Descriptor=="")==0
  # Clean off additional descriptor phrases (large, egg-shaped, etc)
  isorhiza$Descriptor[grep("atrichous",isorhiza$Descriptor)]<-"atrichous"
  isorhiza$Descriptor[grep("holotrichous",isorhiza$Descriptor)]<-"holotrichous"
    # Check - should be 6 values (a-,apo-,basi-,holo-,homo-,mero-)
    View(unique(isorhiza$Descriptor))
  
  # Save
  write.table(isorhiza,"Supplemental-materials/isorhizaM4/isorhiza_datatable.txt",row.names=F,quote=F,sep="\t")
  
  # Now subset as binary to dataframes
  atrichous<-unique(subset(isorhiza[c("Species", "Cnidocyte","Descriptor")], Descriptor == "atrichous"))
  trichous<-unique(subset(isorhiza[c("Species", "Cnidocyte","Descriptor")], Descriptor != "atrichous"))
  
  # Create vectors for species belonging in categories 1,2,3
  spp_both<-intersect(atrichous$Species,trichous$Species) #3
  spp_atrich<-setdiff(as.character(atrichous$Species),spp_both) #1
  spp_trich<-setdiff(as.character(trichous$Species),spp_both) #2
  
  #343 species contain descriptive info of isohiza
  #111 will be labeled 0
  #total on tree should be 454
  
  # Mcat 0 is all species that have cnidocyte descriptions that do not include isorhiza
    # isorhiza with no trich description are essentially missing data and are pruned from tree
  iso_no_desc<-setdiff(isorhiza.all$Species,isorhiza$Species) # union(spp_atrich,spp_trich) is all species with isorhiza with trich decriptions
  
  # Prune these species from the dataset and tree
  isorhiza.tree<-drop.tip(read.nexus("Supplemental-materials/data/cnidaria-final_477.tre"),iso_no_desc)
  # Save a copy of tree
  writeNexus(isorhiza.tree, "Supplemental-materials/isorhizaM4/isorhizaM4-final_451.tre")
  

# CREATE M4 DATASET
  # Make an empty df, populate states$Species with all species on 451 tree
  spp451<-isorhiza.tree$tip.label
  isorhizaM4 <- setNames(data.frame(matrix(ncol=2, nrow=length(spp451))), c("Species", "State"))
  isorhizaM4[,1]<-spp451
  
  # Populate
  isorhizaM4$State<-ifelse(isorhizaM4$Species %in% spp_both,3,ifelse(isorhizaM4$Species %in% spp_trich,2,ifelse(isorhizaM4$Species %in% spp_atrich,1,0)))
  
  # Save
  write.table(isorhizaM4,"Supplemental-materials/isorhizaM4/isorhizaM4.txt",row.names=F,quote=F,sep="\t")
  
