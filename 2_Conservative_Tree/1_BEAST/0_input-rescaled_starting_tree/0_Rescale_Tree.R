# R version 4.2.1

library(RRphylo)
library(phytools)


#Read the topology
phy <- read.nexus("0_Initial_Tree.tree")

# Find node IDs

# export tree with node labels
pdf("1_nodelabels.pdf", width=20, height=100)
plot(phy, cex=0.5)
nodelabels(col= "red", cex=0.5, frame="none")
dev.off()

# Root --> NODE 909
# Bilateria (crown @ ~550mya) -> NODE 1815
# Medusozoa (stem @ 529mya) -> NODE 1343
# Anthozoa (crown @ 530mya) -> NODE 913
# Stylaster (stem @ 66mya) --> NODE 1572
# Poritidae (stem @ 72mya) --> NODE 1258
# Limnomedusa (stem @ 242mya) --> NODE 1413
# Scleractinia (stem @ 259.9mya) --> NODE 1192
# Cubozoa (crown @ 307mya)  --> NODE 1376
# Scyphozoa (stem @ 521mya)  --> NODE 1346
# Zoantharia (stem @ 513mya)  --> NODE 1284
# Leptothecata (crown@ 485.4mya)  --> NODE 1641

# Scale tree
	# Add 3My to crown groups
nodeAges<-c(635,555,529,533,66,72,242,259.9,310,521,513,486)
names(nodeAges)<-c(909,1815,1343,913,1572,1258,1413,1192,1376,1346,1284,1641)
scaleTree(phy,node.ages=nodeAges)->treeS2
writeNexus(treeS2, "2_rescaled_starting_tree.tree")