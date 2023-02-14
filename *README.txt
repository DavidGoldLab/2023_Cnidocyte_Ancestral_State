
1) Summary of process          
2) directory structure
3) packages necessary
	3.1) local R packages (data cleanup)
	3.2) hpc R packages
	3.3) downloading R packages on a cluster
	3.4) pointing bash session toward R packages
4) running r and bash scripts




 
Example bash files are meant to be made for each cnidocyte so they can be run in parallel. At the top of each script, fill in "type" and other starting variables for each specific cnidocyte and your results.
 
Directories are currently pointed in a generic fashion, please update each script to locate your input and output folders when running.

Required software:
 	* R v3.5.1 https://cran.r-project.org/bin/
 	* starttree - available on our ...
 	* IQTree v1.6.12 https://github.com/Cibiv/IQ-TREE/releases/tag/v1.6.12
 	* BEAST v1.10.4 
 	* 


 
SUMMARY:
 
#####################
## 0) Data cleanup ##
#####################
 	

 	
 	Input:
 	-starting treefile (tree_cnid_635.tre)
 	-multiple sequence alignments 5 genes, fasta format (alignments_18/*.fasta)
* 	-starttree input file, txt format ()
	
	0.1)
 	
 
#########################################
## 1) Fossil-calibrated Molecular Clock##
#########################################
 
 	Generate a molecular clock that is calibrated on chronological time using BEAST v1.10.4 (https://beast.community/2018-11-14.10.4_released.html). Branch lengths are estimated from a 5-gene dataset, and the tree is calibrated on a chronological scale using fossil data.
 
 	Input: 
 	-starttree initialized tree (pruned_tree_input.nex)
 	-concatenated genetic data (cat-genes-pruned.nex)
 	-fossil calibration lists, monophyletic species lists for each clade (fossil_calibrations/*.txt)
 	-IQTree model output decisions
* 	-BEAST settings file
 	
 	Notes: fossil calibration lists are defined as crown or stem fossils based on individual record; noted in the filename for each list file.
 	
 	BEAUti (BEAST GUI) is launched within the data_cleanup.Rmd from step (0). Settings for the xml file are input within each BEAUti tab. For help with settings, see https://beast.community/ How-To-Guides tab on left.
 
 
####################################### 
## 2) Ancestral State Reconstruction ##
#######################################
 
 	Conducted in two parts using the Rpackage corHMM (https://cran.r-project.org/web/packages/corHMM/index.html). This is conducted in two steps to ensure that consensus values can be reached for each parameter (due to the size of the models; number of species and model parameters) Model Testing determines optimal values of Qmatrix parameters under each model for cnidocyte dataset (ER or ARD models of substitution; 1 rate category or 2 rate category hidden markov models. see https://doi.org/10.1111/2041-210X.13534 for more information). Best-fit model is the chosen based on AIC values, and this Qmat data is input as a fixed value into the stochastic character mapping process, which then simulates logical histories based on the model of substitution and outputs a likely set of histories.

 	
 	2.0) DATA PREP
 	
 		Input:
 		-time-calibrated tree (cnidaria-final_908.con.tre)
 		-character data (cnidome_v3.txt)
 	
 		2.0.1) Create a time-calibrated pruned tree containing only the species with character data (477 species)
 		2.0.2) Create M2 data files for each cnidocyte where presence and absence of each type are coded as a binary (0/1)
 	
 
 	2.1) MODEL TESTING
 	
 		*Run on linux hpc, or on Rstudio without slurm wrapper
 	
 		Input:
 		-M2 datasets (character-data_477/TYPE_v3.txt)
 		-time-calibrated tree pruned to species with character data (cnidaria-final_477.tre)
 	
 	2.2) STOCHASTIC MAPPING
 	
 		*Run on linux hpc, or on Rstudio without slurm wrapper
 	
 		Input:
 		-Model testing output (model-object/cor.ER.TYPE.RData)
 		-character data (cnidome_v3.txt)
 		-time-calibrated tree pruned to species with character data (cnidaria-final_477.tre)
 		-tree processing functions from thacklr github (in script)
 	
 	
 
 