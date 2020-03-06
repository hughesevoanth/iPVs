# iPVs
## identification of principal variables

### Authors: David Hughes 
##### Date started: 6th March 2020

### About:

This repo contains an R package called iPVs. The purpose of iPVs is to aid users with with the identification of representative (or independent) principal variables in an inter-correlated data set. The count of representative principal variables may also be used as a count of the number of independent variables (Me) in your data set. 

### Installation instructions of the iPVs package

	1. insure that the devtools library is installed on your local machine
		 > ifelse("devtools" %in% rownames(installed.packages()), 
		 NA, install.packages("devtools") )
		 
	2. install iPVs
		> devtools::install_github("hughesevoanth/iPVs")

___
		
## Describing the independent principal variable pipeline (iPVs)

This suite of scripts can help you identify a set of independent features or variables in a dataset with correlated features, such as that found in microbiomics or metabolomics. It is based on Spearman's rank correlations and a tree cutting step to identify clusters of features that are then used to identify a tagging variable using principal variable analysis (PVA). 

Once you have installed the moosefun package as described above, open an R session and follow the example below for a relatively straight forward run through of the iPVs pipeline using a single wrapper function. 

#### iPVs()

	## load some needed packages
	library(moosefun)
	library(Hmisc)
	library(tidyverse)
	
	## load your dataset (mine is a flat text file tab-delimited)
	n = "pathtomydataset/mydataset.txt"
	mydata = read.table(n, header = TRUE, sep = "\t", as.is = TRUE)
	
	## perform any QC that you would like
	## I will filter on feature missingness
	fmis = apply(mydata, 2, function(x){  sum(is.na(x))/length(x) })
	w = which(fmis > 0.2)
	wdata = mydata[, -w]
	
	## NO data transformation is required
	## as we will be using the non-parametric 
	## Spearman's rank correlation. 
	
	## A quick and easy wrapper function to do everything for you is:
	mypvs = iPVs(wdata)
	
	## Done !
	
#### results: "mypvs" will contain a list of 4 objects:
	
	1. iPV_table -- a data.frame containing:
		- the PVs identifiers
		- the variance explained by the PV for its cluster of features|variables
		- a single character containing all feature IDs in the PVs cluster (as found in PV_cluster_members), written as "ID1:ID2:ID3:"
	2. PV_cluster_members -- a list providing all feature|variable ids in each cluster, supplied as a vector.
	3. PVresults -- a list providing all of the PVA results for each cluster in PV_cluster_members
	4. workingdata -- a list returning your:
			- provided dataset
			- its correlation matrix
			- its distance matrix
			- the initial tree of the complete dataset
			- eigenvalues derived from the correlation matrix
			- pca derived from the correlation matrix
			- varexp in PCA
			- and estimates of Me, the effective number of markers. 

#### plot your tree with some color coding for the iPVs
	
	## load a needed R package
	library(dendextend)

	## extract the IDs for your PVs
	pv_ids = as.character(mypvs$iPV_table$PVs )

	## define your tree as a dendrogram
	dend = mypvs$workingdata$tree %>% as.dendrogram

	## create a vector of colors to color your tree labels
	n = labels(dend)
	pcol = rep("black", length(n))
	w = which(n %in% pv_ids ); pcol[w] = "medium blue"

	## redefine elements of dendrogram
	dend = dend %>% 
	set("labels_cex", 0.5) %>% 
  	set("labels_col", pcol) %>% 
  	set("branches_k_color",  value = pcol)

	## plot the dendrogram
	dend %>% plot(main = paste0( "-- Principle Variables --" ) )
	abline(h = 0.5, col = "red", lwd = 3)

#### alternative - long hand

If you would perfer to do everything step by step then feel free to print the function iPVs() to screen in an R session like 
	
	> iPVs 

and use its code to run through the steps of the pipeline. 


---
### Some  other functions in the iPVs suite include


#### treecut.pvs( tree, variabledata, cutheight )
- which will extract principal variables for a single step tree cut, given a specified cut height. 
- you must provide an hclust object, your raw variable data set (variables in cols) and the height at which to cut your tree. It is recommended to generate your tree with hclust( mydistmat, method = "complete"), when your distance matrix is 1 - abs(correlation coefficient).

#### treecut.sumstats(tree, cormat, cutheight)
- will calculate a few summary statistics for any desired tree cut height. 
- specifically, for any identified cluster with an N > 1 it will provide:
	- the number of features in the cluster
	- the minimum rho observed among members of the cluster
	- the mean rho observed among members of the cluster
	- the max rho observed among members of the cluster
	- the standard deviation in rho among members of the cluster, assuming N > 2

#### treecut.iterator.4minrho( tree, cormat, cutheights)
- will calculate the minimum rho observed (1) among members of identified clusters and then (2) among all of those minimums, given a vector (series) of cut heights. In effect it is identifying the longest branch in any cluster at a given cut height and returning the dissimularity value. 

