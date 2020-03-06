# Random useful functions for Moose's work

### Authors: David Hughes 
##### Date started: 31st Jan. 2020

### About:

This repo contains a variety of helpful functions that I come back to time and time again.  Rather than adding them into a variety of different packages I will try to maintain this catch-all that I and my colleagues can source when necessary. 

### Installation instructions of the moosefun package

	1. insure that the devtools library is installed on your local machine
		 > ifelse("devtools" %in% rownames(installed.packages()), 
		 NA, 
		 install.packages("devtools"))
		 
	2. install moosefun
		> devtools::install_github("hughesevoanth/moosefun")

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


---

## DESCRIPTION OF OTHER FUNCTIONS


### ztransform()
- taken from the GenABEL package to perform z-transformations. This function is necessary to run my edited version of the rank normal transformation function.

### rntransform()
- taken from the GenABEL package to perform rank normal transformations, but **edited to randomly split tied values**. 

### moose_biplot()
-  a function to plot a biplot or a PCA with a loadings plot on top of it. 
-  however, the uniqueness here is that we are not plotting the loading from variables used in the construction of the PCA. 
- rather, we are passing a novel set of variables | traits | phenotypes that will be correlated to the PC-axis (1 and 2) and then plotted. 
- as an example:
	1. generate a prcomp object:
		> pca = prcomp( iris[, 1:4] )
	
	2. or you can generate a probabilistic pca (pcaRes) object
		> pca = ppca( as.matrix( iris[, 1:4] ), nPCs = 4)
	
	3. run the function
		> moose_biplot(PCA = pca, dataframe_of_phenotypes = iris[, 1:4], 
             plot_top_N_phenotypes = 3, 
             grouping1 = iris$Species, grouping1NAME = "species",
             grouping2 = iris$Species, grouping2NAME =  "species",
             scalearrows = FALSE )

	- the dataframe_of_phenotypes can be any matrix of quantitative trait with the same number of row as passed to the prcomp() or ppca() functions.
	- *grouping1* dictates the color scheme and the ellipses to be drawn **Currently limited to 9 groups**
	- *grouping2* dictates the plot shapes to be used. **Currently limited to 5 groups**
	- *scalearrows* allows you to scale the largest correlated trait to a rho of 1, and all other arrows in correspoinding manner. This can be done to aid in visualization. Note:  if *scalearrows* is set to TRUE, the relative length of the arrows remain informative.
	