##!/usr/bin/env Rscript
args<-commandArgs(TRUE)

#With the read.tree function from the ape package, read trees from external file and save as a list called trees.
library(ape)

trees <- ape::read.tree(args[1])

#(optional) Read or get gene names somewhere (same order as the trees) and save it as a vector called names.

names <- readLines(args[2])

#Run phylter on your trees (see details below for possible options).
if (!requireNamespace("phylter", quietly = TRUE))
  install.packages("phylter")
library("phylter")

results <- phylter(trees, gene.names = names)
#results <- phylter(trees)

#To get the list of outliers detected by phylter, simply type:

#results$Final$Outliers


#In addition, many functions allow looking at the outliers detected and comparing before and after phyltering.
# Get a summary: nb of outliers, gain in concordance, etc.
#summary(results)
# Show the number of species in each gene, and how many per gene are outliers
#plot(results, "genes") 
# Show the number of genes where each species is found, and how many are outliers
#plot(results, "species") 
# Compare before and after genes x species matrices, highlighting missing data and outliers 
# identified (not efficient for large datasets)
#plot2WR(results) 
# Plot the dispersion of data before and after outlier removal. One dot represents one 
# gene x species association
#plotDispersion(results) 
# Plot the genes x genes matrix showing pairwise correlation between genes
#plotRV(results) 
# Plot optimization scores during optimization
#plotopti(results) 

#Save the results of the analysis to an external file, for example to perform cleaning on raw alignments or pruning gene trees based on the results from phylter.
write.phylter(results, file = args[3])
