
library(adephylo)
library(ape)
library(ctc)      # installation via bioconductor
library(phangorn)
library(Rfast)
library(vegan)


#----------------------------------------------------------------------------
# FUNCTION TO COMBINE THE PREDICTIONS OF THREE CNNS,
# AND TURN THEM INTO A TREE FILE

# REQUIRES FIVE FILES IN d:
# the averaged predictions of the three CNNs:
# -> "family_predictions_averaged_by_family.csv"
# -> "order_predictions_averaged_by_family.csv"
# -> "clade_predictions_averaged_by_family.csv"

# the class names of the family-CNN:
# -> "classnames_family.txt"

# a list of the orders and subclasses to which the families belong:
# -> "taxonomy.dat"
#----------------------------------------------------------------------------
make_ensemble_tree <- function(d)
{ 
  #__ load prediction files
  preds_family <- read.csv(paste(d, "family_predictions_averaged_by_family.csv", sep = ""),row.names=1)
  preds_order <- read.csv(paste(d, "order_predictions_averaged_by_family.csv", sep = ""),row.names=1)
  preds_clade <- read.csv(paste(d, "clade_predictions_averaged_by_family.csv", sep = ""),row.names=1)
  
  #__ load class names of the family-CNN
  class_names <- readLines(paste(d,"classnames_family.txt", sep = ""))
  no_classes <- length(class_names)
  
  #__ get taxonomy data for links between clades & orders & families
  taxa <- read.csv(paste(d,"taxonomy.dat",sep = ""))
  
  # define weights of order and subclass predictions
  order_frac = 1
  clade_frac = 1
  
  # make look-up table for family-subclass
  clade_names <- colnames(preds_clade)
  family_clade <- unique(taxa[c("family","clade")])
  all_fams_of_clade <- lapply(clade_names, function(z) family_clade[which(family_clade$clade==z),])
  fam_of_clade <- lapply(all_fams_of_clade, function(z) z$family[which(z$family %in% class_names)])
  
  # make look-up table for family-order
  order_names <- colnames(preds_order)
  family_order <- unique(taxa[c("family","order")])
  all_fams_of_order <- lapply(order_names, function(z) family_order[which(family_order$order==z),])
  fam_of_order <- lapply(all_fams_of_order, function(z) z$family[which(z$family %in% class_names)])
  
  #__ combine predictions
  # predictions of the family-CNN as basis for the combined predictions
  combined_preds <- preds_family
  
  # loop through the rows with the predictions
  for (i in 1:nrow(combined_preds)){
    # get all predictions of current row
    cp <- combined_preds[i,]
    
    #__ add predictions of each family of each order
    for (j in 1:length(fam_of_order)){
      # get indices of families of current order
      y <- .subset2(fam_of_order,j)
      
      # update predictions of current image with the predictions of these families
      cp[y] <- cp[y] + (.subset2(preds_order,j)[i] * order_frac)
    }
    
    #__ add predictions of each family of each clade (as above)
    for (j in 1:length(fam_of_clade)){
      y <- .subset2(fam_of_clade,j)
      cp[y] <- cp[y] + (.subset2(preds_clade,j)[i] * clade_frac)
    }
    
    # normalize to sum(combined_preds) = 1
    combined_preds[i,] <- cp / sum(cp)
  }
  
  #__ make tree from combined predictions
  # turn predictions into distance matrix
  dist.mat <- vegdist(combined_preds)
  
  # hierarchical clustering of the distance matrix
  hc <- hclust(dist.mat)
  
  # turn cluster into tree
  tree <- as.phylo(hc)
  
  return(tree)
}
