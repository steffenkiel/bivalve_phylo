
library(adephylo)
library(ape)

#---------------------------------------------------------------------------
# TREE EVALUATION FUNCTIONS
#  to calculate the 'tree coherency index' and the 'node distance'-indices
#  for a given tree

# THESE FUNCTIONS REQUIRE:
# a file of type 'tree'; the tip names must:
#  (1) include those given in the "classnames_family.txt" file,
#  (2) all be included in the "taxonomy.dat" file.

# the class names of the family-CNN:
# -> "classnames_family.txt"

# list of the orders and subclasses to which the families belong:
# -> "taxonomy.dat"

# the two helper functions 'non_nested_subtrees' and 'is_subtree'
#---------------------------------------------------------------------------

#-------------------------------------------------
# FUNCTION: tree_coherency_index
#  returns the tree coherency index for 'tree'
#-------------------------------------------------
tree_coherency_index <- function(tree, taxa, class_names)
{
  # make df of unique "family"/"clade" values
  taxa <- unique(taxa[c("family","clade")])
  row.names(taxa) <- NULL
  
  # change predicted family names to their clade (=subclass) names
  res <- lapply(tree$tip.label, function(z) taxa$clade[which(taxa$family == z)])
  tree$tip.label <- unlist(res)
  
  # get subtrees
  st <- subtrees(tree)
  
  # get indices of the monophyletic subtrees
  no_clades <- unlist(lapply(st, function(z) length(unique(z$tip.label))))
  mono_subs <- which(no_clades==1)
  
  # number of tips of each subtree
  no_tips <- as.numeric(unlist(lapply(st, function(z) z$Ntip)))
  
  # number of non-nested, monophyletic subtrees
  nnst <- non_nested_subtrees(st)
  
  # tree coherency index
  no_tips_nnst <- sum(unlist(lapply(nnst, function(z) z$Ntip)))
  tree_coherency_index <- (no_tips_nnst / length(nnst)) * (no_tips_nnst / length(tree$tip.label))
  
  return(tree_coherency_index)
}


#------------------------------------------------------------------------------
# FUNCTION: node_distance_training
#  returns the node distances of each family in 'class_names' in the 'tree'
# INPUT:
#  tree <- the tree to be tested
#  class_names <- class names of the family-CNN (the training families)
#  taxa <- df with order and subclass for each family
# OUTPUT:
#  a df with information for each training family:
#    class_names <- its name
#    clade <- its subclass
#    order <- its order
#    clade_distance <- node distance to the nearest member of its subclass
#    order_distance <- node distance to the nearest member of its order
#------------------------------------------------------------------------------
node_distance_training <- function(tree, class_names, taxa)
{
  #___ make df with families in class_names and their orders and clades
  df <- data.frame(class_names)
  df$clade <- unlist(lapply(df$class_names, function(z) unique(taxa$clade[which(taxa$family==z)])))
  df$order <- unlist(lapply(df$class_names, function(z) unique(taxa$order[which(taxa$family==z)])))
  
  # add columns for results
  df$clade_distance <- NA
  df$order_distance <- NA
  
  # get subtrees
  st <- subtrees(tree)
  
  # loop through all training families
  for (i in 1:length(class_names)){
    
    # get all families of the subclass (clade) of the current family (excl itself)
    clade_fams <- unique(taxa$family[which(taxa$clade==.subset(df$clade,i))])
    clade_fams <- clade_fams[clade_fams %in% class_names]
    clade_fams <- clade_fams[clade_fams != .subset(class_names,i)]
    
    #__ get neighbor distance etc. if it's not a single-family-clade
    if (length(clade_fams) > 0){
      #__ get all subtrees with these and the current family
      wrk_sts <- vector()   # vector for indices of subtrees
      for (j in 1:length(clade_fams)){
        # make pair of current family and the i-th of its clade
        test_pair <- c(.subset(class_names,i), .subset(clade_fams,j))
        
        # get indices of the subtrees with this pair
        wrk_sts <- c(wrk_sts, which(unlist(lapply(st, function(z) all(test_pair %in% z$tip.label)))))
      }
      # reduce to unique
      wrk_sts <- unique(wrk_sts)
      
      # if subtrees have been found:
      if (length(wrk_sts) > 0){
        #___ get smallest of those subtrees
        nn <- which.min(unlist(lapply(st[wrk_sts], function(z) z$Nnode)))
        
        # extract subtree
        smallest_st <- st[[wrk_sts[nn]]]
        
        # identify clade members among its tip names
        tlsst <- smallest_st$tip.label
        kin <- tlsst[tlsst %in% clade_fams]
        
        # get node distances between current family and their kin
        x <- as.matrix(distTips(smallest_st, method = "nNodes"))
        
        # extract smallest value
        df$clade_distance[i] <- min(unlist(lapply(kin, function(z) x[.subset(class_names,i), z])))
      } else {
        # if no subtree with kins exists, use total number of nodes as distance
        df$clade_distance[i] <- tree$Nnode
      }
    }
    
    # AND NOW THE SAME FOR ORDER
    # get all trained families of order of the present family
    order_fams <- unique(taxa$family[which(taxa$order==df$order[i])])
    order_fams <- order_fams[order_fams %in% class_names]
    order_fams <- order_fams[order_fams != .subset(class_names,i)]
    
    #__ get neighbor distance etc. if it's not a single-family-order
    if (length(order_fams) > 0){
      #__ get all subtrees with these and the current family
      wrk_sts <- vector()   # vector for indices of subtrees
      for (j in 1:length(order_fams)){
        # make pair of current family and the i-th of its order
        test_pair <- c(.subset(class_names,i), .subset(order_fams,j))
        
        # get indices of the subtrees with this pair
        wrk_sts <- c(wrk_sts, which(unlist(lapply(st, function(z) all(test_pair %in% z$tip.label)))))
      }
      # reduce to unique
      wrk_sts <- unique(wrk_sts)
      
      # if subtrees have been found:
      if (length(wrk_sts) > 0){
        # get shortest of those subtrees
        nn <- which.min(unlist(lapply(st[wrk_sts], function(z) z$Nnode)))
        smallest_st <- st[[wrk_sts[nn]]]

        # identify order members among its tip names
        tlsst <- smallest_st$tip.label
        kin <- tlsst[tlsst %in% order_fams]
        
        # get node distances between current family and their kin
        x <- as.matrix(distTips(smallest_st, method = "nNodes"))
        
        # extract smallest value
        df$order_distance[i] <- min(unlist(lapply(kin, function(z) x[.subset(class_names,i), z])))
      } else {
        # if no subtree with kins exists, use total number of nodes as distance
        df$order_distance[i] <- tree$Nnode
      }
    }
  }
  return(df)
}


#------------------------------------------------------------------------------
# FUNCTION: node_distance_target
#  returns the node distances of each target family in the 'tree'
# INPUT:
#  tree <- the tree to be tested
#  class_names <- class names of the family-CNN (the training families)
#  taxa <- df with order and subclass for each family
#  target_fams <- the names of the target families
# OUTPUT:
#  a df with information for each target family:
#    class_names <- its name
#    clade <- its subclass
#    order <- its order
#    clade_distance <- node distance to the nearest member of its subclass
#    order_distance <- node distance to the nearest member of its order
#------------------------------------------------------------------------------
node_distance_target <- function(tree, class_names, taxa, target_fams)
{
  #___ make df with target_fams and their orders and clades
  df <- data.frame(target_fams)
  df$clade <- unlist(lapply(df$target_fams, function(z) unique(taxa$clade[which(taxa$family==z)])))
  df$order <- unlist(lapply(df$target_fams, function(z) unique(taxa$order[which(taxa$family==z)])))
  
  # add columns for results
  df$clade_distance <- NA
  df$order_distance <- NA
  
  st <- subtrees(tree)
  
  for (i in 1:length(target_fams)){
    
    # get all trained families of the clade of the current not_trained_family
    clade_fams <- unique(taxa$family[which(taxa$clade==.subset(df$clade,i))])
    clade_fams <- clade_fams[clade_fams %in% class_names]
    
    #__ get neighbor distance etc. if there are families in clade_fams
    if (length(clade_fams) > 0){
      #__ get all subtrees with these and the current family
      wrk_sts <- vector()   # vector for indices of subtrees
      for (j in 1:length(clade_fams)){
        # make pair of current family and the i-th of its clade
        test_pair <- c(.subset(target_fams,i), .subset(clade_fams,j))
        
        # get indices of the subtrees with this pair
        wrk_sts <- c(wrk_sts, which(unlist(lapply(st, function(z) all(test_pair %in% z$tip.label)))))
      }
      # reduce to unique
      wrk_sts <- unique(wrk_sts)
      
      #___ get smallest of those subtrees
      nn <- which.min(unlist(lapply(st[wrk_sts], function(z) z$Nnode)))
      
      # extract subtree
      smallest_st <- st[[wrk_sts[nn]]]
      
      # identify clade members among its tip names
      tlsst <- smallest_st$tip.label
      kin <- tlsst[tlsst %in% clade_fams]
      
      # get node distances between current family and their kin
      x <- as.matrix(distTips(smallest_st, method = "nNodes"))
      
      # extract smallest value
      df$clade_distance[i] <- min(unlist(lapply(kin, function(z) x[target_fams[i], z])))
    }
    
    # AND NOW THE SAME FOR ORDER
    # get all trained families of order of the present family
    order_fams <- unique(taxa$family[which(taxa$order==df$order[i])])
    order_fams <- order_fams[order_fams %in% class_names]
    
    #__ get neighbor distance etc. if there are families in order_fams
    if (length(order_fams) > 0){
      #__ get all subtrees with these and the current family
      wrk_sts <- vector()   # vector for indices of subtrees
      for (j in 1:length(order_fams)){
        # make pair of current family and the i-th of its order
        test_pair <- c(.subset(target_fams,i), .subset(order_fams,j))
        
        # get indices of the subtrees with this pair
        wrk_sts <- c(wrk_sts, which(unlist(lapply(st, function(z) all(test_pair %in% z$tip.label)))))
      }
      # reduce to unique
      wrk_sts <- unique(wrk_sts)
      
      # get shortest of those subtrees
      nn <- which.min(unlist(lapply(st[wrk_sts], function(z) z$Nnode)))
      smallest_st <- st[[wrk_sts[nn]]]

      # identify order members among its tip names
      tlsst <- smallest_st$tip.label
      kin <- tlsst[tlsst %in% order_fams]
      
      # get node distances between current family and their kin
      x <- as.matrix(distTips(smallest_st, method = "nNodes"))
      
      # extract smallest value
      df$order_distance[i] <- min(unlist(lapply(kin, function(z) x[target_fams[i], z])))
    }
  }
  return(df)
}


#___________________________________________
#            helper functions

#------------------------------------------------------------
# takes a list of subtrees (from ape's subtree function)
# and returns a list of non-nested, monophyletic subtrees
#------------------------------------------------------------
non_nested_subtrees <- function(st){
  # identify monophyletic subtrees
  no_clades <- unlist(lapply(st, function(z) length(unique(z$tip.label))))
  mono_subs <- as.numeric(which(no_clades==1))
  
  # get subset of monophyletic subtrees
  mono_st <- lapply(mono_subs,function(z) .subset2(st,z))
  mono_index <- seq(length(mono_st))
  
  # monophyletic subtrees that are nested within any of the others
  nested_trees <- unlist(lapply(mono_index, function(z) is_subtree(mono_st,z)))
  
  # get indices of the non-nested trees
  non_nested <- which(!nested_trees)
  
  # finally: a list of the non-nested trees
  out <- lapply(non_nested,function(z) .subset2(mono_st,z))
  
  return(out)
}


#----------------------------------------------------------------------
# check if the i-th element of the tree list is a subtree
# i.e., if its node label(s) is/are elements of any of the other trees
#----------------------------------------------------------------------
is_subtree <- function(st,i){
  num_trees <- length(st)
  if (i > num_trees) stop(print("Index out of bounds."))
  
  # make vector with the indices of the trees to search
  tree_index <- seq(num_trees)
  tree_index <- tree_index[!tree_index==i]
  
  # get node label(s) of the tree to test
  nl <- .subset2(st,i)$node.label
  
  # get node label(s) of all other trees
  other_nl <- lapply(tree_index, function(z) .subset2(st,z)$node.label)
  
  # is the node label(s) a subset of the node lables of the other trees?
  res <- unlist(lapply(other_nl, function(z) all(nl %in% z)))
  
  return(TRUE %in% res)
}


