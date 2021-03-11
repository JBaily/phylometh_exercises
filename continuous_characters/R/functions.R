CleanData <- function(phy, data) {
  original_taxa <- row.names(phy)
  if_diffs <- name.check(phy, data)
  
  text_1 <- paste0("The species ", if_diffs$tree_not_data, " has been removed
                   as they are not in the data.")
  print(text_1)
  struc_file <- file("results/taxa_removed")
  writeLines(text_1, struc_file)
  close(struc_file)
  
  treedata(phy, data, sort = TRUE, warnings = TRUE)
  
}  

data_prune <- function(phy, data){
  
  if_diffs <- name.check(phy, data)
  to_remove <- if_diffs$data_not_tree

  if(identical(to_remove, character(0))==TRUE){
    
    print("There are no species in the data that are not in the tree")
    data
    
  }
  else{
    
    print("There are species in the data that are not in the tree")
    
    all <-rbind(to_remove,data) 
    all[!duplicated(all,fromLast = FALSE)&!duplicated(all,fromLast = TRUE),]
    all
  }
  
}

VisualizeData <- function(phy, data, file_1, file_2, file_3) {
  
  plot(phy, type = "fan", show.tip.label=TRUE, edge.width = 0.1)
  newdata <- data
  print(newdata)
  
  plot_tree(phy, file_1, TRUE)
  write.csv(newdata, file = file_2)

  
  phytools::contMap(phy, data, plot=TRUE, legend=TRUE)
  print("Visualize success")
  
  pdf(file=file_3)
  phytools::contMap(phy, data, plot=TRUE, legend=TRUE)
  dev.off()
}

name_changer <- function(data, var){
  
  data <- tibble::rownames_to_column(data, "SPECIES")
  colnames(data)[2] <- var
  print(data)
  
}

first_to_row <- function(data, var){
  
  temp <- data.frame(data[,-1], row.names=data[,1])
  colnames(temp)[1] <- var
  print(temp)
}

set_species_names <- function(charac_cont, species){
  
  traits <- charac_cont
  traits <- as.numeric(traits)
  print(traits)
  names(traits) <- species
  print(traits)
  
}

name_changer_2 <- function(data, names, type){
  
  new <- data
  row.names(new) <- names
  colnames(new)[1] <- c(type)
  return(new)
}

swap_orient <- function(data){
  
  temp <- as.data.frame(t(data))
  
}

plot_tree <- function(tree, file, label) {
  pdf(file=file)
  plot(tree, type = "fan", show.tip.label=label, edge.width = 0.1)
  dev.off()
}

evo_rate <- function(tree, data, model_1) {
 
  BM1 <- geiger::fitContinuous(tree, data, model = model_1)
  print(BM1)
}