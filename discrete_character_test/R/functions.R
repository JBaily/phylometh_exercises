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

VisualizeData <- function(phy, data, file_1, file_2) {
  
  plot(phy, type = "fan", show.tip.label=TRUE, edge.width = 0.1)
  newdata <- data
  print(newdata)
  
  plot_tree(phy, file_1, TRUE)
  write.csv(newdata, file = file_2)
  
}

plot_tree <- function(tree, file, label) {
  pdf(file=file)
  plot(tree, type = "fan", show.tip.label=label, edge.width = 0.1)
  dev.off()
}

plot_sim_tree <- function(simm, label) {
  filename = paste0("results/SimmMap_",label)
  pdf(file=filename)
  phytools::plotSimmap(simm[[1]], fsize = 0.5)
  dev.off()
}

print_plotAnc <- function(phyDat, recon, type) {
  
  for(i in 1:6){
    
    phangorn::plotAnc(phyDat, recon, i)
    
    filename = paste0("results/",type,"_",i)
    pdf(file=filename)
    phangorn::plotAnc(phyDat, recon, i)
    dev.off()
  } 

}

print_rates <- function(rates, variable) {
  
  filename = paste0("results/transition_rates_",variable)
  sink(file=filename)
  print(rates)
  sink(file=NULL)
  print("Done writing to file")
  
}

replace_NA <- function(blah) {
  
  model = blah
  temp_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
  temp_matrix[[1,2]] <- model[[1,2]]
  print(temp_matrix[[1,2]])
  temp_matrix[[2,1]] <- model[[2,1]]
  print(temp_matrix[[2,1]])
  temp_matrix[[1,1]] <- 0
  temp_matrix[[2,2]] <- 0
  print(temp_matrix)
}

remove_cols <- function(table, to_remove){
  
  table %>%
  select(remove[[1]])
  mutate(
    remove[1,] = NULL
  )
  print(table)
}