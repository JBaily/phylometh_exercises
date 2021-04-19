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
  plot(tree, type = "phylogram", show.tip.label=label, edge.width = 0.1)
  dev.off()
}

image_save <- function(image, file){
  
  jpeg(filename = file)
  image
  dev.off()
}


compare.multi.trees <- function(trees_1, trees_2, max.depth, y_start, name_1, name_2){
  
  plot(x=c(0, -1*max.depth), y=c(y_start, ape::Ntip(trees_2[[1]])), 
       log="y", type="n", bty="n", xlab="Time", ylab="N")
  
  colors=c(rgb(1,0,0,0.5), rgb(0, 0, 0, 0.5))
  list.of.both <- list(trees_1, trees_2)
  for (i in sequence(2)) {
    tree.list <- list.of.both[[i]]
    for (j in sequence(length(tree.list))) {
      ape::ltt.lines(tree.list[[j]], col=colors[[i]])
    }
  }
  legend("topleft", legend=c(name_1, name_2), fill=colors)
  
}

sim.converted <- function(sim.dat){
  
  temp <- sim.dat
  
  temp[temp[,2]==3,2] = 1
  temp[temp[,2]==4,2] = 2
  temp[,2] = temp[,2] - 1
  
  return(temp)
  
}

edit.rates <- function(rates, to.change, value){
  
  temp <- rates 
  change <- to.change
  temp[change] <- value
  
  return(temp)
  
  
}

make.hisse.list <- function(){
  
  
  hisse.list = list()
  load("data/testrecon1.rda")
  hisse.list[[1]] = pp.recon
  load("data/testrecon2.rda")
  hisse.list[[2]] = pp.recon
  load("data/testrecon3.rda")
  hisse.list[[3]] = pp.recon
  
  return(hisse.list)
  
}

set.node.labels <- function(tree, check, places, new.node){
  
  temp <- tree
  
  if (check == TRUE){
    
    for (i in 1:length(new.node)){
      
      temp$node.label[[i]] <- places[[i]]
      temp$node.state[[i]] <- new.node[[i]]
      
    }
    
  }
  
  else {
    
    blah <- new.node
    temp$node.label <- places
    temp$node.state <- new.node
    
  }
  
  return(temp)
  
}







