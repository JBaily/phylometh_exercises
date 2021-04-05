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

save_me <- function(thing, name){
  
  temp <- thing
  #name <- deparse(quote(thing))
  filename <- paste0("results/", name)
  print(filename)
  sink(filename)
  print(temp)
  sink()
  
}

save_jpeg <- function(pic, name){
  
  filename = paste0("results/", name)
  jpeg(file=filename)
  plot(pic)
  dev.off()
  
}

change_tips <- function(tree, new_labels){
  
  temp <- tree
  temp$node.label <- new_labels
  return(temp)
  
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

is.better <- function(rate1, name1, delta1, rate2, name2, delta2){
  
  if(rate1 > rate2){
    
    return(paste0(name2, " is better than ", name1, ",with a support of a delta",
                  " AIC of ", delta1))
    
  }
  
  else{
    
    return(paste0(name1, " is better than ", name2, ",with a support of a delta",
                  " AIC of ", delta2))
    
  }
  
}

evo_rate <- function(tree, data, model_1) {
 
  BM1 <- geiger::fitContinuous(tree, data, model = model_1, ncores=1)
  
  filename = paste0("results/model_",model_1)
  sink(filename)
  print(BM1)
  sink()
  
  return(BM1)
}

ouwie_data <- function(cont_data, regime_data){
  
  names <- name_changer(cont_data, "species")
  regime = swap_orient(regime_data)
  temp2 <- data.frame("species"= names[,1], "regime"=regime[,1],
                     "continuous"=names[,2])
  print(temp2)
  return(temp2)
  
}

list_names <- function(models, values){
  
  for(i in 1:length(models)){
    names(values)[i] <- models[[i]]
  }
  
  return(values)
}

all_OUwie <- function(models, phylo, data_form){
  
  #mod <- models[[1]]
  blah2 <- OUwie::OUwie(phy=phylo, data=data_form, model=models[[1]],
                simmap.tree=FALSE, diagn=FALSE)
  blah <- list(blah2)
  
  for(i in 2:length(models)){
    
    product <- OUwie(phy=phylo, data=data_form, model=models[[i]],
                     simmap.tree=FALSE, diagn=FALSE)
    
    filename = paste0("results/OUwie_model_",models[[i]])
    sink(filename)
    print(product)
    sink()
    
    blah[[length(blah)+1]] <- product
  }
  for(i in 1:length(models)){
    names(blah)[i] <- models[[i]]
  }
  
  print(blah)
  return(blah)
}

mass_likelihood <- function(tree, data, model, alpha, likelihood, best_mod){
for (i in sequence(length(alpha))) {
  likelihood[i] <-  OUwie.fixed(tree, data, model=model, 
                alpha=rep(alpha[i],2), 
                sigma.sq=best_mod$solution[2,],
                theta=best_mod$theta[,1])$loglik
}
  return(likelihood)
}

plot_alpha <- function(x, y, xlim, ylim, best){
  
  plot(x=x, y=y, xlab="alpha values", ylab="likelihood values", 
             type="l", bty="n", xlim=xlim, ylim=ylim)
  points(x=best$solution[1,1], y=best$loglik, pch=16, col="red")
  text(x=best$solution[1,1], y=best$loglik, "unconstrained best", 
       pos=4, col="red")
  
  abline(h= best$loglik - 2, lty="dotted") #Two log-likelihood

  jpeg(filename="results/alpha_vs_likelihood")
  plot(x=x, y=y, xlab="alpha values", ylab="likelihood values", 
       type="l", bty="n", xlim=xlim, ylim=ylim)
  points(x=best$solution[1,1], y=best$loglik, pch=16, col="red")
  text(x=best$solution[1,1], y=best$loglik, "unconstrained best", 
       pos=4, col="red")
  abline(h= best$loglik - 2, lty="dotted")
  dev.off()
} 

mass_likelihood_theta <- function(tree, data, model, reps, best, theta1, theta2, 
                                  likelihood){
  for (i in sequence(reps)) {
    likelihood[i] <-  OUwie.fixed(tree, data, model=model, 
                                  alpha=best$solution[1,], 
                                  sigma.sq=best$solution[2,],
                                  theta=c(theta1[i],
                                          theta2[i]))$loglik
  }
  print(likelihood)
  return(likelihood)
  
}

inter_graph <- function(points_inter, theta1, theta2, best, Reg, trait){
  
  contour(points_inter, xlim=range(c(theta1, theta2)), 
          ylim=range(c(theta1, theta2)), xlab="Theta 1", ylab="Theta 2", 
          levels=c(2,5,10), add=FALSE,lwd=1, bty="n", asp=1)
  
  points(x=best$theta[1,1], y=best$theta[2,1], col="red", pch=16)
  points(x=trait[which(Reg=="tentacles")],y=rep(min(c(theta1, theta2)), 
                                                length(which(Reg=="tentacles"))),
         pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 1, plotted along x axis
  
  points(y=trait[which(Reg=="no_tentacles")],x=rep(min(c(theta1, theta2)), 
                                                   length(which(Reg=="no_tentacles"))), 
         pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 2, plotted along y axis
  
  jpeg(filename="results/interpolated_graph")
  contour(points_inter, xlim=range(c(theta1, theta2)), 
          ylim=range(c(theta1, theta2)), xlab="Theta 1", ylab="Theta 2", 
          levels=c(2,5,10), add=FALSE,lwd=1, bty="n", asp=1)
  points(x=best$theta[1,1], y=best$theta[2,1], col="red", pch=16)
  points(x=trait[which(Reg=="tentacles")],y=rep(min(c(theta1, theta2)), 
                                                length(which(Reg=="tentacles"))),
         pch=18, col=rgb(0,0,0,.3)) 
  points(y=trait[which(Reg=="no_tentacles")],x=rep(min(c(theta1, theta2)), 
                                                   length(which(Reg=="no_tentacles"))), 
         pch=18, col=rgb(0,0,0,.3))
  dev.off()

}

name_it <- function(Vector, Names){
  
  temp <- Vector
  names(temp) <- Names  
  return(temp)
  print(temp)
}

save_me_general <- function(thing, name){
  
  filename <- paste0("results/",name)
  sink(file=filename)
  print(thing)
  sink()
}


