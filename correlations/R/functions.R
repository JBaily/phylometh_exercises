plot_tree <- function(tree, file, label) {
  pdf(file=file)
  plot(tree, type = "fan", show.tip.label=label, edge.width = 0.1)
  dev.off()
}

set_rows_data <- function(data, names){
  
  temp <- data
  rownames(temp) <- names
  print(temp)
  
}

set_names_data <- function(data, names){
  
  temp <- data
  names(temp) <- names
  print(temp)
  
  file_name <- paste0("data/",data,"_dataset")
  pdf(file=file_name)
  temp
  dev.off()
  
}

print_info <- function(thing, name){
  
  file_name <- paste0("results/",name,".txt")
  sink(file_name)
  print(thing)
  sink()
  
}

new_table <- function(x_values, y_values){
 
  x <- x_values[,1] 
  y <- y_values[,1]
  
  temp <- data.frame(x,y)
}

col_nam_frame <- function(name, place, frame){
  
  temp <- frame
  names(temp)[place] <- name
  return(temp)
  
}

name_null <- function(PIC){
  
  temp <- PIC
  names(temp) <- NULL
  print(temp)
  
}

name_changer_2 <- function(data, type){
  
  new <- as.data.frame(data)
  new[,2] <- new[,1]
  new[,1]<- rownames(new)
  rownames(new) <- NULL
  return(new)
}

get_dis_mat_p1 <- function(dis_mat, tree, file_1, file_2, file_3){
  
  dis_mat[which(grepl("Hylobates", dis_mat[,1])),2]<-1
  trait1<-dis_mat[,2]
  names(trait1)<-dis_mat[,1]
  print(trait1)
  
  sink(file_1)
  print(trait1)
  sink()
  
  tree_1 <- ape::multi2di(tree)
  #return(list(trait1, tree_1))
  plotSimmap(make.simmap(tree_1, trait1), pts=FALSE, fsize=0.8)
  
  sink(file_2)
  make.simmap(tree_1, trait1)
  sink()
  pdf(file=file_3)
  plotSimmap(make.simmap(tree_1, trait1), pts=FALSE, fsize=0.8)
  dev.off()
  
  rate.mat.er<-corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=2, 
                              model="ER")
  sink(file_2)
  print(rate.mat.er)
  sink()
  
  return(list("trait"=trait1, "tree"=tree_1, "rate"=rate.mat.er))
}
 
plot_with_line <- function(x, y, line, file_name){
  
  plot(x,y)
  abline(line, col="blue")
  
  pdf(file=file_name)
  plot(x,y)
  abline(line, col="blue")
  dev.off()
}

get_corHMM <- function(tree, trait, rate, type) {
  
  temp<-corHMM(tree,trait[,c(1,2)],
                rate.cat=1,rate.mat=rate,node.states=type)
  print(temp)
  
  file_name = paste0("results/corHMM_",type,"_",names(trait)[2])
  sink(file_name)
  print(temp)
  sink()
  
}