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

get_dis_mat_p1 <- function(dis_mat, tree, file_1, file_2, file_3, model, nstates){
  
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
  
  rate.mat <- corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, 
                                      nstates=nstates, model=model)
  
  meep = list("trait"=trait1, "tree"=tree_1, "rate"=rate.mat)
  return(meep)
}
 
get_dis_mat_p2 <- function(ntraits, nstates, model){

  
  rate.mat <- corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=ntraits, 
                                      model = model, nstates = nstates)
  print(rate.mat)
  return(rate.mat)
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
  
  return(temp)
  
}

save_me_general <- function(thing, name){
  
  temp <- thing
  filename = paste0("results/",name)
  sink(file = filename)
  print(temp)
  sink()
  
}

assign_four <- function(state, tree, trait){
  
  for(i in sequence(Ntip(tree))) {
    if(trait[i,2]==0 && trait[i,3]==0) {
      state[i]<-0
    }
    if(trait[i,2]==0 && trait[i,3]==1) {
      state[i]<-1
    }
    if(trait[i,2]==1 && trait[i,3]==0) {
      state[i]<-2
    }
    if(trait[i,2]==1 && trait[i,3]==1) {
      state[i]<-3
    }
  }
  
  return(state)
}

remove_4 <- function(rate){
  
  temp <- rate
  temp <-temp[-4,]
  temp <- temp[,-4]
  
  return(temp)
}

edit_value <- function(rate, new_value, row, column){
  
  temp <- rate
  temp[row, column] <- new_value
  return(temp)
  
}

gtr_conv <- function(rate){
  
  gtr.4state <- rate
  
  gtr.4state <- corHMM:::rate.par.eq(rate, c(1,4))
  gtr.4state <- corHMM:::rate.par.eq(rate, c(2,6))
  gtr.4state <- corHMM:::rate.par.eq(rate, c(3,8))
  gtr.4state <- corHMM:::rate.par.eq(rate, c(4,6))
  gtr.4state <- corHMM:::rate.par.eq(rate, c(5,7))
  gtr.4state <- corHMM:::rate.par.eq(rate, c(6,7))

  return(gtr.4state)
  
}

add_zeroes <- function(matrix){
  temp <- matrix 
  for(i in 1:dim(temp)[1]){
    
    temp[i,i] <- 0
    
  }
  return(temp)
}

rate_pagel <- function(rate, drop_vector){
  
  temp <- corHMM:::rate.par.drop(rate, drop.par=drop_vector)
  
  print(temp)
  return(temp)
  
}

rate_change_terms <- function(rate, rows, columns, values){
  
  temp <- rate
  for(i in 1:length(rows)) {
    
    temp[rows[i], columns[i]] <- values[i]
    
  }
  
  print(temp)
  return(temp)
}