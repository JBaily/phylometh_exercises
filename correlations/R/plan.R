plan <- drake_plan(
  
  #Continuous characters
  tree.primates = read.tree(text="((((Homo:0.21,Pongo:0.21):0.28,Macaca:0.49)
                             :0.13,Ateles:0.62):0.38,Galago:1.00);"), 
  
  X = c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968),
  Y = c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259),
  
  X_temp = set_names_data(X, c("Homo", "Pongo", "Macaca", "Ateles", "Galago")),
  Y_temp = set_names_data(Y, c("Homo", "Pongo", "Macaca", "Ateles", "Galago")),
  
  pic.X = pic(X, tree.primates),
  pic.Y = pic(Y, tree.primates), 
  
  X_pos = t(as.data.frame(lapply(name_null(pic.X), abs))),
  Y_pos = t(as.data.frame(lapply(name_null(pic.Y), abs))),
  
  X_cor = as.data.frame(set_rows_data(X_pos, NULL)),
  Y_cor = as.data.frame(set_rows_data(Y_pos, NULL)),
  
  X_lin_reg = lm(Y_cor[,1] ~ 0 + ., X_cor),
  print_info(X_lin_reg, "linear_regression"),
  
  plot_table = as.data.frame(new_table(as.data.frame(X_cor), 
                                       as.data.frame(Y_cor))),

  plot_1 = plot_with_line(plot_table$x, plot_table$y, X_lin_reg, 
                          "results/cont_graph"),

  
  #Discrete characters - single
  #load("~/GitHub_files/phylometh_exercises/correlations/data/primates.rda"),
  list_prim = ls(primates),
  print(list_prim),
  print(primates),
  
  require(phytools),
  
  dis_mat_temp = get_dis_mat_p1(primates$trait, primates$tree, 
                                file_1="results/trait_1.txt", 
                                file_2="results/simmap_info.txt",
                                file_3="results/trait1_simmap"),
  
  #This matrix refers to the rate at which the traits change. 
  # In this instance, it means that the rate of going from state 1 to state 2 
  # and vice-versa are identical (1:1). 
  
  testing = as.data.frame(name_changer_2(dis_mat_temp$trait, "Trait1")),
  test_named = col_nam_frame("Trait1", 2, testing),
  pp.er = get_corHMM(dis_mat_temp$tree, test_named, 
                     dis_mat_temp$rate, type="marginal"),
  
  
  )