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
                                file_3="results/trait1_simmap", 
                                model="ER", nstates=2),
  
  #This matrix refers to the rate at which the traits change. 
  # In this instance, it means that the rate of going from state 1 to state 2 
  # and vice-versa are identical. 
  
  testing = as.data.frame(name_changer_2(dis_mat_temp$trait, "Trait1")),
  test_named = col_nam_frame("Trait1", 2, testing),
  pp.er = get_corHMM(dis_mat_temp$tree, test_named, 
                     dis_mat_temp$rate, type="marginal"),
  print(pp.er),
  save_me_general(pp.er, "pp_er_rates"),
  
  #This matrix refers to the optimized transition rates that underlie the 
  #evolution of this binary character, using the general rate matrix estimated
  #by the ER model. 
  
  dis_mat_temp_ARD = get_dis_mat_p1(primates$trait, primates$tree, 
                                file_1="results/trait_1_ard.txt", 
                                file_2="results/simmap_info_ard.txt",
                                file_3="results/trait1_simmap_ard", 
                                model="ARD", nstates=2),
  pp.ard = get_corHMM(dis_mat_temp_ARD$tree, test_named, 
                     dis_mat_temp_ARD$rate, type="marginal"),
  print(pp.ard),
  save_me_general(pp.ard, "pp_ard_rates"),
  
  #This matrix refers to the optimized transition rates that underlie the
  #evolution of this binary charaacter, using the general rate matrix estimated
  #by the ARD model. 
  
  #In terms of which one is better, the ARD model is SLIGHTLY better, with a 
  #SLIGHTLY higher loglik value (ARD:-23.4031 > ER: -23.41535).
  
  dis_mat_ER_4state = get_dis_mat_p1(primates$trait, primates$tree, 
                                    file_1="results/trait_1_er_4.txt", 
                                    file_2="results/simmap_info_er_4.txt",
                                    file_3="results/trait1_simmap_er_4", 
                                    model="ER", nstates=4),
 
  rate.mat_ER_4state = remove_4(dis_mat_ER_4state$rate),
  print(rate.mat_ER_4state),
  
  fourstate.trait = rep(NA,Ntip(primates$tree)),
  fourstate.trait.assigned = assign_four(fourstate.trait, primates$tree, 
                                         primates$trait),
  fourstate.data = data.frame(Genus_sp=primates$trait[,1], 
                              T1=fourstate.trait.assigned),
  save_me_general(fourstate.data, "fourstate_data"),
  
  four_ray_er = rayDISC(primates$tree, fourstate.data, ntraits=1, model="ER", 
                node.states="marginal"), 
  ray_ard_four = rayDISC(primates$tree, fourstate.data, ntraits=1, 
                        rate.mat=rate.mat_ER_4state, node.states="marginal",
                        model = "ARD"),
  save_me_general(four_ray_er, "Four_state_ER_rates_rayDISC"),
  save_me_general(ray_ard_four, "Four_state_ARD_rates_rayDISC"),
  
  dis_mat_ARD_4state = get_dis_mat_p1(primates$trait, primates$tree, 
                                     file_1="results/trait_1_er_4.txt", 
                                     file_2="results/simmap_info_er_4.txt",
                                     file_3="results/trait1_simmap_er_4", 
                                     model="ARD", nstates=4),
  
  gtr.4state = gtr_conv(dis_mat_ARD_4state$rate),
  gtr.4state.edit = edit_value(remove_4(gtr.4state), 3, 2, 3),
  print(gtr.4state.edit),
  
  
  four_ray_gtr_ard = rayDISC(primates$tree, fourstate.data, ntraits=1, 
                             rate.mat= gtr.4state.edit, model="ARD",
                             node.states="marginal"),
  
  save_me_general(four_ray_gtr_ard, "four_rayDISC_ARD_gtr"),
  
  dis_mat_pag94 = get_dis_mat_p2(2, nstates=2, model="ARD"),
  rate.mat.pag94 = rate_pagel(dis_mat_ARD_4state$rate, c(3,5,8,10)),
  
  ##Route 1###
  #Note:  I will be using the "4" (actually 3) state data from the gtr part. 
  
  #Construct a model to test if state 1 can never be lost
  #Four state data where "1" = 0 = 00, "2" = 1 = 01, "3" = 2 = 10, 
  # and "4" = 3 = 11. 
  #Three state data where "1" = 0 = 00, "2" = 1 = 01, and "3" = 2 = 11.
  #I am going to assume that by "state 1", this sub-question means the 1 
  #in either trait, so X1 or 1X. 
  #Model should look like this:
  #   1   2   3   4      or this   1  2   3 , if rayDISC is being difficult.
  # 1 NA  4   7   10            1 NA  3   5
  # 2 NA  NA  NA  11            2 NA  NA  6
  # 3 NA  NA  NA  12            3 NA  NA  NA
  # 4 NA  NA  NA  NA
  
  #Remove c(1,2,3,5,6,8,9), or c(1,2,4)
  #First version
  
  rate.1.never.4D = rate_pagel(dis_mat_ARD_4state$rate, c(1,2,3,5,6,8,9)),
  
  #Second version
  rate_3D = remove_4(dis_mat_ARD_4state$rate),
  rate.2.temp.3D = rate_change_terms(rate_3D, rows = c(2,3,3), 
                                      columns = c(1,1,2), c(NA,NA,NA)),
  rate.2.never.3D = rate_change_terms(rate.2.temp.3D, rows = c(1,1,2), 
                                      columns = c(2,3,3), c(1,2,3)),
  
  #Experiment with the effects of frequencies at the root.
  #I will use the rayDISC function with all of the same inputs as 
  # "four_ray_gtr_ard", but with different values for the "root.p"
  #variable. I will try "yang", "maddfitz", in addition to c(0.33,0.33,0.33),
  # c(0.5, 0.25, 0.25), c(0.25, 0.5. 0.25), c(1,0,0), c(0,1,0), c(0,0,1).  
  
  root.freq.33x3 = rayDISC(primates$tree, fourstate.data, ntraits=1, 
                             rate.mat= gtr.4state.edit, model="ARD",
                             node.states="marginal", 
                           root.p = c(0.33,0.33,0.33)),
  
  root.freq.5.25.25 = rayDISC(primates$tree, fourstate.data, ntraits=1, 
                           rate.mat= gtr.4state.edit, model="ARD",
                           node.states="marginal", 
                           root.p = c(0.5, 0.25, 0.25)),
  
  root.freq.25.5.25 = rayDISC(primates$tree, fourstate.data, ntraits=1, 
                           rate.mat= gtr.4state.edit, model="ARD",
                           node.states="marginal", 
                           root.p = c(0.25, 0.5, 0.25)),
  
  root.freq.1.0.0 = rayDISC(primates$tree, fourstate.data, ntraits=1, 
                           rate.mat= gtr.4state.edit, model="ARD",
                           node.states="marginal", root.p = c(1,0,0)),
  
  root.freq.0.1.0 = rayDISC(primates$tree, fourstate.data, ntraits=1, 
                           rate.mat= gtr.4state.edit, model="ARD",
                           node.states="marginal", 
                           root.p = c(0,1,0)),
  
  root.freq.0.0.1 = rayDISC(primates$tree, fourstate.data, ntraits=1, 
                            rate.mat= gtr.4state.edit, model="ARD",
                            node.states="marginal", 
                            root.p = c(0,0,1)),
  
  root.freq.yang = rayDISC(primates$tree, fourstate.data, ntraits=1, 
                            rate.mat= gtr.4state.edit, model="ARD",
                            node.states="marginal", 
                            root.p = "yang"),
  
  root.freq.maddfitz = rayDISC(primates$tree, fourstate.data, ntraits=1, 
                            rate.mat= gtr.4state.edit, model="ARD",
                            node.states="marginal", 
                            root.p = "maddfitz"),
  
  log_lik_value.33x3 = paste0("The log likelihood of the 33x3 root model",
                              " is ", root.freq.33x3$loglik),
  log_lik_value.5.25.25 = paste0("The log likelihood of the 5.25.25 root model",
                              " is ", root.freq.5.25.25$loglik),
  log_lik_value.25.5.25 = paste0("The log likelihood of the 25.5.25 root model",
                              " is ", root.freq.25.5.25$loglik),
  log_lik_value.1.0.0 = paste0("The log likelihood of the 1.0.0 root model",
                              " is ", root.freq.1.0.0$loglik),
  log_lik_value.0.1.0 = paste0("The log likelihood of the 0.1.0 root model",
                              " is ", root.freq.0.1.0$loglik),
  log_lik_value.0.0.1 = paste0("The log likelihood of the 0.0.1 root model",
                              " is ", root.freq.0.0.1$loglik),
  log_lik_value.yang = paste0("The log likelihood of the yang root model",
                               " is ", root.freq.yang$loglik),
  log_lik_value.maddfitz = paste0("The log likelihood of the maddfitz root model",
                               " is ", root.freq.maddfitz$loglik),
  
  print(log_lik_value.33x3),
  print(log_lik_value.5.25.25),
  print(log_lik_value.25.5.25),
  print(log_lik_value.1.0.0),
  print(log_lik_value.0.1.0),
  print(log_lik_value.0.0.1),
  print(log_lik_value.yang),
  print(log_lik_value.maddfitz),
  
  #It appears that fixing the root probabilities, at least with this tree
  #and dataset, does change the overall log likelihood of a model by a 
  #substantial degree. Some of the models are better than others though, 
  #with order from best to worst being 0.1.0 > 0.0.1 > maddfitz > 1.0.0 
  # > 25.5.25 > 5.25.25 > 33x3 > yang. 
  # I'm honestly a tad surprised the model prefers fixing 01 or 11 at the 
  # root instead of 00, as we tend to think of 00 being a preferred root
  # state, but upon further reflection I have no idea what this trait 
  # data actually represents, so I have no grounds for asserting that 00
  # is the most likely root state here. 
  # It is curious that all of them have higher log likelihoods than the 
  # regular "four_ray_gtr_ard" model, which I have no good explanation for.
  
  
  
  #Create and use a model to see if transitions from 00 go to 
  # 11 only via 01.
  
  # I will be using the 3 state model instead of the 4 state model to avoid
  # the issues that occur since the 10 state is not observed in the data. 
  # So, in this case, I am assuming that the limitation is only imposed
  # going 00 -> 11, not the reverse. 
  #   1   2   3
  # 1 NA  3   NA
  # 2 1   NA  5
  # 3 2   4   NA
  
  #Making model
  rate_test = remove_4(dis_mat_ARD_4state$rate),
  rate.temp.test = rate_change_terms(rate_test, rows = c(1), 
                                     columns = c(3), c(NA)),
  rate.never.test = rate_change_terms(rate.temp.test, rows = c(1,3,2), 
                                      columns = c(2,2,3), c(3,4,5)),
  
  #Testing model 
  fixed_model_test = rayDISC(primates$tree, fourstate.data, ntraits=1, 
                             rate.mat= rate.never.test, model="ARD",
                             node.states="marginal"),
  print(paste0("The log likelihood of the usual model is ",
               four_ray_gtr_ard$loglik," and the log likelihood of the ",
               "fixed transition path model is ",fixed_model_test$loglik))
  
  #So, while the usual model is still better in this particular case, I
  #would still say that it's within the realm of possibility that the trait
  #goes through the fixed transition path from 00 to 11 as the difference
  #between the two models is only ~0.2 loglik units. 
  
  )