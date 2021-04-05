plan <- drake_plan(
  
  #Getting a tree and some data -- lets try something simple like cats.
  
  tree = read.tree("data/ot_1134.tre.txt"),
  print(tree),
  tree_image = plot_tree(tree, file=file_out("results/initial_tree"), FALSE),
  cont.data = read.csv(file="data/taxa_binary.csv", sep = "",
                           stringsAsFactors=FALSE, row.names = 1, header= TRUE),
  
  #Cleaning the data -- making sure that taxa match and there is no missing data
  cleaned.cont = CleanData(tree, cont.data), 
  new_tree = cleaned.cont$phy,
  new_data = data_prune(tree, cont.data),
  
  two_times = name_changer_2(cleaned.cont$data, cleaned.cont$phy$tip.label, 
                             "mantle_length"),
  phy_temp = as.data.frame(two_times),
  
  phytools_ver = name_changer(phy_temp, "mantle_length"),
  SPECIES = phytools_ver$SPECIES,
  mantle = phytools_ver$mantle_length,
  phy_ver_name = set_species_names(mantle, SPECIES),
  
  #Look at the data
  plot_1 = VisualizeData(new_tree, 
                        phy_ver_name, "results/new_tree", 
                        "results/new_data", "results/contMap"),
  
  #Rates of evolution
  
  phy_conv = as.data.frame(phy_ver_name),
  phy_conv_swap = name_changer(phy_conv, "mantle_length"),
  phy_row_nam = first_to_row(phy_conv_swap, "mantle_length"),
  rate_BM1 = evo_rate(new_tree, phy_row_nam, model="BM"),
  
  print(paste("The rate of evolution is", rate_BM1$opt$sigsq, 
  "in units of centimeters per million years")),
  
  rate_OU1 = evo_rate(new_tree, phy_row_nam, model="OU"),
  print(paste("The rate of evolution is", rate_OU1$opt$sigsq, 
              "in units of meters per hundred million years")),
  par(mfcol = (c(1,2))),
  plot(new_tree, show.tip.label=FALSE),
  ou.tree = rescale(new_tree, model="OU", rate_OU1$opt$alpha),
  plot(ou.tree),
  
  save_jpeg(ou.tree, "OU_tree"),

  
  #Compare Trees
  AIC.BM1 = rate_BM1$opt$aic,
  AIC.OU1 = rate_OU1$opt$aic,
  AICc.BM1 = rate_BM1$opt$aicc,
  AICc.OU1 = rate_OU1$opt$aicc,
  delta.AIC.BM1 = AICc.BM1 - AIC.BM1,
  delta.AIC.OU1 = AICc.OU1 - AIC.OU1,
  
  #The OU model is superior, as the AIC is slightly lower. However, it should
  #be noted that the delta AIC of the BM model is lower, indicating slightly 
  #higher support for that model. If I had more data, I'd feel more confident 
  #using the OU model, but as of right now, they both seem comparable in terms
  #of utility. 
  
  print(paste("The AIC of BM1 is",AIC.BM1, "and the AIC of OU1 is", AIC.OU1,
              ".  Therefore,", is.better(AIC.OU1, name1="OU",delta.AIC.BM1, 
                                         AIC.BM1, name2="BM",delta.AIC.OU1))),
  
  #OUwie runs
  
  single.data = read.csv(file="data/taxa_binary_discrete.csv", sep = "",
                       stringsAsFactors=FALSE, row.names = 1, header= TRUE),
  cleaned.single = CleanData(tree, single.data), 
  new_data_single = data_prune(tree, single.data),
  two_times_single = name_changer_2(cleaned.single$data, 
                                    cleaned.single$phy$tip.label, "Tentacles"),
  phy_temp_single = swap_orient(as.data.frame(two_times_single)),
  one.discrete.char = phy_temp_single,
  reconstruction_info = ace(one.discrete.char, new_tree, type="discrete",
                             method="ML", CI=TRUE),
  save_me(reconstruction_info, "reconstruction_info"),
  best.states = colnames(reconstruction_info$lik.anc)[apply(
    reconstruction_info$lik.anc, 1, which.max)],
  
  labeled.tree = change_tips(new_tree, best.states), 
  ouwie_continuous = ouwie_data(phy_row_nam, one.discrete.char), 
  nodeBased.OUMV = OUwie(labeled.tree, ouwie_continuous,model="OUMV", 
                          simmap.tree=FALSE, diagn=FALSE),
  print(nodeBased.OUMV),
  
  #The first row of numbers refer to how well the given model fits the data.
  #(lnL, AIC, AICc, BIC) and model and ntax (number of taxa) are what they 
  #say on the tin. 
  
  #The numbers under the rate reference how much of the data we see are driven
  #by either Brownian motion ("random" chance) or natural selection/some other
  #active selective force. In this case, the no_tentacles trait seems to be 
  #more heavily influenced by active forces, whereas the tentacles trait is 
  #acquired more through random drift/chance/etc. I interpret this as the alpha
  #value corresponds to the selective force condition and the sigma value 
  #corresponds to the rate of the Brownian motion walk aka its strength. The
  #larger the value for each of these traits, the more influential they are. 
  
  #The optima group relate to the rate of evolution of each trait across the
  #phylogeny. Estimate what it sounds like, and se is "standard error". 
  #Tentacles is being gained at a more rapid rate than no_tentacles is being 
  #lost, but the SE for each of those isn't stellar. It's better in tentacles, 
  #as it can't switch "direction" within it, but still not amazing in my opinon.
  
  #Half-life is what it says on the print out--"(another way of reporting alpha)
  
  #Run all OUwie models
  models = c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"),
  results = all_OUwie(models, phy=labeled.tree, data=ouwie_continuous),
  AICc.values = list(results$BM1$AICc, results$BMS$AICc, results$OU1$AICc,
                     results$OUM$AICc, results$OUMV$AICc, results$OUMA$AICc,
                     results$OUMVA$AICc),
  AICc.values_2 = list_names(models, AICc.values),
  AICc.values_3 = list_names(models, c(AICc.values_2$BM1, AICc.values_2$BMS,
                                       AICc.values_2$OU1, AICc.values_2$OUM, 
                                       AICc.values_2$OUMV, AICc.values_2$OUMA,
                                       AICc.values_2$OUMVA)),
  print(AICc.values_3),
  AICc.values_4 = AICc.values_3 - min(AICc.values_3),
  print(AICc.values_4),
  best = results[[which.min(AICc.values)]],

  
  #Calculate point likelihood values
  alpha.values = seq(from= 0.0000001 , to= 2 , length.out=50),
  likelihood.values = rep(NA, length(alpha.values)),
  print(likelihood.values),
  likelihood.values_2 = mass_likelihood(labeled.tree, ouwie_continuous, 
                                         "OUMV", alpha.values, 
                                        likelihood.values, best),
  print(likelihood.values_2),
  
  plot_alpha_1 = plot_alpha(x= alpha.values , y= likelihood.values_2, 
                          xlim=c(-1,2), ylim=c(-37,130), best=best),
  
  #Look at theta parameters
  require("akima"),
  nreps=400,
  #I had to alter 5*best$theta[1,2] here, as it is naturally 0 and that messes 
  #with the interpolation down the line. 
  theta1.points = c(best$theta[1,1],
                    rnorm(nreps-1,best$theta[1,1], 5*0.0001)), 
                    #center on optimal value, have extra variance
  theta2.points = c(best$theta[2,1], 
                    rnorm(nreps-1, best$theta[2,1], 5*best$theta[2,2])), 
                    #center on optimal value, have extra variance
  likelihood.values_theta = rep(NA,nreps),
  likelihood.values_theta_2 = mass_likelihood_theta(labeled.tree, 
                                                    ouwie_continuous, "OUMV", 
                                                    nreps, best, theta1.points, 
                                                    theta2.points, 
                                                    likelihood.values_theta),
  
  likelihood.differences = (-(likelihood.values_theta_2-
                                max(likelihood.values_theta_2))),
  
  #interpolate
  interpolated.points = interp(x=theta1.points, y=theta2.points, 
                               z= likelihood.differences, linear=FALSE, 
                               extrap=TRUE, xo=seq(min(theta1.points), max(theta1.points), length = 400), 
                               yo=seq(min(theta2.points), max(theta2.points), length = 400)),
  interpolated.graph = inter_graph(interpolated.points, theta1.points, 
                                   theta2.points, best, ouwie_continuous$regime,
                                   ouwie_continuous$continuous),
  
  #discrete
  
  trait.ordered = data.frame(ouwie_continuous[,2], ouwie_continuous[,2],
                             row.names=ouwie_continuous[,1]),
  trait.ordered_2 = trait.ordered[labeled.tree$tip.label,],
  z_trait = trait.ordered[,1],
  z_trait_2 = as.vector(t(as.data.frame(z_trait))),
  z_trait_3 = name_it(z_trait_2, labeled.tree$tip.label),
  
  tree.mapped = make.simmap(labeled.tree, z_trait_3 , model="ER", nsim=1),
  save_me_general(tree.mapped$Q, "simmap_tentacles"),
  leg = c("black","red"),
  leg_name = c(1,2),
  named_leg = name_it(leg, leg_name),
  tree.splot = plotSimmap(tree.mapped,named_leg,pts=FALSE,ftype="off", lwd=1),
  
  #ouwie.swap = 
  simmapBased = OUwie(tree.mapped,ouwie_continuous,model="OUMV", 
                      simmap.tree=TRUE, diagn=FALSE),
  print(simmapBased),
  print(best),
  
  
  #How this compares to best model from above? 
  #The best model is much better, with a lnL of ~125 as compared to the lnL of 
  #the simmap model where the lnL is -10.95. Even comparing the OUMV result
  #from the ouwie with continuous data, the simmap falls short (sim: -10.9538,
  # cont: -10.80084), though only by a little bit.  
  
  #In this case, I think they're comparable, in the sense that you can look at 
  #both and choose one for use in further analysis, as the OUMV loglik is 
  #similar regardless of whether or not it was made from mantle length data
  #or tentacle presence data. 
  
  #***Side-note***$#
  #I'm slightly skeptical as to whether the OUMA model is as spectacular as 
  #the loglik wants me to think it is--and very well could be an artifact of 
  #how few data points are in this particular data set, but I'm going with it
  #as the "best" model for the purposes of this exercise. I just can't think of
  #any errors I could have made during the coding process that would have 
  #caused any issues with that command (especially given how similar the loglik
  #was in between the OUMV model developed by the command and the simmap with
  #totally different character data). I could very well be interpreting it wrong
  #if you (Brian) notice anything wrong with my conceptual understanding of the
  #model outputs, please let me know! 
  )
