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
  
  pdf(file="results/OU_tree"),
  plot(new_tree, show.tip.label=FALSE),
  plot(ou.tree),
  dev.off(),
  
  #Compare Trees
  AIC.BM1 = rate_BM1$opt$aicc,
  AIC.OU1 = rate_OU1$opt$aicc,
  delta.AIC.BM1 = AIC.BM1 - rate_BM1$opt$aic,
  delta.AIC.OU1 = AIC.OU1 - rate_OU1$opt$aic,
  
  print(paste("The delta AIC of BM1 is",delta.AIC.BM1,
              "and the delta AIC of OU1 is", delta.AIC.OU1,
              ".  Therefore,", is.better(delta.AIC.OU1, name1="OU", 
                                         delta.AIC.BM1, name2="BM"))),
  
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
  #labeled.tree <- ________________
  #nodeBased.OUMV <- OUwie(tree, cleaned.continuous,model="OUMV", 
  #                        simmap.tree=FALSE, diagn=FALSE)
  #print(nodeBased.OUMV)
  
  #Run all OUwie models
  #models <- c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA")
  #results <- lapply(models, RunSingleOUwieModel, phy=tree, data=trait)
  #AICc.values<-sapply(results, "[[", "AICc")
  #names(AICc.values)<-models
  #AICc.values<-AICc.values-min(AICc.values)
  #print(AICc.values) #The best model is the one with smallest AICc score
  #best<-results[[which.min(AICc.values)]] #store for later
  #print(best) #prints info on best model
  
  #Calculate point likelihood values
  #alpha.values<-seq(from= _______________ , 
  #to= _______________ , length.out=50)
  #likelihood.values <- rep(NA, length(alpha.values))
  #for (iteration in sequence(length(alpha.values))) {
  #  likelihood.values[iteration] <- OUwie.fixed(tree, trait, model="OUMV", alpha=rep(alpha.values[iteration],2), sigma.sq=best$solution[2,], theta=best$theta[,1])$loglik
  #}
  #plot(x= _______________ , y= _______________, xlab="_______________", ylab="_______________", type="l", bty="n")
  #points(x=best$solution[1,1], y=best$loglik, pch=16, col="red")
  #text(x=best$solution[1,1], y=best$loglik, "unconstrained best", pos=4, col="red")
  #abline(h=_______________, lty="dotted") #Two log-likelihood
  
  #Look at theta parameters
  #require("akima")
  #nreps<-400
  #theta1.points<-c(best$theta[1,1], rnorm(nreps-1, best$theta[1,1], 5*best$theta[1,2])) #center on optimal value, have extra variance
  #theta2.points<-c(best$theta[2,1], rnorm(nreps-1, best$theta[2,1], 5*best$theta[2,2])) #center on optimal value, have extra variance
  #likelihood.values<-rep(NA,nreps)
  
  #for (iteration in sequence(nreps)) {
  #  likelihood.values[iteration] <- OUwie.fixed(tree, trait, model="OUMV", alpha=best$solution[1,], sigma.sq=best$solution[2,], theta=c(theta1.points[iteration], theta2.points[iteration]))$loglik
  #}
  
  #likelihood.differences<-(-(likelihood.values-max(likelihood.values)))
  
  #interpolate
  #interpolated.points<-interp(x=theta1.points, y=theta2.points, z= likelihood.differences, linear=FALSE, extrap=TRUE, xo=seq(min(theta1.points), max(theta1.points), length = 400), yo=seq(min(theta2.points), max(theta2.points), length = 400)
  #contour(interpolated.points, xlim=range(c(theta1.points, theta2.points)),ylim=range(c(theta1.points, theta2.points)), xlab="Theta 1", ylab="Theta 2", levels=c(2,5,10),add=FALSE,lwd=1, bty="n", asp=1)
  #points(x=best$theta[1,1], y=best$theta[2,1], col="red", pch=16)
  #points(x=trait$X[which(trait$Reg==1)],y=rep(min(c(theta1.points, theta2.points)), length(which(trait$Reg==1))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 1, plotted along x axis
  #points(y=trait$X[which(trait$Reg==2)],x=rep(min(c(theta1.points, theta2.points)), length(which(trait$Reg==2))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 2, plotted along y axis

  #discrete
  #trait.ordered<-data.frame(trait[,2], trait[,2],row.names=trait[,1])
  #trait.ordered<- trait.ordered[tree$tip.label,]
  #z<-trait.ordered[,1]
  #names(z)<-rownames(trait.ordered)
  #tree.mapped<-make.simmap(tree,z,model="ER",nsim=1)
  #leg<-c("black","red")
  #names(leg)<-c(1,2)
  #plotSimmap(tree.mapped,leg,pts=FALSE,ftype="off", lwd=1)
  
  #simmapBased<-OUwie(tree.mapped,trait,model="OUMV", simmap.tree=TRUE, diagn=FALSE)
  #print(simmapBased)
  #print(best)
  
  
  #How compares to best model from above? 
  
  )
