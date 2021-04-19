plan <- drake_plan(
  
  my.tree = TreeSim::sim.bd.taxa(n=300, numbsim=1, lambda=0.1, mu=0)[[1]], 
  my.tree.plot = plot_tree(my.tree, file= "results/intial_tree", label=FALSE), 
  lin.time = ape::ltt.plot(my.tree),
  log.time.y = ape::ltt.plot(my.tree, log="y"),
  yule.trees = TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=0.1, 
                                     mu=0, complete=FALSE),
  lin.time.multi = ape::mltt.plot(yule.trees, legend = FALSE, log = "y"),
  
  #Birth-death trees
  bd.trees = TreeSim::sim.bd.taxa(n=300, numbsim=10, 
                                  lambda=1, mu=.9, complete=FALSE),
  bd.multi = ape::mltt.plot(bd.trees, log="y", legend=FALSE),
  
  #Comparison
  
  depth.range = range(unlist(lapply(yule.trees,ape::branching.times)), 
                      unlist(lapply(bd.trees,ape::branching.times))),
  max.depth = sum(abs(depth.range)),
  
  comp.graph = compare.multi.trees(bd.trees, yule.trees, max.depth, y_start = 1,
                                   "Birth-death", "Yule"), 
  comp.graph.zoom = compare.multi.trees(bd.trees, yule.trees, 5, y_start = 200,
                                   "Birth-death", "Yule"), 
  
  #Test trees 
  #Speciation >> Extinction
  trees.1 = TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=2, mu=0.1, 
                                 complete=FALSE),
  test.1 = ape::mltt.plot(trees.1, log="y", legend=FALSE),
  
  #Lambda-mu difference constant
  trees.2 = TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=1, mu=0.9, 
                                 complete=FALSE),
  test.2 = ape::mltt.plot(trees.2, log="y", legend=FALSE),
  trees.3 = TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=0.5, mu=0.4,
                                 complete=FALSE),
  test.3 = ape::mltt.plot(trees.3, log="y", legend=FALSE),
  
  #Lambda-mu sum constant
  trees.4 = TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=1, mu=0.1),
  test.4 = ape::mltt.plot(trees.4, log="y", legend=FALSE),
  trees.5 = TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=0.6, mu=0.5),
  test.5 = ape::mltt.plot(trees.5, log="y", legend=FALSE),
  
  #When the speciation rate >> extinction rate, the time it takes for the trees
  #to coalesce to a more or linear log increase it much shorter, i.e. over a few
  #time units as compared to ~50 time units. 
  
  #When the lambda-mu difference is constant, but the actual values vary, 
  #the higher combo (1 and 0.9) leads to the trees taking longer to start 
  #properly diversifying (-28 time units vs -38 time units), but having a 
  #a slightly faster rate of diversification, resulting in the slope being 
  #steeper than that of the lower combo (0.5 and 0.4).
  
  #When the lambda-mu sum is constant, but the actual values are not, some 
  #fairly interesting things happen. When there is a very large difference 
  #between lambda and mu, the trees take very little time to reach the same
  #number of taxa (less than 8 time units) in comparison to the lower but 
  #similar combination (almost up to 60 time units). Furthermore, the 
  #combination with the large difference results in largely uniform trees, with
  #the overall time required to diversify by the same amount varying within a
  #range of only 4 time units and the trees have about the same shape. In 
  #comparison, the trees simulated with the similar, but smaller, lambda and mu
  #are erratically shaped and can differ in the amount of time to converge on 
  #the same level of diversification by 30 time units. 
  
  
  #Tree + Trait models of diversification 
  #0A, 1A, 0B, 1B
 
  #phy = tree.musse(c(.1,.1,.1,.2,.03,.03,.03,.03,.01,.01,0,.01,0,.01,.01,0,.01,0,.01,.01), max.taxa=50, 
  #                 x0=1, include.extinct=FALSE), Will not run in Drake--
  #                                               defaulting to console versoin.
  phy = test2,
  sim.dat.true = data.frame(names(phy$tip.state), phy$tip.state),
  sim.dat = sim.dat.true,
  
  sim.dat.convert = sim.converted(sim.dat),
  plot2 = plot(phy), 
  
  transform.var = knitr::kable(cbind(sim.dat.convert, 
                                     true.char=sim.dat.true$phy.tip.state)),
  
  #reparametrization 
  
  #BISSE equivalent
  turnover.anc.bisse = c(1,1,0,0),
  eps.anc.bisse = c(1,1,0,0),
  
  #Separate rates for 0A and 1A
  turnover.anc.sep = c(1,2,0,0),
  
  #Full HISSE w/Yule no extinction
  turnover.anc.hisse = c(1,2,3,4),
  eps.anc.hisse = c(0,0,0,0),
  
  
  #set up transition matrix
  
  #default
  trans.rates = TransMatMaker.old(hidden.states=TRUE),
  print(trans.rates),
  
  #remove dual transitions between observed and hidden
  trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10)),
  print(trans.rates.nodual),
  
  #set rates 1 and 6 as equal
  trans.rates.nodual.equal16 = ParEqual(trans.rates.nodual, c(1,6)),
  print(trans.rates.nodual.equal16),
  
  #set all as equal 
  trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, 
                                         c(1,2,1,3,1,4,1,5,1,6,1,7,1,8)),
  print(trans.rates.nodual.allequal),
  
  #matrix need to run BiSSE in HiSSE
  trans.rates.bisse = TransMatMaker.old(hidden.states=FALSE),
  print(trans.rates.bisse),
  
  pp = hisse.old(phy, sim.dat, f=c(1,1), hidden.states=TRUE, 
                 turnover.anc=turnover.anc.hisse, eps.anc=eps.anc.hisse, 
                 trans.rate=trans.rates.nodual.allequal),
  
  #Common Mistake - when hidden state is associated with ONLY one observed
  #state. You must manually remove transitions from the nonexistent pair (0B)
  #from the transition rate matrix. 
  turnover.an.mis = c(1,2,0,3),
  eps.anc.mis = c(1,2,0,3),
  
  trans.rates.no0b = TransMatMaker.old(hidden.states=TRUE),
  trans.rates.nodual.no0B = ParDrop(trans.rates.no0b, c(2,3,5,7,8,9,10,12)),
  print(trans.rates.nodual.no0B),
  
  
  
  #CID-2 model
  
  turnover.anc.cid2 = c(1,1,2,2),
  eps.anc.cid2 = c(1,1,2,2),
  
    #Usual HiSSE
  trans.rates.cid2 = TransMatMaker.old(hidden.states=TRUE),
  trans.rates.nodual.cid2 = ParDrop(trans.rates.cid2, c(3,5,8,10)),
    #All equal
  trans.rates.nodual.allequal.cid2 = ParEqual(trans.rates.nodual.cid2, 
                                         c(1,2,1,3,1,4,1,5,1,6,1,7,1,8)),
  print(trans.rates.nodual.allequal.cid2),
  
    #Three rates: A <--> B, 0->1, and 1->0. 
  
  trans.threerates = trans.rates.nodual.cid2,
  print(trans.threerates),
  # Set all transitions from 0->1 to be governed by a single rate:
  to.change.1 = cbind(c(1,3), c(2,4)),
  trans.threerates.1 = edit.rates(trans.threerates, to.change.1, 1),
  print(trans.threerates.1),
  # Now set all transitions from 1->0 to be governed by a single rate:
  to.change.2 = cbind(c(2,4), c(1,3)),
  trans.threerates.2 = edit.rates(trans.threerates.1, to.change.2, 2),
  print(trans.threerates.2),
  # Finally, set all transitions between the hidden state to be a single rate (essentially giving
  # you an estimate of the rate by which shifts in diversification occur:
  to.change.3 = cbind(c(1,3,2,4), c(3,1,4,2)),
  trans.threerates.3 = edit.rates(trans.threerates.2, to.change.3, 3),
  print(trans.threerates.3),
  
  pp.cid2 = hisse.old(phy, sim.dat, f=c(1,1), hidden.states=TRUE, 
                 turnover.anc = turnover.anc.cid2, eps.anc = eps.anc.cid2, 
                 trans.rate=trans.threerates.3), 
  
  #Plotting HiSSE reconstructions
  print(class(pp.recon)),
  print(pp.recon),
  
  #red<-->blue = rate, white<-->black for state
  #note that the default puts blue at the lowest and red at the highest, this
  #does NOT take into account the actual differences between the two, so 
  #pay attention to the numbers attached to blue and red. 
  
  hisse.plot1 = plot.hisse.states(pp.recon, rate.param="net.div", 
                                  show.tip.label=FALSE),
  
  hisse.plot1.acc.colors = plot.hisse.states(pp.recon, rate.param="net.div", 
                                             show.tip.label=FALSE, 
                                             rate.range=c(0,0.072)),
  
  #Modeling-averaging approach 
  print(pp.recon$aic),
  
  #If null, do this: 
  #pp.recon = MarginRecon(phy, sim.dat, f=c(1,1), hidden.states=TRUE, 
  #pars=pp$solution, aic=pp$aic, n.cores=2)
  
  hisse.results.list = make.hisse.list(),
  
  #second plot is based on the null 2-state model and the other assumes four
  #free turnover rates. 
  
  plot.hisse.3x.1 = plot.hisse.states(hisse.results.list[[1]], 
                                      rate.param="net.div", 
                                      show.tip.label=FALSE, 
                                      rate.range=c(0,0.072)),
  
  plot.hisse.3x.2 = plot.hisse.states(hisse.results.list[[2]], 
                                      rate.param="net.div", 
                                      show.tip.label=FALSE, 
                                      rate.range=c(0,0.072)),
  
  plot.hisse.3x.3 = plot.hisse.states(hisse.results.list[[3]], 
                                      rate.param="net.div", 
                                      show.tip.label=FALSE),
  
  #My data -- since my data are not discrete and don't really lend themselves
  #to diversification analysis, I will just use the standard cephalopod data
  #I got from googling (tentacles [squids and cuttlefishes] vs. no tentacles
  #[octopuses]) and the tree I got from the Open Tree of Life. 
  
  #I will just attempt to make a HiSSE plot, with three rates-- one between 
  #hidden states and one going each way between the observed states. 
  
  tree = read.tree("data/ot_1134.tre_try.txt"),
  print(tree),
  tree_image = plot_tree(tree, file=file_out("results/initial_tree"), FALSE),
  discrete.data = read.csv(file="data/taxa_binary.csv", sep = "",
                           stringsAsFactors=FALSE, row.names = 1, header= TRUE),
  cleaned.discrete = CleanData(tree, discrete.data),
  phy.tree = cleaned.discrete$phy,
  tree.node.labels = set.node.labels(phy.tree, check = FALSE, 
                                     places = rep("node", 16),
                                     new.node = rep(1,16)),
  
  tree.node.labels.ten = set.node.labels(phy.tree, check = FALSE, 
                                     places = c("node1", "node2", "node3", 
                                                "node4", "node5", "node6",
                                                "node7", "node8", "node9",
                                                "node10", "node11", "node12",
                                                "node13", "node14", "node15",
                                                "node16"),
                                     new.node = c(1,1,1,1,1,1,1,1,1,1,
                                                  0,0,0,0,0,0)),
  
  species.names = cleaned.discrete$phy$tip.label, 
  traits = cleaned.discrete$data[,1],
  hisse.data = data.frame(species.names,traits),
  print(hisse.data),
  
  hisse.transitions = trans.threerates.3, 
  turnover.ceph = c(1,1,2,2),
  eps.ceph = c(1,1,2,2), 
  
  hisse.ceph = hisse.old(tree.node.labels.ten, hisse.data, f=c(1,1), 
                         hidden.states=TRUE, turnover.anc = turnover.ceph, 
                         eps.anc = eps.ceph, trans.rate=hisse.transitions),
  
  ceph.recon = MarginRecon.old(tree.node.labels.ten, hisse.data, f=c(1,1), 
                           hidden.states=TRUE, pars=hisse.ceph$solution,
                           AIC=hisse.ceph$AIC, n.cores=1),
  
  hisse.ceph.plot = plot.hisse.states(ceph.recon, rate.param="net.div", 
                                      show.tip.label=FALSE),
  
  #Upon looking at the HiSSE plot using the three state transition matrix and
  #the tentacles or no tentacles data I have, there appears to be no change
  #in diversity across the clades with or without tentacles. The "hidden"
  #states map exactly to tentacles or no tentacles, which either makes me 
  #think I did this incorrectly, or there are no "hidden" traits influencing
  #diversification in this particular situation. 
)
