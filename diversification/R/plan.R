plan <- drake_plan(
  
  my.tree = TreeSim::sim.bd.taxa(n=300, numbsim=1, lambda=0.1, mu=0)[[1]], 
  my.tree.plot = plot_tree(my.tree, file= "results/intial_tree", label=FALSE), 
  lin.time = ape::ltt.plot(my.tree),
  lin.time.y = ape::ltt.plot(my.tree, log="y"),
  yule.trees = TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=0.1, 
                                     mu=0, complete=FALSE),
  lin.time.multi = ape::mltt.plot(yule.trees)
)
