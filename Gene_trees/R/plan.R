plan <- drake_plan(
  phy_1 = get_study_tree("ot_485", "tree1"),
  plot(phy_1, cex=0.3),
  
  phy_2 = drop.random(phy_1, Ntip(phy_1) - 10),
  plot(phy_2),
  axisPhylo(),
  
  gene.tree = phybase::sim.coaltree.phylo(phy_2, pop.size=1e-12),
  plot(gene.tree),
  
  plot(cophylo(phy_2, gene.tree, cbind(sort(phy_2$tip.label), 
                                       sort(gene.tree$tip.label))), cex=0.3),
  newtree = rcoal(7),
  comp_2 = compare.random_trees(newtree),
  
  tip.rows = which(newtree$edge[,2]<=Ntip(newtree)),
  newtree_2 = newtree,
  compare.even_trees(newtree_2,tip.rows),
  
  newtree_2.clado = compute.brlen(newtree_2),
  
  gene_tree = phybase::sim.coaltree.phylo(newtree_2),
  gene_tree.clado = compute.brlen(gene_tree),
  
  final = plot(cophylo(newtree_2.clado, gene_tree.clado, 
               cbind(sort(newtree_2.clado$tip.label),
                sort(newtree_2.clado$tip.label))))
)
