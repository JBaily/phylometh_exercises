compare.random_trees <- function(speciestree) {
  
  speciestree$edge.length = speciestree$edge.length / 
    (10*max(branching.times(speciestree)))
  
  gene.tree_2 = phybase::sim.coaltree.phylo(speciestree)
  
  plot(cophylo(speciestree, gene.tree_2, cbind(sort(speciestree$tip.label), 
                                                sort(gene.tree_2$tip.label))))
}

compare.even_trees <- function(speciestree, tip_rows) {
  
  speciestree$edge.length[tip_rows] <- 100 + speciestree$edge.length[tip_rows]
  
  gene.tree3 <- phybase::sim.coaltree.phylo(speciestree)
  
  plot(cophylo(speciestree, gene.tree3, cbind(sort(speciestree$tip.label), 
                                                sort(gene.tree3$tip.label))))
}