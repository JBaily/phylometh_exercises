plan <- drake_plan(
  taxon_names = read.csv(file=file_in("data/taxa.csv"), stringsAsFactors=FALSE),
  tree.id = rotl::tnrs_match_names(taxon_names$Taxon[1])$ott_id,
  tree = rotl::tol_subtree(ott_id=tree.id, label_format = "name"),
  tree_resolution = paste0("It has ", ape::Nnode(tree), " nodes resolved of ",(ape::Ntip(tree)-1), " possible"),
  tree_print = plot_tree(tree, file=file_out("results/tree.pdf")),
  #phy = ape::rcoal(20+round(stats::runif(1,1,20))),
  #even = is_even(phy),
  #plotted = plot_tree(phy, file_out("results/tree.pdf")),
  #save_even = save(phy, even, file=file_out("results/even.rda"))
)