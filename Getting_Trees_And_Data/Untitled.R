# Steps used for week 2 exercise
# Most are from Phylometh website

GetTreeFromOpenTree <- function() {
  library(rotl)
  library(ape)
  
  #getting phylogeny of the Tardigrada phylum from Open Tree and making it 
  #a tree
  
  tardigrada.id <- rotl::tnrs_match_names("eutardigrada")
  tardigrada.tree <- rotl::tol_subtree(ott_ids = tardigrada.id$ott.id)
  
  #plotting
  
  ape::plot.phylo("Tardigrada", type=fan, cex=0.2)
}

print(paste("The Tardigrada tree has ", ape::Ntip(tardigrada.tree), 
            "terminals and ", Nnode(tardigrada.tree), 
            " internal nodes out of ", ape::Ntip(tardigrada.tree)-2, 
            " possible, which means it is ", 
            round(100*(ape::Nnode(tardigrada.tree)-1)/
                    (ape::Ntip(tardigrada.tree)-3), 2), "% resolved", sep=""))



