file <- file.path("/Users/lab/GitHub_files/phylometh_exercises/Dating_trees/Tree_2_annotated.tree")
nexus <- read.nexus(file)
nexus
phylo <- read.beast(file)
phylo
lol <- write.tree(nexus)
write.table(lol, file = "export_txt", sep="\t")
