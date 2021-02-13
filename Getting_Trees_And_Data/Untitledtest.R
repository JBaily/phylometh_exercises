##testing with instructions from 
#https://cran.r-project.org/web/packages/rotl/vignettes/rotl.html

taxa <- c("Hyla", "Salmo", "Diadema", "Nautilus")
resolved_names <- tnrs_match_names(taxa)
print(resolved_names)
#resolved_names <- tnrs_match_names(taxa, context_name = "Animals")
#print(resolved_names)

my_tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id)

plot(my_tree, no.margin = TRUE)
