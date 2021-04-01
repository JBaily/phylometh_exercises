plan <- drake_plan(
  #Getting a tree and some data -- lets try something simple like cats.

  tree = read.tree("data/ot_1134.tre_try.txt"),
  print(tree),
  tree_image = plot_tree(tree, file=file_out("results/initial_tree"), FALSE),
  discrete.data = read.csv(file="data/taxa_try_2.csv", sep = "",
                           stringsAsFactors=FALSE, row.names = 1, header= TRUE),
  
  #Cleaning the data -- making sure that taxa match and there is no missing data
  cleaned.discrete = CleanData(tree, discrete.data), 
  
  #Look at the data
  plot_1 = VisualizeData(cleaned.discrete$phy, 
                         cleaned.discrete$data, "results/new_tree", 
                         "results/new_data"),
  
  #Plot parsimony estimates
  cleaned.discrete.phyDat = phangorn::phyDat(cleaned.discrete$data, 
                                             type="USER", 
                                             levels = c("no_venom","venom",
                                                        "tentacles","pen",
                                                        "cuttlebone",
                                                        "no_tentacles",
                                                        "no_internal_shell",
                                                        "triangular_fins",
                                                        "slit_eyes","long_fins",
                                                        "circular_eyes",
                                                        "rectangular_eyes",
                                                        "no_fins","blue_rings",
                                                        "no_blue_rings")), 
  
  anc.p = phangorn::ancestral.pars(cleaned.discrete$phy, cleaned.discrete.phyDat),
  plot_3_series = print_plotAnc(cleaned.discrete$phy, anc.p, type="parsimony_anc_recon"),
  
  #Plot likelihood estimates
  anc.ml = phangorn::ancestral.pml(phangorn::pml(cleaned.discrete$phy, 
                                                 cleaned.discrete.phyDat), 
                                   type="ml"),
  plot_4_series = print_plotAnc(cleaned.discrete$phy, anc.ml, type="ML_anc_recon"),
  
  #How does this differ from parsimony? 
    #ML differs from parsimony in that it gives the a pie chart depicting the 
    #likelihood of each ancestral node being any given state. 
  
  #Why does it differ from parsimony?
    #Because parsimony treats the characters only contribute steps towards the 
    #overall tree length, which results in an "either or" situation for the
    #ancestral nodes, as the model will pick whichever one minimizes the overall
    #number of steps. ML treats them as values that contribute to one overall 
    #score, so likelihoods in the forms of percents that can be represented as 
    #pieces of a pie chart can be generated for each node. 
  
  #What does uncertainty mean?
   #In this context, I think that uncertainty refers to how sure we are that the
   #tree we are using is even correct in the first place. It could also refer to 
   #our limitations in accounting for all of the possible ancestral state 
   #reconstructions, aka the further we go back in the tree the more difficult
   #it is to say that there is a 20% of state a, 40% of state b, and 40% of state
   #c at this deep ancestral node. I believe this type of uncertainty stems from 
   #algorithmic limitations. 
  
  #Use the corHMM package to see if
  # 1.) Can estimate transition rates between states? Do it. 
  
 # converted_phangorn_ml = corHMM::ConvertPhangornReconstructions(anc.ml),
  HMM_phylo = cleaned.discrete$phy,
  data_conversion = as.data.frame(cleaned.discrete$data),
  testing_data = tibble::rownames_to_column(data_conversion, "Species"),
  matrix_data = as.matrix(testing_data),
  to_remove = c("Tentacles"),
  HMM_data = remove_cols(matrix_data,to_remove),
   
 
  legend_cor = corHMM::getStateMat4Dat(HMM_data),
  
  #transition_rates = corHMM::corHMM(HMM_phylo, HMM_data,
  #                                  rate.cat = 7, rate.mat = legend_cor$rate.mat,
  #                                  model="ARD", fixed.nodes = TRUE),
  #print_rates(transition_rates, "tentacles"),
  
  # 2.) How could you examine if transition rates are equal?
  #  You just print your rates variable and look at the values across the diagonal
  # in the case of tentacles, going from state 1 to 2 is .307 and the opposite is
  # 0.223, aka it is more likely to go from no tentacles to tentacles within if 
  # you are a cephalopod. 
  
  # 3.) Think about the Lewis (2001) MKV model. 
  #     Are your traits all variable? 
  #     - Yes, none of the traits are constants or autapomorphies. They are all
  #       variable. 
  #     Will using this make sense for your data? 
  #     - I assume so. Considering the points the Lewis paper made concerning
  #       that the Mkv model does tend to yield better log-likelihood units than
  #       parsimonious trees. 
  #     Try using it. Do results change?
  #     - No, from what I can tell. Of course, I'm likely doing it incorrectly,
  #       haha. 
  
  #anc.mkv = phangorn::ancestral.pml(phangorn::pml(cleaned.discrete$phy, 
  #                                               cleaned.discrete.phyDat, 
  #                                               Mkv = TRUE), type="ml"),
  #plot_5_series = print_plotAnc(cleaned.discrete$phy, anc.ml, type="Mkv_anc_recon"),
  
  # 4.) How could you test the order of state evolution?
  # you would use the makeSimmap function in corHMM, This will show you "when"
  # certain character states would likely have showed up on the tree based upon 
  # the transition rates. 
 
  #markov_plot = plotMKmodel(transition_rates, display="square", 
   #                         color=c("blue","red")),
  #model = transition_rates$solution,
  #replaced = replace_NA(transition_rates$solution),
  #history_char = makeSimmap(transition_rates$phy, transition_rates$data, 
   #                         replaced, rate.cat = 1, nSim = 1),
  #plot_gen(transition_rates$phy, transition_rates$data, replaced)
  #simmap_plot = plot_sim_tree(history_char, "tentacles")
)
