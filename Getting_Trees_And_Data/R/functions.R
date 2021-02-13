plot_tree <- function(phy, file) {
  pdf(file=file, width=20, height=20)
  plot(phy)
  plot(phy, type="fan")
  plot(phy, type="fan", show.tip.label=FALSE, edge.width=0.1)
  dev.off()
}