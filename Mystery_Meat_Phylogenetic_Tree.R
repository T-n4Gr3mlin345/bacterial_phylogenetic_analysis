#Installing Packages
if (!require(remotes)) {
  install.packages("remotes")
}

remotes::install_github("evo-palaeo/treesurgeon")
library(treesurgeon)

# Import aligned protein FASTA file
aligned <- ape::read.FASTA("~/Documents/Bioinformatics MSc/Bioinformatics/Mystery Meat/aligned_taxa_seq.fasta", type = "AA")

#Convert alignnment into a matrix
aligned_matrix <- as.matrix(aligned)

#Calculate amino acid distances
aa_dist <- ape::dist.aa(aligned_matrix)

# Estimate a neighbour-joining tree
nj_tree <- ape::nj(aa_dist)

#Plot the unrooted tree
plot(
  nj_tree,
  cex = 0.6,
  no.margin = TRUE,
  type = "unrooted"
)

# outgroup
OG <- which(nj_tree$tip.label == "Anaerolineales_bacterium_MIQ4461479.1")
# Root the tree using the outgroup
rooted_NJtree <- root_node(nj_tree, OG)

# Plot the rooted tree (phylogram)
plot(rooted_NJtree, cex = 0.6, no.margin = TRUE)

# Plot the rooted tree (cladogram)
par(mar = c(2, 2, 3, 10))

plot(
  rooted_NJtree,
  type = "phylogram",
  use.edge.length = FALSE,
  cex = 0.75,
  font = 1,
  edge.width = 1.2,
  label.offset = 0.003,
  main = "Rooted Neighbour-Joining Cladogram"
)

box()
