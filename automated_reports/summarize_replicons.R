#######################################
# REPLICON TREE VISUALISATION

# Input: NEWICK file created either by ClonalFrameML (.newick) or by IQtree (.treefile)

# Output: JPEG file 

# Clear R environment
rm(list = ls())
library(data.table)
library(readxl)
library(ape)
library(phytools)

## GLOBAL PARAMETERS ##
file_extension="newick"
out_prefix="GTR_tree"
replicon="Cluster 11 B. cereus"
metadata_file="metadata.xlsx"
cluster_id=11
mlst_species="bcereus"

## PROGRAM ##
cat("Searching tree file.\n")
{
  setwd(sprintf("phylogeny/cluster_%d/%s/phylogenetic_tree/clonalframeml", cluster_id, mlst_species))
  treefiles <- dir(pattern = sprintf("%s$", file_extension), recursive = TRUE, full.names = TRUE)
  # Reading tree
  # t <- read.nexus(treefiles)
  t <- read.tree(treefiles)
  # Set root to midpoint
  t <- midpoint.root(t)
}

cat("Renaming tips of tree.\n")
{
  # Rename tips
  metadata <- data.table(read_excel(metadata_file, sheet = 1, col_types = "text"))
  time_to_change_df <- sapply(t$tip.label, function(x) metadata$GLIMS[metadata$SAMPLE_ID==x])
  t$tip.label <- time_to_change_df
}

cat("Output tree file and color clonal samples.")
{
  # Function to color clonal samples together
  colorClonal <- function(id1, id2, coloring){
    bottom <- which(t$tip.label == id1)
    top <- which(t$tip.label == id2)
    return(rect(pp$xx[bottom] - 0.01, pp$yy[bottom] - 0.05, pp$xx[bottom] + 0.02, pp$yy[top] + 0.05,
         col = coloring, border = NA))
    
  }
  # Get tree in JPEG file
  jpeg(file = sprintf("%s.jpeg", out_prefix), 2000, 3000, quality = 98, pointsize = 48)
  layout(matrix(1:1, 1))
  plot.phylo(t, type = "phylogram", main = replicon,
             align.tip.label = TRUE, label.offset = 0, direction = "l", edge.width = 2, cex.main = 2, plot = FALSE)
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  par(new = TRUE)
  plot.phylo(t, type = "phylogram", main = replicon,
             align.tip.label = TRUE, label.offset = 0, direction = "l", edge.width = 2, cex.main = 2)
  axisPhylo(lwd = 5, cex.axis = 1)
  mtext("Substitutions/site", side = 1, line = 3, cex = 1)
  # legend(x = 'bottomright', inset=c(0,0), border = "black", cex = 1,
  #        fill = c("lightsalmon"), 
  #        legend = c("Cluster A"))
  dev.off()
}

# cat("To display circular plots...\n")
# {
#    jpeg(file = "circular_plots.jpeg", 2000, 2000, quality = 98, pointsize = 48)
#    layout(matrix(1:1, 1))
#    plot(t, "f", FALSE, cex = 0.5)
#    dev.off()
# }




