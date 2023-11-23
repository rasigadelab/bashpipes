#######################################
# REPLICON TREE VISUALISATION

# Input: NEWICK file created either by ClonalFrameML (.newick) or by IQtree (.treefile)

# Output: JPEG file 

# Clear R environment
rm(list = ls())
library(data.table)
library(readxl)
# library(stringr)
library(ape)
library(phytools)
# library(ggtree)
# library(treeio)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("gtable")

## GLOBAL PARAMETERS ##
file_extension="newick"
out_prefix="GTR_tree"
replicon="Cluster 3 S.marcescens"
metadata_file="~/Projets/5.Suivi/Metadata_eq_rasigade.xlsx"
cluster_id=3
mlst_species="chromosome"

## PROGRAM ##
cat("Searching tree file.\n")
{
  setwd(sprintf("~/Projets/2.Coding_projects/3.R/1.Automated_reports/phylogeny/cluster_%d/%s/phylogenetic_tree/clonalframeml", cluster_id, mlst_species))
  treefiles <- dir(pattern = sprintf("%s$", file_extension), recursive = TRUE, full.names = TRUE)
  # Lecture de l'arbre
  # t <- read.nexus(treefiles)
  t <- read.tree(treefiles)
  #Set root to midpoint
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
    return(rect(pp$xx[bottom] - 0.1, pp$yy[bottom] - 0.5, pp$xx[bottom] + 2, pp$yy[top] + 0.5,
         col = coloring, border = NA))
    
  }
  # Get tree in JPEG file
  jpeg(file = sprintf("%s.jpeg", out_prefix), 3000, 6000, quality = 98, pointsize = 48)
  layout(matrix(1:1, 1))
  plot.phylo(t, type = "phylogram", main = replicon,
             align.tip.label = TRUE, label.offset = 0, direction = "l", edge.width = 2, cex.main = 3.5, plot = FALSE)
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  colorClonal("022156756801-01", "022176281806-02", "lightsalmon")
  # colorClonal("022106460701-01", "022186972005-01", "lightgoldenrod1")
  # colorClonal("022081308401-01", "022081308401-01", "lightblue")
  # colorClonal("022184231101-02", "022184231101-02", "lightblue")
  par(new = TRUE)
  plot.phylo(t, type = "phylogram", main = replicon,
             align.tip.label = TRUE, label.offset = 0, direction = "l", edge.width = 2, cex.main = 3.5)
  axisPhylo(lwd = 5, cex.axis = 2)
  mtext("Substitutions/site", side = 1, line = 3, cex = 2)
  legend(x = 'bottomright', inset=c(0,0), border = "black", cex = 2,
         fill = c("lightsalmon"), 
         legend = c("Cluster A"))
  dev.off()
}

# cat("To display circular plots...\n")
# {
#    jpeg(file = "~/Projets/2.Coding_projects/3.R/11.Summarize_replicon/circular_plots.jpeg", 2000, 2000, quality = 98, pointsize = 48)
#    layout(matrix(1:1, 1))
#    plot(t, "f", FALSE, cex = 0.5)
#    dev.off()
# }




