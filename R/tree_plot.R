# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("treeio")
# BiocManager::install("ggtree")

# browseVignettes("ggtree")
# browseVignettes("treeio")

library("treeio")
library("ggtree")

setwd("birds")


beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)
beast_tree

# vals <- get.data(beast_tree)$rate
vals <- rnorm(151, 10, 10)
ggtree(beast_tree, aes(color=vals)) +
  theme(legend.position="right") + 
  geom_text(aes(label=))

# ggtree(beast_tree, aes(color=rate)) +
#   # scale_color_continuous(low='darkgreen', high='red') +
#   theme(legend.position="right")

# beast_tree2 <- rescale_tree(beast_tree, branch.length = 'rate')
# ggtree(beast_tree2) + theme_tree2()


tree <- read.iqtree("./data/interim/iqtree_runs/drun1/anc.treefile")
ggtree(tree, aes(color='edge.length'))





read.ne