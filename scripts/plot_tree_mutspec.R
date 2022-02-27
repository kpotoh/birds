# devtools::install_github('r-lib/systemfonts')
# BiocManager::install("svglite", dependencies = TRUE)


rm(list=ls(all=T))
setwd("birds")


library(ggtree)
library(ggplot2)
library(treeio)
library(dplyr)
library(ape)


bog_tree = read.tree("./data/interim/iqtree_runs/drun1/anc.treefile")

sbs <- "T>G"
mutspec_file <- paste('./data/processed/sbs_on_tree/', sbs, '_edge_mutspec.tsv', sep="")
grad_val = read.table(mutspec_file, sep="\t", header = T)
names(grad_val) = c('RefNode', 'label', 'MutSpec')

test = full_join(as_tibble(bog_tree), tibble(grad_val[,c(2,3)]), by = 'label')
der = as.treedata(test)

out_image_file <- paste('./figures/', sbs, '_tree.svg', sep = '')
svglite::svglite(out_image_file, width = 20, height = 20)
ggtree(der, aes(color = der@data$MutSpec)) + geom_tiplab(size = 4, colour = "darkgray") +
  scale_color_continuous(low="green", high="red") +
  theme(legend.position="bottom")+
  labs(col=paste(sbs, 'share'))

dev.off()






############## for all sbs
for (n1 in c('A', 'C', 'G', 'T')) {
  for (n2 in c('A', 'C', 'G', 'T')) {
    if (n1 == n2) {
      next
    }
    sbs <- paste(n1, n2, sep = ">")
    print(sbs)
    mutspec_file <- paste('./data/processed/sbs_on_tree/', sbs, '_edge_mutspec.tsv', sep="")
    grad_val = read.table(mutspec_file, sep="\t", header = T)
    names(grad_val) = c('RefNode', 'label', 'MutSpec')
    
    test = full_join(as_tibble(bog_tree), tibble(grad_val[,c(2,3)]), by = 'label')
    der = as.treedata(test)
    
    out_image_file <- paste('./figures/', sbs, '_tree.svg', sep = '')
    svglite::svglite(out_image_file, width = 20, height = 20)
    ggtree(der, aes(color = der@data$MutSpec)) + geom_tiplab(size = 4, colour = "darkgray") +
      scale_color_continuous(low="green", high="red") +
      theme(legend.position="bottom")+
      labs(col=paste(sbs, 'share'))
    
    dev.off()
    
  }
}

