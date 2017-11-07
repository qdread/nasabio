# Create consensus tree from Ericson 1000 trees.
# Use to make a single phylogenetic tree for calculating BBS phylogenetic diversity.
# QDR - NASA project - 12 Sep 2017

# Load Brian O'Meara's code for getting branch lengths for a consensus tree
source('~/GitHub/nasabio/data_extraction/consensusbranchlength.r')

# Load 1000 ericson trees
erictree <- read.nexus('C:/Users/Q/Dropbox/projects/nasabiodiv/ericson1000.tre')

# Get consensus topology
eric_cons <- consensus(erictree)

# Use list of individual trees and consensus topology to get consensus tree with branch lengths.
eric_cons_brlen <- consensusBrlen(focalTree = eric_cons, sourceTreeList = erictree)


# Brian's code does not seem to work so let's use phytools instead
library(phytools)

eric_cons_edges <- consensus.edges(trees = erictree, method = 'least.squares')
eric_cons_edges <- reroot(eric_cons_edges, node.number = 1, position = 10)
save(eric_cons_edges, file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/ericson_cons_tree.r')
write.tree(eric_cons_edges, file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/ericson_cons.tre')
