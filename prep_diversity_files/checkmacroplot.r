# Check to see if the trees with the new distance values are in the old dataset

fia <- read.csv('/mnt/research/nasabio/data/fia/treedata10nov/finley_trees_continental_US_most_recent_evaluations_nov8_2017.csv', stringsAsFactors = FALSE)
treedists <- read.csv('/mnt/research/nasabio/data/fia/finley_pnw_tree_distances_to_subp_center_nov26_2018.csv')

table(treedists$TREE_CN %in% fia$TREE_CN)

macro_trees <- fia$TREE_CN[fia$TPA_UNADJ %in% min(fia$TPA_UNADJ)]

table(macro_trees %in% treedists$TREE_CN)