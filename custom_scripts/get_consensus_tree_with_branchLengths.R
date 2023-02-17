#!/usr/bin/env Rscript

# This script takes a consensus tree and a set of trees and generates a tree with the same topology but branch lengths

# define environment
library(argparser, quietly=TRUE)
library(phytools, quietly=TRUE)

# print the traceback on exit
#options(error=function()traceback(2))
options(warn=1)

# parse cmd line args
argp = arg_parser("Adds branch lengths to consensus tree")

argp = add_argument(argp, "--input_consensus_treefile", help="The input")
argp = add_argument(argp, "--all_trees_file", help="The file with all the trees")
argp = add_argument(argp, "--output_consensus_treefile", help="The outfile of the tree with the added branches")

opt = parse_args(argp)

# load trees
all_trees_object = ape::read.tree(opt$all_trees_file)
input_consensus_tree = ape::read.tree(opt$input_consensus_treefile)

# make the input tree dicothomous
input_consensus_tree = multi2di(input_consensus_tree, random=FALSE) 

# get the branch lengths
consensus_tree_withBL = consensus.edges(all_trees_object, method="mean.edge", consensus.tree=input_consensus_tree)

# write tree
ape::write.tree(consensus_tree_withBL, file=opt$output_consensus_treefile)

