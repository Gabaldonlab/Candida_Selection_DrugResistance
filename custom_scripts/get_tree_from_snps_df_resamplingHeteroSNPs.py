#!/usr/bin/env python

# This is a pipeline to be ran on Candida_mine_env. It takes a df with diploid SNPs (already filtered for not having positions with low coverage or INDELs) and it runs iqtree to generate a tree from an alignment where the heterozygous SNPs are taken (or not) randomly.

##### DEFINE ENVIRONMENT #######

# module imports
import argparse, os
import pandas as pd
import numpy as np
from argparse import RawTextHelpFormatter
import copy as cp
import pickle
import string
import shutil 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
import sys
from shutil import copyfile
import time

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# import functions
import Cmine_functions as fun

###################################

########## PARSE CMD ARGS ########## 


description = """
Runs a phylogenetic tree reconstruction taking randomly heterozygous SNPs.
"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

parser.add_argument("-o", "--outdir", dest="outdir", action="store", required=True, help="outdir")
parser.add_argument("-r", "--ref", dest="ref", required=True, help="Reference genome. Has to end with .fasta.")
parser.add_argument("-thr", "--threads", dest="threads", default=16, type=int, help="Number of threads, Default: 16")
parser.add_argument("--snps_df_file", dest="snps_df_file", action="store", type=str, required=True, help="The file with SNPs")
parser.add_argument("--samples_file", dest="samples_file", action="store", type=str, required=True, help="The file with the necessary samples")

opt = parser.parse_args()

####################################

####### process inputs #####

print("running tree reconstruction into %s with %i threads"%(opt.outdir, opt.threads))
fun.make_folder(opt.outdir)

# get the samples
sorted_samples = sorted(pd.read_csv(opt.samples_file, sep="\t").sampleID.apply(int))

# define the final file
final_file = "%s/tree_was_generated.txt"%opt.outdir

# define chromosomes with some SNPs
chroms_with_SNPs = set(fun.load_object(opt.snps_df_file)["#CHROM"])
chr_to_len = fun.get_chr_to_len(opt.ref)
all_chroms = set(chr_to_len)

strange_chroms = chroms_with_SNPs.difference(all_chroms)
if len(strange_chroms)>0: raise ValueError("There are some strange chrom names with SNPs: %s"%strange_chroms)

chroms_without_SNPs =  all_chroms.difference(chroms_with_SNPs)
if len(chroms_without_SNPs)>0: print("There are some chroms without SNPs: %s. These are their lengths:\n%s"%(chroms_without_SNPs, "\n".join(["%s: %i"%(c, chr_to_len[c]) for c in chroms_without_SNPs])))

# define the refrence genome as upper, onl
reference_genome = "%s/reference_genome.chroms_with_SNPs.fasta"%opt.outdir
all_chroms_seqRecords = [SeqRecord(Seq(str(seq.seq).upper()), id=seq.id, description="", name="") for seq in SeqIO.parse(opt.ref, "fasta") if seq.id in chroms_with_SNPs]
SeqIO.write(all_chroms_seqRecords, reference_genome, "fasta")

############################

##### GENERATE MULTIFASTA ########

# generate multifasta with the alternative genome and ranomly chosen SNPs. It only contains positions where there was some SNP
multifasta = "%s/multifasta_several_samples_randomlyChosenHetSNPs.fasta"%opt.outdir
fun.generate_multifasta_from_snps_df_file_FastaAlternateReferenceMaker(opt.snps_df_file, reference_genome, multifasta, sorted_samples, threads=opt.threads, expected_GTs={"1/1", "0/1"}, pickRandomHetSNPs=True)

# generate a fasta only with the variable sites
multifasta_onlyVariableSites = fun.get_multifasta_onlyVariableSites_snpSites(multifasta)


"""
# slow way

# load snps df
print("loading snps")
snps_df = fun.load_object(opt.snps_df_file)
snps_df["sampleID"] = snps_df.sampleID.apply(int)
snps_df = snps_df.set_index("sampleID", drop=False)

# generate multifasta with the alternative genomes, with randomly chosen SNPs (only positions with SNPs)
multifasta = "%s/multifasta_several_samples_randomlyChosenHetSNPs.fasta"%opt.outdir
positions_df_file = fun.generate_multifasta_from_snps_df(snps_df, multifasta, sorted_samples, threads=1, generate_one_aln_each_chrom=True, pickRandomHetSNPs=True) # threads=1 was for a problem with the parallelization

# generate a concatenated alignment with only constant positions
multifasta_onlyVariableSites = "%s.onlyVariableSites.fasta"%multifasta
fun.get_multifasta_onlyVariableSites(positions_df_file, multifasta_onlyVariableSites, sorted_samples, replace=False)
"""
#####################################

######### BUILD TREE ###########

# get unrooted tree form the aln
outfileprefix= "%s/iqtree_unroted"%opt.outdir
iqtree_file = "%s.iqtree"%outfileprefix
if fun.file_is_empty(iqtree_file):

    # remove previous files. This is not necessary if I want to run from previous checkpoints
    """
    for f in os.listdir(opt.outdir):
        if f.startswith(fun.get_file(outfileprefix)): fun.remove_file("%s/%s"%(opt.outdir, f))
    """

    # keep cmd
    print("running iqtree")
    fun.run_cmd("iqtree -s '%s' -pre %s --mem 25G -m GTR+F+ASC+G4 -T AUTO -ntmax %i"%(multifasta_onlyVariableSites, outfileprefix, opt.threads), env="Candida_mine_env") # GTR+ASC is reccommended by iqtree for SNP data. GTR+F+ASC+G4 is like the GTRGAMMA model (used in a recent C. tropicalis paper for resampling) of iqtree but with ascertainty correction. In the homozygous SNPs-tree I found that the TVM model was more frequent, although the GTR was closest. I decide to use this GTR+F+ASC+G4 model as it should be good enough

################################

# clean files
for f in os.listdir(opt.outdir):
    if f.startswith("reference_genome."): fun.remove_file("%s/%s"%(opt.outdir, f))

# write the final file
open(final_file, "w").write("The tree was generated!!\n")
print("tree correctly generated")