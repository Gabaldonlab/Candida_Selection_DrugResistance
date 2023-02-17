#!/usr/bin/env python

# This is a pipeline to be ran on Candida_mine_env. It calculates the number of SNPs per window between one sample ID and the others

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
It calculates the number of SNPs per window between one sample ID and the others
"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

parser.add_argument("--outdir", dest="outdir", action="store", required=True, help="outdir")
parser.add_argument("--outfile", dest="outfile", action="store", required=True, help="An outfile with chromosome, start, end, query_sampleID, subject_sampleID and number_vars, with no header ")

parser.add_argument("--threads", dest="threads", default=16, type=int, help="Number of threads, Default: 16")
parser.add_argument("--df_variants_file", dest="df_variants_file", action="store", type=str, required=True, help="A .py file with a df with the variants")
parser.add_argument("--query_sampleID", dest="query_sampleID", action="store", type=str, required=True, help="A string with the query sampleID")
parser.add_argument("--filtered_windows_bed", dest="filtered_windows_bed", action="store", type=str, required=True, help="A .bed file with the target windows on which to calculate the variant numbers.")
parser.add_argument("--samples_withNoVars", dest="samples_withNoVars", action="store", type=str, required=True, help="A comma sepparated string with sampleIDs that have no variants")


opt = parser.parse_args()

####################################

print("running pairwise var calculation into %s"%opt.outdir)

# make the outdir
fun.make_folder(opt.outdir)

# load vars df
print("loading %s"%opt.df_variants_file)
df_vars = fun.load_object(opt.df_variants_file)
df_vars["sampleID"] = df_vars.sampleID.apply(str)
df_vars["numeric_variantID"] =  df_vars.numeric_variantID.apply(int)

# define all the expected samples
if opt.samples_withNoVars!="no_samples": samples_withNoVars = set(opt.samples_withNoVars.split(","))
else: samples_withNoVars = set()

all_expected_samples = set(df_vars.sampleID).union(samples_withNoVars)

# map each sample to a set of the different variants as compared to opt.query_sampleID
print("getting sampleID_to_differentVars")
sampleID_to_differentVars = fun.get_sampleID_to_set_different_vars_compared_to_query_sampleID(df_vars[["sampleID", "numeric_variantID"]], opt.query_sampleID, opt.threads, all_expected_samples)

# check that the query sample has 0 different vars
if len(sampleID_to_differentVars[opt.query_sampleID])!=0: raise ValueError("There are some different vars with the query sample: %s"%sampleID_to_differentVars[opt.query_sampleID])

# get the outfile. a dataset with the pariwise distances for all samples
print("generating pairwise distances")
fun.generate_pairwise_variantDistances_outfile(df_vars[["numeric_variantID", "#CHROM", "POS"]].drop_duplicates(), opt.filtered_windows_bed, opt.query_sampleID, sampleID_to_differentVars, opt.threads, opt.outfile, opt.outdir, all_expected_samples)

# clean
print("the script finished correctly")
