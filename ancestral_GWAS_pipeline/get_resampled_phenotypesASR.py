#!/usr/bin/env python

##### DEFINE ENVIRONMENT #######

# general module imports
import argparse, os, time, sys
from argparse import RawTextHelpFormatter
from ete3 import Tree
import pandas as pd
import multiprocessing as multiproc

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# import functions
import ancestral_GWAS_functions as fun

################################

###### ARGS ######

description = """
Inputs a table with the phenotypes (--phenotypes) and a filename to save the pastml results of 10,000 resamples (--resampled_phenotypes_pastml_out). It runs pastml 10,000 times on a random permutation of the phenotypes.
"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

# mandatory args
parser.add_argument("--phenotypes", dest="phenotypes", required=True, help="The table with the phenotypes. It should be a .tab file with a header and 'sampleID' and 'phenotype' (binary 1/0).")
parser.add_argument("--treefile", dest="treefile", required=True, help="The tree file as generated in the get_AST_mutations.py")
parser.add_argument("--resampled_phenotypes_pastml_out", dest="resampled_phenotypes_pastml_out", required=True, help="The filename that has the pastml ran on resampled phenotypes.") 
parser.add_argument("--tmpdir", dest="tmpdir", required=True, help="The output directory, where all tmp files will be written.")
parser.add_argument("--pastml_prediction_method", dest="pastml_prediction_method", required=True, help="The prediction method passed to pastml to reconstruct the ancestral phenotypes.")

# optional args
parser.add_argument("--replace", dest="replace", action="store_true", default=False, help="Re-run all the steps by deleting the output directory.")
parser.add_argument("--threads", dest="threads", type=int, default=4, help="The number of threads.")

# parse all the args
opt = parser.parse_args()

##################


####### CODE ########

# debug
if not fun.file_is_empty(opt.resampled_phenotypes_pastml_out): 
	print("WARNING: %s already exists! Exiting..."%opt.resampled_phenotypes_pastml_out)
	sys.exit(0)

# prepare input files
fun.make_folder(opt.tmpdir)

# prepare the inputs of a function that will resample the phenotypes and run pastml on it
nresamples = 10000
inputs_fn = [(opt.phenotypes, I+1, opt.tmpdir, opt.pastml_prediction_method, nresamples, opt.treefile) for I in range(nresamples)]

fun.print_with_runtime("Running pastml jobs for resampled phenotypes in %i threads for %i resamples"%(opt.threads, nresamples))
with multiproc.Pool(opt.threads) as pool:
    df_pastml_results = pd.concat(pool.starmap(fun.get_resampled_phenotypes_df_pastml_run_one_resample, inputs_fn, chunksize=1))
    pool.close()
    pool.terminate()

#####################

# clean
fun.delete_folder(opt.tmpdir)
fun.save_object(df_pastml_results, opt.resampled_phenotypes_pastml_out)

fun.print_with_runtime("SUCCESS! The generating of 10000 resampled phenotypes and pastml run worked")
