#!/usr/bin/env python

# Runs pastml for one character renaming the tree nodes so that they have the first and last leaf names. Note that this will remove the data table

# imports 
import os, sys
from ete3 import Tree
import pandas as pd
import numpy as np

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# import functions
import ancestral_GWAS_functions as fun

# parse the arguments
treefile, data_table, workdir, outfile, html_file, pastml_prediction_method = sys.argv[1:]

# runs pastml as a function
fun.run_pastml_py_as_a_function(treefile, data_table, workdir, outfile, html_file, pastml_prediction_method, None)

