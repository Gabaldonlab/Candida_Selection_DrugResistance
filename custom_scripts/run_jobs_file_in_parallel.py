#!/usr/bin/env python

# This script runs any jobs file in parallel for some cmd

# define env
import multiprocessing as multiproc
import os, sys


# parse cmd line args
jobs_file, env, threads = sys.argv[1:]
threads = int(threads)
print("running %s with conda env %s on %i threads"%(jobs_file, env, threads))

# functions
def run_cmd_simple(cmd):

    """This function runs a cmd"""
    # run
    out_stat = os.system(cmd) 
    if out_stat!=0: raise ValueError("\n%s\n did not finish correctly. Out status: %i"%(cmd, out_stat)) 

# run
inputs_fn = [(l.strip(),) for l in open(jobs_file, "r").readlines()]
print("Running %i jobs..."%len(inputs_fn))

with multiproc.Pool(threads) as pool:
    pool.starmap(run_cmd_simple, inputs_fn)
    pool.close()
    pool.terminate()

print("running worked well")
