#!/usr/bin/env python

# runs the get_df_diversity_piN_piS_all_steps from Cmine_functions on Candida_mine_env from 

# modules
import os, sys

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    run_in_cluster = False 
    threads = 4   
else:
    run_in_cluster = True    
    ParentDir = "/gpfs/projects/bsc40/mschikora"
    threads = 48

# import functions
CurDir = "%s/CandidaMine_data_generation/v1"%ParentDir
sys.path.insert(0, CurDir)
import Cmine_functions as fun

# parse argv
ProcessedDataDir, species, type_var, type_vars_appearance, type_vars_SimpleRepeats, df_piN_piS_file, PlotsDir, sampleID_to_clade_file = sys.argv[1:]

# load files
df_piN_piS = fun.load_object(df_piN_piS_file)
sampleID_to_clade = fun.load_object(sampleID_to_clade_file)

# run function
fun.get_df_diversity_piN_piS_all_steps(ProcessedDataDir, species, type_var, type_vars_appearance, type_vars_SimpleRepeats, df_piN_piS, PlotsDir, sampleID_to_clade, threads)
                        
# log
print("get_df_diversity_piN_piS_all_steps was ran sucessfully for %s %s %s %s"%(species, type_var, type_vars_appearance, type_vars_SimpleRepeats))
