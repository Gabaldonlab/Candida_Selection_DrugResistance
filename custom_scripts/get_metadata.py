#!/usr/bin/env python


# This should be run with the Candida_mine_env
# This script is useful to download all the metadata from candidamine.
##### DEFINE ENVIRONMENT #######

# module imports
import os
import sys
from ete3 import NCBITaxa
from Bio import SeqIO
import numpy as np
import multiprocessing as multiproc
import pandas as pd

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    run_in_cluster = False    
    threads = 4
else:
    run_in_cluster = True    
    ParentDir = "/gpfs/projects/bsc40/mschikora"
    threads = 48

# define the dir where all perSVade code is
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir
sys.path.insert(0, perSVade_dir)

# import functions
import Cmine_functions as cmine_fun

# define dirs
CurDir = "%s/CandidaMine_data_generation/v1"%ParentDir
sys.path.insert(0, CurDir)
DataDir = "%s/data"%CurDir
PlotsDir = "%s/plots"%CurDir
manually_curated_data = "%s/manually_curated_data"%CurDir

# define paths
perSVade_py = "%s/perSVade.py"%perSVade_dir

####################################

print("WARNING: This should be only run when there are all variant calling files created")
# define the interesting taxIDs
interesting_taxIDs = [5476, 5480, 273372, 273371, 5482, 5478, 498019] # all except 273372 (metapsilosis)

# Candida albicans C. parapsilosis, metapsilosis, orthopsilosis C. tropicalis C. glabrata C. auris

# load the coverage 
 
# get the taxIDs under these taxIDs
ncbi = NCBITaxa()

# init the joined metadata df
all_df_metadata = pd.DataFrame()

for taxID in interesting_taxIDs:
    print(taxID)

    # get the dir
    sciName = "_".join(ncbi.get_taxid_translator([taxID])[taxID].split()).replace("[","").replace("]","")
    taxID_dir = "%s/%s_%i"%(DataDir, sciName, taxID); cmine_fun.make_folder(taxID_dir)
    smallVars_dir = "%s/varCall_output"%taxID_dir

    # define the metadata_dir
    metadata_dir = "%s/metadata"%taxID_dir; cmine_fun.make_folder(metadata_dir)

    #  define the SRRs
    srrs_with_varcall = {f for f in os.listdir(smallVars_dir) if os.path.isdir("%s/%s/smallVars_CNV_output"%(smallVars_dir, f))}
    print("There are %i srrs"%len(srrs_with_varcall))

    # get the sra data
    sra_df = cmine_fun.load_object("%s/obtainingSRAdata/sra.SRA_runInfo_df.py"%(taxID_dir))
    sra_df = sra_df[sra_df.Run.isin(srrs_with_varcall)]

    # get the metadata
    metadata_df_file = "%s/obtainingSRAdata/sra_metadata_manually_curated.tab"%taxID_dir

    if cmine_fun.file_is_empty(metadata_df_file) or True:

        # get the metadata
        df_metadata = cmine_fun.get_sra_metadata_df(sra_df, metadata_dir, manually_curated_data, taxID)

        # save
        cmine_fun.save_df_as_tab(df_metadata, metadata_df_file)
        
    # load 
    df_metadata = cmine_fun.get_tab_as_df_or_empty_df(metadata_df_file)

    # add fields
    df_metadata["species_name"] = sciName
    df_metadata["taxID"] = taxID

    # keep 
    all_df_metadata = all_df_metadata.append(df_metadata)

# save
cmine_fun.save_df_as_tab(all_df_metadata, "%s/sra_metadata_manually_curated_all_species.tab"%DataDir)


