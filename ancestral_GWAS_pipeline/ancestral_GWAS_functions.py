#!/usr/bin/env python

########## ENV ###########

# import modules
import sys, os, shutil, pickle, psutil, subprocess
import pandas as pd
from ete3 import Tree, NodeStyle, TreeStyle, RectFace, CircleFace, TextFace
import numpy as np
import random
import itertools
# dask.dataframe as dask_df
import time
from datetime import date
import multiprocessing as multiproc
import copy as cp
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats import multitest
from matplotlib.lines import Line2D
from scipy.stats import chi2_contingency 



# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])
EnvName = EnvDir.split("/")[-1]
CondaDir =  "/".join(sys.executable.split("/")[0:-4])

# executables
PASTML = "%s/bin/pastml"%EnvDir
RUN_PASTML_PY = "%s/run_pastml.py"%CWD
RUN_GWAS_PY = "%s/run_GWAS.py"%CWD
GET_RESAMPLED_PHENOTYPES_ASR_PY = "%s/get_resampled_phenotypesASR.py"%CWD

##########################

# define variables
#all_pval_methods = ["simulation", "phenotypes", "RelToBranchLen", "RelToBranchLen_wReplace"]
fields_ASR_df = ["phenotype", "phenotype_transition", "genotype_transition", "genotype_appearance", "node_name", "group_name"]

# define functions
def get_fullpath(x):

    """Takes a path and substitutes it bu the full path"""

    # normal
    if x.startswith("/"): return x

    # a ./    
    elif x.startswith("./"): return "%s/%s"%(os.getcwd(), "/".join(x.split("/")[1:]))

    # others (including ../)
    else: return "%s/%s"%(os.getcwd(), x)

def remove_file(f):

    if os.path.isfile(f): 

        try: run_cmd("rm %s > /dev/null 2>&1"%f)
        except: pass

def delete_folder(f):

    if os.path.isdir(f): shutil.rmtree(f)


def make_folder(f):

    if not os.path.isdir(f): os.mkdir(f)


def delete_file_or_folder(f):

    """Takes a path and removes it"""

    if os.path.isdir(f): shutil.rmtree(f)
    if os.path.isfile(f): os.unlink(f)

def file_is_empty(path): 
    
    """ask if a file is empty or does not exist """

    if not os.path.isfile(path):
        return_val = True
    elif os.stat(path).st_size==0:
        return_val = True
    else:
        return_val = False
            
    return return_val


def save_df_as_tab(df, file):

    """Takes a df and saves it as tab"""

    file_tmp = "%s.tmp"%file
    df.to_csv(file_tmp, sep="\t", index=False, header=True)
    os.rename(file_tmp, file)


def get_date_and_time_for_print():

    """Gets the date of today"""

    current_day = date.today().strftime("%d/%m/%Y")
    current_time = time.strftime("%H:%M:%S", time.localtime())

    return "[%s, %s]"%(current_day, current_time)


def print_with_runtime(*x):

    """prints with runtime info"""

    print(get_date_and_time_for_print(), *x)


def get_tab_as_df_or_empty_df(file):

    """Gets df from file or empty df"""

    nlines = len([l for l in open(file, "r").readlines() if len(l)>1])

    if nlines==0: return pd.DataFrame()
    else: return pd.read_csv(file, sep="\t")

def save_object(obj, filename):
    
    """ This is for saving python objects """

    filename_tmp = "%s.tmp"%filename
    remove_file(filename_tmp)
    
    with open(filename_tmp, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

    os.rename(filename_tmp, filename)

def load_object(filename):
    
    """ This is for loading python objects  in a fast way"""
    
    return pickle.load(open(filename,"rb"))


def make_flat_listOflists(LoL):

    return list(itertools.chain.from_iterable(LoL))


def get_correct_tree(treefile, min_support=0):

    """Gets filtering out branches with low support"""

    # load
    tree = Tree(treefile)

    # collapse low support bracnhes
    for n in tree.traverse():

        # for internal nodes, if the support is below the minimum, set each of the children as children of the parent
        if n.is_leaf() is False and len(n.get_ancestors())>0 and n.support<min_support: n.delete()

    return tree


def get_correct_tree_midpointRooted(treefile, min_support=0):

    """Gets a tree with midpointRooting and filtering out branches with low support"""

    # load
    tree = Tree(treefile)

    # get the midpoint rooted tree
    midpoint =  tree.get_midpoint_outgroup()
    tree.set_outgroup(midpoint)

    # collapse low support bracnhes
    for n in tree.traverse():

        # for internal nodes, if the support is below the minimum, set each of the children as children of the parent
        if n.is_leaf() is False and len(n.get_ancestors())>0 and n.support<min_support: n.delete()

    return tree


def get_vars_df_only_subset_samples_for_GWAS(variants_table, interesting_samples, outdir, varID_field):

    """Gets a df only with the variants for GWAS. It only keeps the varID field and the sampleID field"""

    # define the file
    vars_df_file = "%s/variants_df.py"%outdir

    if file_is_empty(vars_df_file):
        print("getting variants df")

        # load vars
        print("loading vars")
        ending_vars = variants_table.split(".")[-1]
        if ending_vars=="tab": vars_df = get_tab_as_df_or_empty_df(variants_table)
        elif ending_vars=="py": vars_df = load_object(variants_table)
        else: raise ValueError("--variants_table is incorrect. Make sure that it ends with .tab or .py.")

        # log
        print("making extra calculations")

        # keep only vars that are in interesting_samples
        vars_df["sampleID"] = vars_df.sampleID.apply(str)
        vars_df = vars_df[(vars_df.sampleID.isin(interesting_samples))][["sampleID", varID_field]].drop_duplicates().sort_values(by=["sampleID", varID_field])

        # define the samples with vars
        samples_with_vars = set(vars_df.sampleID)
        samples_without_vars = interesting_samples.difference(samples_with_vars)
        if (len(samples_without_vars)/len(interesting_samples))>0.1: print("WARNING: There are %i/%i samples with no variants"%(len(samples_without_vars), len(interesting_samples)))
        
        # remove variants that are in all samples
        var_to_fractionSamplesWithVar = vars_df.groupby(varID_field).apply(len)/len(interesting_samples)
        variable_variants = set(var_to_fractionSamplesWithVar[var_to_fractionSamplesWithVar<1].index)
        print("There are %i/%i vars that are not in all samples"%(len(variable_variants), len(set(vars_df[varID_field]))))
        vars_df = vars_df[vars_df[varID_field].isin(variable_variants)]

        # save
        print("saving")
        save_object(vars_df, vars_df_file)

    return load_object(vars_df_file).reset_index(drop=True)

def get_uniqueVals_df(df): return set.union(*[set(df[col]) for col in df.columns])


def get_pastml_cmds_chunk_mutations(df_chunk, outdir, treefile, pastml_prediction_method):

    """Takes a chunk of a df with mutations and gets the cmds for pastml"""

    # get the cmds for each mutation (row)
    all_cmds = list(df_chunk.apply(get_pastml_cmd_one_mutation, outdir=outdir, treefile=treefile, pastml_prediction_method=pastml_prediction_method, axis=1))
    
    return all_cmds

def get_pastml_cmd_one_mutation(r, outdir, treefile, pastml_prediction_method):

    """Generates a cmd for one mutation that runs pastml. It takes a series where the index is the sampleID and the values are the presence absence"""

    # r.name is the mutation

    # define the outdir
    workdir_mut = "%s/%s_workdir"%(outdir, r.name)
    outfile = "%s/%s_pastml_reconstruction.out"%(outdir, r.name)
    html_file = "%s/%s_pastml_reconstruction.html"%(outdir, r.name)

    # only define cmds if not done
    if file_is_empty(outfile): 
    
        # define the mutations table. This will be removed after pastm run
        mutations_table = "%s/%s_mutation_presence.tab"%(outdir, r.name)
        if file_is_empty(mutations_table): 
            df = pd.DataFrame({"sampleID":r.index, r.name:r.values})
            save_df_as_tab(df, mutations_table)

        # return the cmd, which will remove the mutations_table
        return "%s %s %s %s %s %s %s"%(RUN_PASTML_PY, treefile, mutations_table, workdir_mut, outfile, html_file, pastml_prediction_method)

    else: return None


def check_that_tree_has_no_politomies(tree):

    """Raises an error if tree has politomies"""

    for n in tree.traverse():
        if n.is_leaf() is False and len(n.get_children())!=2: raise ValueError("node %s has !=2 children"%(n.name))


def get_tree_with_internalNodeNames(tree):

    """Adds internal node names"""

    for I, n in enumerate(tree.traverse()): 
        if not n.is_leaf(): n.name = "node_%s_to_%s"%(n.get_leaf_names()[0], n.get_leaf_names()[-1])

    return tree


def run_cmd_simple(cmd):

    """This function runs a cmd"""
    # run
    out_stat = os.system(cmd) 
    if out_stat!=0: raise ValueError("\n%s\n did not finish correctly. Out status: %i"%(cmd, out_stat)) 


def run_cmd(cmd, env="ancestral_GWAS_env"):

    """This function runs a cmd with a given env"""

    # define the cmds
    SOURCE_CONDA_CMD = "source %s/etc/profile.d/conda.sh"%CondaDir
    cmd_prefix = "%s && conda activate %s &&"%(SOURCE_CONDA_CMD, env)

    # define the running
    cmd_to_run = "%s %s"%(cmd_prefix, cmd)

    # run
    out_stat = os.system(cmd_to_run) 
    if out_stat!=0: raise ValueError("\n%s\n did not finish correctly. Out status: %i"%(cmd_to_run, out_stat))

def soft_link_files(origin, target):

    """This function takes an origin file and makes it accessible through a link (target)"""

    if file_is_empty(target):

        # rename as full paths
        origin = get_fullpath(origin)
        target = get_fullpath(target)

        # check that the origin exists
        if file_is_empty(origin): raise ValueError("The origin %s should exist"%origin)

        # remove previous link
        try: run_cmd("rm %s > /dev/null 2>&1"%target)
        except: pass

        soft_linking_std = "%s.softlinking.std"%(target)
        #print("softlinking. The std is in %s"%soft_linking_std)
        run_cmd("ln -s %s %s > %s 2>&1"%(origin, target, soft_linking_std))
        remove_file(soft_linking_std)

    # check that it worked
    if file_is_empty(target): raise ValueError("The target %s should exist"%target)

def get_node_state_all(r):

    """Requires all the values to be either 0 or 1 to be classified as 0 or 1."""

    # define all states
    all_states = set(r.values)

    # only return a state if all the nodes are as in all states
    if len(all_states)==1: return next(iter(all_states))
    else: return np.nan


def get_node_state_all_ignoreNaN(r):

    """The same as get_node_state_all, but ignoring NaNs"""

    # define all states
    all_states = set(r[~pd.isna(r)].values)

    # only return a state if all the nodes are as in all states
    if len(all_states)==1: return next(iter(all_states))
    else: return np.nan


def get_series_ASR_states_one_pastml_output(ID, output_pastml, ASR_methods, classification_method):

    """Gets a series were the index is the 'node' in output pastml and the value is either 1, 0 or NaN. """

    # load pastml df
    if type(output_pastml)==str: df_pastml = get_tab_as_df_or_empty_df(output_pastml)
    else: df_pastml = output_pastml

    # keep only the methods related to ASR_methods
    interesting_columns = list(map(lambda m: "%s_%s"%(ID, m), ASR_methods))
    df_pastml = df_pastml.set_index("node")[interesting_columns].applymap(float)

    # get the states of each node
    classification_method_to_fn = {"all" : get_node_state_all, "all_ignoreNaN":get_node_state_all_ignoreNaN}

    states_series = df_pastml.apply(classification_method_to_fn[classification_method], axis=1)

    return states_series

def get_ASR_states_integrate_methods(names, list_outputs_pastml, ASR_methods_str, classification_method, min_support, node_to_support):

    """Defines each node in list_outputs_pastml (named as names) (all should be the same node as 1, 0 or NaN). It considers only the ASR methods in ASR_methods_str. The criteria used is classification_method. Nodes that have low support (<min_support) will get a state of NaN"""

    # get as list
    ASR_methods = ASR_methods_str.split(",")

    # make sure that the classification method is correct
    allowed_classification_methods = {"all", "all_ignoreNaN"}
    if classification_method not in allowed_classification_methods: raise ValueError("classification_method '%s' is invalid. It should be in %s"%(classification_method, allowed_classification_methods))

    # get the series with the states
    list_states_series = list(map(lambda I: get_series_ASR_states_one_pastml_output(names[I], list_outputs_pastml[I], ASR_methods, classification_method), range(len(names))))

    # check that the indices are the same
    all_indices = set(map(lambda x: tuple(sorted(x.index)), list_states_series))
    if len(all_indices)!=1: raise ValueError("not all the indices are unique")

    # get df
    states_df = pd.DataFrame(dict(zip(names, list_states_series)))

    # modify the states df so that nodes with low support get nan
    def get_r_states_according_to_support(r):
        if node_to_support[r.name]>=min_support: return r
        else: return pd.Series(np.nan, index=r.index)
    states_df = states_df.apply(get_r_states_according_to_support, axis=1)

    return states_df.sort_index()


def generate_resampled_phenotype_states_df(resampled_phenotypes_pastml_out, resampled_phenotype_states_df_file, nresamples, ASR_methods_phenotypes, classification_method_phenotypes, threads, correct_treefile, min_support, node_to_support):

    """Generates a df file that has the states of the phenotypes for a different number of resamples"""

    if file_is_empty(resampled_phenotype_states_df_file):
        print("generating %s"%resampled_phenotype_states_df_file)

        ########## GENERATE LONG DF WITH ALL THE DF STATES #############

        # get the df (and filter by resamples)
        resampled_phenotypes_pastml_df = load_object(resampled_phenotypes_pastml_out)
        print("keeping important resamples")
        resampled_phenotypes_pastml_df = resampled_phenotypes_pastml_df[resampled_phenotypes_pastml_df.resample_I.isin(set(range(1, nresamples+1)))]

        # generate the phenotype states df for each resample
        print("get resampled states  in parallel")
        inputs_fn = [(["phenotype"], [df_I], ASR_methods_phenotypes, classification_method_phenotypes, min_support, node_to_support) for rI, df_I in resampled_phenotypes_pastml_df.groupby("resample_I")]

        def get_df_resampled_states_one_resample(x):
            rI, df_I = x
            df_I["resample_I"] = rI+1
            df_I["node"] = df_I.index
            return df_I

        with multiproc.Pool(threads) as pool:
            df_resampled_states =  pd.concat(list(map(get_df_resampled_states_one_resample,  enumerate(pool.starmap(get_ASR_states_integrate_methods, inputs_fn)))))
            pool.close()
            pool.terminate()

        # old way: 
        #df_resampled_states = resampled_phenotypes_pastml_df.groupby("resample_I").apply(lambda df_I: get_ASR_states_integrate_methods(["phenotype"], [df_I], ASR_methods_phenotypes, classification_method_phenotypes, min_support, node_to_support))
        #df_resampled_states["resample_I"] = df_resampled_states.index.get_level_values(0)

        # add fields
        df_resampled_states["resample_I"] = "resample_" + df_resampled_states["resample_I"].apply(str)
        df_resampled_states = df_resampled_states.reset_index(drop=True)

        # checks
        print("checking")
        if len(set(df_resampled_states.groupby("resample_I").apply(len)))!=1: raise ValueError("not all resamples have the same number of nodes")
        if len(df_resampled_states[["node", "resample_I"]].drop_duplicates())!=len(df_resampled_states): raise ValueError("the combination of node and resample_I should be the same")
        rID_to_nodes = df_resampled_states.groupby("resample_I").apply(lambda df_I: tuple(df_I.node))
        if len(set(rID_to_nodes))!=1: raise ValueError("the set of nodes is not the same in all resamples")

        ################################################################

        ########## GET THE PHENOTYPE PRESENCE DF ###########
        print("getting phenotype presence")
        # this is for phyC

        # generate a df that has each column one resample and each index one state
        df_phenotype = df_resampled_states.pivot(index="node", columns="resample_I", values="phenotype")

        # checks
        df_phenotype_checks = df_phenotype[~df_phenotype.apply(lambda r: all(pd.isna(r)), axis=1)]
        nodes_same_pheno_all_resamples = set(df_phenotype_checks[df_phenotype_checks.applymap(str).apply(set, axis=1).apply(len)==1].index)
        if len(nodes_same_pheno_all_resamples)>0: print("WARNING: There are some nodes (%i/%i) with the same phenotype across all resamplings: %s"%(len(nodes_same_pheno_all_resamples), len(df_phenotype_checks), sorted(nodes_same_pheno_all_resamples)))

        if sum(pd.isna(df_resampled_states.phenotype))!=sum(df_phenotype.apply(pd.isna, axis=1).apply(sum)): raise ValueError("some combinations are not existing in the pivot")

        ####################################################

        ############ GET THE PHENOTYPE TRANSITIONS DF #############
        # this is for syncronous
        print("getting phenotype transitions")

        # map each node to the parent, to get the node transitions
        tree = get_tree_with_internalNodeNames(Tree(correct_treefile))
        tree.name = "root"
        node_to_parent = {n.name : n.get_ancestors()[0].name for n in tree.traverse() if n.name!="root"}

        df_phenotype_transition = df_resampled_states.groupby("resample_I").apply(lambda df_I: get_series_transition(df_I.set_index("node")[["phenotype"]], node_to_parent, "transition")).transpose()

        # checks
        if sorted(df_phenotype.index)!=sorted(df_phenotype_transition.index): raise ValueError("the indices should be the same")
        if sorted(df_phenotype.columns)!=sorted(df_phenotype_transition.columns): raise ValueError("the columns should be the same")

        ###########################################################
        
        # save both
        print("saving")        
        save_object({"phenotype":df_phenotype, "phenotype_transition":df_phenotype_transition}, resampled_phenotype_states_df_file)


def get_series_transition(node_states_df, node_to_parent, type_transition):

    """Gets a series with the node states and returns 1, 0 or NaN if the node is a transition (i.e. it is different) as compared to the parent. type_transition can also be is appearance.  """

    # check tha the node_states_df is one
    if len(node_states_df.columns)!=1: raise ValueError("The node_states_df should have 1 column.")

    # get the series
    col_name = node_states_df.columns[0]
    node_to_state = dict(node_states_df[col_name])

    # define a function that takes the state and returns whether it is a transition
    def get_is_transition(r):

        # the root is nan
        if r.name=="root": return np.nan

        # the others
        parent_state = node_to_state[node_to_parent[r.name]]
        current_state = r[col_name]

        if pd.isna(parent_state) or pd.isna(current_state): return np.nan
        elif parent_state!=current_state: return 1.0
        else: return 0.0

    # define a similar function to get_is_transition, but only concerning appearance
    def get_is_appearance(r):

        # the root is nan
        if r.name=="root": return np.nan

        # the others
        parent_state = node_to_state[node_to_parent[r.name]]
        current_state = r[col_name]

        if pd.isna(parent_state) or pd.isna(current_state): return np.nan
        if parent_state==0.0 and current_state==1.0: return 1.0
        else: return 0.0

    # get series
    type_transition_to_fn = {"transition":get_is_transition, "appearance":get_is_appearance}
    return node_states_df.apply(type_transition_to_fn[type_transition], axis=1)

def get_genotype_transition_from_all_transitions_r(r):

    """Gets a row with several genotype transitions (for different mutations) and returns whether it is a genotype transition"""

    # get all vals
    all_transitions = set(r.values)

    # return as a function
    if 1.0 in all_transitions: return 1.0
    elif all_transitions=={0.0}: return 0.0 # only return 0.0 if there is no ambiguous mutation
    else: return np.nan

def get_node_info_df(genotype_states_df, phenotype_states_df, treefile, sorted_mutations):

    """Gets the df were each row is one node and there are columns for each mutation or phenotype indicating the state. This function creates a df that has information about each node."""

    # print the time
    #start_time = time.time()
    #print("gettng node info df")

    # load the tree with the internal names, and change the root
    tree = get_tree_with_internalNodeNames(Tree(treefile))
    tree.name = "root"

    # check that the names are all the same
    all_nodes = sorted(genotype_states_df.index)
    if {n.name for n in tree.traverse()}!=set(genotype_states_df.index): raise ValueError("Not all nodes are considered")


    # map each node to the parent
    node_to_parent = {n.name : n.get_ancestors()[0].name for n in tree.traverse() if n.name!="root"}

    # init the info df with the node distance
    info_df = pd.DataFrame({"edge_len" : dict(zip(all_nodes , map(lambda n: (tree&n).dist, all_nodes) ))})
    if any(pd.isna(info_df.edge_len)): raise ValueError("There can't be NaNs in edge length")

    # add the phenotype and phenotype transition, meaning that the parent phenotype is different than the curr
    if phenotype_states_df is not None:

        info_df["phenotype"] = phenotype_states_df.phenotype
        info_df["phenotype_transition"] = get_series_transition(phenotype_states_df, node_to_parent, "transition")

    # add the genotype transition
    list_genotypes_transitions = list(map(lambda m: get_series_transition(genotype_states_df[[m]], node_to_parent, "transition"), sorted_mutations))
    genotypes_transitions_df = pd.DataFrame(dict(zip(sorted_mutations, list_genotypes_transitions)))
    info_df["genotype_transition"] = genotypes_transitions_df.apply(get_genotype_transition_from_all_transitions_r, axis=1)

    # add the genotype appearance
    list_genotypes_appearances = list(map(lambda m: get_series_transition(genotype_states_df[[m]], node_to_parent, "appearance"), sorted_mutations))
    genotypes_appearances_df = pd.DataFrame(dict(zip(sorted_mutations, list_genotypes_appearances)))
    info_df["genotype_appearance"] = genotypes_appearances_df.apply(get_genotype_transition_from_all_transitions_r, axis=1)

    # add the node
    info_df["node_name"] = info_df.index

    # check
    if sorted(info_df.index)!=sorted(genotype_states_df.index): raise ValueError("the indices should be the same")
    if not all((pd.isna(genotypes_transitions_df)==pd.isna(genotypes_appearances_df)).apply(all)): raise ValueError("The nan patterns of genotypes_transitions_df and genotypes_appearances_df should be the same")

    # add the mutations
    info_df = info_df.merge(genotype_states_df[sorted_mutations].rename(columns=dict(zip(sorted_mutations, list(map(lambda m: "ASR_%s"%m, sorted_mutations))))).loc[info_df.index], left_index=True, right_index=True, validate="one_to_one")
    #for m in sorted_mutations: info_df["ASR_%s"%m] = genotype_states_df[m]

    # define the mutations to the nodes where the transition is nan
    try: mut_to_NaNtransition_nodes = dict(genotypes_transitions_df.transpose().apply(lambda r: set(r[pd.isna(r)].index), axis=1))
    except:
        print(genotypes_transitions_df)
        raise ValueError("calculation of mut_to_NaNtransition_nodes did not work. See genotypes_transitions_df above.")

    #elapsed_time = float(time.time() - start_time)
    #print("get_node_info_df took %.4f s"%(elapsed_time))

    # return
    return info_df, mut_to_NaNtransition_nodes


def get_GT_transition_rates_and_root_state_one_mutation(mut, ASR_mutations_dir, ASR_methods_mutations, classification_method_mutations, treefile):

    """Calculates, for each mutation the genotype appearance / disappearance rates and the root GT state"""

    # warn
    this_would_need_to_be_adapted_to_the_min_support_filtering

    # get a df with the actual states and transitions
    genotype_states_df = get_ASR_states_integrate_methods([mut], ["%s/%s_pastml_reconstruction.out"%(ASR_mutations_dir, mut)], ASR_methods_mutations, classification_method_mutations)
    node_info_df, mut_to_NaNtransition_nodes = get_node_info_df(genotype_states_df, None, treefile, [mut])

    # checks
    if any(pd.isna(node_info_df.genotype_transition)!=pd.isna(node_info_df.genotype_appearance)): raise ValueError("the GT appearance and transition should be nans together")

    # keep only nodes where there is a transition or an appearance
    node_info_df = node_info_df[~(pd.isna(node_info_df.genotype_transition)) & ~(pd.isna(node_info_df.genotype_appearance))]

    # debug empty dfs
    if len(node_info_df)==0:

        gt_appearance_rate = 0.0
        gt_disappearance_rate = 0.0

    else:

        # define the rates of appearance and transition per branch length unit
        gt_appearance_rate = sum(node_info_df.genotype_appearance)/sum(node_info_df.edge_len)
        gt_transition_rate = sum(node_info_df.genotype_transition)/sum(node_info_df.edge_len)
        gt_disappearance_rate = gt_transition_rate - gt_appearance_rate

        # check rates
        for rate in [gt_disappearance_rate, gt_appearance_rate]:
            if pd.isna(rate): raise ValueError("The rate should not be nan")
            if rate<0: raise ValueError("The rate can't be <0")
            if rate==np.inf:  raise ValueError("The rates should not be inf")

    # define the root GT state as a 0 or 1, important for the simulations
    def get_nan_to_0_or_int(x):
        if pd.isna(x): return 0
        else: return int(x)

    root_GT_state = get_nan_to_0_or_int(genotype_states_df.loc["root", mut])

    # clean
    del genotype_states_df
    del node_info_df
    del mut_to_NaNtransition_nodes

    # return the three as a tuple
    return (mut, gt_appearance_rate, gt_disappearance_rate, root_GT_state)

def get_nr_df_mutation_rates(df_mutation_rates, outdir, threads):

    """Gets the non redundant mutation rates and also generates a df that has the representative mutations"""

    # define files
    df_mutation_rates_nr_file = "%s/df_mutation_rates_nr.py"%outdir
    df_mutation_mapping_file = "%s/mutation_to_representative_mutation_GTsimulations_resampling.py"%outdir

    if file_is_empty(df_mutation_rates_nr_file) or file_is_empty(df_mutation_mapping_file):
        print("getting non-redundant files")

        # map each mutation to an ID
        df_mutation_rates["profile"] = df_mutation_rates.apply(lambda r: "%s_%s_%i"%(r.gt_appearance_rate, r.gt_disappearance_rate, r.root_GT_state), axis=1)
        def get_mutations_from_df_p(df_p): return df_p.mutation.values
        profile_to_muts = df_mutation_rates.groupby("profile").apply(get_mutations_from_df_p).apply(sorted)

        # log
        #print("There are %i NR groups of %i mutations with the same mutation rates"%(len(profile_to_muts), len(df_mutation_rates)))

        # get the mutation_to_representative_mutation_df 
        print("concatenating in %i threads"%threads)

        inputs_fn = list(map(lambda x: (x,), profile_to_muts.values))
        with multiproc.Pool(threads) as pool:

            df_mutation_mapping = pd.concat(pool.starmap(get_mapping_df_from_list, inputs_fn)).sort_values(by=["mutation", "rep_mutation"])
            pool.close()
            pool.terminate()


        # check
        if len(set(df_mutation_mapping.rep_mutation))!=len(profile_to_muts): raise ValueError("The length of the data makes no sense")

        # save the mapping
        save_object(df_mutation_mapping, df_mutation_mapping_file)

        # get the df with non-redundant mutations and save
        nr_mutations = set(df_mutation_mapping.rep_mutation)
        df_mutation_rates_nr = df_mutation_rates[df_mutation_rates.mutation.isin(nr_mutations)][["mutation", "gt_appearance_rate", "gt_disappearance_rate", "root_GT_state"]].rename(columns={"mutation":"rep_mutation"})
        save_object(df_mutation_rates_nr, df_mutation_rates_nr_file)

    return load_object(df_mutation_rates_nr_file).reset_index(drop=True)

def get_simulated_GT_node_info_one_mutation(mut, mutation_simulations_dir, treefile, gt_appearance_rate, gt_disappearance_rate, root_GT_state, ASR_methods_mutations, classification_method_mutations):

    """Generates a node info df (which has genotype_transition and genotype appearance) for one mutation, and saves it to mutation_simulations_dir"""

    prefix = "%s/%s_%s_%s"%(mutation_simulations_dir, mut, ASR_methods_mutations, classification_method_mutations)
    genotype_f_to_file = {genotype_f : "%s_%s_square_df.py"%(prefix, genotype_f) for genotype_f in ["genotype_transition", "genotype_appearance"]}

    if any(list(map(file_is_empty, genotype_f_to_file.values()))):

        # load the tree with the internal names, and change the root
        tree = get_tree_with_internalNodeNames(Tree(treefile))
        tree.name = "root"

        # get the node to dist
        node_to_dist = {node.name : node.dist for node in tree.traverse()}

        # map each node to the parents
        node_to_parent = {n.name : n.get_ancestors()[0].name for n in tree.traverse("levelorder") if n.name!="root"} # the root nodes are traversed first, then the children

        # map each node to the probability of appearance of disappearance of the mutation
        node_to_p_appearance = {n : node_to_dist[n]*gt_appearance_rate for n in node_to_parent}
        node_to_p_disappearance = {n : node_to_dist[n]*gt_disappearance_rate for n in node_to_parent}

        # prepare objects that make the running of the resampling faster
        parent_GTstate_to_n_to_prob1 = {0 : {n : p_app for n, p_app in node_to_p_appearance.items()}, # GT appearance
                                        1 : {n : (1-p_disapp) for n, p_disapp in node_to_p_disappearance.items()}} # GT disappearance

        # define a function that returns 1 according to some probability or else 0
        def get_bernoulli_1(p): return 1 if random.random() <= p else 0

        # generate a df with simulated GT states
        def get_simulated_GTstates_series_dicts(I):

            # define a dict that keeps the GT state of all nodes and also the parents
            node_to_gtState = {"root":root_GT_state}
            node_to_parent_gtState = {"root":0}

            # go through each node and assign GT states and transitions 
            for n, parent in node_to_parent.items():

                # define the gt state of the parent
                parent_GTstate = node_to_gtState[node_to_parent[n]]
                node_to_parent_gtState[n] = parent_GTstate

                # define the state of the current GT according to the model
                node_to_gtState[n] = get_bernoulli_1(parent_GTstate_to_n_to_prob1[parent_GTstate][n])

            return (node_to_gtState, node_to_parent_gtState)

        nresamples = 10000
        list_get_simulated_GTstates_series_dicts = list(map(get_simulated_GTstates_series_dicts, range(nresamples)))

        sim_GT_states_df = pd.DataFrame(list(map(lambda x: x[0], list_get_simulated_GTstates_series_dicts))).transpose()
        sim_GT_states_df_parents = pd.DataFrame(list(map(lambda x: x[1], list_get_simulated_GTstates_series_dicts))).transpose().loc[sim_GT_states_df.index, sim_GT_states_df.columns]

        # create one square df for each field
        for genotype_f, file_df in genotype_f_to_file.items():

            # creeate the square df with True or False for each state
            if genotype_f=="genotype_transition": square_df =  (sim_GT_states_df!=sim_GT_states_df_parents) # current is different than parent
            elif genotype_f=="genotype_appearance": square_df =  (sim_GT_states_df>sim_GT_states_df_parents) # current is 1 and parent is 0
            else: raise ValueError("invalid %s"%genotype_f)

            # save
            save_object(square_df, file_df)

        # clean
        del sim_GT_states_df
        del sim_GT_states_df_parents
        del list_get_simulated_GTstates_series_dicts

          
def get_nodes_GenoAndPheno_from_node_info_df(node_info_df, genotype_f, phenotype_f):

    """Returns the number of nodes were you have both a genotype and a phenotype change"""

    return sum((node_info_df[genotype_f]==1.0) & (node_info_df[phenotype_f]==1.0))


def get_plot_pvalue(observed_nodes_GenoAndPheno, resampled_nodes_GenoAndPheno, genotype_f, phenotype_f, pval_method, nresamples, pval, filename_plot):

    """Plots the pvalue into filename_plot"""

    fig = plt.figure(figsize=(3, 3))
    ax = sns.distplot(resampled_nodes_GenoAndPheno,  hist_kws={"color":"gray", "linewidth":1, "alpha":.3}, kde=False, hist=True)
    plt.axvline(observed_nodes_GenoAndPheno, color="black", linewidth=.7,  linestyle="--")
    ax.set_ylabel("# resamples")
    ax.set_xlabel("# nodes w/ %s & %s"%(genotype_f, phenotype_f))
    ax.set_title("empiric distribution (%s vs %s)\npval method: %s, # resamples=%i\np=%s"%(genotype_f, phenotype_f, pval_method, nresamples, pval))
    print("saving %s"%filename_plot)
    fig.savefig(filename_plot, bbox_inches="tight")
    plt.close(fig)

def get_pvalue_resampling_genoAndPheno_null_distribution(df, nresamples, genotype_f, phenotype_f, pval_method):

    """Gets empiric distribution in a hogwash-like manner"""

    # define the initial index
    initial_node_sorting = list(df.index) # at the beginning the index has the nodes as index

    # rest the index
    df_filt = df.reset_index(drop=True)

    # add the edge length
    df_filt["relative_edge_len"] = df_filt["edge_len"] / sum(df_filt.edge_len)
    #df_filt.loc[df_filt.index, "relative_edge_len"] = df_filt.loc[df_filt.index, "edge_len"] / sum(df_filt.edge_len)
    nnodes = len(df_filt)

    # sort the df by edge length
    df_filt = df_filt.sort_values(by="edge_len", ascending=False)

    # define variables needed for reshuffling
    genotypes_array = np.array(cp.deepcopy(df_filt[genotype_f].values))
    phenotypes_array = np.array(cp.deepcopy(df_filt[phenotype_f].values))

    # define the probabilities of each node
    if pval_method in {"RelToBranchLen", "RelToBranchLen_wReplace"}: probabilities_choice = np.array(cp.deepcopy(df_filt.relative_edge_len.values))
    elif pval_method in {"notRelToBranchLen", "notRelToBranchLen_wReplace"}: probabilities_choice = None
    else: raise ValueError("the %s is not correct"%pval_method)

    # define whether to sample with replacement
    if pval_method in {"RelToBranchLen", "notRelToBranchLen"}: replacement = False
    elif pval_method in {"RelToBranchLen_wReplace", "notRelToBranchLen_wReplace"}: replacement = True
    else: raise ValueError("the %s is not correct"%pval_method)

    # create the null df, where each row is 
    def get_get_resampled_genotypes(I): return np.random.choice(genotypes_array, size=nnodes, replace=replacement, p=probabilities_choice)
    df_null_genotypes = pd.DataFrame(list(map(get_get_resampled_genotypes, range(nresamples))))
    df_null_genotypes.columns = df_filt.node_name

    # define a df with the null phenotypes
    df_null_phenotypes = pd.DataFrame(list(map(lambda x: phenotypes_array, range(nresamples))))
    df_null_phenotypes.columns = df_filt.node_name

    return df_null_genotypes[initial_node_sorting], df_null_phenotypes[initial_node_sorting]


def get_square_df_resampled_gts_any_new_mutation(current_df, df_m):

    """Adds the results of one mutation (in df_m) to current_df. If there is any 1.0, it should stay, if there is a mix of 0s and np.nan it should stay."""

    # integrate the results by vectorization
    def get_consensus_transitions_for_one_resample(I): return pd.DataFrame({"current":current_df[I], "new":df_m[I]}).apply(get_genotype_transition_from_all_transitions_r, axis=1)

    consensus_df = pd.DataFrame(dict(zip(current_df.columns, map(get_consensus_transitions_for_one_resample, current_df.columns))))[current_df.columns]

    return consensus_df

def get_simulated_GT_node_info_df_one_mutation(mut, genotype_states_df_m, treefile, genotype_f, nresamples, sorted_interesting_nodes):

    """Generates a node info df (which has genotype_transition and genotype appearance) for one mutation, and saves it to mutation_simulations_dir"""

    # get the node info df for this mutation
    node_info_df, mut_to_NaNtransition_nodes = get_node_info_df(genotype_states_df_m, None, treefile, [mut])

    # define interesting nodes (those where permutation will be performed)
    interesting_nodes = set(sorted_interesting_nodes)

    # keep only nodes of interest
    node_info_df = node_info_df.loc[sorted_interesting_nodes]
    node_info_df = node_info_df[~(pd.isna(node_info_df.genotype_transition)) & ~(pd.isna(node_info_df.genotype_appearance))]

    if len(node_info_df)==0:
        gt_appearance_rate = 0.0
        gt_disappearance_rate = 0.0

    else:

        # define the rates of appearance and transition and per branch length unit
        total_edge_len = sum(node_info_df.edge_len)

        gt_appearance_rate = sum(node_info_df.genotype_appearance==1.0)/total_edge_len
        gt_transition_rate = sum(node_info_df.genotype_transition==1.0)/total_edge_len
        gt_disappearance_rate = gt_transition_rate - gt_appearance_rate

    # check rates
    for rate in [gt_disappearance_rate, gt_appearance_rate]:
        if pd.isna(rate): raise ValueError("The rate should not be nan")
        if rate<0: raise ValueError("The rate can't be <0")
        if rate==np.inf:  raise ValueError("The rates should not be inf")

    # define the root GT state as a 0 or 1, important for the simulations
    def get_nan_to_0_or_int(x):
        if pd.isna(x): return 0
        else: return int(x)

    root_GT_state = get_nan_to_0_or_int(genotype_states_df_m.loc["root", mut])

    # load the tree with the internal names, and change the root
    tree = get_tree_with_internalNodeNames(Tree(treefile))
    tree.name = "root"

    # get the node to dist
    node_to_dist = {node.name : node.dist for node in tree.traverse()}

    # map each node to the parents
    node_to_parent = {n.name : n.get_ancestors()[0].name for n in tree.traverse("levelorder") if n.name!="root"} # the root nodes are traversed first, then the children

    # map each node to the probability of appearance of disappearance of the mutation
    node_to_p_appearance = {n : node_to_dist[n]*gt_appearance_rate for n in node_to_parent}
    node_to_p_disappearance = {n : node_to_dist[n]*gt_disappearance_rate for n in node_to_parent}

    # checks
    if any([p>=1 for p in node_to_p_appearance.values()]): raise ValueError("the p can't be >=1 in node_to_p_appearance")
    if any([p>=1 for p in node_to_p_disappearance.values()]): raise ValueError("the p can't be >=1 in node_to_p_disappearance")

    # prepare objects that make the running of the resampling faster
    parent_GTstate_to_n_to_prob1 = {0 : {n : p_app for n, p_app in node_to_p_appearance.items()}, # GT appearance
                                    1 : {n : (1-p_disapp) for n, p_disapp in node_to_p_disappearance.items()}} # GT disappearance

    # define a function that returns 1 according to some probability or else 0
    def get_bernoulli_1(p): return 1 if random.random() <= p else 0

    # generate a df with simulated GT states
    def get_simulated_GTstates_series_dicts(I):

        # define a dict that keeps the GT state of all nodes and also the parents
        node_to_gtState = {"root":root_GT_state}
        node_to_parent_gtState = {"root":0}

        # go through each node and assign GT states and transitions 
        for n, parent in node_to_parent.items():

            # define the parent
            parent_n = node_to_parent[n]

            # define the gt state of the parent
            parent_GTstate = node_to_gtState[parent_n]
            node_to_parent_gtState[n] = parent_GTstate

            # define the state of the current GT according to the model
            if n in interesting_nodes:  node_to_gtState[n] = get_bernoulli_1(parent_GTstate_to_n_to_prob1[parent_GTstate][n])

            # for nodes that are not considered, propagate the parent tree
            else: node_to_gtState[n] = parent_GTstate

        return (node_to_gtState, node_to_parent_gtState)

    list_get_simulated_GTstates_series_dicts = list(map(get_simulated_GTstates_series_dicts, range(nresamples)))

    sim_GT_states_df = pd.DataFrame(list(map(lambda x: x[0], list_get_simulated_GTstates_series_dicts))).transpose().loc[sorted_interesting_nodes]
    sim_GT_states_df_parents = pd.DataFrame(list(map(lambda x: x[1], list_get_simulated_GTstates_series_dicts))).transpose().loc[sim_GT_states_df.index, sim_GT_states_df.columns]


    # creeate the square df with True or False for each state
    if genotype_f=="genotype_transition": square_df =  (sim_GT_states_df!=sim_GT_states_df_parents) # current is different than parent
    elif genotype_f=="genotype_appearance": square_df =  (sim_GT_states_df>sim_GT_states_df_parents) # current is 1 and parent is 0
    else: raise ValueError("invalid %s"%genotype_f)

    # convert to floats
    square_df = fast_applymap(square_df, float)

    # checks and return
    if sorted(square_df.index)!=sorted(sorted_interesting_nodes): raise ValueError("the square df should have sorted_interesting_nodes as indices")

    # print the number of GTs and median GTs
    resample_to_nGTs = (square_df==1.0).apply(sum, axis=0)
    median_GTs_resamples = np.median(resample_to_nGTs)
    observed_nGTs = sum(node_info_df[genotype_f]==1.0)
    print("median_GTs_resamples: %s. Observed GT: %s"%(median_GTs_resamples, observed_nGTs))

    return square_df.loc[sorted_interesting_nodes]


def get_list_df_resampled_genotypes_simulations(sorted_mutations, genotype_states_df, node_info_df_filt, treefile, genotype_f, phenotype_f, nresamples):

    """gets a set of mutations and their states in different nodes and returns a list with resampled genotypes"""

    print("getting list_df_resampled_genotypes_simulations")

    # check that the nodes make sense
    if any(pd.isna(node_info_df_filt[genotype_f])) or any(pd.isna(node_info_df_filt[phenotype_f])): raise ValueError("there are nans in node_info_df_filt")

    # define all the nodes
    all_nodes = set(node_info_df_filt.index)
    sorted_nodes = list(node_info_df_filt.index)

    # get the list of resampled dfs, one for each mutation
    list_df_resampled_genotypes = list(map(lambda m: get_simulated_GT_node_info_df_one_mutation(m, genotype_states_df[[m]], treefile, genotype_f, nresamples, sorted_nodes), sorted_mutations))

    # check that the columns and indices are all the same
    all_cols_set = set(map(lambda d: tuple(d.columns), list_df_resampled_genotypes))
    all_indices_set = set(map(lambda d: tuple(d.index), list_df_resampled_genotypes))
    if all_cols_set!={tuple(range(nresamples))}: raise ValueError("all the cols should be a tuple from 0 to 10000")
    if all_indices_set!={tuple(sorted_nodes)}: raise ValueError("the indices should be those of node_info_df")

    print("list_df_resampled_genotypes_simulations got")


    return list_df_resampled_genotypes

def get_pvalue_resampling_simulation_null_distribution(df, nresamples, genotype_f, phenotype_f, treefile, sorted_mutations, genotype_states_df):

    """Calculates an empiric pvalue using simulations of the GT values"""

    this_is_all_wrong_you_shoul_also_consider_linkage

    # genetaye a list of resampled genotypes in a simulated manner, one fore each mutation
    list_df_resampled_genotypes = get_list_df_resampled_genotypes_simulations(sorted_mutations, genotype_states_df, df, treefile, genotype_f, phenotype_f, nresamples)

    # keep reduced df of the node info
    df = df[['node_name', phenotype_f]]

    # define the square df that has the genotype_f for each transition. Each row will be a node, each column a resample, and the value a float that indicates if there is some genotype_f. It can also be NaN
    square_df_resampled_gts_any = list_df_resampled_genotypes[0] # init with the first mutation

    if len(list_df_resampled_genotypes)>1: # add additional mutations
        for df_res in list_df_resampled_genotypes[1:]: square_df_resampled_gts_any = get_square_df_resampled_gts_any_new_mutation(square_df_resampled_gts_any, df_res)

    # checks
    if list(square_df_resampled_gts_any.index)!=list(df.index): raise ValueError("the df and the square df sould have the same index")

    # define the phenotypes array
    phenotypes_array =  np.array(cp.deepcopy(df[phenotype_f].values))

    # define a df with the null genotypes, based on resamplings
    df_null_genotypes = square_df_resampled_gts_any.transpose()
    df_null_genotypes.index = list(range(nresamples))

    # define a df with the null phenotypes
    df_null_phenotypes = pd.DataFrame(list(map(lambda x: phenotypes_array, range(nresamples))))
    df_null_phenotypes.columns = df.index

    return df_null_genotypes, df_null_phenotypes

def get_pvalue_resampling_phenotypes_null_distribution(df, nresamples, genotype_f, phenotype_f, pval_method, df_phenotypes_resampled, nodes_short_branches):

    """Takes resampled phenotypes (from nresamples_to_phenotype_f_to_dfResampling) and obtains an emprirical distribution"""

    # get the df_phenotypes_resampled with the nodes of the df (already filtered to only have data for genos and phenos that have some data)
    df_phenotypes_resampled = df_phenotypes_resampled.loc[df.index]

    # check that they are unique
    if len(df_phenotypes_resampled)!=len(set(df_phenotypes_resampled.index)): raise ValueError("df_phenotypes_resampled lacks a unique index")
    if len(df)!=len(set(df.index)): raise ValueError("df lacks a unique index")

    # check that the indices make sense
    if sorted(df.index)!=sorted(df_phenotypes_resampled.index): raise ValueError("the node names are not the same")

    # get the df where there are no nans
    df_filt = df[~pd.isna(df[genotype_f])]

    # keep the df with resampled phenotypes with thes same nodes
    df_phenotypes_resampled = df_phenotypes_resampled.loc[df_filt.index]
    if list(df_phenotypes_resampled.index)!=list(df_filt.index): raise ValueError("the indices should be the same")

    # define the genotypes array
    genotypes_array =  np.array(cp.deepcopy(df_filt[genotype_f].values))

    # define a df with the null genotypes (always the same)
    df_null_genotypes = pd.DataFrame(list(map(lambda x: genotypes_array, range(nresamples))))
    df_null_genotypes.columns = df_filt.index

    # define the null phenotypes, as defined in df_phenotypes_resampled
    df_null_phenotypes = df_phenotypes_resampled.transpose()
    df_null_phenotypes.index = list(range(nresamples))

    return df_null_genotypes, df_null_phenotypes


def get_grouping_df_mut_and_rep_mut(grouping_table, repMuts_df, vars_IDmapping_df, varID_field, interesting_groups):

    """
    Gets a grouping df with representative mutations
    """
    print("get grouping_df")

    # load the df of the groups
    if grouping_table.endswith(".tab"): grouping_df = get_tab_as_df_or_empty_df(grouping_table)
    elif grouping_table.endswith(".py"): grouping_df = load_object(grouping_table)

    # keep interesting groups
    grouping_df = grouping_df[grouping_df.group.isin(interesting_groups)]
    
    # keep only mutually containing dfs
    vars_IDmapping_df = vars_IDmapping_df[vars_IDmapping_df[varID_field].isin(set(grouping_df[varID_field]))]
    grouping_df = grouping_df[grouping_df[varID_field].isin(set(vars_IDmapping_df[varID_field]))]

    # add the numeric mutation in the grouping_df
    var_to_numericMutation = dict(vars_IDmapping_df.set_index(varID_field)["numeric_mutation"])
    grouping_df["numeric_mutation"] = grouping_df[varID_field].map(var_to_numericMutation); check_no_nans_series(grouping_df["numeric_mutation"])

    # add the representative mut
    mut_to_repMut = dict(repMuts_df.set_index("mutation").rep_mutation)
    grouping_df["rep_numeric_mutation"] = grouping_df.numeric_mutation.map(mut_to_repMut); check_no_nans_series(grouping_df["rep_numeric_mutation"])

    # set as index the group 
    grouping_df = grouping_df.set_index("group")

    return grouping_df


def get_GWAS_group_to_mutations(sorted_mutations, grouping_table, varID_field, outdir_ASR_mutations):

    """Maps each group to a list of interesting mutations."""

    if grouping_table is None: group_to_mutations = dict(zip(sorted_mutations, map(lambda x: [x], sorted_mutations))) # each group is just one mutation
    else:

        # load the df of the groups
        print("grouping mutations")
        if grouping_table.endswith(".tab"): grouping_df = get_tab_as_df_or_empty_df(grouping_table)
        elif grouping_table.endswith(".py"): grouping_df = load_object(grouping_table)

        # if there are no groups, return
        if len(grouping_df)==0: raise ValueError("the grouping file %s has no groups"%(grouping_table))

        # load the representative_muts_df
        repMuts_df = get_tab_as_df_or_empty_df("%s/mutation_to_representative_mutation.tab"%outdir_ASR_mutations)
        print("There are %i mutations and %i representative mutations"%(len(set(repMuts_df.mutation)), len(set(repMuts_df.rep_mutation))))

        # load the variants ID mapping and keep only the ones that are in some grouping
        vars_IDmapping_df = get_tab_as_df_or_empty_df("%s/variants_IDmapping.tab"%outdir_ASR_mutations)
        if set(repMuts_df.mutation)!=set(vars_IDmapping_df.numeric_mutation): raise ValueError("The varianst in vars_IDmapping_df should be the same as in repMuts_df")

        vars_IDmapping_df = vars_IDmapping_df[vars_IDmapping_df[varID_field].isin(set(grouping_df[varID_field]))]

        # keep only the variants in grouping_df that are also in  variants_IDmapping. Maybe some of the grouping were not considered because they don't pass the filter of being in all samples
        grouping_df = grouping_df[grouping_df[varID_field].isin(set(vars_IDmapping_df[varID_field]))]

        # add the numeric mutation in the grouping_df
        var_to_numericMutation = dict(vars_IDmapping_df.set_index(varID_field)["numeric_mutation"])
        grouping_df["numeric_mutation"] = grouping_df[varID_field].map(var_to_numericMutation)
        if any(pd.isna(grouping_df["numeric_mutation"])): raise ValueError("There can't be NaNs in numeric_mutation. This means that there are mutations in the grouping that are not in the vars_IDmapping_df")

        # keep only some fields
        grouping_df = grouping_df[["group", "numeric_mutation"]]

        # change the numeric_mutation from grouping_df so that it includes the rep_mut
        mut_to_repMut = dict(repMuts_df.set_index("mutation").rep_mutation)
        grouping_df["numeric_mutation"] = grouping_df.numeric_mutation.map(mut_to_repMut)
        if any(pd.isna(grouping_df["numeric_mutation"])): raise ValueError("There can't be NaNs in numeric_mutation. This means that there are mutations in the grouping that don't have a representative mutation")

        # check that the mutations are in sorted_mutations
        unexpected_mutations = set(grouping_df.numeric_mutation).difference(set(sorted_mutations))
        if len(unexpected_mutations)>0: raise ValueError("There are some mutations for which you did not run the ASR: %s (%i mutations)"%(unexpected_mutations, len(unexpected_mutations)))

        # map each group to the mutations
        grouping_df = grouping_df[["group", "numeric_mutation"]].drop_duplicates()
        group_to_mutations = dict(grouping_df.groupby("group").apply(lambda df_g: sorted(df_g.numeric_mutation)))

    return group_to_mutations


def get_files_significant_GWAS_results(r, significant_groups_ASR_dir, significant_groups_plots_dir, treefile, numericVar_to_var, grouping_df, GWAS_results_dir, mut_to_repMut, varID_field, skip_tree_rendering):

    """Gets the GWAS results for one row that is significant and saves files """

    # define the group_name without '/'
    corrected_group_name = r.group_name.replace("/", "_")

    ######## GET THE DF WITH CHANGED MUTATION NAMES ########

    # keep the file with the ASR
    asr_df_file = "%s/%s_ASR.tab"%(significant_groups_ASR_dir, corrected_group_name)
    if file_is_empty(asr_df_file):

        # init the fields to remove
        fields_to_remove = []

        # load the df with the ASR of this group
        if grouping_df is None: 

            # define the representative mutation, which is in GWAS_results_dir
            real_variant = r.group_name
            rep_mut = mut_to_repMut[r.numeric_mutation]

            # load ASR of the representative mut
            asr_df = get_tab_as_df_or_empty_df("%s/%s_ASR.tab"%(GWAS_results_dir, rep_mut))

            # redefine the field to include the non-rep mutation
            numeric_mutation_f = "ASR_%s"%rep_mut
            asr_df["ASR_%s"%real_variant] = asr_df[numeric_mutation_f]
            fields_to_remove.append(numeric_mutation_f)

        else:

            # load ASR of the group (which has the representative numeric mutations)
            asr_df = get_tab_as_df_or_empty_df("%s/%s_ASR.tab"%(GWAS_results_dir, r.group_name))

            # go through each mutation associated to this group and create a column of it based on the representative field
            for g_name, r_group in grouping_df.loc[{r.group_name}].iterrows():

                # define the field of asr_df that has this mutation, and flag for removal
                numeric_mutation_f = "ASR_%s"%(r_group.rep_numeric_mutation)
                fields_to_remove.append(numeric_mutation_f)

                # create a new column with the variant name
                asr_df["ASR_%s"%(r_group[varID_field])] = asr_df[numeric_mutation_f]

        # remove unnecessary fields
        asr_df = asr_df.drop(columns=fields_to_remove)

        # write 
        save_df_as_tab(asr_df, asr_df_file)

    asr_df = get_tab_as_df_or_empty_df(asr_df_file).set_index("node_name")

    ########################################################


    ######## GET A TREE WITH THE MUTATION INFO ##########

    # define file, and only run if not available
    tree_plot = "%s/%s_%s_results.pdf"%(significant_groups_plots_dir, corrected_group_name, r.method)
    if file_is_empty(tree_plot) and skip_tree_rendering is False:        
        print("Generating tree plot with the data")

        # load the tree with the internal names, and change the root
        tree = get_tree_with_internalNodeNames(Tree(treefile))
        tree.name = "root"

        # define graphics
        phenotype_type_to_color = {"1.0":"red", "0.0":"blue", "nan":"gray"}
        genotype_type_to_color = {"1.0":"red", "0.0":"blue", "nan":"gray"}
        width_rect = 20
        radius_circle = 10

        for n in tree.traverse():

            # define the genotype and phenotype types
            pheno = str(asr_df.loc[n.name, r.phenotype_f])
            geno = str(asr_df.loc[n.name, r.genotype_f])

            # define the color of the branch
            if pheno=="1.0" and geno=="1.0": branch_color = "red"
            elif pheno=="0.0" and geno=="0.0": branch_color  = "blue"
            elif pheno!=geno: branch_color = "purple"
            elif pheno=="nan" or geno=="nan": branch_color = "gray"
            else: raise ValueError("geno/pheno %s/%s combination not permited"%(geno, pheno))

            # set the lines
            nst = NodeStyle()
            nst["hz_line_width"] = 6
            nst["vt_line_width"] = 6
            nst["size"] = 0
            nst["hz_line_color"] = branch_color
            nst["vt_line_color"] = branch_color
            n.set_style(nst)

            # add a square for the phenotype
            n.add_face(RectFace(width_rect, width_rect, fgcolor="black", bgcolor=phenotype_type_to_color[pheno]), column=0)

            # add a circle for the genotype
            n.add_face(CircleFace(radius_circle, genotype_type_to_color[geno], style='sphere'), column=0)


        # get the treestyle
        ts = TreeStyle()
        ts.show_branch_length = False
        ts.show_branch_support = False
        ts.show_leaf_name = False
        ts.draw_guiding_lines = True
        ts.guiding_lines_type = 2 # 0=solid, 1=dashed, 2=dotted.
        ts.guiding_lines_color = "black" 
        ts.legend_position = 3 # bottom-left

        # add the legend for each box
        for pheno, color in phenotype_type_to_color.items():
            ts.legend.add_face(RectFace(width_rect, width_rect, fgcolor="black", bgcolor=color), column=0)
            ts.legend.add_face(TextFace(" %s%s"%({"1.0":"", "0.0":"no ", "nan":"unkown "}[pheno], r.phenotype_f), bold=True, fsize=10), column=1)

        for geno, color in genotype_type_to_color.items():
            ts.legend.add_face(CircleFace(radius_circle, color, style='sphere'), column=0)
            ts.legend.add_face(TextFace(" %s%s"%({"1.0":"", "0.0":"no ", "nan":"unkown "}[geno], r.genotype_f), bold=True, fsize=10), column=1)

        # add the title with the metadata of how this analysis was done
        for f in fields_GWAS_stats_df:

            ts.title.add_face(TextFace("%s: %s"%(f, r[f]), fsize=10), column=0)

        # write
        print("rendering %s"%tree_plot)
        tree_plot_tmp = "%s.tmp.pdf"%tree_plot
        tree.render(file_name=tree_plot_tmp, tree_style=ts)
        os.rename(tree_plot_tmp, tree_plot)

    #####################################################


def rsync_file(origin, dest):

    """syncs one file to the other"""

    if file_is_empty(dest):

        dest_tmp = "%s.tmp"%dest
        run_cmd("rsync --copy-links %s %s"%(origin, dest_tmp))
        os.rename(dest_tmp, dest)

def get_all_cmds_ASR_mutations(square_vars_df, sorted_mutations, threads_generate_cmds, outdir_ASR, correct_treefile, pastml_prediction_method):

    """Gets all the cmds for ASR mutations in parallel, timing the operation"""

    # transpose df and keep only the desired mutations
    square_vars_df_t = square_vars_df.transpose().loc[sorted_mutations]

    # init time
    start_time = time.time()

    # redefine the threads
    threads_generate_cmds = min([threads_generate_cmds, multiproc.cpu_count()])

    # make outdir
    make_folder(outdir_ASR)

    # add a 'chunkID' to square_vars_df_t to process in parallel each chunk
    chunk_size = int(len(square_vars_df_t)/threads_generate_cmds) + 1
    square_vars_df_t["chunkID"] = make_flat_listOflists([[I+1]*chunk_size for I in range(threads_generate_cmds)])[0:len(square_vars_df_t)]

    # get cmds in parallel
    inputs_fn = list(map(lambda x: (x[1], outdir_ASR, correct_treefile, pastml_prediction_method), square_vars_df_t.groupby("chunkID")))
    with multiproc.Pool(threads_generate_cmds) as pool:

        all_cmds = make_flat_listOflists(pool.starmap(get_pastml_cmds_chunk_mutations, inputs_fn))
        pool.close()
        pool.terminate()
        pool.join()


    # remove the nones
    all_cmds = [cmd for cmd in all_cmds if cmd is not None]

    # calculate the time
    elapsed_time = float(time.time() - start_time)

    return all_cmds, elapsed_time

def get_mapping_df_from_list(muts):

    """Takes a list of mutations and returns a df"""

    df = pd.DataFrame({"mutation":muts})
    df["rep_mutation"] = muts[0]

    return df

def get_nr_square_vars_df(square_vars_df, outdir, threads):

    """Gets the square vars df and gets the non redundant ones, keeping a file that maps both of them"""

    # define the mapping file
    mutation_mapping_file = "%s/mutation_to_representative_mutation.tab"%outdir
    if file_is_empty(mutation_mapping_file):
        print("Getting only non-redundant mutations")

        # check df
        if any((square_vars_df==0).apply(all, axis=0)): raise ValueError("There are mutations that are in no samples")
        if any((square_vars_df==1).apply(all, axis=0)): raise ValueError("There are mutations that are all samples")

        # map each mutation to a profile
        mut_to_profile = square_vars_df.applymap(str).agg("".join, axis=0)

        # map each profile to the mutations 
        print("grouping profiles")
        mutations_df = pd.DataFrame({"profile":mut_to_profile, "mutation":mut_to_profile.index}).reset_index()
        def get_mutations_from_df_p(df_p): return df_p.mutation.values
        profile_to_muts = mutations_df.groupby("profile").apply(get_mutations_from_df_p).apply(sorted)

        # get the mutation_to_representative_mutation_df 
        print("concatenating in %i threads"%threads)

        inputs_fn = list(map(lambda x: (x,), profile_to_muts.values))
        with multiproc.Pool(threads) as pool:

            mutation_mapping_df = pd.concat(pool.starmap(get_mapping_df_from_list, inputs_fn)).sort_values(by=["mutation", "rep_mutation"])
            pool.close()
            pool.terminate()

        save_df_as_tab(mutation_mapping_df, mutation_mapping_file)

    # load df
    mutation_mapping_df = get_tab_as_df_or_empty_df(mutation_mapping_file)

    # define the non-redundant mutations from this mappinh the 
    nr_mutations = sorted(set(mutation_mapping_df.rep_mutation))
    nr_square_vars_df = square_vars_df[nr_mutations]
    print("There are %i/%i NR mutations"%(len(nr_mutations), len(square_vars_df.columns)))

    return nr_square_vars_df

def print_available_mem():

    """prints the available mem"""

    available_mem = (psutil.virtual_memory().available)/1e9
    print("There are %.3f Gb of mem available"%available_mem)


def get_x_as_bool(x):
    
    """Returns bool or x"""

    if pd.isna(x): return False
    else: return bool(x)


def fast_applymap(df, fun_apply):

    """Runs applymap on df in a way that is fast"""

    return pd.DataFrame(list(map(lambda idx: df.loc[idx].apply(fun_apply), df.index)))



def calculate_epsilon(nodes_GenoAndPheno, nodes_withGeno, nodes_withPheno):

    """Calculates epsilon"""

    numerator_epsilon = 2*nodes_GenoAndPheno
    denominator_epsilon = nodes_withGeno + nodes_withPheno
    if denominator_epsilon==0: epsilon = np.nan
    else: epsilon = numerator_epsilon / denominator_epsilon

    return epsilon

def calculate_chi_square(nodes_GenoAndPheno, nodes_GenoAndNoPheno, nodes_noGenoAndPheno, nodes_noGenoAndNoPheno):

    """Calculates a chi square statistic for the contingency table 

    nodes_GenoAndPheno, nodes_GenoAndNoPheno
    nodes_noGenoAndPheno, nodes_noGenoAndNoPheno
    
    If correction=True, the Yates correction is applied. This correction reduced the chi square'd statistic, which can be overly conservative for small sample sizes
    """

    # return 0 for cases with not enough variation
    if nodes_GenoAndPheno==0 or nodes_noGenoAndNoPheno==0: return 0.0

    # make the chi square test
    contingency_table = [[nodes_GenoAndPheno, nodes_GenoAndNoPheno], [nodes_noGenoAndPheno, nodes_noGenoAndNoPheno]]
    try: chi_val, p_val, dof, expected =  chi2_contingency(contingency_table, correction=False)
    except: raise ValueError("chi squared failed. This is the table:\n%s\n"%contingency_table)

    return chi_val

def calculate_OR(nodes_GenoAndPheno, nodes_GenoAndNoPheno, nodes_noGenoAndPheno, nodes_noGenoAndNoPheno):

    """Calculates a OR statistic for the contingency table 

    nodes_GenoAndPheno, nodes_GenoAndNoPheno
    nodes_noGenoAndPheno, nodes_noGenoAndNoPheno
    
    """

    # get OR
    contingency_table = [[nodes_GenoAndPheno, nodes_GenoAndNoPheno], [nodes_noGenoAndPheno, nodes_noGenoAndNoPheno]]
    numerator_OR = nodes_GenoAndPheno*nodes_noGenoAndNoPheno
    denominator_OR = nodes_GenoAndNoPheno*nodes_noGenoAndPheno

    if denominator_OR==0 and numerator_OR==0: return 0
    elif denominator_OR==0: return np.inf
    else: return numerator_OR/denominator_OR
    

def get_dict_stats_one_method_gwas(node_info_df, genotype_f, phenotype_f, method_name, ASR_methods_phenotypes, ASR_methods_mutations, classification_method_phenotypes, classification_method_mutations, name, sorted_mutations, outdir_ASR_mutations, nresamples_to_filePhenotypesResampling, nodes_short_branches, treefile, genotype_states_df, phenotypes_1000resamples_stats_prefix):

    """Gets a dict with all the GWAS statistics of one method"""

    # keep dfs
    node_info_df = cp.deepcopy(node_info_df)

    # filter to keep only interesting nodes

    # calculate the number of nodes that have each genotype of phenotype change
    nodes_withGeno = sum(node_info_df[genotype_f]==1.0)
    nodes_withoutGeno = sum(node_info_df[genotype_f]==0.0)
    nodes_unkownGeno = sum(pd.isna(node_info_df[genotype_f]))

    nodes_withPheno = sum(node_info_df[phenotype_f]==1.0)
    nodes_withoutPheno = sum(node_info_df[phenotype_f]==0.0)
    nodes_unkownPheno = sum(pd.isna(node_info_df[phenotype_f]))

    # calculate the number of nodes that have both a genotype and a phenotype change, similar to the N in hogwash.
    nodes_GenoAndPheno = get_nodes_GenoAndPheno_from_node_info_df(node_info_df, genotype_f, phenotype_f)

    # calculate other types of nodes (maybe useful for OR calculations)
    nodes_noGenoAndPheno = sum((node_info_df[genotype_f]==0.0) & (node_info_df[phenotype_f]==1.0))
    nodes_GenoAndNoPheno = sum((node_info_df[genotype_f]==1.0) & (node_info_df[phenotype_f]==0.0))
    nodes_noGenoAndNoPheno = sum((node_info_df[genotype_f]==0.0) & (node_info_df[phenotype_f]==0.0))

    # define a filtered df, from which to calculate some stats
    node_info_df_filt = node_info_df[~(pd.isna(node_info_df[genotype_f])) & ~(pd.isna(node_info_df[phenotype_f]))]

    # calculate some stats
    epsilon = calculate_epsilon(nodes_GenoAndPheno, sum(node_info_df_filt[genotype_f]==1.0), sum(node_info_df_filt[phenotype_f]==1.0))
    observed_OR = calculate_OR(nodes_GenoAndPheno, nodes_GenoAndNoPheno, nodes_noGenoAndPheno, nodes_noGenoAndNoPheno)

    # calculate chi square
    observed_chi_square = calculate_chi_square(nodes_GenoAndPheno, nodes_GenoAndNoPheno, nodes_noGenoAndPheno, nodes_noGenoAndNoPheno)

    # init a dict for this method name
    dict_stats = {"method":method_name, "phenotype_f":phenotype_f, "genotype_f":genotype_f, "nodes_withGeno":nodes_withGeno, "nodes_withoutGeno":nodes_withoutGeno, "nodes_unkownGeno":nodes_unkownGeno, "nodes_withPheno":nodes_withPheno, "nodes_withoutPheno":nodes_withoutPheno, "nodes_unkownPheno":nodes_unkownPheno, "nodes_GenoAndPheno":nodes_GenoAndPheno, "nodes_noGenoAndPheno":nodes_noGenoAndPheno, "nodes_GenoAndNoPheno":nodes_GenoAndNoPheno, "nodes_noGenoAndNoPheno":nodes_noGenoAndNoPheno, "epsilon":epsilon, "group_name":name, "ASR_methods_phenotypes":ASR_methods_phenotypes, "ASR_methods_mutations":ASR_methods_mutations, "classification_method_phenotypes":classification_method_phenotypes, "classification_method_mutations":classification_method_mutations, "OR":observed_OR, "chi_square":observed_chi_square}

    # add the pvals, calculated in different ways
    for pval_method in all_pval_methods:

        # only calculate p value for convergent intems, where there is more genotype appearance / transition in nodes that have the phenotype than nodes that don't have it
        if (nodes_GenoAndPheno>=2) and (nodes_noGenoAndNoPheno>0) and (observed_OR>1):
            
            # check that the epsilon makes sense
            if pd.isna(epsilon) or epsilon<0 or epsilon>1: raise ValueError("epsilon %s is invalid"%epsilon)

            # define a filtered node_info_df, which should be used to calculate the resamples
            if len(node_info_df_filt)<2: raise ValueError("there should be at least 2 rows in node_info_df_filt")
        
            # calculate p values for different numbers of resamples
            for nresamples in [1000, 10000]:

                # get the null distribution
                if pval_method in {"RelToBranchLen", "notRelToBranchLen", "RelToBranchLen_wReplace", "notRelToBranchLen_wReplace"}: df_null_genotypes, df_null_phenotypes = get_pvalue_resampling_genoAndPheno_null_distribution(node_info_df_filt[[genotype_f, phenotype_f, "edge_len", "node_name"]], nresamples, genotype_f, phenotype_f, pval_method)

                elif pval_method=="phenotypes":  df_null_genotypes, df_null_phenotypes = get_pvalue_resampling_phenotypes_null_distribution(node_info_df_filt[[genotype_f]], nresamples, genotype_f, phenotype_f, pval_method, load_object(nresamples_to_filePhenotypesResampling[nresamples])[phenotype_f], nodes_short_branches)

                #elif pval_method=="simulation":  df_null_genotypes, df_null_phenotypes = get_pvalue_resampling_simulation_null_distribution(node_info_df_filt, nresamples, genotype_f, phenotype_f, treefile, sorted_mutations, genotype_states_df)

                else: raise ValueError("pval method %s is invalid"%pval_method)

                # checks
                if len(df_null_genotypes)!=nresamples: raise ValueError("df_null_genotypes.index should be as nresamples")
                if len(df_null_genotypes.columns)!=len(node_info_df_filt): raise ValueError("df_null_genotypes.columns should be as node_info_df_filt")
                if list(df_null_genotypes.index)!=list(df_null_phenotypes.index): raise ValueError("geno and pheno dfs should have the same indices")
                if list(df_null_genotypes.columns)!=list(df_null_phenotypes.columns): raise ValueError("geno and pheno dfs should have the same columns")
                if list(df_null_genotypes.columns)!=list(node_info_df_filt.index): 
                    print(pval_method, "\n", df_null_genotypes, "\n", df_null_phenotypes, "\n", node_info_df_filt, "\n\n\n\n\n")
                    raise ValueError("the columns should be the nodes")

                # Define a df with the resample's values
                df_resamples = pd.DataFrame(index=df_null_genotypes.index)
                df_resamples["nodes_withGeno"] = ((df_null_genotypes==1.0)).apply(sum, axis=1)
                df_resamples["nodes_withPheno"] = ((df_null_phenotypes==1.0)).apply(sum, axis=1)
                df_resamples["nodes_GenoAndPheno"] = ((df_null_genotypes==1.0) & (df_null_phenotypes==1.0)).apply(sum, axis=1)
                df_resamples["nodes_GenoAndNoPheno"] = ((df_null_genotypes==1.0) & (df_null_phenotypes==0.0)).apply(sum, axis=1)
                df_resamples["nodes_noGenoAndPheno"] = ((df_null_genotypes==0.0) & (df_null_phenotypes==1.0)).apply(sum, axis=1)
                df_resamples["nodes_noGenoAndNoPheno"] = ((df_null_genotypes==0.0) & (df_null_phenotypes==0.0)).apply(sum, axis=1)

                # add stats
                df_resamples["chi_square"] = df_resamples.apply(lambda r: calculate_chi_square(r.nodes_GenoAndPheno, r.nodes_GenoAndNoPheno, r.nodes_noGenoAndPheno, r.nodes_noGenoAndNoPheno), axis=1)
                df_resamples["OR"] = df_resamples.apply(lambda r: calculate_OR(r.nodes_GenoAndPheno, r.nodes_GenoAndNoPheno, r.nodes_noGenoAndPheno, r.nodes_noGenoAndNoPheno), axis=1)
                df_resamples["epsilon"] = df_resamples.apply(lambda r: calculate_epsilon(r.nodes_GenoAndPheno, r.nodes_withGeno, r.nodes_withPheno), axis=1)

                # debug
                for f in ["OR", "chi_square", "epsilon"]:
                    if any(pd.isna(df_resamples[f] )): raise ValueError("there can't be nans in %s"%f)
                    if any(df_resamples[f]<0): raise ValueError("there can't be negative values in %s"%f)

                if pd.isna(observed_chi_square): raise ValueError("the chi square can't be nan")
                if pd.isna(observed_OR): raise ValueError("the OR square can't be nan")
                if pd.isna(epsilon): raise ValueError("the OR epsilon can't be nan")

                # calculate p vals
                pval_chi_square = sum(df_resamples.chi_square>=observed_chi_square) / nresamples
                pval_OR = sum(df_resamples.OR>=observed_OR) / nresamples
                pval_GenoAndPheno = sum(df_resamples.nodes_GenoAndPheno>=nodes_GenoAndPheno) / nresamples
                pval_epsilon = sum(df_resamples.epsilon>=epsilon) / nresamples



                # save the df_resamples for some specific cases
                if pval_method=="phenotypes" and nresamples==1000: 
                    #df_resamples["group_name"] = name # this is unneccessary because this df is only for maxT calculation
                    df_resamples["resampleI"] = df_resamples.index
                    phenotypes_1000resamples_stats_file = "%s_%s.tab"%(phenotypes_1000resamples_stats_prefix, method_name)
                    phenotypes_1000resamples_stats_file_tmp = "%s.tmp"%phenotypes_1000resamples_stats_file

                    df_resamples[["chi_square", "OR", "epsilon", "resampleI"]].to_csv(phenotypes_1000resamples_stats_file_tmp, sep="\t", header=False, index=False)
                    os.rename(phenotypes_1000resamples_stats_file_tmp, phenotypes_1000resamples_stats_file)

                # break if the 2 sided pvalue is above some threshold
                min_pval = (1/nresamples)*100
                if all([p>min_pval for p in [pval_chi_square, pval_OR, pval_GenoAndPheno, pval_epsilon]]): break

        else:
            #print("WARNING: Group %s has nnodes_geno_and_pheno_known=%s and fraction_nodes_GenoAndPheno=%s"%(name, nnodes_geno_and_pheno_known, fraction_nodes_GenoAndPheno))
            pval_chi_square = -1
            pval_OR = -1 
            pval_GenoAndPheno = -1
            pval_epsilon = -1

        # check
        if any([pd.isna(p) for p in [pval_chi_square, pval_OR, pval_GenoAndPheno, pval_epsilon]]): raise ValueError("the pval can't be nan")

        # keep
        dict_stats["pval_chi_square_%s"%pval_method] = pval_chi_square
        dict_stats["pval_OR_%s"%pval_method] = pval_OR
        dict_stats["pval_GenoAndPheno_%s"%pval_method] = pval_GenoAndPheno
        dict_stats["pval_epsilon_%s"%pval_method] = pval_epsilon

    return dict_stats


def run_GWAS_as_function(outdir, ASR_mutations_dir, phenotype_states_df_file, tree, name, sorted_mutations, ASR_methods_mutations, classification_method_mutations, ASR_methods_phenotypes, classification_method_phenotypes, max_pval_plot, pct_progress, nresamples_to_filePhenotypesResampling, outdir_ASR_mutations, min_support, node_to_support, interesting_gwas_methods):

    """This is inspired by the scratch/run_GWAS.py script, to be ran as a function"""

    # init the start time
    start_time = time.time()

    # log
    #if pct_progress.endswith("0"): print("Already ran %s%s of GWAS jobs"%(pct_progress, "%"))

    # define final files
    gwas_stats_file = "%s/%s_GWAS_stats.tab"%(outdir, name) # a table with the stats of this group
    ASR_file = "%s/%s_ASR.tab"%(outdir, name) # a table with info about each node
    gwas_stats_file_tmp = "%s.tmp"%gwas_stats_file
    ASR_file_tmp = "%s.tmp"%ASR_file

    # define the 1000 resamples file. These will be dataframes with some stats about each data point
    phenotypes_1000resamples_stats_prefix = "%s/%s_phenotypes_1000resamples_stats"%(outdir, name)

    # remove the tmp files
    for f in [gwas_stats_file_tmp, ASR_file_tmp]: remove_file(f)

    # if the final files are done, skip
    if all([not file_is_empty(f) for f in [gwas_stats_file, ASR_file]]): raise ValueError("The final files were already generated. This function should not be ran if the final files are already generated")

    # check that there are some muts
    if len(sorted_mutations)==0: raise ValueError("There have to be some mutations in group '%s'"%name)

    # get the genotype states for each mutation and node
    ASR_mutations_list = list(map(lambda m: "%s/%s_pastml_reconstruction.out"%(ASR_mutations_dir, m), sorted_mutations))
    genotype_states_df = get_ASR_states_integrate_methods(sorted_mutations, ASR_mutations_list, ASR_methods_mutations, classification_method_mutations, min_support, node_to_support)

    # load the phenotype states, which was generated as get_ASR_states_integrate_methods
    phenotype_states_df = load_object(phenotype_states_df_file)

    # check that the indices (nodes) are the same
    if list(genotype_states_df.index)!=list(phenotype_states_df.index): raise ValueError("The nodes in genotypes and phenotypes tables are not the same ones")

    # add the type of node, and save as ASR_file_tmp
    node_info_df, mut_to_NaNtransition_nodes = get_node_info_df(genotype_states_df, phenotype_states_df, tree, sorted_mutations)
    node_info_df["group_name"] = name
    node_info_df[fields_ASR_df].to_csv(ASR_file_tmp, sep="\t", header=False, index=False)

    # filter out long branches (they are longer than 25% the total tree length)
    total_tree_len = sum(node_info_df.edge_len)
    idx_short_branches = node_info_df.edge_len<=(total_tree_len*0.25)
    #print("There are %i/%i short enough branches"%(sum(idx_short_branches), len(idx_short_branches)))
    node_info_df = node_info_df[idx_short_branches]

    # remove the root
    node_info_df["node"] = node_info_df.index
    node_info_df = node_info_df[node_info_df.node!="root"]
    if len(node_info_df)==0: raise ValueError("There can't 0 nodes") 

    # define the nodes with short branches
    nodes_short_branches = set(idx_short_branches[idx_short_branches==False].index)

    # get the stats for each method
    dict_stats = {}
    for method_name, phenotype_f, genotype_f in [("phyC", "phenotype", "genotype_appearance"), ("synchronous", "phenotype_transition", "genotype_transition")]:
        if method_name not in interesting_gwas_methods: continue # skip some methods

        # get the dict stats for this method
        dict_stats[method_name] = get_dict_stats_one_method_gwas(node_info_df, genotype_f, phenotype_f, method_name, ASR_methods_phenotypes, ASR_methods_mutations, classification_method_phenotypes, classification_method_mutations, name, sorted_mutations, outdir_ASR_mutations, nresamples_to_filePhenotypesResampling, nodes_short_branches, tree, genotype_states_df, phenotypes_1000resamples_stats_prefix)

    # get the stats as df and save
    df_stats = pd.DataFrame(dict_stats).transpose()[fields_GWAS_stats_df]
    df_stats[fields_GWAS_stats_df].to_csv(gwas_stats_file_tmp, sep="\t", header=False, index=False)

    # get the total elapsed time
    elapsed_time = float(time.time() - start_time)

    # clean
    os.rename(ASR_file_tmp, ASR_file)
    os.rename(gwas_stats_file_tmp, gwas_stats_file)

    #print("GWAS for group %s worked well in %.4f s"%(name, elapsed_time))

    # clean
    del node_info_df
    del genotype_states_df
    del phenotype_states_df
    del ASR_mutations_list
    del df_stats

def run_pastml_py_as_a_function(treefile, data_table, workdir, outfile, html_file, pastml_prediction_method, pct_progress):

    """Runs pastml as a function"""

    # record time
    start_time = time.time()

    ########## RUN PASTML ##########

    # print log if specified
    if pct_progress is not None and pct_progress.endswith("0"): print("%s%s pastml jobs completed"%(pct_progress, "%"))
    # nd pct_progress.endswith("0")

    # clean
    outfile_tmp = "%s.tmp"%outfile
    for f in [workdir, outfile, html_file, outfile_tmp]: delete_file_or_folder(f)

    # make the workdir
    make_folder(workdir)

    # re-define the tree to include node names, but this does not include support
    internalNodes_treefile = "%s/tree_internalNodes.nw"%workdir
    tree = get_tree_with_internalNodeNames(Tree(treefile))
    tree.write(outfile=internalNodes_treefile, format=1)

    # check that the tree has no politomies
    check_that_tree_has_no_politomies(tree)

    # get time
    #print("until tree writing took %.4fs"%(time.time()-start_time))


    # run pastml with the re-defined tree
    run_cmd_simple("%s --tree %s --data %s --id_index 0 --out_data %s --work_dir %s --html %s --prediction_method %s --threads 1"%(PASTML, internalNodes_treefile, data_table, outfile_tmp, workdir, html_file, pastml_prediction_method)) 

    # get time
    #print("until pastml run took %.4fs"%(time.time()-start_time))

    ################################

    ########## ONE LINE FOR EACH NODE ##########

    # pastml writes the nodes with uncertainty in different rows, which complicates the analysis. The following step reconfigures the df to include one row per node, and uncertain classifications as NaN

    # get all nodes sorted
    sorted_nodes = ["root"] + list(map(lambda n: n.name, tree.traverse()))
    sorted_nodes = [x for x in sorted_nodes if x!=(tree.name)] # this is to remove the root

    # load as df and check
    df_pastml = get_tab_as_df_or_empty_df(outfile_tmp)
    if set(df_pastml["node"])!=set(sorted_nodes): raise ValueError("The nodes are not the same")

    # redefine so that each node is one row
    def get_consensus_or_nan(series):

        # get the non nans
        unique_vals = set(series[~pd.isna(series)].values)

        # depending on the number of unique vals, return one or the other
        if len(unique_vals)==1: return next(iter(unique_vals))
        elif len(unique_vals)==2: return np.nan
        else: raise ValueError("invalid series: %s"%series)

    def get_r_from_df_node(df_n):
        
        # if there is only one row it means that it is clear
        if len(df_n)==1: return df_n.iloc[0]

        # if there are 2 rows it means that some datasets are uncertain
        elif len(df_n)==2: 

            # define a series to return with the values transformed according to get_consensus_or_nan
            series_return = df_n[[c for c in df_n.columns if c!="node"]].apply(get_consensus_or_nan, axis=0)

            # add the node
            series_return["node"] = df_n.name

            return series_return[df_n.columns]

        else: raise ValueError("invalid df_n: %s"%df_n)

    df_pastml = df_pastml.groupby("node").apply(get_r_from_df_node).set_index("node", drop=True).applymap(float).loc[sorted_nodes]
    df_pastml["node"] = df_pastml.index

    # check 
    if len(df_pastml)!=len(sorted_nodes): raise ValueError("not all nodes are unique")

    # rewrite the df so that the output includes <mutation>_<method> even if there is only one method 
    allowed_predictionMethods = {'MPPA', 'MAP', 'JOINT', 'DOWNPASS', 'ACCTRAN', 'DELTRAN', 'ALL', 'ML', 'MP'}
    if pastml_prediction_method in {'MPPA', 'MAP', 'JOINT', 'DOWNPASS', 'ACCTRAN', 'DELTRAN'}:
        for c in df_pastml.columns:
            if c!="node": 

                df_pastml["%s_%s"%(c, pastml_prediction_method)] = df_pastml[c] 
                df_pastml.pop(c)

    elif pastml_prediction_method not in {'ALL', 'ML', 'MP'}: raise ValueError("pastml_prediction_method %s should be in %s"%(pastml_prediction_method, allowed_predictionMethods))

    # save
    save_df_as_tab(df_pastml, outfile_tmp)

    ############################################

    # clean
    delete_folder(workdir)
    remove_file(data_table)

    # log
    os.rename(outfile_tmp, outfile)

    # print time
    #print("complete pastml run took %.4fs"%(time.time()-start_time))

def check_no_nans_series(x):

    """Raise value error if nans"""

    if any(pd.isna(x)): raise ValueError("There can't be nans in series %s"%x)

def set_nan_to_1(x):

    """Converts nans to 1's"""

    if pd.isna(x): return 1.0
    else: return x


def get_df_cmds_gwas_with_rep_group_name(df_cmds):

    """Adds the representative group name"""

    # checks
    if len(df_cmds)!=len(set(df_cmds.groupName)): raise ValueError("the df_cmds does not have unique groups")

    # add the mutations as a list and log the NR ones
    df_cmds["mutations_list_string"] = df_cmds.mutations_list.str.join("_")
    print("There are %i/%i NR groups of mutations"%(len(set(df_cmds.mutations_list_string)), len(df_cmds)))

    # map each rep group to the groups that it includes
    print("getting mutString_to_groups")
    mutString_to_groups = df_cmds.reset_index(drop=True).groupby("mutations_list_string").apply(lambda df_m: sorted(set(df_m.groupName)))

    # get a df that has the group and the rep group
    print("getting df_rep_groups")

    def get_group_mapping_df_from_list_groups(groups):

        df = pd.DataFrame({"group":groups})
        df["rep_group"] = groups[0]

        return df
    df_rep_groups = pd.concat(mutString_to_groups.apply(get_group_mapping_df_from_list_groups).values).reset_index(drop=True)

    # add to df_cmds
    print("adding rep_groupName to df")
    g_to_repG = dict(df_rep_groups.set_index("group").rep_group)
    df_cmds["rep_groupName"] = df_cmds.groupName.apply(lambda g: g_to_repG[g])

    # checks
    print("checking...")
    if len(set(df_cmds.rep_groupName))!=len(set(df_cmds.mutations_list_string)): raise ValueError("something went wrong with the rep_group definition")

    return df_cmds, df_rep_groups

def chunks(l, n):
    
    """Yield successive n-sized chunks from a list l"""
    
    for i in range(0, len(l), n):
        yield l[i:i + n]


def soft_link_files_fast(origin, target):

    """This function takes an origin file and makes it accessible through a link (target)"""

    if file_is_empty(target):

        # check that the origin exists
        if file_is_empty(origin): raise ValueError("The origin %s should exist"%origin)

        # remove previous link
        try: run_cmd("rm %s > /dev/null 2>&1"%target)
        except: pass

        soft_linking_std = "%s.softlinking.std"%(target)
        run_cmd_simple("ln -s %s %s > %s 2>&1"%(origin, target, soft_linking_std))
        remove_file(soft_linking_std)

    # check that it worked
    if file_is_empty(target): raise ValueError("The target %s should exist"%target)


def copy_gwas_results_redundant_files(outdir, rep_group, group, pct_progress, max_pval_keepASR):

    """Copies the results of rep_group to group in outdir, changing the names."""

    # print progress
    if pct_progress.endswith("0"): print("Already copied %s%s of redundant GWAS jobs"%(pct_progress, "%"))

    # define the files and check that they exist or not
    origin_gwas_stats_file = "%s/%s_GWAS_stats.tab"%(outdir, rep_group) 
    origin_ASR_file = "%s/%s_ASR.tab"%(outdir, rep_group)
    dest_gwas_stats_file = "%s/%s_GWAS_stats.tab"%(outdir, group) 
    dest_ASR_file = "%s/%s_ASR.tab"%(outdir, group)

    # check that the final file does not exist
    if not file_is_empty(dest_gwas_stats_file): raise ValueError("%s should not exist"%dest_gwas_stats_file)

    # load the original ASR df and rename the representative group by the actual one
    asr_df = pd.read_csv(origin_ASR_file, sep="\t", header=None, names=fields_ASR_df)
    if set(asr_df.group_name)!={rep_group}: raise ValueError("the group name should be as the rep_group in ASR df")
    asr_df["group_name"] = group

    # load the original gwas stats df and rename the representative group by the actual one
    gwas_stats_df = pd.read_csv(origin_gwas_stats_file, sep="\t", header=None, names=fields_GWAS_stats_df)
    
    if set(gwas_stats_df.group_name)!={rep_group}: 
        print(gwas_stats_df.group_name, set(gwas_stats_df.group_name), rep_group)
        print(origin_gwas_stats_file)
        raise ValueError("the group name should be as the rep_group")

    gwas_stats_df["group_name"] = group

    # define a function to save a df as tab with no header
    def save_df_as_tab_no_header(df, file):
        file_tmp = "%s.tmp"%file
        df.to_csv(file_tmp, sep="\t", header=False, index=False)
        os.rename(file_tmp, file)

    # save dfs
    save_df_as_tab_no_header(asr_df[fields_ASR_df], dest_ASR_file)
    save_df_as_tab_no_header(gwas_stats_df[fields_GWAS_stats_df], dest_gwas_stats_file)

    # clean
    del asr_df
    del gwas_stats_df


def get_resampled_phenotypes_df_pastml_run_one_resample(phenotypes_file, I, tmpdir, pastml_prediction_method, njobs, treefile):

    """Runs pastml on one resample of the phenotypes"""

    # define the final file
    outfile_pastml = "%s/outfile_pastml_%i.out"%(tmpdir, I)
    if file_is_empty(outfile_pastml):

        # load the df with the phenotypes and check
        phenotypes_df = get_tab_as_df_or_empty_df(phenotypes_file)[["sampleID", "phenotype"]]
        phenotypes_df["sampleID"] = phenotypes_df.sampleID.apply(str)
        strange_phenotypes = set(phenotypes_df.phenotype).difference({0, 1})
        if len(strange_phenotypes)>0: raise ValueError("There are strange phenotypes: %s"%strange_phenotypes)

        # reshuffle the phenotypes and save
        reshuffled_phenotypes_file = "%s/reshuffled_phenotypes_%i.tab"%(tmpdir, I)
        if file_is_empty(reshuffled_phenotypes_file):
            phenotypes_df["phenotype"] = np.random.choice(phenotypes_df.phenotype.values, size=len(phenotypes_df), replace=False, p=None)
            save_df_as_tab(phenotypes_df[["sampleID", "phenotype"]], reshuffled_phenotypes_file)

        # run pastml
        html_file = "%s/outfile_pastml_%i.html"%(tmpdir, I)
        workdir = "%s/workdir_pastml_%i"%(tmpdir, I)
        pct_progress = "%.3f"%((I/njobs)*100)

        if file_is_empty(outfile_pastml): run_pastml_py_as_a_function(treefile, reshuffled_phenotypes_file, workdir, outfile_pastml, html_file, pastml_prediction_method, pct_progress)

        # clean
        remove_file(reshuffled_phenotypes_file)

    # return the loaded df
    df_pastml = get_tab_as_df_or_empty_df(outfile_pastml)
    df_pastml["resample_I"] = I
    return df_pastml

def get_repeated_mutation_gwas_or_asr_df(df, repMut_to_n_mutations, repMut_to_mutations, numericVar_to_var, fields_to_return):

    """Gets a df of gwas, where each row is one group, and returns the duplicated df with all mutations"""

    # checks
    if len(df)!=len(set(df.group_name)): raise ValueError("group name shoul be unique")

    # set the index as the group namne
    df = df.set_index("group_name", drop=False).sort_index()

    # duplicate all the rows, so that each mutation gets n rows according to repMut_to_n_mutations
    df = df.loc[df.index.repeat(repMut_to_n_mutations)].reset_index(drop=True).sort_values(by="group_name")

    # change the group_name, so that it has each mutation, not only the rep_ones
    df["rep_group_name"] = df.group_name
    df["group_name"] = make_flat_listOflists(repMut_to_mutations.values)

    # checks
    if len(df)!=sum(repMut_to_n_mutations.values): raise ValueError("the gwas df should have a row for each mutation")
    if len(df)!=len(set(df.group_name)): raise ValueError("non-unique group names")

    new_repMut_to_mutations = df.set_index("rep_group_name").group_name.groupby("rep_group_name").apply(np.array)
    if list(new_repMut_to_mutations.index)!=list(repMut_to_mutations.index): raise ValueError("The index should be the same")
    if list(new_repMut_to_mutations.apply(tuple))!=list(repMut_to_mutations.apply(tuple)): raise ValueError("the mapping of new mutations did not work properly")

    # redefine the 'group_name' to include the actual variants, keeping the numeric mutation
    df["numeric_mutation"] = df.group_name
    df["group_name"] = df.group_name.map(numericVar_to_var); check_no_nans_series(df.group_name)

    return df[fields_to_return]

def get_asr_df_with_all_mutations(asr_df, repMut_to_n_mutations, repMut_to_mutations, numericVar_to_var, threads, run_without_parallelization):

    """Does something similar to the gwas df but for the asr_df"""

    print("duplicating asr df across all mutations")
    start_time = time.time()

    # checks
    if len(asr_df)!=(len(set(asr_df.node_name)) * len(set(asr_df.group_name))): raise ValueError("There should be one row for each combination of node and group")

    # init all nodes
    all_sorted_nodes = sorted(set(asr_df.node_name))

    # run in parallel the get_repeated_mutation_gwas_or_asr_df
    inputs_fn = [(asr_df[asr_df.node_name==node], repMut_to_n_mutations, repMut_to_mutations, numericVar_to_var, fields_ASR_df) for node in all_sorted_nodes]

    if run_without_parallelization is False:
        
        # in parallel
        with multiproc.Pool(threads) as pool:
            asr_df_final = pd.concat(pool.starmap(get_repeated_mutation_gwas_or_asr_df, inputs_fn))[fields_ASR_df].reset_index(drop=True)
            pool.close()
            pool.terminate()

    else: asr_df_final = pd.concat(map(lambda x: get_repeated_mutation_gwas_or_asr_df(x[0], x[1], x[2], x[3], x[4]), inputs_fn))[fields_ASR_df].reset_index(drop=True)


    # check
    if len(asr_df_final)!=(len(set(asr_df_final.group_name))*len(all_sorted_nodes)): raise ValueError("There should be 2 entries for each group name")
    print("get_asr_df_with_all_mutations took %s seconds"%(time.time()-start_time))

    return asr_df_final


def get_gwas_df_with_all_mutations(gwas_df, repMut_to_n_mutations, repMut_to_mutations, numericVar_to_var, interesting_gwas_methods):

    """Gets the gwas_df or asr_df and returns it with all mutations, not only the representative ones"""

    print("duplicating gwas df across all mutations")

    # define the interested fields
    interesting_fields = list(gwas_df.keys())

    # init the merged df, with data for each method
    merged_gwas_df = pd.DataFrame()

    for method in sorted(set(gwas_df.method)):

        # get the gwas df for this method
        gwas_df_method = gwas_df[gwas_df.method==method]

        # get the duplicated df
        gwas_df_method = get_repeated_mutation_gwas_or_asr_df(gwas_df_method, repMut_to_n_mutations, repMut_to_mutations, numericVar_to_var, interesting_fields)

        # keep
        merged_gwas_df = merged_gwas_df.append(gwas_df_method[interesting_fields]).reset_index(drop=True)

    # check
    if len(merged_gwas_df)!=(len(set(merged_gwas_df.group_name))*len(interesting_gwas_methods)): raise ValueError("There should be 1 or 2 entries for each group name")
    #print("There are %i mutations and the df has %i rows"%(len(set(gwas_df.group_name)), len(gwas_df)))

    return merged_gwas_df


def get_gwas_df_with_fdr_corrected_pvals(gwas_df):

    """Gets the gwas df with the fdr corrected pvals"""
    print("adding FDR-corrected p vals")

    merged_gwas_df = pd.DataFrame()
    for method in sorted(set(gwas_df.method)):

        # get the gwas df for this method
        gwas_df_method = gwas_df[gwas_df.method==method]
        if len(gwas_df_method)!=len(set(gwas_df_method.group_name)): raise ValueError("the gname should be unique")

        # go through each pval field
        for pval_f in pval_fields:

            # get the dfs where it makes sense to calculate pvals. Note that uncalculated p vals are set to -1
            check_no_nans_series(gwas_df_method[pval_f])
            gwas_df_method_validPvals = gwas_df_method[gwas_df_method[pval_f]>=0]

            # map each group name to a corrected pval
            if len(gwas_df_method_validPvals)>0:

                gwas_df_method_validPvals["corrected_pval"] = multitest.fdrcorrection(gwas_df_method_validPvals[pval_f].values)[1]
                gname_to_corrected_pval = dict(gwas_df_method_validPvals.set_index("group_name").corrected_pval)

            else: gname_to_corrected_pval = {}

            # add the others as one
            gname_to_corrected_pval = {**gname_to_corrected_pval, **{gname : 1.0 for gname in set(gwas_df_method[gwas_df_method[pval_f]<0].group_name)}}

            # add to df
            corrected_pval_f = "%s_fdr"%pval_f
            gwas_df_method[corrected_pval_f] = gwas_df_method.group_name.map(gname_to_corrected_pval); check_no_nans_series(gwas_df_method[corrected_pval_f])

        # keep
        merged_gwas_df = merged_gwas_df.append(gwas_df_method).reset_index(drop=True)

    gwas_df = merged_gwas_df

    return gwas_df

def rgb_to_hex(rgb):

    # Helper function to convert colour as RGB tuple to hex string
    rgb = tuple([int(255*val) for val in rgb])
    hex_val = '#%02x%02x%02x'%(rgb[0], rgb[1], rgb[2])

    if len(hex_val)!=7: raise ValueError("%s is not valid"%hex_val)

    return hex_val

def get_value_to_color(values, palette="mako", n=100, type_color="rgb", center=None):

    """TAkes an array and returns the color that each array has. Checj http://seaborn.pydata.org/tutorial/color_palettes.html"""

    # get the colors
    colors = sns.color_palette(palette, n)

    # change the colors
    if type_color=="rgb": colors = colors
    elif type_color=="hex": colors = [rgb_to_hex(c) for c in colors]
    else: raise ValueError("%s is not valid"%palette)

    # if they are strings
    if type(list(values)[0])==str:

        palette_dict = dict(zip(values, colors))
        value_to_color = palette_dict

    # if they are numbers
    else:

        # map eaqually distant numbers to colors
        if center==None:
            min_palette = min(values)
            max_palette = max(values)
        else: 
            max_deviation = max([abs(fn(values)-center) for fn in [min, max]])
            min_palette = center - max_deviation
            max_palette = center + max_deviation

        all_values_palette = list(np.linspace(min_palette, max_palette, n))
        palette_dict = dict(zip(all_values_palette, colors))

        # get value to color
        value_to_color = {v : palette_dict[find_nearest(all_values_palette, v)] for v in values}

    return value_to_color, palette_dict

def generate_pval_distributions(gwas_df, plots_dir, pval_fields_plot):

    """Plots the distributions of p values"""
    print("plotting p value distributions")

    # define the methods
    gwas_methods = ["synchronous", "phyC"]

    # debug
    if len(gwas_df)==0: 
        print("WARNING: There is no single group with nodes_GenoAndPheno>=2")
        return

    # make folder
    make_folder(plots_dir)

    # get colors
    pval_f_to_color = get_value_to_color(pval_fields_plot, n=len(pval_fields_plot), type_color="hex", palette="tab10")[0]

    # get the fig parms
    nrows = len(gwas_methods)
    ncols = len(pval_fields_plot)
    fig = plt.figure(figsize=(ncols*4.5, nrows*4.5)); I=1

    for Ir, gwas_method in enumerate(gwas_methods):
        for Ic, pval_f in enumerate(pval_fields_plot):
            df_plot = gwas_df[(gwas_df.method==gwas_method)]

            # init subplot
            ax = plt.subplot(nrows, ncols, I); I+=1

            # define the hist
            ax = sns.distplot(df_plot[pval_f], hist_kws={"color":pval_f_to_color[pval_f], "linewidth":2, "alpha":1.0}, hist=True, kde=False, rug=False)

            # define the lims
            #ax.set_xlim([0, 1])

            # set the axes
            if Ir==1: ax.set_xlabel("p-val")
            else: 
                ax.set_xticklabels([])
                ax.set_xlabel("")

            if Ic==0: ax.set_ylabel("%s\n# groups"%(gwas_method))
            else: ax.set_ylabel("")

            if len(df_plot)==0:
                ax.set_ylim([0,0]) 
                ax.set_yticks([])
                ax.set_yticklabels([])

            if Ir==0: ax.set_title(pval_f.split("pval_")[1] + "\n%i groups"%len(df_plot))
            else: ax.set_title( "%i groups"%len(df_plot))

    # save and close
    filename_plot = "%s/raw_pval_hist.pdf"%(plots_dir)
    print("saving", filename_plot)
    fig.savefig(filename_plot, bbox_inches='tight')
    plt.close(fig)

def generate_qq_plot(gwas_df, plots_dir, pval_fields_plot):

    """Generates QQ plot for one GWAS result into plotsdir"""
    print("plotting QQ plot")

    # define the methods
    gwas_methods = ["synchronous", "phyC"]

    # debug
    if len(gwas_df)==0: 
        print("WARNING: There is no single group with nodes_GenoAndPheno>=2")
        return

    # make plots dir
    make_folder(plots_dir)

    # define the pseudocount of the pval
    pseudocount_pval = 1/10000 # note that this is 10x more than the resamples

    # define the color
    pval_f_to_color = get_value_to_color(pval_fields_plot, n=len(pval_fields_plot), type_color="hex", palette="tab10")[0]

    # define the lims
    lims = [-0.1, -np.log10(pseudocount_pval) + 0.1]

    # get the fig parms
    nrows = 1
    ncols = len(gwas_methods)
    fig = plt.figure(figsize=(ncols*2.2, nrows*2.2)); I=1

    for Ic, gwas_method in enumerate(gwas_methods):

        # get the df and check
        df_plot = gwas_df[(gwas_df.method==gwas_method)]
        if len(df_plot)!=len(set(df_plot.group_name)): raise ValueError("there should be unique group names")

        # init subplot
        ax = plt.subplot(nrows, ncols, I); I+=1

        if len(df_plot)>0:

            # create a long df with the expected and observed pvals
            df_long = pd.DataFrame()
            for f in pval_fields_plot:
                df = pd.DataFrame({"observed p-value" : df_plot[f].values + pseudocount_pval}).sort_values(by="observed p-value")
                df["type p-value"] = f
                df["expected p-value"] = np.linspace(1/len(df), 1, len(df))
                df["minus_log10pval_observed"] = -np.log10(df["observed p-value"])
                df["minus_log10pval_expected"] = -np.log10(df["expected p-value"])
                df_long = df_long.append(df)

            # plot
            #ax = sns.scatterplot(x="minus_log10pval_expected", y="minus_log10pval_observed", data=df_long, hue="type p-value", edgecolor="none", alpha=.7, style="type p-value", palette="tab10")
            ax = sns.lineplot(x="minus_log10pval_expected", y="minus_log10pval_observed", data=df_long, hue="type p-value", linewidth=2, palette=pval_f_to_color, markers=True)

            # add line
            plt.plot([0, -np.log10(pseudocount_pval)],  [0, -np.log10(pseudocount_pval)], color="gray", linewidth=.7, linestyle="--")

            alpha_line = -np.log10(0.05 + pseudocount_pval)
            plt.axvline(alpha_line, color="gray", linewidth=.7, linestyle="--")
            plt.axhline(alpha_line, color="gray", linewidth=.7, linestyle="--")

        # set legend
        if Ic==(ncols-1): 

            def get_legend_element(color, label): return Line2D([0], [0], color=color, lw=4, label=label, alpha=1.0) 
            legend_elements = [get_legend_element("white", "type p-value")] + [get_legend_element(color, pval_f) for pval_f, color in pval_f_to_color.items()]
            ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(1.0, 0.3))

        elif len(df_plot)>0: ax.get_legend().remove()

        # axes
        ax.set_xlabel("-log(expected p -val)")

        if Ic!=0:
            ax.set_ylabel("")
            ax.set_yticklabels([])

        else: ax.set_ylabel("-log(observed p -val)")

        # title
        ax.set_title("%s\n(%i/%i groups)"%(gwas_method, len(df_plot), sum(gwas_df.method==gwas_method)))

        # get lims
        ax.set_xlim(lims)
        ax.set_ylim(lims)

    plt.subplots_adjust(wspace=0.02, hspace=0.02)

    # save and close
    filename_plot = "%s/QQ_plot.pdf"%(plots_dir)
    print("saving", filename_plot)
    fig.savefig(filename_plot, bbox_inches='tight')
    plt.close(fig)

def get_bytes_as_gb(size_in_bytes):

    """transforms bytes to Gb (from https://thispointer.com/python-get-file-size-in-kb-mb-or-gb-human-readable-format/)"""

    return size_in_bytes/(1024*1024*1024)


def get_set_groups_done_gwas_outdir(suffix, outdir_GWAS): return set(map(lambda x: x.split(suffix)[0], str(subprocess.check_output("cd %s && ls *%s"%(outdir_GWAS, suffix), shell=True)).lstrip("b'").rstrip("'").rstrip("\\n").split("\\n")))


def get_gwas_df_with_maxT_pvals(gwas_df, outdir_GWAS, all_representative_groups, outdir):

    """Adds maxT pvalues to gwas df"""

    print("adding calculating maxT pvalues")
    start_time = time.time()
    final_gwas_df_file = "%s/gwas_df_w_maxTpvals.py"%outdir
    if file_is_empty(final_gwas_df_file):

        # add for each gwas method
        final_gwas_df = pd.DataFrame()
        for gwas_method in sorted(set(gwas_df.method)):

            # keep df of this method
            gwas_df_m = gwas_df[gwas_df.method==gwas_method]

            # define the statistics
            list_statistics = ["chi_square", "OR", "epsilon"]

            # define the groups w expected resampling
            expected_groups_w_resampling = set(gwas_df_m[(gwas_df_m.group_name.isin(all_representative_groups)) & ((gwas_df_m[[f for f in pval_fields if f.endswith('_phenotypes')]]!=-1).apply(any, axis=1))].group_name)

            # if there are such groups, add the maxTs
            if len(expected_groups_w_resampling)>0:

                # validate that you have the expected groups with resampling
                groups_w_resampling = get_set_groups_done_gwas_outdir("_phenotypes_1000resamples_stats_%s.tab"%(gwas_method), outdir_GWAS)
                if groups_w_resampling!=expected_groups_w_resampling: raise ValueError("not all groupings got the resampling file")

                # create the integrated phenotypes file
                print("getting integrated resamples df...")
                integrated_resamples_file = "%s/integrated_phenotype_resamples_df_%s.tab"%(outdir, gwas_method)
                integrated_resamples_file_tmp = "%s.tmp"%integrated_resamples_file

                if file_is_empty(integrated_resamples_file):
                    run_cmd("cat %s/*_phenotypes_1000resamples_stats_%s.tab > %s"%(outdir_GWAS, gwas_method, integrated_resamples_file_tmp))
                    os.rename(integrated_resamples_file_tmp, integrated_resamples_file)

                df_resamples_all = pd.read_csv(integrated_resamples_file, sep="\t", header=None, names=(list_statistics +["resampleI"]))
                for f in df_resamples_all.keys(): check_no_nans_series(df_resamples_all[f])
                if any(df_resamples_all.epsilon)<0 or any(df_resamples_all.epsilon)>1: raise ValueError("epsilon should be btw 0 and 1")
                nresamples = 1000
                if len(df_resamples_all)!=(nresamples*len(expected_groups_w_resampling)): raise ValueError("there should be one I for each group")

                # create maxT null distribution
                print("adding p values...")
                for stat_f in list_statistics:

                    print("getting null distrubution")
                    maxT_null_distribution = df_resamples_all[["resampleI", stat_f]].groupby("resampleI").apply(lambda df: max(df[stat_f])).values

                    print("adding ps")
                    def get_empiric_pval_maxT(x): return sum(maxT_null_distribution>=x)/nresamples
                    gwas_df_m["pval_%s_maxT"%(stat_f)] = gwas_df_m[stat_f].apply(get_empiric_pval_maxT)
                    print("There are %i groups with p<0.05"%(sum(gwas_df_m["pval_%s_maxT"%(stat_f)]<0.05)))
                    check_no_nans_series(gwas_df_m["pval_%s_maxT"%(stat_f)])

            # if not, set to 1.0
            else:
                print("WARNING: There are no resampled pvalues, setting all to 1.0")
                for stat_f in list_statistics: gwas_df_m["pval_%s_maxT"%(stat_f)] = 1.0

            # keep
            final_gwas_df = final_gwas_df.append(gwas_df_m)

        # saving
        save_object(final_gwas_df, final_gwas_df_file)

    final_gwas_df = load_object(final_gwas_df_file)
    print("maxT pval calculations took %.4fs"%(time.time()-start_time))

    return final_gwas_df