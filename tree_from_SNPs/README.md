# Overview

This is a repository with code to get a tree from variant calling data. This is a reproducible implementation of the method used in the paper "Genome-wide signatures of recent selection and drug resistance across Candida opportunistic pathogens".

# Running

Create the conda environment with all dependencies:

`conda env create --file ./tree_from_SNPs_env.yml --name tree_from_SNPs_env` # you may also use mamba instead of conda

Activate environment with:

`conda activate tree_from_SNPs_env`

Run script with: 

`<path to tree_from_SNPs>/scripts/get_tree.py --paths_table <paths.csv>  --mode <haploid / diploid_homozygous / diploid> --reference_genome <genome.fasta> --min_coverage_pos <min coverage>`

This would be an example of `<paths.csv>` to run for three samples in haploid mode (filtering out positions with heterozygous SNPs):

```
sampleID	sorted_bam	vcf_haploid	vcf_diploid
s1	./aligned_reads_s1.bam.sorted	./small_variants_haploid_s1.vcf	./small_variants_diploid_s1.vcf
s2	./aligned_reads_s2.bam.sorted	./small_variants_haploid_s2.vcf	./small_variants_diploid_s2.vcf
s3	./aligned_reads_s3.bam.sorted	./small_variants_haploid_s3.vcf	./small_variants_diploid_s3.vcf
```
Relevant optional arguments are:

- `--batch_mode`: To run in a batch mode for HPC clusters.
- `--reference_genome <ref genome path>`: Path to the reference genome, only necessary for `--mode diploid`.
- `n_resampled_trees <number>`: Number of resampled trees to generate the consensus tree, only necessary for `--mode diploid`.


To better understand how to use the arguments, type `<path to tree_from_SNPs>/scripts/get_tree.py -h`

# Running in BSC

This script is already installed in the MN, which can be used by BSC users.

Activate conda env with:

`source /gpfs/projects/bsc40/mschikora/anaconda3/etc/profile.d/conda.sh && conda activate tree_from_SNPs_env`

The scripts are in `/gpfs/projects/bsc40/mschikora/scripts/tree_from_SNPs`, you can run them from there.