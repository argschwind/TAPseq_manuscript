#!/usr/bin/env python

# extract DGE data from UMI observation counts obtained from the output of the Drop-seq tool
# GatherMolecularBarcodeDistributionByGene. this program also allows for downsampling of the total
# number of genic reads and filtering for chimeric reads based on an approach from
# https://github.com/asncd/schimera

# import modules and functions
from __future__ import print_function
import sys, argparse, os, gzip
import pandas as pd
import numpy as np

# define functions ---------------------------------------------------------------------------------

# helper function to print to stderr
def eprint(*args, **kwargs):
  print(*args, file = sys.stderr, **kwargs)
  
# filter for cell barcodes on whitelist
def filter_cell_barcodes(umi_obs, whitelist):
  
  # get umi observations on whitelist
  bc_filter = umi_obs['Cell Barcode'].isin(whitelist)
  filt_umi_obs = umi_obs[bc_filter]
  
  # get number of filtered out reads and barcodes
  n_filt_bc = len(umi_obs[~bc_filter]['Cell Barcode'].unique().tolist())
  n_filt_reads = umi_obs[~bc_filter]['Num_Obs'].sum()
  perc_filt_reads = n_filt_reads / umi_obs['Num_Obs'].sum() * 100
  
  # report filtered barcodes and reads to stderr
  eprint(' Removed ' + str(n_filt_reads) + ' reads (' + str(round(perc_filt_reads, 3)) + '%) ' +
         'from ' + str(n_filt_bc) + ' cell barcodes')
  
  return(filt_umi_obs)
  
# sample a specified number of total genic reads from the umi counts
def sample_reads(umi_obs, n, seed = None):
  
  # repeat each bc-umi-gene combination by the number of observations to recreate pool of reads
  reads = umi_obs.reindex(umi_obs.index.repeat(umi_obs.Num_Obs))
  
  # randomly sample n reads
  sample_reads = reads.sample(n = n, replace = False, random_state = seed)
  
  # count the number of times each read (bc-umi-gene combination) is found in the downsampled data
  read_counts = sample_reads.groupby(['Cell Barcode', 'Gene', 'Molecular_Barcode']).size()
  
  # reset index to "ungroup" data frame and rename counts column
  read_counts = read_counts.reset_index().rename(columns = {0:'Num_Obs'})
  
  return(read_counts)

# calculate transcripts per transcripts from umi observations (output from  Drop-seq tool:
# GatherMolecularBarcodeDistributionByGene)
def calculate_tpt(umi_obs):

  # calculate the total number of reads per cell-umi barcode combination
  reads = umi_obs.groupby(['Cell Barcode', 'Molecular_Barcode'])['Num_Obs'].sum()

  # convert to DataFrame with barcodes as columns and rename Num_obs column
  reads = pd.DataFrame(reads).reset_index()
  reads = reads.rename(columns = {'Num_Obs': 'total_reads'})

  # add total reads to umi_counts
  umi_obs = pd.merge(umi_obs, reads, how = 'left', on = ['Cell Barcode', 'Molecular_Barcode'])

  # calculate the fraction of total reads for each gene-cell-umi combination, where total reads is
  # the sum of reads for that cell-umi combination. this metric is also called transcripts per
  # transcript (tpt).
  umi_obs['tpt'] = np.divide(umi_obs['Num_Obs'], umi_obs['total_reads'])
  umi_tpt = umi_obs.drop(columns = 'total_reads')

  return(umi_tpt)

# compute tpt histogram from calculate_tpt() output
def compute_tpt_histogram(umi_tpt):
  
  # calculate number of bc-umi-gene combinations per tpt value
  tpt_hist = pd.value_counts(umi_tpt.tpt).to_frame().reset_index()
  
  # rename columns and sort according to tpt value
  tpt_hist.columns = ['tpt','transcripts']
  tpt_hist = tpt_hist.sort_values('tpt', ascending = False)
  
  return(tpt_hist)

# filter umi observations based on minimum tbt value
def filter_tpt(umi_tpt, tpt_threshold):

  # retain observations with tbt >= tbt_threshold
  filter_logical = umi_tpt['tpt'] >= tpt_threshold
  umi_filt = umi_tpt[filter_logical]

  eprint('  Removed ' + str(np.round(100 * (1 - np.mean(filter_logical)), 3)) + '% of transcripts')

  return(umi_filt)
  
# calculate dge summary based of umi observations. calculates number of genic reads, molecules and
# detected genes per cell
def dge_summary(umi_obs):
  
  # calculate number of genic reads, detected transcripts and genes per cell barcode
  agg_funs = {'Num_Obs': ['sum', 'count'], 'Gene' : 'nunique'}
  stats = umi_obs.groupby('Cell Barcode').agg(agg_funs)
  
  # sort according to number of genic reads
  stats = stats.sort_values([('Num_Obs', 'sum')], ascending = False)
  
  # make cell barcode a regular column and set new column names
  stats = stats.reset_index()
  stats.columns = ['cell_barcode', 'genic_reads', 'transcripts', 'genes']
  
  return(stats)
  
# convert umi observations 
def create_dge_matrix(umi_obs):
  
  # count number of umis per cell barcode and gene
  umi_counts = umi_obs.groupby(['Cell Barcode', 'Gene']).size()
  
  # convert to wide format ('unstack') with missing values set to 0 to create dge expression matrix
  dge_matrix = umi_counts.unstack(level = 0, fill_value = 0)
  
  # rename 'Gene' index to 'GENE' to mirror Drop-seq tools DGE output
  dge_matrix.index.rename('GENE', inplace = True)
  
  return(dge_matrix)

# execute if run as main program -------------------------------------------------------------------

# parse command line arguments
if __name__ == '__main__':
  
  # function to check if an argument contains None or a positive integer
  def none_or_int(arg):
    if arg == 'None':
      return(None)
    elif arg.isdigit():
      return(int(arg))
    else:
      msg = "invalid None or int value: %r" % arg
      raise argparse.ArgumentTypeError(msg)

  parser = argparse.ArgumentParser(description = 'Extract DGE data from UMI observations. This \
    program allows for downsampling of total genic reads and filtering for potential chimeric \
    reads. It returns a text file containing the gene expression matrix.')

  parser.add_argument('-i', '--inputfile', help = 'File containing UMI observations as produced by \
    GatherMolecularBarcodeDistributionByGene from Drop-seq tools.', required = True)
  parser.add_argument('-o', '--outputfile', help = 'Output file for filtered digital gene \
    expression matrix. TPT histogram and DGE summary files will be written into the same directory',
    required = True)
  parser.add_argument('-w', '--whitelist', help = 'File containing possible true cell barcodes \
    without header and one cell barcode per line. If provided only cell barcodes on this \
    whitelist will be considered. Default = None, which turns cell barcode filtering off.',
    default = None)
  parser.add_argument('--tpt_threshold', help = 'Minimum transcript per transcript value for \
    gene-cell-umi combinations to pass chimeric reads filtering. Default = 0, which does not \
    remove any reads', type = float, default = 0.0)
  parser.add_argument('--sample', help = 'Number of genic reads (int) to be drawn from input for \
    downsampling. Default = None, which turns sampling off.', type = none_or_int, default = 0)
  parser.add_argument('--seed', help = 'Seed (int) for the random number generator used when \
    downsampling. Default = None', type = none_or_int, default = None)

  args = parser.parse_args()

  # read umi observations input file
  eprint('Reading input file...')
  umi_obs = pd.read_csv(args.inputfile, sep = '\t')
  
  # filter for cell barcodes on whitelist (if whitelist file is provided)
  if args.whitelist:
    eprint('Filtering for cell barcodes on provided whitelist...')
    if args.whitelist.endswith('.gz'):
      with gzip.open(args.whitelist, 'rt') as f:
        whitelist = f.read().splitlines()
    else:
      with open(args.whitelist, 'r') as f:
        whitelist = f.read().splitlines()
    umi_obs = filter_cell_barcodes(umi_obs, whitelist)
  
  # downsample total genic reads if specified
  if args.sample:
    if args.sample <= umi_obs.Num_Obs.sum():
      eprint('Downsampling to ' + str(args.sample) + ' reads...')
      umi_obs = sample_reads(umi_obs, n = args.sample, seed = args.seed)
    else:
      eprint('Desired sampling size larger than number of input reads! Exiting.')
      quit()

  # calculate tpt
  eprint('Calculating TPT...')
  umi_tpt = calculate_tpt(umi_obs)
  
  # compute tpt histogram
  tpt_hist = compute_tpt_histogram(umi_tpt)
  
  # write histogram to text file
  hist_file = os.path.splitext(args.outputfile)[0] + '_tpt_histogram.txt'
  tpt_hist.to_csv(hist_file, sep = '\t', index = False)

  # filter transcripts based on minimum tbt
  eprint('Filtering for chimeric reads (min TPT = ' + str(args.tpt_threshold) + ')...')
  umi_filt = filter_tpt(umi_tpt, tpt_threshold = args.tpt_threshold)
  
  # calculate dge summary statistics
  dge_stats = dge_summary(umi_filt)
  
  # write dge summary to file
  dge_stats_file = os.path.splitext(args.outputfile)[0] + '_summary.txt'
  dge_stats.to_csv(dge_stats_file, sep = '\t', index = False)
  
  # create expression matrix and save to tab deliminted .txt file
  dge = create_dge_matrix(umi_filt)
  eprint('Writing DGE matrix to file...')
  dge.to_csv(args.outputfile, sep = '\t', index = True)
  
  eprint('Done!')
