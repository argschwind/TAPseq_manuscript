#!/usr/bin/env python

# collapse perturbation status data by targets of gRNA vectors

# import modules and functions
from __future__ import print_function
import sys, argparse
import pandas as pd

# define functions ---------------------------------------------------------------------------------

# helper function to print to stderr
def eprint(*args, **kwargs):
  print(*args, file = sys.stderr, **kwargs)
  
# collapse perturbations by target
def collapse_perts(perts, targets):
  
  # add targets to perturbation status and remove VECTOR from data.frame
  perts = pd.merge(targets, perts, how = 'right', on = 'VECTOR').drop(columns = 'VECTOR')
  
  # count number of perturbations per target
  perts = perts.groupby('TARGET').sum()

  return perts

# execute if run as main program -------------------------------------------------------------------

if __name__ == '__main__':

  # parse command line arguments
  parser = argparse.ArgumentParser(description = 'Collapse perturbations by gRNA targets.')

  parser.add_argument('-i', '--inputfile', help = 'File containing perturbation status data to \
    merge. Typically output of perturbation_status.R', required = True)
  parser.add_argument('-o', '--outputfile', help = 'Output file for collapsed perturbation data',
    required = True)
  parser.add_argument('-t', '--targets', help = 'Two column .csv file listing every gRNA vector \
    (column named VECTOR) and their targets (column named TARGET).', required = True)
  parser.add_argument('-c', '--count', help = 'Should perturbations be counted rather than \
    returned as binary matrix (default)', action='store_true')

  args = parser.parse_args()
  
  # load perturbation and vector target files
  eprint('Loading input files...')
  perts = pd.read_csv(args.inputfile, sep = '\t')
  targets = pd.read_csv(args.targets)
  
  # collapse pertubations by vector targets
  eprint('Collapsing perturbations...')
  perts = collapse_perts(perts, targets)
  
  # convert to binary matrix unless specified otherwise
  if not args.count:
    eprint('Converting collapsed perturbations to binary...')
    perts[perts > 0] = 1
  
  eprint('Writing collapsed perturbations to file...')
  perts.to_csv(args.outputfile, sep = '\t', index = True)
  
  eprint('Done!')
