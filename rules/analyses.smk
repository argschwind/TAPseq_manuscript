## perform all downstream analyses

### input, output and shell paths are all relative to the project directory ###

# define required python functions -----------------------------------------------------------------

# create all dge output filenames for a sample id based on the number cells
# specified for that sample under 'dge_ncells' in config file
def dge_outfiles(sample):
  ncells = config["dge_ncells"][sample]  # get number of cells for sample
  files = expand("data/{sample}/dge_top_{ncells}_cells.txt",
             sample = sample, ncells = ncells)
  return files
  
# get fastq read 2 (genome read) file for a given sample
def fastq_read2(sample):
  fastq_dir = config["samples"][sample]
  f2 = glob.glob(fastq_dir + "/*" + sample  + "_2_sequence.txt.gz")
  return f2
  
# rules --------------------------------------------------------------------------------------------
