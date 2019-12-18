## perform nature methods revision analyses

### input, output and shell paths are all relative to the project directory ###

# downsample new data for revisions    
rule downsample_revision_data:
  input:
    dge = [expand("data/wtxmmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
             rpc = config["downsample"]["reads_per_cell"]["wtxmmix"]),
           expand("data/wtxlung/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
             rpc = config["downsample"]["reads_per_cell"]["wtxlung"]),
           expand("data/tapmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
             rpc = config["downsample"]["reads_per_cell"]["tapmix"]),
           expand("data/taplung/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
             rpc = config["downsample"]["reads_per_cell"]["taplung"]),
           expand("data/19s005246/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
             rpc = config["downsample"]["reads_per_cell"]["19s005246"]),
           expand("data/taphumanmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
             rpc = config["downsample"]["reads_per_cell"]["taphumanmix"])]
