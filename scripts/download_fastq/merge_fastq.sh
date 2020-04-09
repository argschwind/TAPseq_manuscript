#!/bin/bash

# some samples were split into multiple sequencing runs. this script  merges all read 1 and 2 fastq
# files into one read 1 and one read 2 fastq file. the i7 barcode of each sample (10x lane) is added
# to the barcode read to avoid cell barcode collisions.

# directory of the script
sd=$(dirname $0)

# input directory passed as argument
indir=$1

# output files passed as arguments
r1_outfile=$2
r2_outfile=$3

# delete potentially existing output files to avoid attaching to them
if [ -f $r1_outfile ]; then
  rm $r1_outfile
fi
  
if [ -f $r2_outfile ]; then
  rm $r2_outfile
fi

# function to simply append all genome reads to one compressed fastq file
process_read() {
  local indir=$1
  local outfile=$2
  local read=$3
  for i in ${indir}/*_${read}.fastq.gz; do
    echo -e "\tprocessing ${i}..."
    cat $i >> $outfile
  done
}

echo "Merging barcode reads:"
process_read $indir $r1_outfile "1"

echo -e "\nMerging genome reads:"
process_read $indir $r2_outfile "2"

echo -e "\nDone!"
