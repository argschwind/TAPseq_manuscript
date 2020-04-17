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

# BARCODE READS ------------------------------------------------------------------------------------

# function to get i7 barcode from barcode table, add it to the beginning of every barcode read and
# append output to compressed fastq file
process_read1 () {
  local indir=$1
  local outfile=$2
  local i7_file=$3
  for i in ${indir}/*_1.fastq.gz; do
    bn=$(basename $i)
    sample=$(echo $bn | perl -pe 's/(\w+)_1.fastq.gz/$1/')
    i7=$(awk -v sample="$sample" '$1 == sample{print $2}' $i7_file)
    echo -e "\tprocessing ${i} (i7 = ${i7})..."
    zcat $i | sed "2~4s/^/${i7}/" | sed '4~4s/^/IIIIIIII/' | gzip >> $outfile
  done
}

# add all genome reads into one fastq file
echo "Merging barcode reads:"
process_read1 $indir $r1_outfile ${sd}/i7_barcodes.txt

# GENOME READS -------------------------------------------------------------------------------------

# function to simply append all genome reads to one compressed fastq file
process_read2() {
  local indir=$1
  local outfile=$2
  for i in ${indir}/*_2.fastq.gz; do
    echo -e "\tprocessing ${i}..."
    cat $i >> $outfile
  done
}

# add all genome reads into one fastq file
echo -e "\nMerging genome reads:"
process_read2 $indir $r2_outfile

echo -e "\nDone!"
