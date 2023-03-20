#!/bin/bash


set -e

ACCESSIONS=$1
OUTDIR=$2

# Construction of the accession list
ACC_LIST=$(cat $ACCESSIONS | cut -d '.' -f 1 | xargs)

# Download
datasets download genome accession $ACC_LIST
mkdir -p $OUTDIR
mv ncbi_dataset.zip $OUTDIR

cd $OUTDIR

# Prepare data
unzip -o ncbi_dataset.zip
mv ncbi_dataset/data/*/*.fna .
rm -r ncbi_dataset.zip ncbi_dataset README.md

# Count kmers
echo "" > db_list.txt
for f in *.fna;
do
	kmc -fm -ci0 -k31 $f ${f}_k31 /tmp;
	echo "$OUTDIR${f}_k31" >> db_list.txt
done

cd -
