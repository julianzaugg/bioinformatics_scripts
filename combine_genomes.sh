#!/bin/bash

combine_genomes () {

cd $1
rm $2/genomes_combined.$3

for file in $1/*.$3; do

filename=$(basename $file)
reduced_name=$(echo $filename | sed "s/.$3//g")
header_top=$(head -n 10 $file | grep ">" | sed "s/>//g")
new_header=">$reduced_name"
echo $new_header >> $2/genomes_combined.$3
grep -hv ">" $file | tr -d "\n"  >> $2/genomes_combined.$3
echo "\n" >> $2/genomes_combined.$3

done

cd $2

}

# $1 Input directory, FULL PATH
# $2 Output directory, FULL PATH
# $3 filename extension for input and output genome files, e.g. fna, fasta, etc.

combine_genomes $1 $2 $3
