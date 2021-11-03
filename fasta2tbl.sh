#!/usr/bin/awk -f

# Taken from https://bioinformatics.stackexchange.com/questions/2649/how-to-convert-fasta-file-to-tab-delimited-file

{
if (substr($1,1,1)==">")
if (NR>1)
printf "\n%s\t", substr($0,2,length($0)-1)
else
printf "%s\t", substr($0,2,length($0)-1)
else
printf "%s", $0
}END{printf "\n"}
