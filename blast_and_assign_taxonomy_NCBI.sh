module load blast/2.9.0.binary

# Location of rankedlineage.dmp file from NCBI
# Download from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
NCBI_TAXID_LINEAGES="rankedlineage.dmp"

# NCBI database.
DATABASE="/srv/db/ncbi/20190520/nt"

# Queries to BLAST
QUERY="MY_QUERIES.fasta"

# Output location
OUTPUT="out"
mkdir -p $OUTPUT

THREADS=20
EVALUE_THRESH="0.0000000001"

# BLAST query sequences
blastn \
-query $QUERY \
-task blastn \
-out $OUTPUT/blast_out.tsv \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxid sscinames scomnames sskingdoms" \
-db $DATABASE \
-num_threads $THREADS \
-evalue $EVALUE_THRESH \
-max_hsps 1

# Just get the top hit for each query ($1), assumes first hit is always the best
# !a[$1]++ = if the entry does not exist, ! will make it true and the default command (print) will be applied. Else false.
awk '! a[$1]++' $OUTPUT/blast_out.tsv > $OUTPUT/blast_out_top.tsv

# Extract the taxids from blast hits. Sort and get unique set.
awk -F "\t" '{print $16}' $OUTPUT/blast_out_top.tsv | sort -u > $OUTPUT/unique_taxids.tsv

# Loop through taxids and get the corresponding taxonomy from the ranked lineage dump from NCBI
# Write these to a mapping file
rm $OUTPUT/taxids_taxonomy.tsv
while read taxid; do
#taxonomy=$(grep "^$taxid\s" $NCBI_TAXID_LINEAGES | sed "s/[\s\t ]//g" | cut -d "|" -f 2- | awk -F "|" '{for (i=NF;i>0;i--){printf $i";"}}' | sed "s/^;//g;s/;$//g")
taxonomy=$(grep "^$taxid\s" $NCBI_TAXID_LINEAGES | sed "s/[\t]//g" | cut -d "|" -f 2- | awk -F "|" '{for (i=NF;i>0;i--){printf $i";"}}' | sed "s/^;//g;s/;$//g")
echo -e "$taxid\t$taxonomy" >> $OUTPUT/taxids_taxonomy.tsv

done < $OUTPUT/unique_taxids.tsv


#blastdbcmd -db $DATABASE -entry_batch unique_accessions.tsv -outfmt "%a\t%L\t%T" > unique_accessions_taxonomies.tsv

# Now combine the top BLAST results for each query (when available) with the corresponding taxid and various taxonomy information
echo -e "query_id\tsubject_id\tpct_identity\taln_length\tn_of_mismatches\tgap_openings\tq_start\tq_end\ts_start\ts_end\te_value\tbit_score\tquery_length\tsubject_length\tsubject_title\tsubject_taxid\tsubject_scientific_name\tsubject_common_name\tsubject_kingdom\tTaxonomy" > $OUTPUT/blast_final.tsv
while read line; do
taxid=$(echo -e "$line" | awk -F "\t" '{print $16}')
taxonomy=$(grep -P "^$taxid\t" $OUTPUT/taxids_taxonomy.tsv | awk -F "\t" '{print $2}')
echo -e "$line\t$taxonomy" >> $OUTPUT/blast_final.tsv
done < $OUTPUT/blast_out_top.tsv