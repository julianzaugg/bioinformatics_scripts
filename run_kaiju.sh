# Convience wrapper for running kaiju on input fastq.gz files

# TODO, improve input argument checking

function run_kaiju() {
# $1 Input fastq.gz file
# $2 kaiju nodes.dmp file
# $3 kaiju .fmi file
# $4 output directory
# $5 output prefix (default name string preceeding .fastq.gz)
# $6 output suffix (default to _kaiju_out.tsv)
# $7 threads (default 5)
# $8 kaiju names file (optional)

# if [ -z "$1" ]; then
#     echo "No fastq.gz file specified \$1, exiting"
#     exit 1
fi
if [[ ! $1 == *.fastq.gz ]]; then
    echo "Input read file must end with fastq.gz \$1, exiting"
    exit 1
fi
# if [ -z "$2" ]; then
#     echo "No kaiju nodes.dmp file specified \$2, exiting"
#     exit 1
# fi
if [[ ! $2 == *.dmp ]]; then
    echo "No kaiju nodes.dmp file specified \$2, exiting"
    exit 1
fi
# if [ -z "$3" ]; then
#     echo "No kaiju .fmi file specified \$3, exiting"
#     exit 1
# fi
if [[ ! $3 == *.fi ]]; then
    echo "No kaiju .fmi file specified \$3, exiting"
    exit 1
fi
if [ -z "$4" ]; then
    echo "No output directory specified \$4, exiting"
    exit 1
fi
if [ -z "$5" ]; then
    echo "No output prefix specified, defaulting to content preceeding .fastq.gz \$5"
    PREFIX=$4
else
    PREFIX=${i%.fastq.gz}
fi
if [ -z "$6" ]; then
    echo "No output suffix specified, defaulting to kiaju_out \$6"
    SUFFIX="kiaju_out"
else
    SUFFIX=$6
fi
THREADS=${7:-5}

if [[ ! -f $4/${PREFIX}_${SUFFIX}.tsv ]]; then
kaiju \
-i $1 \
-t $2 \
-f $3 \
-o $4/${PREFIX}_${SUFFIX}.tsv \
-z $THREADS \
-a greedy

fi 

# If names file provided, generate output for each taxonomy level
if [ ! -z "$8" ]; then
    
    if [[ ! -f $4/${PREFIX}_${SUFFIX}.krona ]]; then
        kaiju2krona -t $2 -n $8 -i $4/${PREFIX}_${SUFFIX}.tsv -o $4/${PREFIX}_${SUFFIX}.krona
    fi
    if [[ ! -f $4/kaiju2table_done ]]; then
        kaiju2table -t $2 -n $8 -p -r species -o $4/${PREFIX}_${SUFFIX}_species.tsv $4/${PREFIX}_${SUFFIX}.tsv
        kaiju2table -t $2 -n $8 -p -r genus -o $4/${PREFIX}_${SUFFIX}_genus.tsv $4/${PREFIX}_${SUFFIX}.tsv
        kaiju2table -t $2 -n $8 -p -r family -o $4/${PREFIX}_${SUFFIX}_family.tsv $4/${PREFIX}_${SUFFIX}.tsv
        kaiju2table -t $2 -n $8 -p -r order -o $4/${PREFIX}_${SUFFIX}_order.tsv $4/${PREFIX}_${SUFFIX}.tsv
        kaiju2table -t $2 -n $8 -p -r class -o $4/${PREFIX}_${SUFFIX}_class.tsv $4/${PREFIX}_${SUFFIX}.tsv
        kaiju2table -t $2 -n $8 -p -r phylum -o $4/${PREFIX}_${SUFFIX}_phylum.tsv $4/${PREFIX}_${SUFFIX}.tsv
        touch $4/kaiju2table_done
    fi
fi

}

export run_kaiju

#### EXAMPLE ####
# READS_INPUT_DIR="my/input_read/dir"
# KAIJU_OUTPUT_DIR="my/output_kaiju/dir"
# KAIJU_VIRAL_DB="kaiju_viral_db"

THREADS=30

cd $READS_INPUT_DIR

for i in *.fastq.gz; do
    reduced_name=${i%.fastq.gz}
    mkdir -p $KAIJU_OUTPUT_DIR/${reduced_name}
    run_kaiju $i $KAIJU_VIRAL_DB/nodes.dmp $KAIJU_VIRAL_DB/kaiju_db_viruses.fmi $KAIJU_OUTPUT_DIR/${reduced_name} $reduced_name "viruses" $THREADS $KAIJU_VIRAL_DB/names.dmp
done


