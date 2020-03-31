module load miniconda3/1.1
conda activate fastspar_0.0.10
module load parallel/20180222

OTU_TABLE="otu_counts.tsv"
OUTDIR="fastpar_results"

THREADS=1
PARALLEL_THREADS=36
BOOTSTRAP_PERMUTATIONS=1000
FASTSPAR_ITERATIONS=50

run_fastspar(){
# $1 Input OTU table
# $2 Output base directory
# $3 Prefix on files
PREFIX=${3:-fastspar}
mkdir -p $2/bootstrap_counts
mkdir -p $2/bootstrap_correlations
mkdir -p $2/bootstrap_covariances

# First calculate the median correlations and covariances
fastspar --otu_table $1 --correlation $2/${PREFIX}_correlation.tsv --covariance $2/${PREFIX}_covariance.tsv #--threads $THREADS

# Create random data sets from original table
fastspar_bootstrap --otu_table $1 --number $BOOTSTRAP_PERMUTATIONS --prefix $2/bootstrap_counts/$PREFIX #--threads $THREADS

# Run FastSpar (SPARCC) on all bootstraped tables
parallel -j $PARALLEL_THREADS fastspar --yes --otu_table {} --correlation $2/bootstrap_correlations/cor_{/} --covariance $2/bootstrap_covariances/cov_{/} --iterations 5 --threads $THREADS ::: $2/bootstrap_counts/*.tsv

# Calculate p-values
fastspar_pvalues --otu_table $1 --correlation $2/${PREFIX}_correlation.tsv --prefix $2/bootstrap_correlations/cor_$PREFIX --permutations $BOOTSTRAP_PERMUTATIONS --outfile $2/${PREFIX}_pvalues.tsv #--threads $THREADS

}

export -f run_fastspar
export THREADS
export PARALLEL_THREADS
export BOOTSTRAP_PERMUTATIONS
export FASTSPAR_ITERATIONS
run_fastspar $OTU_TABLE $OUTDIR "OTU"
