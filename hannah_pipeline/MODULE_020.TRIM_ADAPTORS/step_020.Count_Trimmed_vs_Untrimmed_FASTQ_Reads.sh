project_name=2025.NYNRIN.DOD.ELCHEVA
batch_name=2025_01

raw_fastq_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/${project_name}/FASTQ_FILES/RAW/BATCH_${batch_name}
trim_fastq_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/${project_name}/FASTQ_FILES/FASTQ_FILES.TRIMMED/BATCH_${batch_name}

# Output directory
output_dir=/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/${project_name}/READ_COUNTS_TRIMMED/${batch_name}/

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Output files
reads_counts="$output_dir/reads_counts.txt"
trimmed_reads_counts="$output_dir/trimmed_reads_counts.txt"
combined_counts="$output_dir/combined_reads_counts.txt"

# Function to count reads
count_reads() {
    local dir=$1
    local output=$2

    cd "$dir" || { echo "Error: Directory $dir not found"; exit 1; }

    # Count reads in each file
    for file in *.fastq.gz; do
        num_reads=$(($(wc -l < "$file") / 4))
        echo -e "$file\t$num_reads"
    done > "$output"
}

# Count reads in non-trimmed files
count_reads "$raw_fastq_dir" "$reads_counts"

# Count reads in trimmed files
count_reads "$trim_fastq_dir" "$trimmed_reads_counts"

# Combine the results into a single file
paste "$reads_counts" "$trimmed_reads_counts" > "$combined_counts"

echo "Read counts before and after trimming have been saved to $combined_counts."

