#!/bin/bash
#SBATCH --mail-user=Christelle.ColinLeitzinger@moffitt.org
#SBATCH --mail-type=FAIL,END
#SBATCH --cpus-per-task 2
#SBATCH --ntasks 1
#SBATCH --time 0-200:00:00
#SBATCH --mem=800gb
#SBATCH --output noise_filter-%j.out
ml R/3.6
cd /share/lab_gillis/Christelle/intronic_skewness/pipeline/noise_filter/


# Define the directory containing your input files
# INPUT_DIR="/share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad/small_file_for_testing"
INPUT_DIR="/share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad/v4.1/Non_coding/comment_removed/AB_extract"
# Ensure the directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Directory $INPUT_DIR not found."
    exit 1
fi

# Loop through each file in the directory
for file in "$INPUT_DIR"/*; do
    # Check if it's a file (not a directory)
    if [ -f "$file" ]; then
    	# Get the base filename without path (e.g., 'data.txt')
        filename=$(basename "$file")
        # Run Rscript in the background using &
        # Use nohup to keep it running if you close the terminal
        nohup Rscript new_noise_filter.R "$file" > "/share/lab_gillis/Christelle/intronic_skewness/pipeline/noise_filter/${filename}.log" 2>&1 &
    fi
done

# Wait for all background processes to complete before finishing
wait
echo "All simultaneous jobs have finished."

