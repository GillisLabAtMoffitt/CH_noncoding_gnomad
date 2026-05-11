#!/bin/bash
#SBATCH --mail-user=Christelle.ColinLeitzinger@moffitt.org
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks 1
#SBATCH --time 0-5:00:00
#SBATCH --mem=200gb
#SBATCH --output download_genome-%j.out

cd /share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad


download_dir="/share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad/v4.1/"
mkdir -p "${download_dir}"


chr_list=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

for chr in "${chr_list[@]}"; do
    echo "Downloading chromosome ${chr}..."
    wget -P "${download_dir}" \
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"
done