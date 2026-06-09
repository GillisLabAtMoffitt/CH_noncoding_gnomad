#! /bin/bash
#SBATCH --nodes 1
#SBATCH --cpus-per-task 2
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb

 cd $SLURM_SUBMIT_DIR

ml BEDTools
ml BCFtools


bedtools intersect -v -header \
-a /share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad/v4.1/gnomad.genomes.v4.1.sites.${chr}.vcf.bgz \
-b /share/lab_gillis/Christelle/intronic_skewness/other_data/exons.merged.bed \
> /share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad/v4.1/Non_coding/gnomad.genomes.v4.1.sites.${chr}.non_coding.vcf

#compress the vcf and create index
ml tabixpp/1.1.0-GCC-10.2.0

bgzip -f /share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad/v4.1/Non_coding/gnomad.genomes.v4.1.sites.${chr}.non_coding.vcf
tabix -f /share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad/v4.1/Non_coding/gnomad.genomes.v4.1.sites.${chr}.non_coding.vcf.gz

#bcftools view -h /share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad/v4.1/Non_coding/gnomad.genomes.v4.1.sites.${chr}.non_coding.vcf.gz | tail -1 \
#> /share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad/v4.1/Non_coding/comment_removed/gnomad.genomes.v4.1.sites.${chr}.non_coding.clean.vcf


#bcftools view -H /share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad/v4.1/Non_coding/gnomad.genomes.v4.1.sites.${chr}.non_coding.vcf.gz \
#>> /share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad/v4.1/Non_coding/comment_removed/gnomad.genomes.v4.1.sites.${chr}.non_coding.clean.vcf

#bgzip -f /share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad/v4.1/Non_coding/comment_removed/gnomad.genomes.v4.1.sites.${chr}.non_coding.clean.vcf


