# CH_noncoding_gnomad
Identify noncoding CH variants in gnomAD database

Identification of noncoding CH variants in gnomAD
Chromosome-level genome data from gnomAD v4.1 were downloaded from
https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes using the download_genome.sh script.
Comprehensive gene annotations, including exon coordinates, were obtained from GENCODE Release 49 (GRCh38.p14), 
downloaded from https://www.gencodegenes.org/human/

Noncoding regions were extracted by running the Extract_nonCoding.sh script. 
The resulting chromosome-level noncoding datasets were then used to identify clonal hematopoiesis (CH) 
variants in gnomAD noncoding regions using the pipeline described below.

...