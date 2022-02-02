#!/usr/bin/env bash

wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
    rm snpEff_latest_core.zip
        cd ./snpEff/
        mkdir -p IL3000 genome


# annotation file
cd IL3000
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/Trypanosoma_congolense/all_assembly_versions/GCA_000227395.2_ASM22739v2/GCA_000227395.2_ASM22739v2_genomic.gtf.gz
gunzip GCA_000227395.2_ASM22739v2_genomic.gtf.gz
mv GCA_000227395.2_ASM22739v2_genomic.gtf genes.gtf

# genome
 cd ../genome
 wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/Trypanosoma_congolense/all_assembly_versions/GCA_000227395.2_ASM22739v2/GCA_000227395.2_ASM22739v2_genomic.fna.gz
 gunzip GCA_000227395.2_ASM22739v2_genomic.fna.gz
 mv GCA_000227395.2_ASM22739v2_genomic.fna IL3000.fa

# config
cd ../
nano snpEff.config | echo "# Database for Trypanosoma congolense, IL3000" >> snpEff.config |
echo "IL3000.genome : Trypanosoma congolense " >> snpEff.config |
echo "IL3000.reference :https://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/Trypanosoma_congolense/latest_assembly_versions/GCA_000227395.2_ASM22739v2/GCA_000227395.2_ASM22739v2_genomic.fna.gz" >> snpEff.config

# build
java -jar snpEff.jar build -gtf22 -v IL3000
