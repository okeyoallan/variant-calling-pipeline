#!/usr/bin/env bash

# This script is a snipet of steps to follow to creating a snpeff database for your organism
# in the cases your organism's  database is not pre-built in the snpeff database
# install snpeff binaries
        wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
# unpack snpeff executables
        unzip snpEff_latest_core.zip
# remove the compressed file
        rm snpEff_latest_core.zip
# move into the sneff directory
        cd ./snpEff/
       
# You can copy the lines of code above, run them on the command line interface,
# then run the command 
# java -jar snpEff.jar databases | grep <your organism's name> for example 
# java -jar snpEff.jar databases | grep Trypanosoma 
# to check if your study has a database pre-built
# In the event the output is null, follow carefully how this script was used to build a database for Trypanosoma congolense
# to build one for your organism 


# make a directory 'IL3000' and  'genomes' in './data' directory
        mkdir -p ./data/IL3000 ./data/genomes
        cd ./data/IL3000
	# download Trypanosoma congolense gff file
        wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/Trypanosoma_congolense/latest_assembly_versions/GCA_000227395.2_ASM22739v2/GCA_000227395.2_ASM22739v2_genomic.gff.gz

# unzip then rename genes.gff
        gunzip GCA_000227395.2_ASM22739v2_genomic.gff.gz
        mv GCA_000227395.2_ASM22739v2_genomic.gff  genes.gff

# download gtf file (annotated genes), unzip then rename to genes.gtf
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/Trypanosoma_congolense/latest_assembly_versions/GCA_000227395.2_ASM22739v2/GCA_000227395.2_ASM22739v2_genomic.gtf.gz
gunzip GCA_000227395.2_ASM22739v2_genomic.gtf.gz
mv GCA_000227395.2_ASM22739v2_genomic.gtf genes.gtf

# download proteins, unzip, then rename to proteins.fa
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/Trypanosoma_congolense/latest_assembly_versions/GCA_000227395.2_ASM22739v2/GCA_000227395.2_ASM22739v2_protein.faa.gz
gunzip GCA_000227395.2_ASM22739v2_protein.faa.gz
mv GCA_000227395.2_ASM22739v2_protein.faa proteins.fa

# download CDS, unzip then rename to cds.fa
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/Trypanosoma_congolense/latest_assembly_versions/GCA_000227395.2_ASM22739v2/GCA_000227395.2_ASM22739v2_translated_cds.faa.gz
gunzip GCA_000227395.2_ASM22739v2_translated_cds.faa.gz
mv GCA_000227395.2_ASM22739v2_translated_cds.faa cds.fa

# move into 'genomes' dir and download the reference genome
        cd ../genomes
        
        wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/Trypanosoma_congolense/latest_assembly_versions/GCA_000227395.2_ASM22739v2/GCA_000227395.2_ASM22739v2_genomic.fna.gz

# Unzip then rename to IL3000.fa
        gunzip GCA_000227395.2_ASM22739v2_genomic.fna.gz
        mv GCA_000227395.2_ASM22739v2_genomic.fna  IL3000.fa

# move back into 'snpEff' dir and edit the 'snpEff.config' file
        cd ../../
        
	
# edit 'snpEff.config' to add a genome

        nano snpEff.config | echo "# Database for Trypanosoma congolense" >> snpEff.config | echo "IL3000.genome : Trypanosoma congolense IL3000" >> snpEff.config |
 echo "IL3000.reference :  https://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/Trypanosoma_congolense/latest_assembly_versions/GCA_000227395.2_ASM22739v2/GCA_000227395.2_ASM22739v2_genomic.fna.gz" >> snpEff.config

# build a database for annotating the variants
        java -Xmx20g -jar snpEff.jar build -v IL3000 2>&1 | tee IL3000.build
	
# Annotate variants 
	java -jar snpEff.jar  eff IL3000 ${decomvar}  > ${annota_var}
