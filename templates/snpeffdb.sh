#!/usr/bin/bash
# install snpeff
	wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
# unpack snpeff executables
	unzip snpEff_latest_core.zip 
# remove the compressed file
	rm snpEff_latest_core.zip
# move into the sneff directory 
  	cd ./snpEff/
# make a directory 'IL3000' and  'genomes' in './data' directory
	mkdir -p ./data/IL3000 ./data/genomes
	cd ./data/IL3000 
# download Trypanosoma congolense gff file
	wget https://tritrypdb.org/common/downloads/Current_Release/TcongolenseIL3000/gff/data/TriTrypDB-55_TcongolenseIL3000.gff 

# rename TriTrypDB-55_TcongolenseIL3000.gff to genes.gff

	mv TriTrypDB-55_TcongolenseIL3000.gff  genes.gff
# move into 'genomes' dir and download the reference genome
	cd ../genomes
	wget https://tritrypdb.org/common/downloads/Current_Release/TcongolenseIL3000/fasta/data/TriTrypDB-55_TcongolenseIL3000_Genome.fasta 

# rename TriTrypDB-55_TcongolenseIL3000_Genome.fasta to Genome.fasta

	mv TriTrypDB-55_TcongolenseIL3000_Genome.fasta Genome.fasta

# move back into 'snpEff' dir and edit the 'snpEff.config' file
	cd ../../
# edit 'snpEff.config' to add a genome

	nano snpEff.config | echo "# Database for Trypanosoma congolense" >> snpEff.config | echo "IL3000.genome : Trypanosoma congolense IL3000" >> snpEff.config |
 echo "IL3000.reference : https://tritrypdb.org/common/downloads/Current_Release/TcongolenseIL3000/fasta/data/TriTrypDB-55_TcongolenseIL3000_Genome.fasta" >> snpEff.config

#		      # Database for Trypanosoma congolense 
#		      IL3000.genome : Trypanosoma congolense IL3000 
#		      IL3000.reference : https://tritrypdb.org/common/downloads/Current_Release/TcongolenseIL3000/fasta/data/TriTrypDB-55_TcongolenseIL3000_Genome.fasta

# echo "# Database for Trypanosoma congolense" > file.txt
# echo "IL3000.genome : Trypanosoma congolense IL3000" >> file.txt
# echo "IL3000.reference : https://tritrypdb.org/common/downloads/Current_Release/TcongolenseIL3000/fasta/data/TriTrypDB-55_TcongolenseIL3000_Genome.fasta" >>file.txt

#
# use 'tail' to check added lines as; [optional]

#	tail -n 5 snpEff.config
# build a database for annotating the variants
	java -jar snpEff.jar build -gff3 -v IL3000 
	java ijar snpEff.jar  eff IL3000 ${decomvar}  > ${annota_var}
