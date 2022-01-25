
#!/usr/bin/bash

#Indexing the file
bwa index ${ref_ch}

#Aligning reads
bwa mem -t 8 -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPU:${sample_id}\\tPL:ILLUMINA" -M ${ref_ch} ${read_R1} ${read_R2} > ${aligned}


