#!/bin/bash

# Create dirs
mkdir -p temp/gene_counts/
mkdir -p temp/hisat2_log/
mkdir -p temp/htseq_log/
mkdir -p temp/bam_stat/

# Get all files to be processed
selected_files=$(ls raw/bb875_rawdata/ | grep "fq.gz" | sed -e 's/_..fq.gz//g' | uniq -)

# Process one by one
for f in $selected_files; 
do
  # Check whether processed/processing already and skip if it does
  if [ -f temp/gene_counts/$f'.gene_counts' ] || [ -f temp/$f'.bam' ];
  then
  # echo "temp/gene_counts/$f'.gene_counts' exists already, skipping ...";
  echo "$f processed already, skipping ...";
  continue;
  fi;
  
  echo 'Processing' $f;
  file1=$f'_1.fq.gz'
  file2=$f'_2.fq.gz'
  
  # Align to grch38 using hisat: +/-15'
  echo '  Aligning ...';
  (hisat2 -p 8 --no-mixed --no-discordant --rna-strandness RF --no-unal -x downloads/genomes/igenomes/hg38/genome -1 raw/bb875_rawdata/$file1 -2 raw/bb875_rawdata/$file2 | samtools view -bS - > temp/$f'.bam') >& temp/hisat2_log/$f'_log.txt'
  
  # Quantify to gencode 29 using htseq-count: +/-45'
  echo '  Quantifying ...';
  (samtools view temp/$f'.bam' | htseq-count -m intersection-strict -s reverse -i gene_id - downloads/gencode/gencode.v29.annotation.gtf > temp/gene_counts/$f'.gene_counts') >& /dev/stdout | tail -n200 > temp/htseq_log/$f'_log.txt'
  # samtools view -H raw/AS1h_1.bam # "chr1" format, "1" format for igenomes grch38
  
  # Get bam stats: includes insert sizes, read lengths, ...
  samtools stats temp/$f'.bam' > temp/bam_stat/$f'_stat_log.txt'
  
  # Rm bam file
  echo '  Finishing ...';
  rm temp/$f'.bam'
  echo ''

done