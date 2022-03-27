#For downloading Pair-end SRR files in FASTQ file format as zip files
while read f; do fastq-dump --gzip --split-3 "$f";done<srr_number_pair_end
fastqc *.fastq.gz #Accessing the quality of the downloaded fastq files

#trimming the downloaded fastq files in 15 base pair in both upstream and downstream regions 
for fqfile in `ls *gz | sed 's/_1.fastq.gz//g' | sed 's/_2.fastq.gz//g'|sort -u`; do
   cutadapt -u 15 -U 15 -o ${fqfile}_1_trim.fastq.gz -p ${fqfile}_2_trim.fastq.gz ${fqfile}_1.fastq.gz ${fqfile}_2.fastq.gz
done
fastqc *_trim.fastq.gz #Again accessing the quality of the trimmed fastq files                                                                     

#Downloading transcripts index                                                  
wget http://ftp.ensembl.org/pub/release-104/fasta/saccharomyces_cerevisiae/cdna\
/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz


#getting transcripts  index                                                     
kallisto index -i transcripts.idx Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.\gz

#For quantifying the abundances of Transcripts
for fqfile in `ls *gz | sed 's/_1_trim.fastq.gz//g' | sed 's/_2_trim.fastq.gz//\
g' | sort -u`; do
   kallisto quant -i transcripts.idx -o ${fqfile} -b 100 ${fqfile}_1_trim.fastq\
.gz ${fqfile}_2_trim.fastq.gz
done
