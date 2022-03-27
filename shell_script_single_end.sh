#For downloading Single-end SRR files in FASTQ file format as zip files 
while read f; do fastq-dump --gzip "$f";done<srr_number_single_end
fastqc *.fastq.gz #Accessing the quality of the downloaded fastq files

# trimming the fastq files in 15 base pair from both upstream and downstream regions
for fqfile in `ls *gz | sed 's/.fastq.gz//g' | sort -u`; do
   cutadapt -u 15 -o ${fqfile}_trim.fastq.gz ${fqfile}.fastq.gz
done
fastqc *_trim.fastq.gz  #again checking the quality of the trimmed fastq files                                                                     

#Downloading transcripts index                                                  
wget http://ftp.ensembl.org/pub/release-104/fasta/saccharomyces_cerevisiae/cdna\
     /Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz

#retriving transcripts index of Saccharomyces cerevisae
kallisto index -i transcripts.idx Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.\gz

#quantifying the abundances of transcripts
for fqfile in `ls *gz | sed 's/_trim.fastq.gz//g' | sort -u`; do
   kallisto quant -i transcripts.idx -o ${fqfile} --single -l 200 -s 20 -b 100 ${fqfile}_trim.fastq.gz 
done
