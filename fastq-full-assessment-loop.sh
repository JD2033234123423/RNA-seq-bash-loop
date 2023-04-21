#!/bin/bash
## loop for renaming file in numerical order
n=1
for i in *_1.fastq
do
new=$(printf "%01d_1.fastq" "$n")
mv "${i}" "$new"
let n=n+1
done
n=1
for i in *_2.fastq
do
new=$(printf "%01d_2.fastq" "$n")
mv "${i}" "$new"
let n=n+1
done
## loop for chaning files into can-bc(number)
for i in {1..14}
do
mv ${i}_1.fastq can-bc${i}_1.fastq
mv ${i}_2.fastq can-bc${i}_2.fastq
done
##listing of various sequence trimming and modifiying modules
for i in {1..14} 
do
## fastqc for quality assessment of raw reads
module load fastqc-0.11.9-gcc-8.3.1-m5eb2rt
fastqc can-bc${i}_1.fastq
fastqc can-bc${i}_2.fastq
## trimmomatic for automatic trimming of raw reads and automatic trimming of paired end reads
module load  trimmomatic-0.39-gcc-8.3.1-uaedgca
trimmomatic PE can-bc${i}_1.fastq can-bc${i}_2.fastq -baseout can-${i}-trim.fastq ILLUMINACLIP:/mnt/clusters/hawker/data/classdata/Bioinformatics/REFS/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
module unload  trimmomatic-0.39-gcc-8.3.1-uaedgca
## fastqc for quality assessment of trimmed paired reads
fastqc can-${i}-trim_1P.fastq
fastqc can-${i}-trim_2P.fastq
module unload fastqc-0.11.9-gcc-8.3.1-m5eb2rt
##SPAdes for assembly of trimmed reads
module load SPAdes/
spades.py -t 4 -1 can-${i}-trim_1P.fastq -2 can-${i}-trim_2P.fastq -o ass-${i}
module unload SPAdes
## comparision of mtDNA to database using blast, placing format out as a .tsv
module load blast 
blastx -query ass-${i}/contigs.fasta -db ~/classdata/Bioinformatics/REFS/blastdb/mito_pro -outfmt 6 -num_threads 4 -out blast${i}
module unload blast
## automatic genome annotation of assemblies using mitochondria kingdom and coding table for verterbrates
module load prokka 
prokka --outdir prokka-${i} --kingdom Mitochondria --gcode 2 ass-${i}/contigs.fasta
module purge

done
## quast for quality comparision of the assemblies done using SPAdes
module load quast 
quast.py -o quast_summary ass*/contigs.fasta
module unload quast
## multiqc for assessment overview of all data
module load multiqc/1.9 
multiqc .
module unload multiqc/1.9
## making new directory for fastqs files to be placed in
mkdir fastqs
##movement of fastq files into the directory 
mv *fastq fastqs/

for i in {1..14}
do
##bedtools for abstraction of 16S gene from fasta files for phylogenetic comparision
module load bedtools/2.29.1
grep '16S' prokka-${i}/PROKKA*.gff > ass-${i}-16S.gff
bedtools getfasta -fi ass-${i}/contigs.fasta -bed ass-${i}-16S.gff > ass-${i}-16S.fasta
## awk for movement and naming of long contig files into assembly data names
awk '/^>/ {gsub(/.fa(sta)?$/,"",FILENAME);printf(">%s\n",FILENAME);next;} {print}' ass-${i}-16S.fasta >> all_16S_genes.fasta
module purge
done 
#set up loop for snippy alignment

#Load modules for snippy
module load snippy/v4.6.0
module load samtools-1.10-gcc-8.3.1-zyy7kv4
module load python/3.7.4

for i in {1..14}
do

#snippy command 
snippy --outdir snip-${i} --reference kf857179.fasta --R1 fastqs/can-${i}-trim_1P.fastq --R2 fastqs/can-${i}-trim_2P.fastq --force

done

#calculate alignments 
snippy-core --ref kf857179.fasta --prefix output-alignment snip*

#unload snippy modules
module load samtools-1.10-gcc-8.3.1-zyy7kv4
module load python/3.7.4
module unload snippy/v4.6.0


# snp-sites module
module load snp-sites/2.5.1

#varibel sites between sequences
snp-sites -cb -o genomes.aln output-alignment.full.aln

#unload snp-sites module
module unload snp-sites/2.5.1

#load fasttree module
module load fasttree-2.1.10-gcc-8.3.1-x6l3qac

#create tree with FastTree
FastTreeMP -nt -gtr < genomes.aln > tree.tre