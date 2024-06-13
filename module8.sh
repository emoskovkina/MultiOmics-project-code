docker run --rm -v $PWD:/data ncbi/sra-tools:3.0.1 fasterq-dump --outdir /data SRR13269806
docker run --rm -v $PWD:/data ncbi/sra-tools:3.0.1 fasterq-dump --outdir /data SRR13269807
docker run -v $PWD:$PWD -w=$PWD --rm staphb/bowtie2 bowtie2 -k 1 --threads 10 -x ./Genomes/Genomes/GRCh38_noalt_as/GRCh38_noalt_as -U SRR13269806.fastq -S TAZ_ChIPseq.sam
docker run -v $PWD:$PWD -w=$PWD --rm biocontainers/bowtie2:v2.4.1_cv1 bowtie2 -k 1 --threads 10 -x ./GRCh38_noalt_as/GRCh38_noalt_as -U SRR13269807.fastq -S input.sam
docker run -v $PWD:$PWD -w=$PWD --rm miguelpmachado/sambamba:0.7.1-01 sambamba view -t 10 -S -f bam TAZ_ChIPseq.sam -o temp.bam
docker run -v $PWD:$PWD -w=$PWD --rm miguelpmachado/sambamba:0.7.1-01 sambamba view -t 10 -S -f bam input.sam -o temp2.bam
docker run -v $PWD:$PWD -w=$PWD --rm miguelpmachado/sambamba:0.7.1-01 sambamba sort -t 10 -o TAZ_ChIPseq.bam temp.bam
docker run -v $PWD:$PWD -w=$PWD --rm miguelpmachado/sambamba:0.7.1-01 sambamba sort -t 10 -o input.bam temp2.bam
docker run -v $PWD:$PWD -w=$PWD --rm resolwebio/chipseq:5.1.3 macs2 callpeak -t TAZ_ChIPseq.bam -c input.bam -f BAM -g 2.7e9 -q 0.05 -n TAZ --outdir macs2
docker run -v $PWD:$PWD -w=$PWD --rm biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1 bedtools intersect -v -a macs2/TAZ_summits.bed -b Genomes/Genomes/hg38.blacklist.bed > output.bed
docker run -v $PWD:$PWD -w=$PWD --rm biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1 bedtools slop -i output.bed -g Genomes/Genomes/genomeSize.txt -b 200 > extended.bed
docker run -v $PWD:$PWD -w=$PWD --rm biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1 bedtools getfasta -fi Genomes/Genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -bed extended.bed > out.fasta
docker run -v $PWD:$PWD -w=$PWD --rm nfcore/chipseq:latest findMotifs.pl out.fasta human findMotifsOutput/ -fasta Genomes/Genomes/gencode.v26.annotation_promoters.fasta
