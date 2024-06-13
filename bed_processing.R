#BiocManager::install("ChIPpeakAnno")
#BiocManager::install("EnsDb.Hsapiens.v86")
#BiocManager::install("org.Hs.eg.db")

library("ChIPpeakAnno")
library("GenomicRanges")

library("rtracklayer")
setwd("~/nova/mo/module8/")
bed <- import("extended.bed")
gr1 <- GRanges(seqnames=bed$V1, ranges=IRanges(start=bed$V2, end=bed$V3)) 

#read in the bed file
df1<-read.table("extended.bed", header=FALSE)

#convert the peaks to a GRanges object
gr1 <- GRanges(seqnames=df1$V1, ranges=IRanges(start=df1$V2, end=df1$V3))



##annotation
library("EnsDb.Hsapiens.v86")


## create annotation file from EnsDb (Ensembl) or TxDb (transcript annotation) packages

annoData <- toGRanges(EnsDb.Hsapiens.v86, feature="gene")
annoData[1:2]

# annotate the peak GRanges object with peaks mapped to gene with a -4000 and 500 bp window around the TSS

anno.gr1 <- annotatePeakInBatch(gr1, 
                                AnnotationData=annoData, 
                                output="nearestBiDirectionalPromoters",
                                bindingRegion=c(-4000, 500))

#trim out of bound ranges
anno.gr1 <- trim(anno.gr1)




#annotate with Entrez IDs

library("org.Hs.eg.db")
anno.gr1 <- addGeneIDs(anno.gr1,"org.Hs.eg.db",IDs2Add = "entrez_id")

# list annotated peaks
head(anno.gr1)

diff_expressed = read.table("/home/eka/nova/mo/project/kallisto_output/geneNames_TAZ_S89A_overexpression_cutoff_logFC_1.txt")
#gene_names = read.table("/home/eka/nova/mo/project/kallisto_output/geneNames.txt")
head(diff_expressed)
common_genes <- intersect(diff_expressed$V1, anno.gr1$gene_name)
length(common_genes)
dim(diff_expressed)
write.table(common_genes, file="geneNames_regulated_and_chip_seq_TAZ_S89A.txt", quote=F, row.names=F, col.names=F, sep="\t")
