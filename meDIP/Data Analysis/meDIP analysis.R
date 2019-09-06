#loading required libraries
library(scales)
library(ggplot2)
library(reshape2)
library(ggthemes)
library(csaw)
library(ChIPpeakAnno)
library(RMariaDB)
library(GenomicFeatures)
library(ggsci)
library(clusterProfiler)
library(BSgenome.Drerio.UCSC.danRer10)
#downloading genome files
txdb <- makeTxDbFromUCSC("danRer10", "ensGene")

#Processing Positive Enrichment Data
#downloading and fixing input files
setwd("/Volumes/Personal/Tet2 Bam/MACS Output")
meDIP_1<-read.csv("meDIP2x_TIvsCI_peaks.csv")
meDIP_2<-read.csv("meDIP2x_TIIvsCII_peaks.csv")
meDIP_3<-read.csv("meDIP2x_TIIIvsCIII_peaks.csv")
hmeDIP_1<-read.csv("hmeDIP3x_TIvsCI_peaks.csv")
hmeDIP_2<-read.csv("hmeDIP3x_TIIvsCII_peaks.csv")
hmeDIP_3<-read.csv("hmeDIP3x_TIIIvsCIII_peaks.csv")
colnames(meDIP_1)[7]<-"-log10pvalue"
colnames(meDIP_1)[9]<-"-log10qvalue"
colnames(meDIP_2)[7]<-"-log10pvalue"
colnames(meDIP_2)[9]<-"-log10qvalue"
colnames(meDIP_3)[7]<-"-log10pvalue"
colnames(meDIP_3)[9]<-"-log10qvalue"
colnames(hmeDIP_1)[7]<-"-log10pvalue"
colnames(hmeDIP_1)[9]<-"-log10qvalue"
colnames(hmeDIP_2)[7]<-"-log10pvalue"
colnames(hmeDIP_2)[9]<-"-log10qvalue"
colnames(hmeDIP_3)[7]<-"-log10pvalue"
colnames(hmeDIP_3)[9]<-"-log10qvalue"
#overlapping genomic regions
me_gr1 <- toGRanges(meDIP_1, format="MACS", header=FALSE, skip=3)
me_gr2 <- toGRanges(meDIP_2, format="MACS", header=FALSE, skip=3)
me_gr3 <- toGRanges(meDIP_3, format="MACS", header=FALSE, skip=3)
hme_gr1 <- toGRanges(hmeDIP_1, format="MACS", header=FALSE, skip=3)
hme_gr2 <- toGRanges(hmeDIP_2, format="MACS", header=FALSE, skip=3)
hme_gr3 <- toGRanges(hmeDIP_3, format="MACS", header=FALSE, skip=3) 
ol_me <- findOverlapsOfPeaks(me_gr1, me_gr2, me_gr3)
ol_hme <- findOverlapsOfPeaks(hme_gr1, hme_gr2, hme_gr3)
overlaps_me <- ol_me$peaklist[["me_gr1///me_gr2///me_gr3"]]
overlaps_hme <- ol_hme$peaklist[["hme_gr1///hme_gr2///hme_gr3"]]
#annotation of the peaks into genomic regions
annoData <- toGRanges(txdb, feature="gene")
binOverFeature(overlaps_me, annotationData=annoData,
               radius=5000, nbins=20, FUN=length, errFun=0,
               ylab="count", 
               main="Distribution of aggregated peak numbers around TSS - Methylation")

binOverFeature(overlaps_hme, annotationData=annoData,
               radius=5000, nbins=20, FUN=length, errFun=0,
               ylab="count", 
               main="Distribution of aggregated peak numbers around TSS - Hydroxymethylation")

aCR_me<-assignChromosomeRegion(overlaps_me, nucleotideLevel=FALSE, 
                               precedence=c("Promoters", "immediateDownstream", 
                                            "fiveUTRs", "threeUTRs", 
                                            "Exons", "Introns"), TxDb=txdb)

aCR_hme<-assignChromosomeRegion(overlaps_hme, nucleotideLevel=FALSE, 
                                precedence=c("Promoters", "immediateDownstream", 
                                             "fiveUTRs", "threeUTRs", 
                                             "Exons", "Introns"), TxDb=txdb)
#extracting percentage data
perdata_me<-aCR_me$percentage
perdata_me<-as.data.frame(perdata_me)
perdata_me<- cbind(subjectHits = rownames(perdata_me), perdata_me)
rownames(perdata_me)<-NULL
perdata_hme<-aCR_hme$percentage
perdata_hme<-as.data.frame(perdata_hme)
perdata_me$Sample<-"Methylation"
perdata_hme$Sample<-"Hydroxymethylation"
names(perdata_hme)[2]<-paste("Frequency")
names(perdata_me)[2]<-paste("Frequency")
percent_df<-rbind(perdata_hme, perdata_me)
newlabs <- c("Promoters", "Immediate downstream", "5'UTRs", "3'UTRs", "Exons", "Introns", "Intergenic Region")

#Processing of negative enrichment data
#input files
setwd("/Volumes/Personal/Tet2 Bam/macs2_neg_ER_xls4Vino")
negmeDIP_1<-read.csv("negER_meDIP2x_TIvsCI_peaks.csv")
negmeDIP_2<-read.csv("negER_meDIP2x_TIIvsCII_peaks.csv")
negmeDIP_3<-read.csv("negER_meDIP2x_TIIIvsCIII_peaks.csv")
neghmeDIP_1<-read.csv("negER_hmeDIP2x_TIvsCI_peaks.csv")
neghmeDIP_2<-read.csv("negER_hmeDIP2x_TIIvsCII_peaks.csv")
neghmeDIP_3<-read.csv("negER_hmeDIP2x_TIIIvsCIII_peaks.csv")
colnames(negmeDIP_1)[7]<-"-log10pvalue"
colnames(negmeDIP_1)[9]<-"-log10qvalue"
colnames(negmeDIP_2)[7]<-"-log10pvalue"
colnames(negmeDIP_2)[9]<-"-log10qvalue"
colnames(negmeDIP_3)[7]<-"-log10pvalue"
colnames(negmeDIP_3)[9]<-"-log10qvalue"
colnames(neghmeDIP_1)[7]<-"-log10pvalue"
colnames(neghmeDIP_1)[9]<-"-log10qvalue"
colnames(neghmeDIP_2)[7]<-"-log10pvalue"
colnames(neghmeDIP_2)[9]<-"-log10qvalue"
colnames(neghmeDIP_3)[7]<-"-log10pvalue"
colnames(neghmeDIP_3)[9]<-"-log10qvalue"
#Overlapping peaks
negme_gr1 <- toGRanges(negmeDIP_1, format="MACS", header=FALSE, skip=3)
negme_gr2 <- toGRanges(negmeDIP_2, format="MACS", header=FALSE, skip=3)
negme_gr3 <- toGRanges(negmeDIP_3, format="MACS", header=FALSE, skip=3)
neghme_gr1 <- toGRanges(neghmeDIP_1, format="MACS", header=FALSE, skip=3)
neghme_gr2 <- toGRanges(neghmeDIP_2, format="MACS", header=FALSE, skip=3)
neghme_gr3 <- toGRanges(neghmeDIP_3, format="MACS", header=FALSE, skip=3) 
ol_me_neg <- findOverlapsOfPeaks(negme_gr1, negme_gr2,negme_gr3)
ol_hme_neg <- findOverlapsOfPeaks(neghme_gr1, neghme_gr2, neghme_gr3)
overlaps_me_neg <- ol_me_neg$peaklist[["negme_gr1///negme_gr2///negme_gr3"]]
overlaps_hme_neg <- ol_hme_neg$peaklist[["neghme_gr1///neghme_gr2///neghme_gr3"]]
#Annotating genomic ranges to genomic features
annoData <- toGRanges(txdb, feature="gene")
binOverFeature(overlaps_me_neg, annotationData=annoData,
               radius=5000, nbins=20, FUN=length, errFun=0,
               ylab="count", 
               main="Distribution of aggregated peak numbers around TSS - Methylation")

binOverFeature(overlaps_hme_neg, annotationData=annoData,
               radius=5000, nbins=20, FUN=length, errFun=0,
               ylab="count", 
               main="Distribution of aggregated peak numbers around TSS - Hydroxymethylation")

aCR_me_neg<-assignChromosomeRegion(overlaps_me_neg, nucleotideLevel=FALSE, 
                                   precedence=c("Promoters", "immediateDownstream", 
                                                "fiveUTRs", "threeUTRs", 
                                                "Exons", "Introns"), TxDb=txdb)

aCR_hme_neg<-assignChromosomeRegion(overlaps_hme_neg, nucleotideLevel=FALSE, 
                                    precedence=c("Promoters", "immediateDownstream", 
                                                 "fiveUTRs", "threeUTRs", 
                                                 "Exons", "Introns"), TxDb=txdb)
#Extracting percentage data
perdata_me_neg<-aCR_me_neg$percentage
perdata_me_neg<-as.data.frame(perdata_me_neg)
perdata_me_neg<- cbind(subjectHits = rownames(perdata_me_neg), perdata_me_neg)
rownames(perdata_me_neg)<-NULL
perdata_hme_neg<-aCR_hme_neg$percentage
perdata_hme_neg<-as.data.frame(perdata_hme_neg)
perdata_me_neg$Sample<-"Methylation"
perdata_hme_neg$Sample<-"Hydroxymethylation"
names(perdata_hme_neg)[2]<-paste("Frequency")
names(perdata_me_neg)[2]<-paste("Frequency")
percent_negdf<-rbind(perdata_hme_neg, perdata_me_neg)
newlabs <- c("Promoters", "Immediate downstream", "5'UTRs", "3'UTRs", "Exons", "Introns", "Intergenic Region")

#plotting negative enrichment plot
q<-ggplot(data=percent_negdf, aes(x=subjectHits, y=Frequency, fill = Sample)) + geom_bar(stat="identity", position=position_dodge(1)) + xlab("") + ylab("Frequency") + theme_bw(base_size = 22) + theme ()+ scale_fill_aaas()+coord_flip()+scale_x_discrete(labels= newlabs)
q+theme(legend.position = c(0.6, 0.2), legend.title = element_blank())+theme(legend.background = element_rect(size=0.5, linetype="solid", colour ="black"))

#plotting both positive and negative enrichment dat
percent_df$Process<-"Positive Peaks"
percent_negdf$Process<-"Negative Peaks"
plot_df<-rbind(percent_df, percent_negdf)
plot_df$Sample= factor(plot_df$Sample, levels=c('Methylation','Hydroxymethylation'))
r<-ggplot(data=plot_df, aes(x=subjectHits, y=Frequency, fill = Process)) + geom_bar(stat="identity", position=position_dodge(0.7), width = 0.7) + facet_grid(Sample~.) + xlab("") + ylab("Frequency") + theme_bw(base_size = 22) + theme ()+ scale_fill_aaas()+coord_flip()+scale_x_discrete(labels= newlabs)
r+theme(legend.position = c(0.65, 0.1), legend.title = element_blank())+theme(legend.background = element_rect(size=0.5, linetype="solid", colour ="black"))

#extracting sequence files to perform rsat analysis
genes <- genes(txdb)
prom <- promoters(genes, upstream=2000, downstream=200)
prompk_hme = findOverlaps(overlaps_hme, prom)
prompk_me = findOverlaps(overlaps_me, prom)
prom_hme = genes[subjectHits(prompk_hme)]
prom_me = genes[subjectHits(prompk_me)]
geneids_hme <- unique(names(prom_hme))
geneids_me <- unique(names(prom_me))
promtme = subsetByOverlaps(overlaps_me, prom)
promthme = subsetByOverlaps(overlaps_hme, prom)
seq1<-getAllPeakSequence(promtme, genome = Drerio)
seq2<-getAllPeakSequence(promthme, genome = Drerio)
write2FASTA(seq1, file="promoter_me.fa")
write2FASTA(seq2, file="promoter_hme.fa")
