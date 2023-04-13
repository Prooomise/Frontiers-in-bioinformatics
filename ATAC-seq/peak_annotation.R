library(ChIPseeker)
library(GenomicFeatures)
library(patchwork)

peak <- readPeakFile("./SRR23404203.bed")
spompe <- makeTxDbFromGFF("./Mus_musculus.GRCm39.109.chr.gtf")
# peaks <- list(peak1 = peak1, peak2 = peak2)

peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000), TxDb = spompe)
write.table(peakAnno, file = 'peak.txt',sep = '\t', quote = FALSE, row.names = FALSE)

pdf("./peakannotation.pdf")
plotAnnoBar(peakAnno)
vennpie(peakAnno)
plotAnnoPie(peakAnno)
plotDistToTSS(peakAnno)
dev.off()