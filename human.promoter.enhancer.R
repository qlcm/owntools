library("matrixStats")
library(GenomicRanges)
library(data.table)

gene.multipromoter <- read.table("~/Project/CpGDenLowess/Humanresult/multipromoter/gene.multi.promoter.csv", header = TRUE, sep = "\t") 


tissue.score <- gene.multipromoter[,11:28]
tissue.sd <- rowSds(as.matrix(tissue.score))
gene.multipromoter$tissuesd <- tissue.sd
diff.gene <- gene.multipromoter$geneSymbol[gene.multipromoter$tissuesd >= 0.2]
gene.multipromoter.sel <- gene.multipromoter[gene.multipromoter$geneSymbol %in% diff.gene,]
gene.multipromoter.sel <- gene.multipromoter.sel[order(gene.multipromoter.sel$geneSymbol),]

dt.gene.sel <- data.table(gene.multipromoter.sel)
dt.gene.sel.subset <- dt.gene.sel[,.(promStart = unique(promStart), promEnd = unique(promEnd)), by =  c("geneSymbol", "chrom", "tissuesd")] 
df.gene <- as.data.frame(table(dt.gene.sel.subset$geneSymbol))
df.gene.sel <- df.gene[df.gene$Freq >= 2,]
dt.gene <- dt.gene.sel.subset[dt.gene.sel.subset$geneSymbol %in% df.gene.sel$Var1,]
dt.gene.sd <- dt.gene[dt.gene$tissuesd >= 0.15,]


tissues <- list.files("/home/fsch/Project/CpGDenLowess/Humanresult/tissues/")
gene.list <- list()
for (ts in tissues){
	G.regions <- read.table(paste("/home/fsch/Project/CpGDenLowess/Humanresult/tissues/", ts,"/regionvalue/region/regions.G.csv", sep = ""), sep = "\t", header = TRUE)
	G.regions <- G.regions[G.regions$meanScore < mean(G.regions$meanScore),]
	G.Range <- GRanges(seqnames = G.regions$chrom, ranges = IRanges(G.regions$start, G.regions$end))
	Gene.Range <- GRanges(seqnames = dt.gene.sd$chrom, ranges = IRanges(dt.gene.sd$promStart, dt.gene.sd$promEnd))
	Range.Overlap <- findOverlaps(G.Range, Gene.Range)
	gene.prom.diff <- cbind(G.regions[ Range.Overlap@queryHits,], dt.gene.sd[Range.Overlap@subjectHits,])
	gene.prom.diff$tissues <- ts
	gene.list[[ts]] <- gene.prom.diff
}
gene.diff <- do.call(rbind,gene.list)


##
write.table(gene.diff,"/home/fsch/Project/CpGDenLowess/Humanresult/multipromoter/diff.gene.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
######################
rm.genes <- c("TCRA", "Vmn2r47", "Snord116", "Gm15093", "AK139882", "Vmn1r100", "ENSMUSG00000076836", "TCR-alpha chain", "AK016172", "Gm4836")
gene.names <- unique(gene.diff$geneSymbol)[!unique(gene.diff$geneSymbol) %in% rm.genes]

gene.s <- gene.multipromoter[gene.multipromoter$geneSymbol %in% gene.names,]
gene.s <- gene.s[order(gene.s$geneSymbol),]
dt.gene.s <- data.table(gene.s)
dt.gene.s <- data.table(gene.s, key = c("geneSymbol","promStart"))                                                                                                           
dt.gene.s.unique <- unique(dt.gene.s)

dt.gene.s.unique$promId <- dt.gene.s.unique[,.(id = 1:length(chrom)), by = geneSymbol]$id
dt.gene.s.unique$placenta <- NULL
tissue.matrix <- as.matrix(dt.gene.s.unique[,18:33, with = FALSE])
tissue.matrix[tissue.matrix >= 0.4] <- 1
tissue.matrix[tissue.matrix < 0.4] <- 0
df.tissue <- as.data.frame(tissue.matrix)
df.tissue$promId <- dt.gene.s.unique$promId
df.tissue$geneSymbol <- dt.gene.s.unique$geneSymbol
df.tissue$geneSymbolPromoter <- paste("P",df.tissue$promId, sep = "")

write.table(df.tissue, "/home/fsch/Project/CpGDenLowess/result2/multipromoter/TissuediffPromNew.csv", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


################ enhancers 
~/Project/CpGDenLowess/codes/RegionGetValues.Rscript -m ~/Rendata/ -l "BSgenome.Mmusculus.UCSC.mm9" -n "Mmusculus" -t "mm9" -o ./ mouse.enhancers.regions
##
library("matrixStats")
library(GenomicRanges)
library(data.table)


mouse.enhancer <- read.table("~/Project/CpGDenLowess/enhancer/FANTOM5/mm9/mouse.enhancers.bed", sep = "\t")
bonemarrow.gap <- read.table("/home/fsch/Project/CpGDenLowess/mm9result/tissues/bonemarrow/regionValue/value/Gap/region.tissues.score.txt", header = T, sep = "\t")
bonemarrow.overlap <- read.table("/home/fsch/Project/CpGDenLowess/mm9result/tissues/bonemarrow/regionValue/value/overlap/region.tissues.score.txt", header = T, sep = "\t")

bonemarrow.gap <- bonemarrow.gap[!is.na(bonemarrow.gap$bonemarrow_score) ,]
bonemarrow.gap.sel <- bonemarrow.gap[bonemarrow.gap$bonemarrow_score < 0.5,]
bonemarrow.gap.sel <- bonemarrow.gap.sel[bonemarrow.gap.sel$meanDensity < 0.15,]

bonemarrow.overlap <- bonemarrow.overlap[!is.na(bonemarrow.overlap$bonemarrow_score) ,]
bonemarrow.overlap.sel <- bonemarrow.overlap[bonemarrow.overlap$bonemarrow_score > 0.5,]
bonemarrow.overlap.sel <- bonemarrow.overlap.sel[bonemarrow.overlap.sel$meanDensity >= 0.15,]


E.Range <- GRanges(seqnames = mouse.enhancer$V1, ranges = IRanges(mouse.enhancer$V2, mouse.enhancer$V3))
Gap.Range <- GRanges(seqnames = bonemarrow.gap.sel$chrom, ranges = IRanges(bonemarrow.gap.sel$start, bonemarrow.gap.sel$end))
Range.G.Overlap <- findOverlaps(E.Range, Gap.Range)

Overlap.Range <- GRanges(seqnames = bonemarrow.overlap.sel$chrom, ranges = IRanges(bonemarrow.overlap.sel$start, bonemarrow.overlap.sel$end))
Range.O.Overlap <- findOverlaps(E.Range, Overlap.Range)

# > 22648 44459 length(unique(Range.O.Overlap@subjectHits)) -> 1470 length(unique(Range.O.Overlap@queryHits)) enhancer -> 1669
# > 13163 44459 length(unique(Range.G.Overlap@subjectHits)) -> 796 length(unique(Range.G.Overlap@queryHits)) -> 986

Gap.enhancer <- cbind(mouse.enhancer[ Range.G.Overlap@queryHits,], bonemarrow.gap.sel[Range.G.Overlap@subjectHits,])
Overlap.enhancer <- cbind(mouse.enhancer[ Range.O.Overlap@queryHits,], bonemarrow.overlap.sel[Range.O.Overlap@subjectHits,])

Overlap.enhancer.diff <- Overlap.enhancer[Overlap.enhancer$diffScore > 0.4,]
# > length(unique(Overlap.enhancer.diff$V4)) -> 723
Gap.enhancer.diff <- Gap.enhancer[Gap.enhancer$diffScore > 0.4,]
# > length(unique(Gap.enhancer.diff$V4)) -> 742

### super enhancer
library("matrixStats")
library(GenomicRanges)
library(data.table)

mouse.enhancer <- read.table("/home/fsch/Project/CpGDenLowess/enhancer/dbSUPER/BoneMarrow.bed", sep = "\t")
bonemarrow.gap <- read.table("/home/fsch/Project/CpGDenLowess/mm9result/tissues/bonemarrow/regionValue/value/Gap/region.tissues.score.txt", header = T, sep = "\t")
bonemarrow.overlap <- read.table("/home/fsch/Project/CpGDenLowess/mm9result/tissues/bonemarrow/regionValue/value/overlap/region.tissues.score.txt", header = T, sep = "\t")

bonemarrow.gap <- bonemarrow.gap[!is.na(bonemarrow.gap$bonemarrow_score) ,]
bonemarrow.gap.sel <- bonemarrow.gap[bonemarrow.gap$bonemarrow_score < 0.5,]
bonemarrow.gap.sel <- bonemarrow.gap.sel[bonemarrow.gap.sel$meanDensity < 0.15,]

bonemarrow.overlap <- bonemarrow.overlap[!is.na(bonemarrow.overlap$bonemarrow_score) ,]
bonemarrow.overlap.sel <- bonemarrow.overlap[bonemarrow.overlap$bonemarrow_score > 0.5,]
bonemarrow.overlap.sel <- bonemarrow.overlap.sel[bonemarrow.overlap.sel$meanDensity >= 0.15,]


E.Range <- GRanges(seqnames = mouse.enhancer$V1, ranges = IRanges(mouse.enhancer$V2, mouse.enhancer$V3))
Gap.Range <- GRanges(seqnames = bonemarrow.gap.sel$chrom, ranges = IRanges(bonemarrow.gap.sel$start, bonemarrow.gap.sel$end))
Range.G.Overlap <- findOverlaps(E.Range, Gap.Range)

Overlap.Range <- GRanges(seqnames = bonemarrow.overlap.sel$chrom, ranges = IRanges(bonemarrow.overlap.sel$start, bonemarrow.overlap.sel$end))
Range.O.Overlap <- findOverlaps(E.Range, Overlap.Range)

# > 22648 350 length(unique(Range.O.Overlap@subjectHits)) -> 229 length(unique(Range.O.Overlap@queryHits)) -> 138
# > 13163 350 length(unique(Range.G.Overlap@subjectHits)) -> 46 length(unique(Range.G.Overlap@queryHits)) -> 41

Gap.enhancer <- cbind(mouse.enhancer[ Range.G.Overlap@queryHits,], bonemarrow.gap.sel[Range.G.Overlap@subjectHits,])
Overlap.enhancer <- cbind(mouse.enhancer[ Range.O.Overlap@queryHits,], bonemarrow.overlap.sel[Range.O.Overlap@subjectHits,])

Overlap.enhancer.diff <- Overlap.enhancer[Overlap.enhancer$diffScore > 0.4,]
# > length(unique(Overlap.enhancer.diff$V4)) -> 67
Gap.enhancer.diff <- Gap.enhancer[Gap.enhancer$diffScore > 0.4,]
# > length(unique(Gap.enhancer.diff$V4)) -> 31

############### bonemarrow h3k4me1 enhancer
library("matrixStats")
library(GenomicRanges)
library(data.table)

bonemarrow.h3k4me1 <- fread("~/Renchipseq/Macs2bedGraph/bonemarrowH3K4me1/bonemarrowH3K4me1sorted.br", sep = "\t")  
colnames(bonemarrow.h3k4me1) <- c("chrom", "start", "end", "value")
bonemarrow.h3k4me1$end  <- bonemarrow.h3k4me1$end - 1
dt.bonemarrow.h3k4me1 <- data.table(bonemarrow.h3k4me1)

dt.bonemarrow.h3k4me1 <- dt.bonemarrow.h3k4me1[dt.bonemarrow.h3k4me1$chrom != "chrM",]
dt.bonemarrow.h3k4me1 <- dt.bonemarrow.h3k4me1[dt.bonemarrow.h3k4me1$chrom != "chrY",]
dt.bonemarrow.h3k4me1$id <- 1:nrow(dt.bonemarrow.h3k4me1)
dt.bonemarrow.h3k4me1 <- dt.bonemarrow.h3k4me1[dt.bonemarrow.h3k4me1$value > 0,]

mouse.enhancer <- read.table("~/Project/CpGDenLowess/enhancer/FANTOM5/mm9/mouse.enhancers.bed", sep = "\t")
dt.mouse.enhancer <- data.table(mouse.enhancer)
dt.mouse.enhancer$id <- 1:nrow(dt.mouse.enhancer)

h3k4me1.list <- list()
chr.names <- as.character(unique(dt.bonemarrow.h3k4me1$chrom))
for (chr.name in chr.names){
	dt.bonemarrow.h3k4me1.exp <- dt.bonemarrow.h3k4me1[dt.bonemarrow.h3k4me1$chrom == chr.name,][,.(pos = start:end), by = c("chrom", "value","id")]
	dt.mouse.enhancer.exp <- dt.mouse.enhancer[dt.mouse.enhancer$V1 == chr.name,][,.(pos = V2:V3), by = c("V1", "id")]
	be.merge <- merge(dt.bonemarrow.h3k4me1.exp, dt.mouse.enhancer.exp, by = "pos")
	h3k4me1.list[[chr.name]] <- be.merge[,.(max(value)), by = "id.y"] 

}

h3k4me1.value <- do.call(rbind, h3k4me1.list)
colnames(h3k4me1.value) <- c("id", "value")
dt.mouse.enhancer.value <- merge(dt.mouse.enhancer, h3k4me1.value, by = "id")

bonemarrow.enhancer <- dt.mouse.enhancer.value[dt.mouse.enhancer.value$value > 5, ]

bonemarrow.gap <- read.table("/home/fsch/Project/CpGDenLowess/mm9result/tissues/bonemarrow/regionValue/value/Gap/region.tissues.score.txt", header = T, sep = "\t")
bonemarrow.overlap <- read.table("/home/fsch/Project/CpGDenLowess/mm9result/tissues/bonemarrow/regionValue/value/overlap/region.tissues.score.txt", header = T, sep = "\t")

bonemarrow.gap <- bonemarrow.gap[!is.na(bonemarrow.gap$bonemarrow_score) ,]
bonemarrow.gap.sel <- bonemarrow.gap[bonemarrow.gap$bonemarrow_score < 0.5,]
bonemarrow.gap.sel <- bonemarrow.gap.sel[bonemarrow.gap.sel$meanDensity < 0.15,]

bonemarrow.overlap <- bonemarrow.overlap[!is.na(bonemarrow.overlap$bonemarrow_score) ,]
bonemarrow.overlap.sel <- bonemarrow.overlap[bonemarrow.overlap$bonemarrow_score > 0.5,]
bonemarrow.overlap.sel <- bonemarrow.overlap.sel[bonemarrow.overlap.sel$meanDensity >= 0.15,]


E.Range <- GRanges(seqnames = bonemarrow.enhancer$V1, ranges = IRanges(bonemarrow.enhancer$V2, bonemarrow.enhancer$V3))
Gap.Range <- GRanges(seqnames = bonemarrow.gap.sel$chrom, ranges = IRanges(bonemarrow.gap.sel$start, bonemarrow.gap.sel$end))
Range.G.Overlap <- findOverlaps(E.Range, Gap.Range)

Overlap.Range <- GRanges(seqnames = bonemarrow.overlap.sel$chrom, ranges = IRanges(bonemarrow.overlap.sel$start, bonemarrow.overlap.sel$end))
Range.O.Overlap <- findOverlaps(E.Range, Overlap.Range)

# > 22648 16493 length(unique(Range.O.Overlap@subjectHits)) -> 562 length(unique(Range.O.Overlap@queryHits)) -> 622
# > 13163 16493 length(unique(Range.G.Overlap@subjectHits)) -> 411 length(unique(Range.G.Overlap@queryHits)) -> 492

Gap.enhancer <- cbind(bonemarrow.enhancer[ Range.G.Overlap@queryHits,], bonemarrow.gap.sel[Range.G.Overlap@subjectHits,])
Overlap.enhancer <- cbind(bonemarrow.enhancer[ Range.O.Overlap@queryHits,], bonemarrow.overlap.sel[Range.O.Overlap@subjectHits,])

Overlap.enhancer.diff <- Overlap.enhancer[Overlap.enhancer$diffScore > 0.4,]
# > length(unique(Overlap.enhancer.diff$V4)) -> 270
Gap.enhancer.diff <- Gap.enhancer[Gap.enhancer$diffScore > 0.4,]
# > length(unique(Gap.enhancer.diff$V4)) -> 387



################################### #############################################
region.tissues <- read.table("~/Project/CpGDenLowess/enhancer/FANTOM5/mm9/region.tissues.score.txt", header = T, sep = "\t")
region.tissues <- region.tissues[!is.na(region.tissues$diffScore),]
region.tissues.sel <- region.tissues[region.tissues$diffScore > 0.4,]
region.tissues.sel <- region.tissues.sel[region.tissues.sel$cgnum > 10,]

