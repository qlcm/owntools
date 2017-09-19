library("data.table")
load("/home/qzzh/Renbing/DMR/getvalue/chrX.Rdata")
chrX <- chrX[order(chrX$start),]
H3k4me1 <- fread("/home/qzzh/Renbing/Chip_seq/ToGB/H3k4me1_treat.bdg")
chrX.H3k4me1 <- H3k4me1[H3k4me1$V1 == "chrX",]
colnames(chrX.H3k4me1) <- c("chrom","start","end","value")
chrX.H3k4me1 <- chrX.H3k4me1[chrX.H3k4me1$value > 0 ,]
chrX.H3k4me1$start <- chrX.H3k4me1$start + 1
chrX.H3k4me1$marker <- c(1:nrow(chrX.H3k4me1))
dt.chrX.H3k4me1 <- chrX.H3k4me1[,.(chrom = chrom, posi = c(start : end), value = value),by = marker]
dt.chrX.H3k4me1$marker <- NULL
chrX.dmr.region <- chrX[,1:3]
chrX.dmr.region$marker <- c(1:nrow(chrX.dmr.region))
dt.chrX.dmr.region <- as.data.table(chrX.dmr.region)
dt.chrX.dmr.region$chrom <- as.character(dt.chrX.dmr.region$chrom)
dt.chrX.dmr.region$start <- as.integer(dt.chrX.dmr.region$start)
dt.chrX.dmr.region$end <- as.integer(dt.chrX.dmr.region$end)
dt.chrX.dmr.region.split <- dt.chrX.dmr.region[,.(chrom = chrom, posi = c(start : end)),by = marker]
dt.dmr.value <- merge(dt.chrX.H3k4me1,dt.chrX.dmr.region.split,by.x = c("posi","chrom"),by.y = c("posi","chrom"))
dt.dmr.chip.value <- dt.dmr.value[,.(H3k4me1 = sum(value)),by = marker]
chrX.add.value <- merge(dt.chrX.dmr.region,dt.dmr.chip.value,by="marker",all = TRUE)
chrX$HCC_H3k4me1 <- chrX.add.value$H3k4me1
tmp <- grep("[^chrX]",ls())
rm(list=ls()[tmp])
gc()
save(chrX,file="/home/qzzh/chrX.Rdata")
