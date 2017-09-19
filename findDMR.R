library("BSgenome.Hsapiens.UCSC.hg38")
library("methyutils")
library("reshape2")
library("data.table")
library(zoo)
args <- commandArgs(T)
USAGE <- function(){
	message("USAGE:\n\tRscript findDMR.R <mtbrV.cg.path,normal cancer> <outpath(dmr result)>\n")
					}
USAGE()
chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
input1 <- args[1]
input2 <- args[2]
output <- args[3]
MultiAddPv <- function(mg.data){
		sn <- (ncol(mg.data)-2) / 2
		mg.data <- mg.data[order(mg.data$chrom,mg.data$posi),]
		m <- melt(mg.data, id=c("chrom", "posi"),variable.name="read", value.name="count")
		m <- m[order(m$chrom,m$posi),]
		cp <- MultiChisqTest(m$count,nrow=sn,ncol=2)
		mg.data$pv<-cp
		return(mg.data)
							}

for (chr in chrs){
     message(chr," is running ",date())
     dna.seq <- Hsapiens[[chr]]
     I1 <- paste(input1,"/",chr,".Rdata",sep="")
     load(I1)
     cg.mtbr$rC <- cg.mtbr$rC_p + cg.mtbr$rC_n
     cg.mtbr$rT <- cg.mtbr$rT_p + cg.mtbr$rT_n
     rC <- integer(length(dna.seq))
     rC[cg.mtbr$posi] <- cg.mtbr$rC
     HMEC.rC <- rC
     rT <- integer(length(dna.seq))
     rT[cg.mtbr$posi] <- cg.mtbr$rT
     HMEC.rT <- rT
     HMEC <- data.frame(chrom=cg.mtbr$chrom,posi=cg.mtbr$posi,HMEC_rC=cg.mtbr$rC,HMEC_rT=cg.mtbr$rT)
     HMEC[is.na(HMEC)] <- 0
     I2 <- paste(input2,"/",chr,".Rdata",sep="")
     load(I2)
     cg.mtbr$rC <- cg.mtbr$rC_p + cg.mtbr$rC_n
     cg.mtbr$rT <- cg.mtbr$rT_p + cg.mtbr$rT_n
     rC <- integer(length(dna.seq))
     rC[cg.mtbr$posi] <- cg.mtbr$rC
     HCC.rC <- rC
     rT <- integer(length(dna.seq))
     rT[cg.mtbr$posi] <- cg.mtbr$rT
     HCC.rT <- rT
     HCC <- data.frame(chrom=cg.mtbr$chrom,posi=cg.mtbr$posi,HCC_rC=cg.mtbr$rC,HCC_rT=cg.mtbr$rT)
     HCC[is.na(HCC)] <- 0
     all <- merge(HMEC,HCC)
     all[is.na(all)] <- 0
     addPv <- MultiAddPv(all)
     addPv[addPv < 0] <- 1
     dmr <- data.frame(chrom = as.character(),start = as.integer(),end = as.integer())
     dl <- addPv$pv >0 & addPv$pv <= 0.05
     region.index  <- regionAsBed(dl,10,20,chr)
     region.postion <- data.frame(chrom = as.character(addPv$chrom[region.index$start]),start = as.integer(addPv$posi[region.index$start]),end = as.integer(addPv$posi[region.index$end]))
     dmr <- data.frame(chrom = c(as.character(dmr$chrom),as.character(region.postion$chrom)),start = c(dmr$start,region.postion$start), end   = c(dmr$end,region.postion$end))
     dmr$length <- dmr$end - dmr$start
     cposition <- GetCcontextPosition(dna.seq)
     dmr$cgNum <- apply(dmr[,c("start","end")],1,function(rg) sum(cposition[rg["start"]:rg["end"]]))
     dmr$HMEC_rC <- apply(dmr[, c("start", "end")], 1, function(rg) sum(HMEC.rC[rg["start"]:rg["end"]]))
     dmr$HMEC_rT <- apply(dmr[, c("start", "end")], 1, function(rg) sum(HMEC.rT[rg["start"]:rg["end"]]))
     dmr$HMEC_coverNum <- dmr$HMEC_rC + dmr$HMEC_rT
     dmr$HMEC_score <- dmr$HMEC_rC / dmr$HMEC_coverNum
     dmr$HCC_rC <- apply(dmr[, c("start", "end")], 1, function(rg) sum(HCC.rC[rg["start"]:rg["end"]]))
     dmr$HCC_rT <- apply(dmr[, c("start", "end")], 1, function(rg) sum(HCC.rT[rg["start"]:rg["end"]]))
     dmr$HCC_coverNum <- dmr$HCC_rC + dmr$HCC_rT
     dmr$HCC_score <- dmr$HCC_rC / dmr$HCC_coverNum
     dmr$diff <- dmr$HMEC_score - dmr$HCC_score
     O <- paste(output,"/",chr,sep="")
     write.table(dmr,file=O,quote=F,col.name=T,row.names=F)
									}
