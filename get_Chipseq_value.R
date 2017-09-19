library("data.table")
library(methyutils)
library("BSgenome.Hsapiens.UCSC.hg38")
args <- commandArgs(T)
USAGE <- function(){
 message("USAGE:\n\tRscript getvalue.R <tissueregion> <chip-seq(filepath)> <outpath(filepath and filename .Rdata)>\n")
}
USAGE()
region <- read.table(args[1],header=T)
region$id <- 1:nrow(region)
Chipdata <- list.files(args[2])
chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
  sort.region <- data.frame()
  for(c in chrs){
	  test <- region[region$chrom == c,]
	  test <- test[order(test$start,test$end),]
	  #test <- as.data.table(test)
	  #split.test <- test[,.(chrom = chrom,posi = (start:end)),by = id]
	  #cons <- fread(paste("/home/qzzh/major_project/cgdensity/phastCons20way/",c,sep=""))
	  #colnames(cons) <- c("chrom","start","end","score")
	  #cons$start <- as.integer(cons$start + 1)
	  #cons$id <- 1:nrow(cons)
	  #cons <- cons[,.(chrom = chrom,posi=c(start:end),score=score),by="id"]
	  #cons$id <- NULL
	  #cons <- cons[cons$posi %in% split.test$posi,] 
	  #mg <- merge(split.test,cons,by=c("chrom","posi"))
	  #mg <- mg[,.(Cons = mean(score)),by=id]
	  #test<- merge(test,mg,by="id")
	  #test$id <- NULL
	  #test <- as.data.frame(test)
	  #dna.seq <- Hsapiens[[c]]
	  #cposition <- GetCcontextPosition(dna.seq)
	  #test$cgNum <- apply(test[,c("start","end")], 1, function(rg)sum(cposition[rg["start"]:rg["end"]]))
	  sort.region <- rbind(sort.region,test)
  }
  region <- sort.region
	
  for (chip in Chipdata){
       addvalue <- data.frame()
       data <- data.frame()
       I <- paste(args[2],"/",chip,sep="")      
       tmp1 <- fread(I)
       message(chip," is going ",date())
	for (chr in chrs){
	  message(chr," is going ",date())
          tmp2 <- tmp1[tmp1$V1 == chr,]
          colnames(tmp2) <- c("chrom","posi","value")
	  tmp2 <- tmp2[tmp2$value > 0,]
          tmp.region <- region[region$chrom == chr,]
	  tmp.region <- tmp.region[,1:3]
          colnames(tmp.region) <-  c("chrom","start","end")
          tmp.region$marker <- c(1:nrow(tmp.region))
          dt.tmp.region <- as.data.table(tmp.region)
          dt.tmp.region$chrom <- as.character(dt.tmp.region$chrom)
          dt.tmp.region$start <- as.integer(dt.tmp.region$start)
          dt.tmp.region$end <- as.integer(dt.tmp.region$end)
          dt.tmp.region.split <- dt.tmp.region[,.(chrom = chrom, posi = c(start : end)),by = marker]
          dt.value <- merge(dt.tmp.region.split,tmp2,by="posi")
          dt.chip.value <- dt.value[,.(score = mean(value)),by = marker]
          tmp.add.value <- merge(tmp.region,dt.chip.value,by="marker",all.x = TRUE)
	  tmp.add.value[is.na(tmp.add.value)] <- 0
          tmp3 <- data.frame(value=tmp.add.value$score)
          colnames(tmp3) <- c(chip)
	  addvalue <- rbind(addvalue,tmp3)
	  addvalue[is.na(addvalue)] <- 0
	                   }
      region <- cbind(region,addvalue)
                       }
save(region,file=args[3])
