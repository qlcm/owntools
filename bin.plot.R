args <- commandArgs(T)
USAGE <- function(){
	         message("USAGE:\n\tRscript plot.R <regionfile> <valuefile> <binnum> <regionsize> <outpath(filepath)> <tissname>\n")
}
USAGE()

library(ggplot2)
library(reshape2)
library(data.table)
region <- fread(args[1])
value <- fread(args[2])
colnames(value) <- c("chrom","start","end","score")
binnum <- args[3]
size <- args[4]
region$length <- as.integer(((region$end - region$start) / 2) * size)
region$start <- region$start - region$length
region$end <- region$end + region$length
region$length <- NULL
region$id <- c(1:nrow(region))
region <- as.data.table(region)
dt.region <- region[,.(chrom = chrom,posi = (start:end),bin.id = c(rep(1:binnum,each = floor((end-start+1)/binnum)),rep(binnum,((end-start+1) - floor((end -start+1) / binnum) * binnum)))),by = id]
region.addvalue <- data.table()
type <- unique(region$chrom)
random.chrom <- sample(type,1000,replace = T)
random.posi <- sample(c(10000:46709983),1000)
random <- data.frame(chrom = random.chrom,start=random.posi - 2000,end = random.posi + 2000)
random$id <- c(1:nrow(random))
random <- as.data.table(random)
dt.random <- random[,.(chrom = chrom,posi = (start:end),bin.id = c(rep(1:binnum,each = floor((end-start+1)/binnum)),rep(binnum,((end-start+1) - floor((end -start+1) / binnum)*binnum)))),by = id]
random.addvalue <- data.table()
for(t in type){
	region.split <- dt.region[dt.region$chrom == t,]
	region.split <- as.data.table(region.split)
	random.split <- dt.random[dt.random$chrom == t,]
	random.split <- as.data.table(random.split)
	value.split <- value[value$chrom == t,]
	value.split$start <- as.integer(value.split$start + 1)
	value.split$length <- value.split$end - value.split$start
	value.split.eq0 <- value.split[value.split$length == 0,]
	value.split.gt0 <- value.split[value.split$length > 0,]
	value.split.gt0$length <- NULL
	value.split.gt0$id <- c(1:nrow(value.split.gt0))
	value.split.gt0 <- as.data.table(value.split.gt0)
	dt.value.split.gt0 <- value.split.gt0[,.(chrom = chrom,posi=c(start:end),score=score),by="id"]
	dt.value.split.gt0$id <- NULL
	value.split.eq0$end <- NULL
	value.split.eq0$length <- NULL
	colnames(value.split.eq0) <- c("chrom","posi","score")
	value.split <- rbind(value.split.eq0,dt.value.split.gt0)
	value.split <- value.split[value.split >0,]
	dt.value.split <- as.data.table(value.split)
	region.addvalue.split <- merge(region.split,dt.value.split,by=c("chrom","posi"))
	region.addvalue <- rbind(region.addvalue,region.addvalue.split)
	random.addvalue.split <- merge(random.split,dt.value.split,by=c("chrom","posi"))
	random.addvalue <- rbind(random.addvalue,random.addvalue.split)
}
region.result <- region.addvalue[,.(region_binscore = mean(score)),by=bin.id]
random.result <- random.addvalue[,.(random_binscore = mean(score)),by=bin.id]
all <- merge(region.result,random.result,by="bin.id")
all.mt <- melt(all,id="bin.id",variable.name="type", value.name="binscore")
tissuename <- args[6]
outpath1 <- paste(args[5],"/",tissuename,".Rdata",sep="")
save(all.mt,file = outpath1)
plot <- ggplot(all,aes(x=bin.id,y=binscore,colour=type)) + stat_smooth(se=FALSE) + labs(title= tissuename)
outpath2 <- paste(args[5],"/",tissuename,".pdf",sep="")
ggsave(plot,file=outpath2)

