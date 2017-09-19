args <- commandArgs(T)
USAGE <- function(){ 
	 message("USAGE:\n\tRscript SNP_marker.R <snpfile>  <outpath(filepath and filename)>\n")
}
USAGE()
library(data.table)
dataname<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrM","chrX","chrY")
raw <- fread(args[1],header=F)
raw <- data.frame(raw)
colnames(raw) <- c("chrom","posi","value")
for (chr in dataname){
wig <- raw[raw$chrom == chr,]	
wig$chrom <- NULL
name <- paste("variableStep chrom=",chr," span=1",sep="")
write.table(name,file=args[2],append=T,quote=F,row.names=F,col.names=F)
write.table(wig,file=args[2],append=T,quote=F,row.names=F,col.names=F)
}






