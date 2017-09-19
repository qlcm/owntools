args <- commandArgs(T)
USAGE <- function(){
	         message("USAGE:\n\tRscript get_region_methy.R <regionfile> <mtbr> <outpath(filepath and filename)>\n")
}
USAGE()
Get_methy_mtbr <- function(cg.mtbr,region,tissuename){
	 colnames(region) <- c("chrom","start","end","id")
	 region$id <- 1:nrow(region)                                          
  	region <- as.data.table(region)
  	methy <- data.frame(chrom = cg.mtbr$chrom,posi = cg.mtbr$posi,methyvalue = with(cg.mtbr,((rC_n + rC_p)/(rC_n + rC_p + rT_n + rT_p))))
	methy[is.na(methy)] <- 0
   	region <- region[,.(chrom=chrom,posi=c(start:end)),by="id"]
    	addvalue <- merge(region,methy,by=c("chrom","posi"))
    	addvalue <- as.data.table(addvalue)
     	result <- addvalue[,.(value=mean(methyvalue)),by="id"]
      	colnames(result) <- c("id",tissuename)
      	result <- as.data.frame(result)  
      return(result)
}
library(data.table)
region <- read.table(args[1],header=T)
colnames(region) <- c("chrom","start","end")
chrom <- unique(region$chrom)
tissue <- list.files(args[2])
all <- data.frame()
for(chr in chrom){
	region.sp <- region[region$chrom == chr,]
	region.sp$id <- 1:nrow(region.sp)
	value <- data.frame(id = 1:nrow(region.sp))
	for (t in tissue){
		message(chr,t,"  is going  ",date())
		f <- paste(args[2],t,"/","mtbr_cg/",chr,".Rdata",sep="")
		load(f)
		result <- Get_methy_mtbr(cg.mtbr,region.sp,t)
		value <- merge(value,result,by="id",all.x=T)
		value[is.na(value)] <- 0
	}
	region.sp <- cbind(region.sp,value)
	all <- rbind(all,region.sp)
}
write.table(all,file=args[3],quote=F,row.names=F,col.names=T)

