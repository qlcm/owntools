RegionAsBed <- function(marker, cf.length=5, tolerance.length=NULL,to_up_len = 2,chrom) {
	#library(zoo)
	r <- rle(marker)
	if(is.null(tolerance.length)){
		end <- with(r,cumsum(lengths)[values & lengths>=cf.length])
		start <- end - with(r,lengths[values & lengths>=cf.length]) + 1
		if(length(start) == 0 || length(end) == 0){
			df <- data.frame(chrom=character(), start=integer(), end=integer())
		} else {
			df <- data.frame(chrom, start, end)
		}
	} else {
		end<-cumsum(r$lengths)
		start<- end - r$lengths + 1
		df.tmp <- data.frame(start = start  , end = end , value = r$values)
		df.tmp$len <- df.tmp$end - df.tmp$start + 1 
		df.tolerance <- df.tmp[!df.tmp$value & df.tmp$len <= tolerance.length,] 
		tolerance_rowname <- as.numeric(rownames(df.tolerance))
		tolerance_data <- data.frame(tolerance_rowname = tolerance_rowname)
		tolerance_data$tolerance_up <- tolerance_rowname - 1 
		tolerance_data$tolerance_down <- tolerance_rowname + 1 
		tolerance_data_sel <- tolerance_data[with(tolerance_data, tolerance_up > 0 & tolerance_down <= nrow(df.tmp)),]		
		tolerance_mark <- tolerance_data_sel$tolerance_rowname[df.tmp$len[tolerance_data_sel$tolerance_up] >= to_up_len & df.tmp$len[tolerance_data_sel$tolerance_down] >= to_up_len]


		dt.to <- data.table(df.tmp[tolerance_mark,])
		if (nrow(dt.to)!=0) {
			dt.to$id <- 1:nrow(dt.to)
			dt.to.index <- dt.to[,.(index = start:end), by = id] 		
			marker[dt.to.index$index] <- TRUE
		}
		
		r <- rle(marker)
        	end <- with(r,cumsum(lengths)[values & lengths>=cf.length])
	        start <- end - with(r,lengths[values & lengths>=cf.length]) + 1
        	if(length(start) == 0 || length(end) == 0){
            		df <- data.frame(chrom=character(), start=integer(), end=integer())
        	} else {
            		df <- data.frame(chrom, start, end)
        	}
    	}	
    return(df)
}
library(data.table)

#args <- commandArgs(T)

#flist <- list.files()

chr<-paste(rep("chr",23),c(1:22,"X"),sep="")
for (i in chr)
{
	outname <- paste(i,".hotspot",sep="")
#	a <- read.table(args[1],header=T)
	a <- read.table(paste("/home/mytong/Project/231WGS/Varscan_zhangyu/Varscan_real.bam/SNVINFOFromVCF.231WGS/hotspot_2D_SP4/hotspot_2D/chr_dis/",i,".dis",sep=""),header=T)
	a$marker <- a$Pvalue < 0.01
	m <- a$marker
#	hotsport <- RegionAsBed(m,5,1,2,args[2])
        hotsport<-RegionAsBed(m,6,1,2,i)
	hotsport$start <- a$Pos[hotsport$start]
	hotsport$end <- a$Pos[hotsport$end]
	write.table(hotsport, file=outname,row.names=F, quote=F, sep="\t")
}

