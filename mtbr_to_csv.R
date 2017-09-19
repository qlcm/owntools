USAGE <- function(){
         message("USAGE:\n\tRscript mtbr_to_csv.R <mtbrfilePath>  <outfilePath>\n")
		message("	The filePath must be include / at the end")
}
USAGE()
library(data.table)
args <-commandArgs(T)
MtbrPath <- args[1]
OutPath <- args[2]
if(!file.exists(MtbrPath)){
  stop("mtbr path \"", MtbrPath ,"\" does not exist.")
}

if(!file.exists(OutPath)){
  dir.create(OutPath)
}

tissues <- list.files(MtbrPath)
for (tissue in tissues){
	message(tissue,"  is start......")
	header <- paste("chrom","posi","rC_n","rC_p","rT_n","rT_p",sep=",")
	write.table(header,file=paste(OutPath,tissue,".mtbr.csv",sep=""),quote=F,append = T,row.name=F,col.names=F)
	splitfile <- list.files(paste(MtbrPath,tissue,"/",tissue,".mtbr.cg",sep=""))
	for (split in splitfile){
		message(split,"  is start.....")
		load(paste(MtbrPath,tissue,"/",tissue,".mtbr.cg/",split,sep=""))
		fwrite(cg.mtbr,file=paste(OutPath,tissue,".mtbr.csv",sep=""),quote=F,append=T,row.name=F,col.names=F,sep=",")
							}
						}



