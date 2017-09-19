dataname<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
filenames<-list.files("/home/qzzh/J.R.Ecker/mtbr")
file_path1="/home/qzzh/J.R.Ecker/mtbr/"
file_path2="/home/qzzh/J.R.Ecker/methydata/"
for (samplename in filenames){
message(samplename," is running",date())
for (z in dataname){
O <- paste(file_path1,samplename,"/",samplename,".mtbr.cg","/",z,".Rdata",sep="")
load(O)
output=data.frame() 
coverage <- cg.mtbr$rC_p + cg.mtbr$rC_n + cg.mtbr$rT_p + cg.mtbr$rT_n
rC <- cg.mtbr$rC_p + cg.mtbr$rC_n
chrom <- cg.mtbr$chrom
posi <- cg.mtbr$posi
score <- rC / coverage
output <- data.frame(chrom,posi,score)
output[is.na(output)] <- 0
Q<-paste(file_path2,samplename,".txt",sep="")
write.table(output,file=Q,append=T,quote=F,row.names=F,col.names=F)
}
message(samplename," is finished",date())
}


