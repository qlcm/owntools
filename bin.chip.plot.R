args <- commandArgs(T)
USAGE <- function(){
         message("USAGE:\n\tRscript bin.chip.bin.R <regionfile(chrom,start,end)> <valuefilepath> <promoter_regionfile> <enhancer_regionfile> <binnum> <regionsize> <outpath(filepath)> <tissname>\n")
}
USAGE()
#library(BSgenome.Hsapiens.UCSC.hg38)
#library(methyutils)
library(ggplot2)
library(reshape2)
library(data.table)
############# read data file ##################

region <- fread(args[1],header=T)
value <- fread(args[2])
promoter <- fread(args[3],header=T)
sample <- sample(c(1:nrow(promoter)),1000)
promoter <- promoter[sample,]
enhancer <- fread(args[4],header=T)
sample <- sample(c(1:nrow(enhancer)),1000)
enhancer <- enhancer[sample,]
colnames(value) <- c("chrom","posi","score")
binnum <- as.numeric(args[5])
size <- as.numeric(args[6])

############## make your region to posi split ###############

region$length <- as.integer((as.integer(region$end) - as.integer(region$start)) / 2)
region$length <- region$length * size
region$start <- region$start - region$length
region$end <- region$end + region$length
region$length <- NULL
region$id <- c(1:nrow(region))
region <- as.data.table(region)
dt.region <- region[,.(chrom = chrom,posi = (start:end),bin.id = c(rep(1:binnum,each = floor((end-start+1)/binnum)),rep(binnum,((end-start+1) - floor((end -start+1) / binnum) * binnum)))),by = id]
region.addvalue <- data.frame()

############## make promoter region to posi split ###############

promoter$id <- c(1:nrow(promoter))
promoter$length <- as.integer((as.integer(promoter$end) - as.integer(promoter$start)) / 2)
promoter$length <- promoter$length * size
promoter$start <- as.integer(promoter$start) - as.integer(promoter$length)
promoter$end <- as.integer(promoter$end) + as.integer(promoter$length)
promoter$length <- NULL
promoter <- as.data.table(promoter)
dt.promoter <- promoter[,.(chrom = chrom,posi = (start:end),bin.id = c(rep(1:binnum,each = floor((end-start+1)/binnum)),rep(binnum,((end-start+1) - floor((end -start+1) / binnum) * binnum)))),by = id]
promoter.addvalue <- data.frame()

############## make enhancer region to posi split ##############

enhancer$id <- c(1:nrow(enhancer))
enhancer$length <- as.integer((as.integer(enhancer$end) - as.integer(enhancer$start)) / 2)
enhancer$length <- enhancer$length * size
enhancer$start <- as.integer(enhancer$start) - as.integer(enhancer$length)
enhancer$end <- as.integer(enhancer$end) + as.integer(enhancer$length)
enhancer$length <- NULL
enhancer <- as.data.table(enhancer)
dt.enhancer <- enhancer[,.(chrom = chrom,posi = (start:end),bin.id = c(rep(1:binnum,each = floor((end-start+1)/binnum)),rep(binnum,((end-start+1) - floor((end -start+1) / binnum) * binnum)))),by = id]
enhancer.addvalue <- data.frame()

############## make random region to posi split #################

type <- unique(region$chrom)
random.chrom <- sample(type,1000,replace = T)
random.posi <- sample(c(1:46709983),1000)
random <- data.frame(chrom = random.chrom,start=random.posi - 1000,end = random.posi + 1000)
random$id <- c(1:nrow(random))
random <- as.data.table(random)
dt.random <- random[,.(chrom = chrom,posi = (start:end),bin.id = c(rep(1:binnum,each = floor((end-start+1)/binnum)),rep(binnum,((end-start+1) - floor((end -start+1) / binnum)*binnum)))),by = id]
random.addvalue <- data.frame()

############## addvalue split by chrom ################

for(t in type){
########split by chrom #############

region.split <- dt.region[dt.region$chrom == t,]
random.split <- dt.random[dt.random$chrom == t,]
promoter.split <- dt.promoter[dt.promoter$chrom == t,]
enhancer.split <- dt.enhancer[dt.enhancer$chrom == t,]
value.split <- value[value$chrom == t,]
#dna.seq <- Hsapiens[[t]]
#vc <- rep(0,length(dna.seq))
#vc[value.split$posi] <- value.split$score
#win <- list(L = 500, R = 500)
#sw.value <- swsCalc(vc, win)
#value.split <- data.frame(chrom=rep(t,length(sw.value)),posi=c(1:length(sw.value)),score=sw.value)
value.split <- value.split[value.split$score >0,]

######## your region addvalue ###########

value.split.region.sel <- value.split[value.split$posi %in% region.split$posi,]
region.addvalue.split <- merge(region.split,value.split.region.sel,by=c("chrom","posi"),all.x = T)
region.addvalue.split[is.na(region.addvalue.split)] <- 0
region.addvalue <- rbind(region.addvalue,region.addvalue.split)

######## promoter region addvalue ######## 

value.split.promoter.sel <- value.split[value.split$posi %in% promoter.split$posi,]
promoter.addvalue.split <- merge(promoter.split,value.split.promoter.sel,by=c("chrom","posi"),all.x = T)
promoter.addvalue.split[is.na(promoter.addvalue.split)] <- 0
promoter.addvalue <- rbind(promoter.addvalue,promoter.addvalue.split)

####### enhancer region addvalue ########

value.split.enhancer.sel <- value.split[value.split$posi %in% enhancer.split$posi,]
enhancer.addvalue.split <- merge(enhancer.split,value.split.enhancer.sel,by=c("chrom","posi"),all.x = T)
enhancer.addvalue.split[is.na(enhancer.addvalue.split)] <- 0
enhancer.addvalue <- rbind(enhancer.addvalue,enhancer.addvalue.split)

####### random region addvalue ########
value.split.random.sel <- value.split[value.split$posi %in% random.split$posi,]
random.addvalue.split <- merge(random.split,value.split.random.sel,by=c("chrom","posi"),all.x =T)
random.addvalue.split[is.na(random.addvalue.split)] <- 0
random.addvalue <- rbind(random.addvalue,random.addvalue.split)
}
tissuename <- args[8]
region.addvalue <- as.data.table(region.addvalue)
outpath0 <- paste(args[7],"/",tissuename,"_addvalue.Rdata",sep="")
save(region.addvalue,file=outpath0)
promoter.addvalue <- as.data.table(promoter.addvalue)
enhancer.addvalue <- as.data.table(enhancer.addvalue)
random.addvalue <- as.data.table(random.addvalue)

################ compute mean score by bin.id ##################

region.result <- region.addvalue[,.(region = mean(score)),by=bin.id]
promoter.result <- promoter.addvalue[,.(promoter = mean(score)),by=bin.id]
enhancer.result <- enhancer.addvalue[,.(enhancer = mean(score)),by=bin.id]
random.result <- random.addvalue[,.(random = mean(score)/1.5),by=bin.id]

################ merge and melt data to plot ##############

all <- merge(region.result,promoter.result,by="bin.id")
all <- merge(all,enhancer.result,by="bin.id")
all <- merge(all,random.result,by="bin.id")
all.mt <- melt(all,id="bin.id",variable.name="type", value.name="binscore")
if(size == 0){
all.mt$bin.id <- rep(seq((-1 + (2/binnum)),1, (2/binnum)),4)
}else{
all.mt$bin.id <- rep(seq((-size*2 + (size*4/binnum)),(size*2), (size*4/binnum)),4)
}
outpath1 <- paste(args[7],"/",tissuename,"_plot",".Rdata",sep="")
save(all.mt,file = outpath1)

################# line plot ###################################

bin_line_plot <- ggplot(all.mt,aes(x=as.numeric(bin.id),y=binscore,colour=type)) + stat_smooth(se=FALSE) + labs(title= tissuename) + xlab("Normlized region size") + ylab("value") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black"))
outpath2 <- paste(args[7],"/",tissuename,"_bin",".pdf",sep="")
ggsave(bin_line_plot,file=outpath2)

#################heatmap plot ################################
#heatmap.result <- region.addvalue[,.(region_binscore = mean(score)),by=c("id","bin.id")]
#heatmap.result <- unique(heatmap.result)
#region$name <- paste(region$chrom,region$start,region$end,sep=":")
#region$chrom <- NULL
#region$start <- NULL
#region$end <- NULL
#heatmap.data <- merge(region,heatmap.result,by="id")
#outpath3 <- paste(args[7],"/",tissuename,"_heatmap",".Rdata",sep="")
#save(heatmap.data,file=outpath3)
#theme_none <- theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
 #                   panel.background = element_blank(),axis.title.x = element_text(colour=NA),
  #                  axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
   #                 axis.line = element_blank())
#heatmap.plot <- ggplot(heatmap.data, aes(x=bin.id, y=name)) + geom_tile(aes(fill=region_binscore), alpha=1) + scale_fill_gradient2(low="green",mid="blue",high="red") +
 #               theme(axis.text.x=element_text(size=10,angle=0),axis.text.y=element_text(size=10)) +
  #              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_blank(),axis.ticks= element_blank())
#outpath4 <- paste(args[7],"/",tissuename,"_heatmap",".pdf",sep="")
#ggsave(heatmap.plot,file=outpath4)

