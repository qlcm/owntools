
library(data.table)
regions <- read.table("~/cgDensity/Convolve/Convolve2/test.regions.csv",header = T)
load("~/Rendata/bonemarrow/mtbr_cg/chr19.Rdata")
cg.mtbr$rC <- cg.mtbr$rC_n + cg.mtbr$rC_p
cg.mtbr$rT <- cg.mtbr$rT_n + cg.mtbr$rT_p
cg.mtbr <- cg.mtbr[(cg.mtbr$rC + cg.mtbr$rT) > 5,]
dt.mtbr <- data.table(cg.mtbr)
dt.regions <- data.table(regions)
dt.regions$id = 1:nrow(dt.regions)
dt.regions.split <- dt.regions[,.(chr = chr, pos = start:end, value = value, class = class), by = id]

chr19.wig <- fread("~/cgDensity/Convolve/Convolve2/chr19.test.wig")
dt.regions.split$value <- chr19.wig$V1[dt.regions.split$pos]

colnames(dt.mtbr) <- c("chr", "pos", "rC_n", "rC_p", "rT_n", "rT_p", "rC", "rT")
dt.mtbr$score <- dt.mtbr$rC / (dt.mtbr$rC + dt.mtbr$rT)
dt.region.values <- merge(dt.mtbr, dt.regions.split, by = c("chr", "pos"))


dt.region.sum.values <- dt.region.values[,.( score = sum(rC) / (sum(rC) + sum(rT))), by = c("id", "chr", "class", "value")]

dt.region.sum.values <- dt.region.values[,.( score = mean(score)), by = c("id", "chr", "class", "value")]

dt.region.cor.values <- dt.region.values[,.( score = cor(score, value, method = "spearman")), by = c("id", "chr", "class")]

dt.region.cor.values.H <- dt.region.cor.values[dt.region.cor.values$class == "H",]
dt.region.cor.values.M <- dt.region.cor.values[dt.region.cor.values$class == "M",]                                                                                            
dt.region.cor.values.L <- dt.region.cor.values[dt.region.cor.values$class == "L",]                                                                                            
hist(dt.region.cor.values.H$score)                                                                                                                                            
hist(dt.region.cor.values.M$score)                                                                                                                                            
hist(dt.region.cor.values.L$score)                                                       






dt.region.sum.values.order <- dt.region.sum.values[order(dt.region.sum.values$id),]
cor(dt.region.sum.values.order$value, dt.region.sum.values.order$score)
dt.region.sum.values.order.H <- dt.region.sum.values.order[dt.region.sum.values.order$class == "H",]                                                                          
cor(dt.region.sum.values.order.H$value, dt.region.sum.values.order.H$score)                                                                                                   
dt.region.sum.values.order.L <- dt.region.sum.values.order[dt.region.sum.values.order$class == "L",]                                                                          
cor(dt.region.sum.values.order.L$value, dt.region.sum.values.order.L$score)                                                                                                   
dt.region.sum.values.order.M <- dt.region.sum.values.order[dt.region.sum.values.order$class == "M",]                                                                          
cor(dt.region.sum.values.order.M$value, dt.region.sum.values.order.M$score)   

plot(dt.region.sum.values.order.H$value,dt.region.sum.values.order.H$score) 

### rescale the CpG Density 

dt.region.values$rescalevalue <- (dt.region.values$value - min(dt.region.values$value)) / (max(dt.region.values$value) - min(dt.region.values$value))
dt.region.values$res <- res
cmpl <- dt.region.values[,.(gap = (sum(res < -0.1) / length(res)), cmpl = (sum(res > -0.1 & res < 0.1) / length(res)), overlap = (sum(res > 0.1) / length(res)), density = median(rescalevalue)), by = c("id", "class")]

meth <- fread("~/cgDensity/Convolve/Convolve2/chr19.meth.some.wig")
meth$pos <- 1:nrow(meth)
meth$chr <- "chr19"
colnames(meth) <- c("methvalue", "pos", "chr")
dt.region.meth <- merge(dt.regions.split, meth, by = c("pos", "chr"))
dt.region.meth.H <- dt.region.meth[dt.region.meth$class == "H"]
dt.region.meth.M <- dt.region.meth[dt.region.meth$class == "M"]
dt.region.meth.L <- dt.region.meth[dt.region.meth$class == "L"]

dt.region.meth.H.value <- dt.region.meth.H[,.(density = max(value), methy = max(methvalue)), by = c("chr", "id", "class")]
dt.region.meth.M.value <- dt.region.meth.M[,.(density = max(value), methy = max(methvalue)), by = c("chr", "id", "class")]
dt.region.meth.L.value <- dt.region.meth.L[,.(density = min(value), methy = min(methvalue)), by = c("chr", "id", "class")]



kde <- read.table("/home/fsch/cgDensity/Convolve/Convolveroots/test2/test/win655/testtime/win2500/chr19.kde")  


library(ggplot2)
colnames(kde) <- c("pos", "cg", "nocg", "valley")
library(reshape2)

kde.melt <- melt(kde,id = "pos")                                              

par()
ggplot(kde.melt,aes(x = pos, y = value,color = variable)) + geom_line() + xlim(0, 0.2)    
