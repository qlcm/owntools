hG <- data.table(df.Gap.CombineRegion.sel)
ld <- hG[,grep("_|id", colnames(hG), value=T), with = FALSE]
m <- melt(ld, id.vars="id")

cs <- colsplit(string=as.character(m$variable), pattern="_", names=c("tissue","valueType"))
mcs <- cbind(m,cs)
dc <- dcast(data=mcs, formula=id + tissue ~ valueType, value.var="value")
####left data
lfd <- hG[,grep("_", colnames(hG), value=T, invert=T), with = FALSE]
####long format, which was got after melting and casting
lmc <- merge(lfd, dc, by="id")
lmc <- data.table(lmc)

