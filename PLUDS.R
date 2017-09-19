pluds <- function(val,binwidth=100, peakwidth=30000,raise.threshold=4,raise.length=5,raise.tolerance=1,fall.threshold=4,fall.length=5,fall.tolerance=1){
  # find peaks in a vector
  #
  # Args:
  #     binwidth:cut the vector into bins
  #     peakwidth: define the max peak width
  #     raise.threshold: value of bin raise 
  #     raise.length: the minimun length of peak rasie region
  #     raise.tolerance: tolerance of the raise region
  # Returns:
  #     a data.frame include chrom, start, end
  # Author:
  #     Jimmy Chen, E-mail: Jimmy_ccj@outlook.com
  
  bin <-val[seq(1,length(val),by=binwidth)]
  bin.val <- diff(bin)
  mark <- bin.val > raise.threshold
  start <- regionAsBed(mark,raise.length,raise.tolerance,chr)
  peak <- data.frame(id=c(1:nrow(start)),
                     chrom=rep(chr,nrow(start)),
                     start=(start$start)*binwidth,
                     end=numeric(nrow(start)),
                     AS=(start$start)*binwidth,
                     AE=start$end*binwidth,
                     DS=numeric(nrow(start)),
                     DE=numeric(nrow(start)),
                     AL=numeric(nrow(start)),
                     DL=numeric(nrow(start)),
                     Rang=numeric(nrow(start)),
                     MinA=numeric(nrow(start)),
                     MaxA=numeric(nrow(start)),
                     MaxD=numeric(nrow(start)),
                     MinD=numeric(nrow(start)))
  
  for(i in 1:nrow(peak)){
    message("Peak id :", peak$id[i] )
    fall.start = peak[i,]$AE
    fall.end   = peak[i,]$AE+peakwidth 
    
    if(fall.end > length(val)){
      fall.end=length(val)
    }
    
    fall.val <-val[fall.start:fall.end]
    bin <- fall.val[seq(1,length(fall.val),by=binwidth)]
    bin.val <- diff(bin)
    mark <- bin.val < -fall.threshold
    end <- regionAsBed(mark,fall.length,fall.tolerance,chr)  
    
    ## if the peak has descending regrion
    if( nrow(end) == 0){
      peak$end[i] <- peak$start[i]   
      message("进度百分比:",round(i/nrow(peak)*100,2),"%")
    }else{
      ## ascending reion information
      x=c(peak[i,]$AS:peak[i,]$AE) 
      ystart=seq(peak[i,]$AS,peak[i,]$AE,by=binwidth)
      yend=seq(peak[i,]$AS+binwidth,peak[i,]$AE,by=binwidth)
      y=data.frame(start=ystart[-length(ystart)],end=yend,val=numeric(length(yend)))
      for( j in 1:nrow(y)){
        y$val[j]=mean(val[y$start[j]:y$end[j]])
      }    
      y=y$val
      x=c(1:length(y))
      lm.sol <- lm(y~1+x)
      
      peak$AL[i] <- as.numeric(coef(lm.sol)[2])
      peak$MaxA[i] <- max(y)
      peak$MinA[i] <- min(y)
      
      ## descending reion information
      peak$DS[i] <- fall.start+end$start[1]*binwidth
      peak$DE[i] <- fall.start+end$end[1]*binwidth
      peak$end[i] <- peak$DE[i]
      
      x=c(peak$DS[i]:peak$DE[i])
      ystart=seq(peak[i,]$DS,peak[i,]$DE,by=binwidth)
      yend=seq(peak[i,]$DS+binwidth,peak[i,]$DE,by=binwidth)
      y=data.frame(start=ystart[-length(ystart)],end=yend,val=numeric(length(yend)))
      
      for( j in 1:nrow(y)){
        y$val[j]=mean(val[y$start[j]:y$end[j]])
      }
      
      y=y$val
      x=c(1:length(y))
      lm.sol <- lm(y~1+x)
      
      peak$DL[i] <- as.numeric(coef(lm.sol)[2])
      peak$MaxD[i] <- max(y)
      peak$MinD[i] <- min(y)
      
      peakVal=val[peak$start[i]:peak$end[i]]
      peak$Rang[i]=max(peakVal)-min(peakVal)
      
      message("进度百分比:",round(i/nrow(peak)*100,2),"%")
    }
  }
  peak <- peak[with(peak,end-start) > 0 ,]
}
