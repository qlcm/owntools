library(RCurl)
library(XML)

download <- function(strURL){
    h <- basicTextGatherer()# �鿴���������ص�ͷ��Ϣ
    txt <- getURL(strURL, headerfunction = h$update,.encoding="gbk") ## �ַ�����ʽ
    htmlParse(txt,asText=T,encoding="gbk")      #ѡ��gbk������ҳ�Ľ���
}

extradress <- function(strURL){
  prefix <- "http://muchong.com/"
  pattern <- "html/[0-9/]+.html"
  links <- getHTMLLinks(strURL)
  needlinks <- gregexpr(pattern,links)
  needlinkslist <- list()
  for (i in which(unlist(needlinks)>0)){
    preadress <- substr(links[i],needlinks[[i]],needlinks[[i]]+attr(needlinks[[i]],'match.length')-1)
    needlinkslist<- c(needlinkslist,list(preadress))
    adresses <- lapply(needlinkslist,function(x)paste(prefix,x,sep=""))
    #adresses <- xpathSApply(adresses, "//a/@href")
  }
  return (adresses)
}

gettopic <- function(doc){
    xmlValue(getNodeSet(doc,'//p')[[2]])
}

greg <- function(pattern,istring){
    gregout <- gregexpr(pattern,istring)
    substr(istring,gregout[[1]],gregout[[1]]+attr(gregout[[1]],'match.length')-1)
}

getinf <- function(topic){
pattern1 <- "��[\u4E00-\u9FA5]+[0-9-]*[\u4E00-\u9FA5]*[��������,;]*[\u4E00-\u9FA5]*[��������,;]*[\u4E00-\u9FA5]*[��������,;]*[\u4E00-\u9FA5]*[��������,;]*[\u4E00-\u9FA5]*(�о���)|(����)"
pattern2 <- "([\u4E00-\u9FA5]*������|[\u4E00-\u9FA5]*�Ŷ�)"  
pattern21 <- "[\u4E00-\u9FA5]*[��������,;]*(����|��ʿ)"
pattern3 <- "[\u4E00-\u9FA5]*[��������,;]*[-A-Za-z0-9_.%]+@[-A-Za-z0-9_.%]+\\.[A-Za-z]+[.A-Za-z]*"
    #ƥ��@163.com�����@abc.edu.cn��������
pattern4 <- "[\u4E00-\u9FA5]+��ʦ"  #ƥ��ĳ��ʦ
pattern5 <- "[\u4E00-\u9FA5]*[��:]*1[3,5,8]{1}[0-9]{1}[0-9]{8}|0[0-9]{2,3}-[0-9]{7,8}(-[0-9]{1,4})?" #ƥ����ϵ�˺ͺ���
pattern6 <- "(��|����)*[\u4E00-\u9FA5]*(���о�|����)Ϊ*[��������,;]*[\u4E00-\u9FA5]*"
pattern7 <- "[\u4E00-\u9FA5]+(��ѧ|ѧԺ|�о�Ժ|�о���)"
pattern8 <-"[-A-Za-z0-9_.%]+@[-A-Za-z0-9_.%]+\\.[A-Za-z]+[.A-Za-z]*" #��ȷƥ������


cate <- greg(pattern1,topic)
proj <- greg(pattern2,topic)
PI <- greg(pattern21,topic)
email <- greg(pattern3,topic)
man <- greg(pattern4,topic)
phone <- greg(pattern5,topic)
direc <- greg(pattern6,topic)
univ <- greg(pattern7,topic)
print(cate)
if (greg("(����|����|ֲ��|ϸ��|ҽѧ|����|ˮ)+",topic) !=""){
    if (man =="" && proj != ""){
        man <- unlist(strsplit(proj,"������")[1])
    }
  
    if (email != ""){
      email <- greg(pattern10,email)
    }
    
    data.frame("���"=cate,"��ѧ"=univ,"����"=proj,"PI"=PI,"��ϵ��"=man,"����"=email,"����"=direc,"�绰"=phone)
}
else{
  return("")
}
}

strURLs="http://muchong.com/html/f430.html"
n=200
dat <- data.frame("URL"="URL","���"="���","��ѧ"="��ѧ","����"="����","PI"="PI","��ϵ��"="��ϵ��","����"="����","����"="����","�绰"="�绰")
strURLs <- c(strURLs,paste(rep("http://muchong.com/html/f430_",n),c(2:n),".html",sep=""))
output <- "a2017.2.18.txt"

for ( strURL in strURLs){
    adresses <- extradress(strURL)
    for (adress in adresses){
      message(adress)
      doc <- download(adress)
      topic <- gettopic(doc)
      inf <- getinf(topic)
      if (inf != ""){
        URL <- data.frame("URL"=adress)
        inf <- cbind(URL,inf)
        dat<- rbind(dat,inf)
      }
    }
}

write.table(dat, file = output, row.names = F, col.names=F,quote = F, sep="\t")  # tab �ָ����ļ�
message("��ɣ�")