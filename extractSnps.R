ub<-1415*100

snps<-c(3,4,39,1193)
fl<-list.files(pattern="^sop_")

N<-length(fl)

for(k in c(1:7,9)){
    dat<-as.matrix(read.table(fl[k],header=FALSE))
    hd<-unlist(strsplit(fl[k],split=".txt"))
    for(j in 1:4){
        xx<-seq(snps[j],ub,1415)
        sdat<-dat[xx,]
        out<-paste(hd,"_",snps[j],".txt",sep="")
        write.table(sdat,file=out,row.names=FALSE,col.names=FALSE,quote=FALSE)
    }
}        

k<-8
dat<-as.matrix(read.table(fl[k],header=FALSE))
hd<-unlist(strsplit(fl[k],split=".txt"))
for(j in 1:4){
    xx<-seq(snps[j],1415,1415)
    sdat<-dat[xx,]
    out<-paste(hd,"_",snps[j],".txt",sep="")
    write.table(sdat,file=out,row.names=FALSE,col.names=FALSE,quote=FALSE)
}
