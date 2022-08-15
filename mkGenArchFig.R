library(data.table)
library(scales)
gu<-read.table("genarch_unc.txt",header=TRUE)

gu[gu[,1] > 0.5,1]<-1-gu[gu[,1] > 0.5,1]

par(mfrow=c(3,1))

sop_drift<-as.matrix(fread("sop_gunc_drift.txt",header=FALSE))

prec_drift<-rep(NA,1415)
for(i in 1:1415){
    prec_drift[i]<-median(1/apply(sop_drift[seq(i,141500,1415),],2,var),na.rm=FALSE)
    }
    

prec_drift[is.infinite(prec_drift)]<-NA


sop_sel<-as.matrix(fread("sop_gunc_weather.txt",header=FALSE))

prec_sel<-rep(NA,1415)
for(i in 1:1415){
    prec_sel[i]<-median(1/apply(sop_sel[seq(i,141500,1415),],2,var),na.rm=FALSE)
    }
    

prec_sel[is.infinite(prec_sel)]<-NA



sop_drsel<-as.matrix(fread("sop_gunc_dr_weather.txt",header=FALSE))

prec_drsel<-rep(NA,1415)
for(i in 1:1415){
    prec_drsel[i]<-median(1/apply(sop_drsel[seq(i,141500,1415),],2,var),na.rm=FALSE)
    }
    

prec_drsel[is.infinite(prec_drsel)]<-NA

cs<-alpha("darkgray",.5)
cl<-1.7;ca<-1.15;cm<-1.6

pdf("sfigGenUncertain.pdf",width=5,height=12)
par(mfrow=c(3,1))
par(mar=c(4.5,5.5,2.5,1))

plot(gu[,1],log10(prec_drift),xlab="Uncertainty",ylab="Precision (log10)",cex.lab=cl,cex.axis=ca,pch=19, col=cs)
title(main="(A) Drift + unc. genetics",cex.main=cm)
o<-cor.test(gu[,1],log10(prec_drift))
mtext(side=3,paste("r = ",round(o$estimate,2),sep=""),line=-1.5,adj=.7)
mtext(side=3,paste("(",round(o$conf.int[1],2),",",round(o$conf.int[2],2),")",sep=""),line=-2.75,adj=.7)

plot(gu[,1],log10(prec_sel),xlab="Uncertainty",ylab="Precision (log10)",cex.lab=cl,cex.axis=ca,pch=19, col=cs)
title(main="(B) Unc. sel. + weather + gen.",cex.main=cm)
o<-cor.test(gu[,1],log10(prec_sel))
mtext(side=3,paste("r = ",round(o$estimate,2),sep=""),line=-1.5,adj=.7)
mtext(side=3,paste("(",round(o$conf.int[1],2),",",round(o$conf.int[2],2),")",sep=""),line=-2.75,adj=.7)

plot(gu[,1],log10(prec_drsel),xlab="Uncertainty",ylab="Precision (log10)",cex.lab=cl,cex.axis=ca,pch=19, col=cs)
title(main="(C) Drift + unc. sel. + w. + gen",cex.main=cm)
o<-cor.test(gu[,1],log10(prec_drsel))
mtext(side=3,paste("r = ",round(o$estimate,2),sep=""),line=-1.5,adj=.7)
mtext(side=3,paste("(",round(o$conf.int[1],2),",",round(o$conf.int[2],2),")",sep=""),line=-2.75,adj=.7)

dev.off()

