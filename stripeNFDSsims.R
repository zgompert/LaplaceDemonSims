ffst<-function(brk=.5,wl=NA,wh=NA,f=NA){
    if(f < brk){w<-wl}
    else{w<-wh}
    return(w)
 }

sdif<-function(w1=NA,w2=NA,f1=NA){
    wbar<-(w1*f1 + w2*(1-f1))
    psf1<-f1 * (w1/wbar)
    psf2<-(1-f1) * (w2/wbar)
    return(c(psf1,psf2))
    }

sdifh2<-function(w1=NA,w2=NA,f1=NA,h2=NA){
    wbar<-(w1*f1 + w2*(1-f1))
    psf1<-f1 * (w1/wbar)
    psf2<-(1-f1) * (w2/wbar)
    Sd<-psf1-f1
    psf1<-(Sd*h2)+f1
    psf2<-1-psf1
    return(c(psf1,psf2))
    }


## simulate uncertainty in function
## wl<-c(.52,.155)
## wh<-c(.155+.05,.52-.05)
## brk<-.85
brk<-runif(100,.7,.9) 
wl<-matrix(NA,nrow=100,ncol=2)
wl[,1]<-rbeta(100,52,48)
wl[,2]<-rbeta(100,62,338)

wh<-matrix(NA,nrow=100,ncol=2)
wh[,1]<-rbeta(100,31+10,169-10)
wh[,2]<-rbeta(100,26-2.5,24+2.5)

Ne<-110 ## from science paper, change in FHA

#586 (.446–.754) ## from Aarons paper
# qnorm(p=c(.5,.025,.975),mean=0.586,sd=.08)
#h2<-rbeta(100,16,4)## true .586, using .8 instead, 
## sd = 0.0873
h2<-rep(.8,100)
vars<-matrix(NA,nrow=5,ncol=100)
library(scales)

## hs = .8
pdf("nfdsSimsh2_80.pdf",width=12,height=8)
cm<-1.5
cl<-1.6
par(mfrow=c(2,3))
par(mar=c(5.5,4.5,2.5,1.5))


### with heritability, true 0.586, using .8 instead
## drift only
p<-rep(.5,100)  
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(0.85,c(0.52,0.155),c(.155+.05,.52-.05),p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
    
}    
vars[1,]<-apply(p,2,var)     
mtext(round(median(1/vars[1,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(A) Drift",cex.main=cm)

## uncertainty in s
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    p[j,i]<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}  
vars[2,]<-apply(p,2,var)    

mtext(round(median(1/vars[2,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(B) Unc. in selection",cex.main=cm)

## uncertainty in s and drift
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}     
vars[3,]<-apply(p,2,var)     
mtext(round(median(1/vars[3,]),1),side=1,adj=.75,line=-2.5,cex=1.5)

title(main="(C) Unc. in selection + drift",cex.main=cm)

brk<-runif(100,.8,.9) 
wl<-matrix(NA,nrow=100,ncol=2)
wl[,1]<-rbeta(100,52*2,48*2)
wl[,2]<-rbeta(100,62*2,338*2)

wh<-matrix(NA,nrow=100,ncol=2)
wh[,1]<-rbeta(100,(31+10)*2,(169-10)*2)
wh[,2]<-rbeta(100,(26-2.5)*2,(24+2.5)*2)

### with heritability , true 0.586, using .8 instead
## uncertainty in s
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    p[j,i]<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}  
vars[4,]<-apply(p,2,var)    
     
mtext(round(median(1/vars[4,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(D) Unc. in selection, 2X",cex.main=cm)

## uncertainty in s and drift
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}     
vars[5,]<-apply(p,2,var)     
mtext(round(median(1/vars[5,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
   
title(main="(E) Unc. in selection + drift, 2X",cex.main=cm)

e_vars<-apply(1/vars,1,quantile,probs=c(.5,.25,.75),na.rm=TRUE)
yub<-max(as.vector(e_vars))*1.05
x<-barplot(e_vars[1,],names.arg=LETTERS[1:5],xlab="Conditions",ylab="Precision in p",cex.lab=cl,ylim=c(0,yub))
segments(x,e_vars[2,],x,e_vars[3,])
title(main="(F) Predictability",cex.main=cm)

dev.off()

e_vars
#        [,1]     [,2]     [,3]     [,4]     [,5]
#50% 272.4602 121.9296 110.8694 242.9725 216.4700
#25% 247.6738 111.4042 101.1779 226.4005 197.1216
#75% 294.8956 136.9264 130.2634 261.1967 237.0451

sqrt(1/e_vars)
#          [,1]       [,2]       [,3]       [,4]       [,5]
#50% 0.06058268 0.09056188 0.09497170 0.06415366 0.06796747
#25% 0.06354187 0.09474345 0.09941622 0.06646015 0.07122508
#75% 0.05823255 0.08545872 0.08761709 0.06187513 0.06495080




##### heritability 1
## simulate uncertainty in function
## wl<-c(.52,.155)
## wh<-c(.155+.05,.52-.05)
## brk<-.85
brk<-runif(100,.7,.9) 
wl<-matrix(NA,nrow=100,ncol=2)
wl[,1]<-rbeta(100,52,48)
wl[,2]<-rbeta(100,62,338)

wh<-matrix(NA,nrow=100,ncol=2)
wh[,1]<-rbeta(100,31+10,169-10)
wh[,2]<-rbeta(100,26-2.5,24+2.5)

Ne<-110 ## from science paper, change in FHA

#586 (.446–.754) ## from Aarons paper
# qnorm(p=c(.5,.025,.975),mean=0.586,sd=.08)
h2<-rep(1,100)
vars<-matrix(NA,nrow=5,ncol=100)
library(scales)

pdf("nfdsSimsh2_100.pdf",width=12,height=8)
cm<-1.5
cl<-1.6
par(mfrow=c(2,3))
par(mar=c(5.5,4.5,2.5,1.5))

## drift only
p<-rep(.5,100)  
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(0.85,c(0.52,0.155),c(.155+.05,.52-.05),p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
    
}    
vars[1,]<-apply(p,2,var)     
mtext(round(median(1/vars[1,]),1),side=1,adj=.75,line=-2.5,cex=1.5)

title(main="(A) Drift",cex.main=cm)

## uncertainty in s
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    p[j,i]<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}  
vars[2,]<-apply(p,2,var)    

mtext(round(median(1/vars[2,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(B) Unc. in selection",cex.main=cm)

## uncertainty in s and drift
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}     
vars[3,]<-apply(p,2,var)     

mtext(round(median(1/vars[3,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(C) Unc. in selection + drift",cex.main=cm)

brk<-runif(100,.8,.9) 
wl<-matrix(NA,nrow=100,ncol=2)
wl[,1]<-rbeta(100,52*2,48*2)
wl[,2]<-rbeta(100,62*2,338*2)

wh<-matrix(NA,nrow=100,ncol=2)
wh[,1]<-rbeta(100,(31+10)*2,(169-10)*2)
wh[,2]<-rbeta(100,(26-2.5)*2,(24+2.5)*2)

### with heritability , true 0.586, using .8 instead
## uncertainty in s
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    p[j,i]<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}  
vars[4,]<-apply(p,2,var)    
     
mtext(round(median(1/vars[4,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(D) Unc. in selection, 2X",cex.main=cm)

## uncertainty in s and drift
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}     
vars[5,]<-apply(p,2,var)     
   
mtext(round(median(1/vars[5,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(E) Unc. in selection + drift, 2X",cex.main=cm)

e_vars<-apply(1/vars,1,quantile,probs=c(.5,.25,.75),na.rm=TRUE)
yub<-max(as.vector(e_vars))*1.05
x<-barplot(e_vars[1,],names.arg=LETTERS[1:5],xlab="Conditions",ylab="Precision in p",cex.lab=cl,ylim=c(0,yub))
segments(x,e_vars[2,],x,e_vars[3,])
title(main="(F) Predictability",cex.main=cm)

dev.off()

e_vars
#        [,1]      [,2]     [,3]     [,4]     [,5]
#50% 190.5820  93.65362 86.86530 173.9576 155.4733
#25% 176.3233  85.48631 81.15326 160.6889 141.8437
#75% 202.7426 103.23617 96.12022 191.7715 171.0895

sqrt(1/e_vars)
#          [,1]       [,2]      [,3]       [,4]       [,5]
#50% 0.07243677 0.10333269 0.1072943 0.07581905 0.08019959
#25% 0.07530870 0.10815627 0.1110061 0.07888731 0.08396436
#75% 0.07023078 0.09842016 0.1019982 0.07221177 0.07645190

##############  weak selection


##### heritability 1
## simulate uncertainty in function
## wl<-c(.52,.155)
## wh<-c(.155+.05,.52-.05)
## brk<-.85
brk<-runif(100,.7,.9) 
wl<-matrix(NA,nrow=100,ncol=2)
wl[,1]<-rbeta(100,52,48)## .52
wl[,2]<-rbeta(100,192,208)## .48

wh<-matrix(NA,nrow=100,ncol=2)
wh[,1]<-rbeta(100,97,104)
wh[,2]<-rbeta(100,28,27)

Ne<-110 ## from science paper, change in FHA

h2<-rep(1,100)
vars<-matrix(NA,nrow=5,ncol=100)
library(scales)

pdf("nfdsSimsh2_100_weak.pdf",width=12,height=8)
cm<-1.5
cl<-1.6
par(mfrow=c(2,3))
par(mar=c(5.5,4.5,2.5,1.5))

## drift only
p<-rep(.5,100)  
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(0.85,c(0.52,0.48),c(.483,.509),p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
    
}    
vars[1,]<-apply(p,2,var)     

mtext(round(median(1/vars[1,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(A) Drift",cex.main=cm)

## uncertainty in s
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    p[j,i]<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}  
vars[2,]<-apply(p,2,var)    

mtext(round(median(1/vars[2,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(B) Unc. in selection",cex.main=cm)

## uncertainty in s and drift
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}     
vars[3,]<-apply(p,2,var)     

mtext(round(median(1/vars[3,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(C) Unc. in selection + drift",cex.main=cm)

brk<-runif(100,.8,.9) 
wl<-matrix(NA,nrow=100,ncol=2)
wl[,1]<-rbeta(100,52*2,48*2)## .52
wl[,2]<-rbeta(100,192*2,208*2)## .48

wh<-matrix(NA,nrow=100,ncol=2)
wh[,1]<-rbeta(100,97*2,104*2)
wh[,2]<-rbeta(100,28*2,27*2)

### with heritability , true 0.586, using .8 instead
## uncertainty in s
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    p[j,i]<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}  
vars[4,]<-apply(p,2,var)    
     
mtext(round(median(1/vars[4,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(D) Unc. in selection, 2X",cex.main=cm)

## uncertainty in s and drift
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}     
vars[5,]<-apply(p,2,var)     
   
mtext(round(median(1/vars[5,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(E) Unc. in selection + drift, 2X",cex.main=cm)

e_vars<-apply(1/vars,1,quantile,probs=c(.5,.25,.75),na.rm=TRUE)
yub<-max(as.vector(e_vars))*1.05
x<-barplot(e_vars[1,],names.arg=LETTERS[1:5],xlab="Conditions",ylab="Precision in p",cex.lab=cl,ylim=c(0,yub))
segments(x,e_vars[2,],x,e_vars[3,])
title(main="(F) Predictability",cex.main=cm)

dev.off()

##############  weak selection, 1%


##### heritability 1
brk<-runif(100,.7,.9) 
wl<-matrix(NA,nrow=100,ncol=2)
wl[,1]<-rbeta(100,50.5,49.5)## 50.5
wl[,2]<-rbeta(100,198,202)## 49.5

wh<-matrix(NA,nrow=100,ncol=2)
wh[,1]<-rbeta(100,99,101)
wh[,2]<-rbeta(100,25.25,24.75)

Ne<-110 ## from science paper, change in FHA

h2<-rep(1,100)
vars<-matrix(NA,nrow=5,ncol=100)
library(scales)

pdf("nfdsSimsh2_100_very_weak.pdf",width=12,height=8)
cm<-1.5
cl<-1.6
par(mfrow=c(2,3))
par(mar=c(5.5,4.5,2.5,1.5))

## drift only
p<-rep(.5,100)  
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(0.85,c(0.505,0.495),c(.495,.505),p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
    
}    
vars[1,]<-apply(p,2,var)     

mtext(round(median(1/vars[1,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(A) Drift",cex.main=cm)

## uncertainty in s
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    p[j,i]<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}  
vars[2,]<-apply(p,2,var)    

mtext(round(median(1/vars[2,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(B) Unc. in selection",cex.main=cm)

## uncertainty in s and drift
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}     
vars[3,]<-apply(p,2,var)     

mtext(round(median(1/vars[3,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(C) Unc. in selection + drift",cex.main=cm)

brk<-runif(100,.8,.9) 
wl<-matrix(NA,nrow=100,ncol=2)
wl[,1]<-rbeta(100,52*2,48*2)## .52
wl[,2]<-rbeta(100,192*2,208*2)## .48

wh<-matrix(NA,nrow=100,ncol=2)
wh[,1]<-rbeta(100,97*2,104*2)
wh[,2]<-rbeta(100,28*2,27*2)

### with heritability , true 0.586, using .8 instead
## uncertainty in s
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    p[j,i]<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}  
vars[4,]<-apply(p,2,var)    
     
mtext(round(median(1/vars[4,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(D) Unc. in selection, 2X",cex.main=cm)

## uncertainty in s and drift
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}     
vars[5,]<-apply(p,2,var)     
   
mtext(round(median(1/vars[5,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(E) Unc. in selection + drift, 2X",cex.main=cm)

e_vars<-apply(1/vars,1,quantile,probs=c(.5,.25,.75),na.rm=TRUE)
yub<-max(as.vector(e_vars))*1.05
x<-barplot(e_vars[1,],names.arg=LETTERS[1:5],xlab="Conditions",ylab="Precision in p",cex.lab=cl,ylim=c(0,yub))
segments(x,e_vars[2,],x,e_vars[3,])
title(main="(F) Predictability",cex.main=cm)

dev.off()

e_vars
#        [,1]     [,2]     [,3]      [,4]     [,5]
#50% 35.29516 6.155319 5.997340 11.328041 10.05917
#25% 33.30428 5.702595 5.526185  9.758039  8.88748
#75% 52.23909 8.291266 8.409779 16.552761 13.82038
sqrt(1/e_vars)
#         [,1]      [,2]      [,3]      [,4]      [,5]
#50% 0.1683226 0.4030647 0.4083388 0.2971136 0.3152963
#25% 0.1732806 0.4187586 0.4253900 0.3201244 0.3354368
#75% 0.1383573 0.3472878 0.3448321 0.2457903 0.2689924

## 2 vs 5 x
ffst<-function(brk=.5,wl=NA,wh=NA,f=NA){
    if(f < brk){w<-wl}
    else{w<-wh}
    return(w)
 }

sdif<-function(w1=NA,w2=NA,f1=NA){
    wbar<-(w1*f1 + w2*(1-f1))
    psf1<-f1 * (w1/wbar)
    psf2<-(1-f1) * (w2/wbar)
    return(c(psf1,psf2))
    }

sdifh2<-function(w1=NA,w2=NA,f1=NA,h2=NA){
    wbar<-(w1*f1 + w2*(1-f1))
    psf1<-f1 * (w1/wbar)
    psf2<-(1-f1) * (w2/wbar)
    Sd<-psf1-f1
    psf1<-(Sd*h2)+f1
    psf2<-1-psf1
    return(c(psf1,psf2))
    }


## simulate uncertainty in function
## wl<-c(.52,.155)
## wh<-c(.155+.05,.52-.05)
## brk<-.85
brk<-runif(100,.8,.9) 
wl<-matrix(NA,nrow=100,ncol=2)
wl[,1]<-rbeta(100,2*52,2*48)
wl[,2]<-rbeta(100,2*62,2*338)

wh<-matrix(NA,nrow=100,ncol=2)
wh[,1]<-rbeta(100,2*(31+10),2*(169-10))
wh[,2]<-rbeta(100,2*(26-2.5),2*(24+2.5))

Ne<-110 ## from science paper, change in FHA

#586 (.446–.754) ## from Aarons paper
# qnorm(p=c(.5,.025,.975),mean=0.586,sd=.08)
#h2<-rbeta(100,16,4)## true .586, using .8 instead, 
## sd = 0.0873
h2<-rep(.8,100)
vars<-matrix(NA,nrow=5,ncol=100)
library(scales)

## hs = .8
pdf("nfdsSimsh2_80_2x5.pdf",width=12,height=8)
cm<-1.5
cl<-1.6
par(mfrow=c(2,3))
par(mar=c(5.5,4.5,2.5,1.5))


### with heritability, true 0.586, using .8 instead
## drift only
p<-rep(.5,100)  
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(0.85,c(0.52,0.155),c(.155+.05,.52-.05),p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
    
}    
vars[1,]<-apply(p,2,var)     
mtext(round(median(1/vars[1,]),1),side=1,adj=.75,line=-2.5,cex=1.5)

title(main="(A) Drift",cex.main=cm)

## uncertainty in s
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    p[j,i]<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}  
vars[2,]<-apply(p,2,var)    

mtext(round(median(1/vars[2,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(B) Unc. in selection, 2X",cex.main=cm)

## uncertainty in s and drift
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}     
vars[3,]<-apply(p,2,var)     

mtext(round(median(1/vars[3,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(C) Unc. in selection + drift, 2X",cex.main=cm)

brk<-runif(100,.8,.9) 
wl<-matrix(NA,nrow=100,ncol=2)
wl[,1]<-rbeta(100,52*5,48*5)
wl[,2]<-rbeta(100,62*5,338*5)

wh<-matrix(NA,nrow=100,ncol=2)
wh[,1]<-rbeta(100,(31+10)*5,(169-10)*5)
wh[,2]<-rbeta(100,(26-2.5)*5,(24+2.5)*5)

### with heritability , true 0.586, using .8 instead
## uncertainty in s
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    p[j,i]<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}  
vars[4,]<-apply(p,2,var)    
     
mtext(round(median(1/vars[4,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(D) Unc. in selection, 5X",cex.main=cm)

## uncertainty in s and drift
p<-rep(.5,100)
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Trait freq. (p)",cex.lab=cl)    

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    if(j < 100){lines(p[j,],col=alpha("gray",.5),lwd=0.6,type='l')}
    else{lines(p[j,],col="black",lwd=1,type='l')}
}     
vars[5,]<-apply(p,2,var)     
   
mtext(round(median(1/vars[5,]),1),side=1,adj=.75,line=-2.5,cex=1.5)
title(main="(E) Unc. in selection + drift, 5X",cex.main=cm)

e_vars<-apply(1/vars,1,quantile,probs=c(.5,.25,.75),na.rm=TRUE)
yub<-max(as.vector(e_vars))*1.05
x<-barplot(e_vars[1,],names.arg=LETTERS[1:5],xlab="Conditions",ylab="Precision in p",cex.lab=cl,ylim=c(0,yub))
segments(x,e_vars[2,],x,e_vars[3,])
title(main="(F) Predictability",cex.main=cm)

dev.off()

e_vars

sqrt(1/e_vars)





