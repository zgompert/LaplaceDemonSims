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


### with heritability, true 0.586, using .8 instead
## drift only
p<-rep(.5,100)  

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(0.85,c(0.52,0.155),c(.155+.05,.52-.05),p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
    
}    
vars[1,]<-apply(p,2,var)     


## uncertainty in s
p<-rep(.5,100)

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    p[j,i]<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    }
}  
vars[2,]<-apply(p,2,var)    


## uncertainty in s and drift
p<-rep(.5,100)

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
}     
vars[3,]<-apply(p,2,var)     


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

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    p[j,i]<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    }
}  
vars[4,]<-apply(p,2,var)    

## uncertainty in s and drift
p<-rep(.5,100)

p<-matrix(.5,nrow=100,ncol=100)
for(j in 1:100){
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[j,i-1])
    epi<-sdifh2(w[1],w[2],p[j,i-1],h2[j])[1]
    p[j,i]<-rbinom(n=1,size=2*Ne,prob=epi)/(2*Ne)
    }
}     
vars[5,]<-apply(p,2,var)     
   
pdf("time_nfdsSimsh2_80.pdf",width=12,height=8)
cm<-1.5
cl<-1.6
par(mfrow=c(2,3))
par(mar=c(5.5,4.5,2.5,1.5))
tit<-c("(A) Drift","(B) Unc. in selection","(C) Unc. in selection + drift","(D) Unc. in selection, 2X",
       "(E) Unc. in selection + drift, 2X")
for(i in 1:5){
	plot(1/vars[i,],type='b',xlab="Time",ylab="Predictability (precision)",cex.lab=cl)   
	title(main=tit[i],cex.main=cm)
}
dev.off()
