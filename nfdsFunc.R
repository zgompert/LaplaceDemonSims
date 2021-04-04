ff<-function(C=NA,s=NA,f=NA){
    w<-C - s * f
    return(w)
 }
 
sdif<-function(w1=NA,w2=NA,f1=NA){
    wbar<-(w1*f1 + w2*(1-f1))
    psf1<-f1 * (w1/wbar)
    psf2<-(1-f1) * (w2/wbar)
    return(c(psf1,psf2))
    }
    
C<-c(.5,.3)
s<-c(.3,.3)

p<-rep(.9,100)
for(i in 2:100){
    w<-ff(C,s,c(p[i-1],1-p[i-1]))
    p[i]<-sdif(w[1],w[2],p[i-1])[1]
    }    
    
plot(p,ylim=c(0,1),type='l')    

ffst<-function(brk=.5,wl=NA,wh=NA,f=NA){
    if(f < brk){w<-wl}
    else{w<-wh}
    return(w)
 }

library(scales)

pdf("stripeUnc.pdf",width=5,height=5)
par(mar=c(5,5,.5,.5))
wl<-c(.52,.155)
wh<-c(.2025,.21) 
wh<-c(.155*1.1,.52*0.9)
brk<-.85

p<-rep(.5,100)
for(i in 2:100){
    w<-ffst(brk,wl,wh,p[i-1])
    p[i]<-sdif(w[1],w[2],p[i-1])[1]
    }    
plot(p,ylim=c(0,1),type='n',xlab="Time",ylab="Stripe freq.",cex.lab=1.5)    
    
brk<-runif(100,.6,.9) 
wl<-matrix(NA,nrow=100,ncol=2)
wl[,1]<-rbeta(100,52,48)
wl[,2]<-rbeta(100,62,338)

wh<-matrix(NA,nrow=100,ncol=2)
wh[,2]<-rbeta(100,26,24)
wh[,1]<-rbeta(100,31,169)

for(j in 1:100){
p<-rep(.5,100)
for(i in 2:100){
    w<-ffst(brk[j],wl[j,],wh[j,],p[i-1])
    p[i]<-sdif(w[1],w[2],p[i-1])[1]
    }    
lines(p,col=alpha("gray",.5),lwd=.5,type='l')    
  }     
    
wl<-c(.52,.155)
wh<-c(.2025,.21) 
wh<-c(.155*1.1,.52*0.9)
brk<-.85

p<-rep(.5,100)
for(i in 2:100){
    w<-ffst(brk,wl,wh,p[i-1])
    p[i]<-sdif(w[1],w[2],p[i-1])[1]
    }    
lines(p)     

dev.off()    
    
    
    
    
    
    
ffst<-function(brk=c(.5,.6),wl=NA,wh=NA,whh=NA,f=NA){
    if(f < brk[1]){w<-wl}
    else if(f<brk[2]){w<-wh}
    else{w<-whh}
    return(w)
 }

wl<-c(.52,.155)
wh<-c(.2025,.21) 
whh<-c(.155*1.1,.52*0.9)
brk<-c(.78,.8)

p<-rep(.5,100)
for(i in 2:100){
    w<-ffst(brk,wl,wh,whh,p[i-1])
    p[i]<-sdif(w[1],w[2],p[i-1])[1]
    }    
plot(p,ylim=c(0,1),type='l')    
        
