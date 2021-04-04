## Hanski metapopulation model, deterministic approximations for the stochastic model

## functions
## calcualte probability of patch occupancy
calcPi<-function(mj=NA,ei=NA,pj=NA){
    sum_m_p<-sum(mj*pj)
    pi<-sum_m_p/(ei+sum_m_p)
    return(pi)
    }

## calculate expected phenotype at time of colonization
calcZiC<-function(mj=NA,pj=NA,zj=NA){
    zic<-sum(mj*zj*pj)/sum(mj*pj)
    return(zic)
    }
    
## calculate expected phenotype more generally
calcZi<-function(mj=NA,pj=NA,zic=NA,Ezi=NA,gamma=NA,sig2=NA,rho=NA,Ai=NA,thetai=NA){
    top<-gamma*sig2*thetai +
        (Ezi + (rho/Ai) * sum(mj*pj)) * zic
    bottom<-Ezi + gamma*sig2 + (rho/Ai) * sum(mj*pj)
    zi<-top/bottom
    return(zi)
    }
    
## calculate intrinsic growth parameter for z*
calcRi<-function(r0=NA,gamma=NA,sig2p=NA,thetai=NA,zi=NA,sel=NA){
    if(sel==1){# hard
        ri<-r0 - (gamma/2)*sig2p - (gamma/2)*(thetai-zi)^2
    }
    else{
        ri<-r0 - (gamma/2)*sig2p
        ri<-rep(ri,length(zi)) ## so it doesn't collapse to a single value 
    }        
    return(ri)
    }

## calculate extinction probability ... need to use actual exponential here, otherwise hanski equation gets off, rate is not the same thing as probability per unit time when rate > 1
calcEi<-function(Ai=NA,ri=NA,v=NA){
    s<-(2*ri/v)
    K<-Ai
    k<-log(K)
    left<-(K^s)/(s*ri)
    right<-1 - ((1+s*k)/(exp(s*k)))
    erate<-1/(left*right)
    ## convert to per gen. extinction prob
    ei<-dexp(1,erate)
    return(ei)
    }
    
## calculate migration/colonization rate from each j to i
calcMigi<-function(C=NA,rj=NA,beta=NA,thetai=NA,zj=NA,Ai=NA,Aj,alpha=NA,dj=NA){
    Kj<-Aj
    mij<-C * Kj * exp(rj) * exp(-1 * beta * abs(thetai - zj)) * 
        Ai * ((alpha^2)/(2*pi)) * exp(-1 * alpha * dj)
    return(mij)
    }




## model parameters and default values

#gamma<-5 ## strength of local stabalizing selection ... set to 0 to kill feedback
#sig2<-1 ## additive genetic variance
#sig2p<-sig2 ## phenotypic variance, assuming h2 = 1
#r1<-1 ## average growth rate

#v<-1 ## variance of growth rate
#C<-0.002 ## dispersal paramter
#alpha<-0.5 ## spatial range of dispersal
#delta<-0.8 ## difference in optimal phenotype; thetai is related, optimal trait value
#f1<-0.5 ## frequency of patch type 1
#rho<-0 ## gene flow
#rho<-40
#beta<-0 ## immigrant selection, set to 0 for no immigrant selection
#beta<-3
#hard<-1 ## logical hard or soft selection

hanski<-function(metapop=NULL,hard=1,gamma=5,sig2=1,sig2p=1,r1=1,v=1,C=0.002,alpha=0.5,delta=0.8,f1=0.5,rho=0,beta=0,Npatch=100,Nsim=100){
    ## set r0 based on provided parameters
    r0<-r1 +(gamma/2)*sig2p 
    
    if(is.null(metapop)==TRUE){## pop was not provided, make one
        metapop<-makePop(Npatch=Npatch,delta=delta,f1=f1)
    }

    ## drop the list, pull out the bits for easier reference
    D<-metapop$D;A<-metapop$A;t1<-metapop$t1;t2<-metapop$t2;theta<-metapop$theta

    ## iterative solve, start all p and z at 0.5 following Hanski

    ## initialize p and z
    P<-rep(0.5,Npatch)
    Z<-P
    Zc<-Z
    M<-matrix(NA,nrow=Npatch,ncol=(Npatch-1))

    for(xx in 1:Nsim){
        ## growth rate
        R<-calcRi(r0=r0,gamma=gamma,sig2p=sig2p,thetai=theta,zi=Z,sel=hard)
        ## extinction
        E<-calcEi(Ai=A,ri=R,v=v)
        for(i in 1:Npatch){
            M[i,]<-calcMigi(C=C,rj=R[-i],beta=beta,thetai=theta[i],zj=Z[-i],A[i],A[-i],alpha=alpha,dj=D[i,-i])
        }

        ## calcualte initial and equil. trait values
        for(i in 1:Npatch){
            Zc[i]<-calcZiC(mj=M[i,],pj=P[-i],zj=Z[-i])
        }
        Ztemp<-Z
        for(i in 1:Npatch){
            Z[i]<-calcZi(mj=M[i,],pj=P[-i],zic=Zc[i],Ezi=E[i],gamma=gamma,sig2=sig2,rho=rho,Ai=A[i],thetai=theta[i])
        }

        ## patch occupancy prob.
        Ptemp<-P
        for(i in 1:Npatch){
            P[i]<-calcPi(mj=M[i,],ei=E[i],pj=P[-i])
        }
        DZ<-dist(rbind(Z,Ztemp))
        DP<-dist(rbind(P,Ptemp))
        cat(xx," ..... E. dist. P = ",DP,"; Z = ",DZ,"\n",sep="")
    }
    RanP<-rbinom(n=Npatch,size=1,prob=P)
    out<-list(Z=Z,P=P,RanP=RanP)
    return(out)
}

## construct area, do this in seperate function to re-use
makePop<-function(Npatch=100,delta=0.8,f1=0.5){
## Ki is set to Ai
    X<-runif(Npatch,0,10) ## 10x10 area
    Y<-runif(Npatch,0,10)
    D<-as.matrix(dist(cbind(X,Y)))## Euclidean, make absolute?
    A<-rlnorm(n=Npatch,meanlog=2,sdlog=sqrt(0.5))
    t1<-0.5-0.5*delta;t2<-t1+delta
    theta<-sample(c(t1,t2),Npatch,replace=TRUE,prob=c(f1,1-f1))
    pop<-list(D=D,A=A,t1=t1,t2=t2,theta=theta,X=X,Y=Y)
    return(pop)
    }

pops<-vector("list",100)
resultsHard<-vector("list",100)
resultsSoft<-vector("list",100)
for(rep in 1:100){
    pops[[rep]]<-makePop(Npatch=100,delta=0.8,f1=0.5)
    resultsHard[[rep]]<-hanski(metapop=pops[[rep]],alpha=.2,gamma=3,sig2=0.01,sig2p=1,rho=40,beta=0,C=0.002,hard=1)
    resultsSoft[[rep]]<-hanski(metapop=pops[[rep]],alpha=.2,gamma=3,sig2=0.01,sig2p=1,rho=40,beta=0,C=0.002,hard=0)
}

mae<-rep(NA,100)
for(rep in 1:100){
    mae[rep]<-mean(abs(resultsSoft[[rep]]$P-resultsHard[[rep]]$P))
}


### trying other values of selection
resultsHardStr<-vector("list",100)
resultsSoftStr<-vector("list",100)
resultsHardWeak<-vector("list",100)
resultsSoftWeak<-vector("list",100)
for(rep in 1:100){
    resultsHardStr[[rep]]<-hanski(metapop=pops[[rep]],alpha=.2,gamma=5,sig2=0.01,sig2p=1,rho=40,beta=0,C=0.002,hard=1)
    resultsSoftStr[[rep]]<-hanski(metapop=pops[[rep]],alpha=.2,gamma=5,sig2=0.01,sig2p=1,rho=40,beta=0,C=0.002,hard=0)
    resultsHardWeak[[rep]]<-hanski(metapop=pops[[rep]],alpha=.2,gamma=1,sig2=0.01,sig2p=1,rho=40,beta=0,C=0.002,hard=1)
    resultsSoftWeak[[rep]]<-hanski(metapop=pops[[rep]],alpha=.2,gamma=1,sig2=0.01,sig2p=1,rho=40,beta=0,C=0.002,hard=0)
}

maeStr<-rep(NA,100)
maeWeak<-rep(NA,100)
for(rep in 1:100){
    maeStr[rep]<-mean(abs(resultsSoftStr[[rep]]$P-resultsHardStr[[rep]]$P))
    maeWeak[rep]<-mean(abs(resultsSoftWeak[[rep]]$P-resultsHardWeak[[rep]]$P))
}

library(scales)
library(RColorBrewer)
cm<-1.6;cl<-1.6;ca<-1.2
pdf("fig_occHardSoft.pdf",width=4,height=12)
cs<-alpha(c("orange","blue"),.7)

par(mfrow=c(3,1))
par(mar=c(4.5,4.5,3,3))
plot(pops[[rep]]$X,pops[[rep]]$Y,pch=19,cex=log(pops[[rep]]$A*.8),axes=FALSE,col=cs[round(pops[[rep]]$theta)+1],xlab="",ylab="")
box()
title("(A) Patch types",cex.main=cm)

plot(resultsSoft[[rep]]$P,resultsHard[[rep]]$P,pch=19,col=cs[round(pops[[rep]]$theta)+1],xlab="Soft selection",ylab="Hard selection",cex.lab=cl,cex.axis=ca)
abline(a=0,b=1)
title("(B) Predicted patch occupancy",cex.main=cm)

mns<-c(mean(maeWeak),mean(mae),mean(maeStr))
sds<-c(sd(maeWeak),sd(mae),sd(maeStr))
x<-barplot(mns,ylim=c(0,0.25),xlab="Selection",ylab="Error",cex.lab=cl,cex.axis=ca,cex.names=1.3,names=c(1,3,5))
segments(x,mns+sds,x,mns-sds)
box()
title("(C) Prediction error",cex.main=cm)

dev.off()




library(scales)
library(RColorBrewer)
cm<-1.5
pdf("fig_occHardSoft.pdf",width=9,height=9)
cs<-alpha(c("orange","blue"),.7)
cprc<-brewer.pal(9,"Greys")
cprob<-rep(cprc[1],100)
cseq<-c(.1,.2,.3,.4,.6,.7,.8,.9)

par(mfrow=c(2,2))
par(mar=c(4.5,4.5,3,3))
plot(pops[[rep]]$X,pops[[rep]]$Y,pch=19,cex=log(pops[[rep]]$A*.8),axes=FALSE,col=cs[round(pops[[rep]]$theta)+1],xlab="",ylab="")
box()
title("(A) Patch types",cex.main=cm)

## soft
cprob<-rep(cprc[1],100)
for(k in 1:8){
    a<-which(resultsSoft[[rep]]$P > cseq[k])
    cprob[a]<-cprc[k+1]
    }
plot(pops[[rep]]$X,pops[[rep]]$Y,pch=19,cex=log(pops[[rep]]$A*.8),axes=FALSE,col=cprob,xlab="",ylab="")
box()
title("(B) Occupancy with soft selection",cex.main=cm)


## Hard
cprob<-rep(cprc[1],100)
for(k in 1:8){
    a<-which(resultsHard[[rep]]$P > cseq[k])
    cprob[a]<-cprc[k+1]
    }
plot(pops[[rep]]$X,pops[[rep]]$Y,pch=19,cex=log(pops[[rep]]$A*.8),axes=FALSE,col=cprob,xlab="",ylab="")
box()
title("(C) Occupancy with hard selection",cex.main=cm)


plot(resultsSoft[[rep]]$P,resultsHard[[rep]]$P,pch=19,col=cs[round(pops[[rep]]$theta)+1],xlab="Soft selection",ylab="Hard selection",cex.lab=1.5,cex.axis=1.1)
abline(a=0,b=1)
title("(D) Predicted occupancy hard vs. soft",cex.main=cm)
dev.off()






