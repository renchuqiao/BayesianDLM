##############################################################
##############################################################
##
## This program is used to cluster data using smoothing 
## spline clustering method. This program 
## is an iteration algorithm using EM
##
## Input:
##	sm is a matrix with each represent the gene expression curve
##	tt is the time points that each sample was sampled
##	clust is the starting point for the clusters.
##		random: sample(rep(#of clusters, replications))
##	nclust is the number of clusters
##	clust.initial: is the initial mean curves
##      zeta: is the log10 based initial variance ratio parameters estimates
##      varht: is the initial variance estimates
##      thres: is the threshold that used for rejection-controlled EM
##      like : is the initial likelihood values 
##
##
## Output: 
##	clust.label: is the label of clusters for each gene
##      clust.initial: is the mean curve of each clust
##	likelihood.rec: is the likelihood value of each iteration
##
##
##  Author: Wenxuan Zhong, Ping Ma
##  Email:  wenxuan@uiuc.edu, pingma@uiuc.edu          
##
#################################################################
#################################################################


#________________________________________________________________
#The next function is used to get posterior probability and the likelihood
#________________________________________________________________
library(mvtnorm)
post=function(xx, nclust,pai,mu,zeta,varht)
  {temp=NULL
   m=length(xx)
   x.which=which(is.na(xx))
   xx[x.which]=0
   ind.p=NULL
      for(i in 1:nclust)
          {
           clust.sigma=matrix(0, nrow=m, ncol=m)
           tempvar=10^zeta[i]*varht[i]
           rho=tempvar/(tempvar+varht[i])
           clust.sigma=(matrix(rho,m,m)+diag(rep((1-rho),m)))*(tempvar+varht[i]) 
           temp[i]=pai[i]*dmvnorm(xx, mu[i,],clust.sigma)
          }
   likelihood=-log(sum(temp))
   print(likelihood)
   ind.pvalue=temp*1000/sum(temp*1000)
   ind.p=which(ind.pvalue==max(ind.pvalue))[1]
   return(list(like=likelihood, clust=ind.p, pvalue=ind.pvalue))
 }

#_________________________________________________________________
# estimate the parameters using EM
#_________________________________________________________________

em.clust=function(sm,tt,clust,nclust,clust.initial,zeta,varht,thres,iter.max){
  iter=1
  like.rec=bic.rec=NULL
  n=nrow(sm)
  m=ncol(sm)
  likelihood.bef=likelihood=0
  
  gamma.prior=NULL
  for(i in 1:nclust)gamma.prior[i]=sum(clust==i)/n
  rm(clust)
  trc=rep(0,nclust)
  mkrk.nominal <- function (levels) 
{
    k <- length(levels)
    if (k < 2) 
        stop("error: factor should have at least two levels")
    code <- 1:k
    names(code) <- as.character(levels)
    env <- list(code = code, table = diag(k))
    fun <- function(x, y, env, outer.prod = FALSE) {
        if (!(is.factor(x) & is.factor(y))) {
            stop("error in rk: inputs are of wrong types")
        }
        x <- as.numeric(env$code[as.character(x)])
        y <- as.numeric(env$code[as.character(y)])
        if (any(is.na(c(x, y)))) {
            stop("error in rk: unknown factor levels")
        }
        if (outer.prod) 
            env$table[x, y]
        else env$table[cbind(x, y)]
    }
    list(fun = fun, env = env)
}
 
  for(iter in 1:iter.max){
         likelihood.bef=likelihood
         like=NULL
         
         tk=matrix(nrow=n,ncol=nclust)
         for(i in 1:n){
         temp=NULL
         temp=post(xx=sm[i,],nclust=nclust,pai=gamma.prior,mu=clust.initial ,zeta=zeta,varht=varht)
         tk[i,]=temp$pvalue
         like[i]=temp$like
         rm(temp)
         }
         likelihood=sum(like)
         print(paste("liklelihood is", likelihood))
         bic=2*likelihood + (sum(trc)+ 4*nclust)*log(n)
         print(paste("bic is ", bic))
         like.rec[iter]=likelihood
         bic.rec[iter]=bic
         
        b.matrix=matrix(0, nrow=n,ncol=nclust)
        my.label=matrix(0, nrow=n,ncol=nclust)
        
         for(kk in 1:nclust){
         em.select=NULL  
         em.select=which(tk[,kk]>=thres)
         label1=NULL
         label1=which(tk[,kk]<thres)
         my.remain=NULL
         my.remain=tk[label1,kk]
         sem.select=NULL
         sem.select=as.vector(tapply(my.remain,1:length(my.remain),
           function(x,p=thres){sample(c(1,0),1,prob=c(x/p,(1-x/p)))}))
         label=NULL
         label=c(em.select,label1[which(sem.select==1)])
         my.label[label,kk]=rep(1,length(label))
         rm(list=c("em.select","label1","my.remain","sem.select"))

         my.weight1=matrix(nrow=m,ncol=length(label))
         for(mm in 1:length(label))my.weight1[,mm]=rep(tk[label[mm],kk],m)      
         temp=NULL
         temp.pred=NULL
         clusti.ini=as.numeric(as.vector(t(sm[label,])))
         tm.ini=rep(tt,length(label))
         my.weight.ini=as.numeric(as.vector(my.weight1))
         my.geneid.ini=as.factor(sort(rep(1:length(label),m)))
         my.which=which(as.numeric(is.na(clusti.ini))==0)
         clusti=clusti.ini[my.which]
         tm=tm.ini[my.which]
         my.weight=my.weight.ini[my.which]
         my.geneid=my.geneid.ini[my.which]
         
         formula=~1|my.geneid
      	 form.wk <- terms.formula(formula)[[2]]
         if (!("|" %in% strsplit(deparse(form.wk), "")[[1]])) 
           stop("error: missing | in grouping formula")
         term.wk <- strsplit(deparse(form.wk), " \\| ")[[1]]
         z2.wk <- eval(parse(text = term.wk[2]))  
         z <- NULL
         lvl.z2 <- levels(z2.wk)
         for (i in lvl.z2) z <- cbind(z, as.numeric(z2.wk == i))
         init <- -1
         env <- length(levels(z2.wk))
         fun <- function(zeta, env) diag(10^(-zeta), env)
         sigma <- list(fun = fun, env = env)    	
         my.ran=list(z = z, sigma = sigma, init = init)
         par <- list(nphi=1,mkphi=mkphi.cubic,mkrk=mkrk.cubic,
            env=c(min(tm),max(tm))+c(-1,1)*(max(tm)-min(tm))*.3)
         temp=ssanova1(clusti~tm,type=list(tm=list("custom",par)), weights=my.weight,
	      random=my.ran, id.basis=1:m, alpha=1)
#              )
         temp.pred=predict(temp,data.frame(tm=tt))
         b.matrix[label,kk]=temp$b
         zeta[kk]=temp$zeta
         varht[kk]=temp$varht            
         clust.initial[kk,]=temp.pred
         trc[kk]=temp$trc
         rm(temp,my.weight1)
         
         gamma.prior[kk]=sum(tk[,kk])/n
        }
  
        b.est=apply(b.matrix,1, function(x) mean(x[which(x!=0)])) 
        
         for(kk in 1:nclust){
         label.ttt=as.numeric(my.label[,kk])
         label=which(label.ttt==1)
         my.weight1=matrix(nrow=m,ncol=length(label))
         for(mm in 1:length(label))my.weight1[,mm]=rep(tk[label[mm],kk],m)      
         temp=NULL
         temp.pred=NULL
         
         my.geneid=as.factor(sort(rep(1:length(label),m)))
         formula=~1|my.geneid
         form.wk <- terms.formula(formula)[[2]]
         if (!("|" %in% strsplit(deparse(form.wk), "")[[1]])) 
           stop("error: missing | in grouping formula")
         term.wk <- strsplit(deparse(form.wk), " \\| ")[[1]]
         z2.wk <- eval(parse(text = term.wk[2]))  
         z <- NULL
         lvl.z2 <- levels(z2.wk)
         for (i in lvl.z2) z <- cbind(z, as.numeric(z2.wk == i))
         init <- -1
         env <- length(levels(z2.wk))
         fun <- function(zeta, env) diag(10^(-zeta), env)
         sigma <- list(fun = fun, env = env)    	
         my.ran=list(z = z, sigma = sigma, init = init)
         
         clusti.ini=as.numeric(as.vector(t(sm[label,])))
         my.zb.offset.ini=as.numeric(my.ran$z%*%b.est[label])
         tm.ini=rep(tt,length(label))
         my.weight.ini=as.numeric(as.vector(my.weight1))
         my.geneid.ini=as.factor(sort(rep(1:length(label),m)))
         my.which=which(as.numeric(is.na(clusti.ini))==0)
         clusti=clusti.ini[my.which]
         tm=tm.ini[my.which]
         my.weight=my.weight.ini[my.which]
         my.geneid=my.geneid.ini[my.which]
         my.zb.offset=my.zb.offset.ini[my.which]
         par <- list(nphi=1,mkphi=mkphi.cubic,mkrk=mkrk.cubic,
            env=c(min(tm),max(tm))+c(-1,1)*(max(tm)-min(tm))*.3)

         temp=ssanova1(clusti~tm, type=list(tm=list("custom",par)),  weights=my.weight,
	      offset=my.zb.offset,id.basis=1:m, alpha=1)
     
         
         temp.pred=predict(temp,data.frame(tm=tt, offset=rep(0,length(tt))))
         clust.initial[kk,]=temp.pred
         trc[kk]=temp$trc
          }          
          
	if(likelihood <= min(like.rec)){
	mean.curve=clust.initial
	var.ratio=zeta
	variance=varht
	pp=gamma.prior	
	minlike=min(like.rec)
        minbic=min(bic.rec)
	}
}    
return(list(mean.curve=mean.curve,var.ratio=var.ratio,variance=variance,pp=pp,likelihood=minlike, bic=minbic))
}

