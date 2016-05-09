##############################################################
# This program implements the SSC algorithm described in     #
# Ma, Castillo-Davis, Zhong, and Liu (2006)                  #
# "A data-driven clustering method for Time Course Gene      #
# Expression Data"  Nucleic Acids Research, 34(4), 1261-1269.# 
#                                                            #
# Author: Wenxuan Zhong, Ping Ma                             #
# Email:  wenxuan@uiuc.edu, pingma@uiuc.edu          #
##############################################################

rm(list=ls())        
#___________________________________________________________________________________
#parameter specification

# "nchain" is how many chains you want to use when performing the RCEM
# "nclust" is the number of clusters (increment this by one until BIC goes up)
# "threshold" is the p-value cut-off for the RCEM algorithm
# replace "data.txt" with your data input file with the following tab-delimited
# format:
#
# genename     time1  time2  time3 ...
# geneA 0.54          0.22  ...
# geneB 0.11   0.76   0.004 ...
# ...  
# Note: missing values must be replaced by blanks, see test data.
nchain = 2
nclust = 118       
my.data = read.table("processed_data_chromosome_1_T.txt", header=T, na.strings =" ", sep="\t")
threshold = 0.1

#end parameter specification
#___________________________________________________________________________________



#_____________________  

my.name=rownames(my.data)
 my.sample=data.matrix(my.data)
n=dim(my.sample)[1]
m=dim(my.sample)[2]
if( n/m < 5) print("The dataset is too small to have so many clusters")
tt=1:m
library(gss)
source("R_code/ssanova1.R")
source("R_code/mcem.R")
like.rec=bic.rec=NULL

my.clust=list()
clust.centers=list()
for(chain in 1:nchain){
print(paste("chain", chain, "is running"), quote = F)
my.sample.vec=as.numeric(my.sample)
my.which=which(is.na(my.sample.vec))
my.sample.vec[my.which]=0
my.sample0=matrix(my.sample.vec, n, m)
clust=kmeans(my.sample,nclust)$cluster
zeta=NULL
varht=NULL
clust.initial=matrix(nrow=nclust,ncol=m)
 
for(kk in 1:nclust){
         label=which(clust==kk)
         temp=NULL
         temp.pred=NULL
         clusti.ini=as.numeric(as.vector(t(my.sample[label,])))
         tm.ini=rep(tt,length(label))
         my.geneid.ini=as.factor(sort(rep(1:length(label),m)))
         my.which=which(as.numeric(is.na(clusti.ini))==0)
         clusti=clusti.ini[my.which]
         tm=tm.ini[my.which]
         my.geneid=my.geneid.ini[my.which]
         formula=~1|my.geneid
         form.wk <- terms.formula(formula)[[2]]
         term.wk <- strsplit(deparse(form.wk), " \\| ")[[1]]
         z2.wk <- eval(parse(text = term.wk[2]))  
         z <- NULL
         lvl.z2 <- levels(z2.wk)
         for (k in lvl.z2) z <- cbind(z, as.numeric(z2.wk == k))
         init <- -1
         env <- length(levels(z2.wk))
         fun <- function(zeta, env) diag(10^(-zeta), env)
         sigma <- list(fun = fun, env = env)    	
         my.ran=list(z = z, sigma = sigma, init = init) 
         par <- list(nphi=1,mkphi=mkphi.cubic,mkrk=mkrk.cubic,
            env=c(min(tm),max(tm))+c(-1,1)*(max(tm)-min(tm))*.00)

         temp=ssanova(clusti~tm,
                random=my.ran,type=list(tm=list("custom",par)),  id.basis=1:m, alpha=1)
         temp.pred=predict(temp,data.frame(tm=tt,offset=rep(0,length(tt))))
         clust.initial[kk,]=temp.pred
         my.geneid=as.factor(sort(rep(1:length(label),m)))
         formula=~1|my.geneid
         form.wk <- terms.formula(formula)[[2]]
         term.wk <- strsplit(deparse(form.wk), " \\| ")[[1]]
         z2.wk <- eval(parse(text = term.wk[2]))  
         z <- NULL
         lvl.z2 <- levels(z2.wk)
         for (k in lvl.z2) z <- cbind(z, as.numeric(z2.wk == k))
         init <- -1
         env <- length(levels(z2.wk))
         fun <- function(zeta, env) diag(10^(-zeta), env)
         sigma <- list(fun = fun, env = env)    	
         my.ran=list(z = z, sigma = sigma, init = init) 
         random.zb=my.ran$z%*%temp$b
         temp.pred=predict(temp,data.frame(tm=tt))
         zeta[kk]=temp$zeta
         varht[kk]=temp$varht
         clust.initial[kk,]=temp.pred
         rm(temp)
         print(paste("cluster",kk, "is calculated"), quote = F)     
 }

clust.sim=NULL
clust.sim=em.clust(my.sample,tt,clust,nclust,clust.initial,zeta,varht,thres=threshold,iter.max=20)# 
clust.centers[[chain]]=clust.sim$mean.curve
my.zeta<-my.varht<-gamma.prior<-NULL
my.zeta=clust.sim$var.ratio
my.varht=clust.sim$variance
gamma.prior=clust.sim$pp
like.rec[chain]=clust.sim$likelihood
bic.rec[chain]=clust.sim$bic
clust.temp=NULL
for(i in 1:n)
	clust.temp[i]=post(my.sample[i,],nclust,gamma.prior,clust.centers[[chain]],my.zeta,my.varht)$clust
my.clust[[chain]]=clust.temp
}

####prepare output

output.clust=NULL


output.clust=my.clust[[which(bic.rec==min(bic.rec))[1]]]
output.center=matrix(nrow=nclust,ncol=m)
output.center=clust.centers[[which(bic.rec==min(bic.rec))[1]]]
output.bic=min(bic.rec)

dump(c("output.clust","output.center","output.bic"),"output.Rdata")


source("R_code/show.R")
gene.name.clust(my.name,nclust,output.clust)
write(my.name,"all_genes.txt")
postscript("raw_curves.ps")
mclust.plot(my.sample,nclust,output.clust)
dev.off()


x.select=my.sample
nclass=nclust
cl=output.clust
pp=dim(x.select)[2]
tt=1:pp
tt.grid=seq(min(tt), max(tt),length=101)
UCL=matrix(nrow=nclass,ncol=length(tt.grid))
LCL=matrix(nrow=nclass,ncol=length(tt.grid))
curv=matrix(nrow=nclass,ncol=length(tt.grid))
library(gss)


for(i in 1:nclass)
  {
    class=x.select[cl==i,]
    nn=dim(class)[1]
    y.ini=as.vector(t(class))
    my.time.ini=rep(tt,nn)
    geneid.ini=as.factor(sort(rep(tt,nn)))
    my.which=which(as.numeric(is.na(y.ini))==0)
    y=y.ini[my.which]
    my.time=my.time.ini[my.which]
    geneid=geneid.ini[my.which]
    formula=~1|geneid
    form.wk <- terms.formula(formula)[[2]]
    term.wk <- strsplit(deparse(form.wk), " \\| ")[[1]]
    z2.wk <- eval(parse(text = term.wk[2]))  
    z <- NULL
    lvl.z2 <- levels(z2.wk)
    for (k in lvl.z2) z <- cbind(z, as.numeric(z2.wk == k))
    init <- -1
    env <- length(levels(z2.wk))
    fun <- function(zeta, env) diag(10^(-zeta), env)
    sigma <- list(fun = fun, env = env)    	
    my.ran=list(z = z, sigma = sigma, init = init)
    par <- list(nphi=1,mkphi=mkphi.cubic,mkrk=mkrk.cubic,
            env=c(min(my.time),max(my.time))+c(-1,1)*(max(my.time)-min(my.time))*.00)
    temp=ssanova(y~my.time, random=my.ran,type=list(my.time=list("custom",par)),  id.basis=rep(tt,2), alpha=1)
    tt.grid=seq(min(tt), max(tt),length=101)
    temp.pred=predict(temp,data.frame(my.time=tt.grid),se.fit=T)
    UCL[i,]=temp.pred$fit+1.96*temp.pred$se.fit
    LCL[i,]=temp.pred$fit-1.96*temp.pred$se.fit
    curv[i,]=temp.pred$fit 
    if( is.infinite(max(UCL[i,])* min(LCL[i,]))) print("Confiddence interval can not be calulated")  
}

postscript("mean_curves.ps")
par(mfrow=c(2,2))
for(i in 1:nclass)
  {
   plot(tt.grid,curv[i,],type="l",xlab="time", ylab=paste("cluster", i ,sep="",  collapse=NULL), ylim=c(min(LCL[i,]), max(UCL[i,])) )
   lines(tt.grid,LCL[i,],col=2)
   lines(tt.grid,UCL[i,],col=2)

        }

dev.off()



##### print cluster number and BIC score ######
print(paste("Clusters =", nclust, sep=" "), quote = F);
print(paste("BIC Score =", round(bic.rec,2), sep=" "), quote = F);


