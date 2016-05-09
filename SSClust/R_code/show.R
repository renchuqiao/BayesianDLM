###
###
# Modules that used to output the gene name for each clust and draw 
# necessary plots.
#
# mclust.plot plot each individual gene expression curves of all clusters 
#  
#	sm is the raw data with each raw represent a gene expression 
# 	value curve.  
#	nclust is the number of clusters that  
#	clust is the cluster membership
#	
#
# gene.name.clust is the function that used to output the gene names for 
# for each cluster.
#	my.name is the gene names for the full data set.
#	nclust is the number of clusters
#	clust is the clust labels obtained from any clustering method
#
# branch.fit is the function to estimate the mean curvec using 
# branching spline.
#	my.clean is the data as sm in mclust.plot
#	clust.ind is the clust label
#	NN is the number of clusters
#	change is the last point that the two branch have common pattern
#	branchl is the number of points that the two branches have different 
#		patterns.
#
# Output: gene name file, raw_clust_plot, center curve and its standard 
#         deviation curve
#
# Author: Wenxuan Zhong, Ping Ma
###
### 

mclust.plot=function(sm,nclust,clust){
pp=dim(sm)[2]
par(mfrow=c(2,2))
for(i in 1:nclust){class=sm[clust==i,]
                    
                   temp.x=1:pp
                   temp.y=class[1,]
                   plot(temp.x,temp.y, type="l", ylim=c(min(class, na.rm=T), max(class, na.rm=T)),
				xlab="time",ylab="gene expression", main=paste("cluster",i,sep="",collapse=NULL))
                   for(j in 1:dim(class)[1]) lines(temp.x, class[j,])
                   
                   }
}



write.name=function(my.names,index,i)
{class=my.name[index==i]
 write(class,paste("cluster",i,".txt",collapse=NULL,sep=""))
}


gene.name.clust=function(my.name,nclust,index){
for(i in 1:nclust)
	write.name(my.name,index,i)
}

mean.fit=function(x.select,nclass,cl){
pp=dim(x.select)[2]
tt=1:pp
UCL=matrix(nrow=nclass,ncol=length(tt))
LCL=matrix(nrow=nclass,ncol=length(tt))
curv=matrix(nrow=nclass,ncol=length(tt))
library(gss)
 

for(i in 1:nclass)
  {
    class=x.select[cl==i,]
    nn=dim(class)[1]
    y=as.vector(t(class))
    my.time=rep(tt,nn)
    geneid=as.factor(sort(rep(tt,nn)))
    formula=~1|geneid
    form.wk <- terms.formula(formula)[[2]]
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
            env=c(min(my.time),max(my.time))+c(-1,1)*(max(my.time)-min(my.time))*.00)

    temp=ssanova(y~my.time,type=list(my.ime=list("custom",par)),  random=my.ran,id.basis=tt,alpha=1)
    temp.pred=predict(temp,data.frame(my.time=tt),se.fit=T)
    UCL[i,]=temp.pred$fit+1.96*temp.pred$se.fit
    LCL[i,]=temp.pred$fit-1.96*temp.pred$se.fit
    curv[i,]=temp.pred$fit
}

postscript(paste("group",nclass,".ps",sep="",collapse=NULL))
par(mfrow=c(2,2))
for(i in 1:nclass)
  {
   plot(tt,curv[i,],type="l")
   lines(tt,LCL[i,],col=2)
   lines(tt,UCL[i,],col=2)

	}

}




branch.fit=function(my.clean,clust.ind,NN,change,branchl){
PP=change+branchl
clust.partial.part1=matrix(nrow=NN,ncol=PP)
clust.partial.sigma.part1=matrix(nrow=NN,ncol=PP)
clust.partial.part2=matrix(nrow=NN,ncol=PP)
clust.partial.sigma.part2=matrix(nrow=NN,ncol=PP)
meas.clean=PP
part1=my.clean[,1:PP]
part2=my.clean[,c(1:change,(change+branchl+1):(change+2*branchl))]
library(gss)
for(i in 1:NN){
                nclust=i
		classi=NULL
		my.clean=NULL
	  	my.clean=part1
		classi=my.clean[clust.ind==nclust,]
		clusti=NULL
		clusti=as.vector(t(classi))
  		time=NULL
		time=rep(1:meas.clean,dim(classi)[1])
		geneid=NULL
  		geneid=as.factor(sort(rep(1:dim(classi)[1],meas.clean)))
                formula=~1|geneid
                form.wk <- terms.formula(formula)[[2]]
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
                clusti=clusti
                time=time
                geneid=geneid 
            par <- list(nphi=1,mkphi=mkphi.cubic,mkrk=mkrk.cubic,
            env=c(min(time),max(time))+c(-1,1)*(max(time)-min(time))*.00)
  		temp=ssanova(clusti~time,type=list(time=list("custom",par)),   random=my.ran,
                		id.basis=1:meas.clean)            
  		temp.pred=predict(temp,data.frame(time=time[1:meas.clean]),se.fit=T)
  		clust.partial.part1[i,]=temp.pred$fit
  		clust.partial.sigma.part1[i,]=temp.pred$se.fit

		classi=NULL
		my.clean=NULL
	  	my.clean=part2
		classi=my.clean[clust.ind==nclust,]
		clusti=NULL
		clusti=as.vector(t(classi))
  		time=NULL
		time=rep(1:meas.clean,dim(classi)[1])
		geneid=NULL
  		geneid=as.factor(sort(rep(1:dim(classi)[1],meas.clean)))
                formula=~1|geneid
                form.wk <- terms.formula(formula)[[2]]
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
                clusti=clusti
                time=time
                geneid=geneid 
              par <- list(nphi=1,mkphi=mkphi.cubic,mkrk=mkrk.cubic,
            env=c(min(time),max(time))+c(-1,1)*(max(time)-min(time))*.00)

  		temp=ssanova(as.vector(t(classi))~time, type=list(time=list("custom",par)),   andom=my.ran,
                		id.basis=1:meas.clean) 
  		temp.pred=predict(temp,data.frame(time=time[1:meas.clean]),
			se.fit=T)
 		clust.partial.part2[i,]=temp.pred$fit
  		clust.partial.sigma.part2[i,]=temp.pred$se.fit
		}
dump(c("clust.partial.part1","clust.partial.sigma.part1",
     "clust.partial.part2","clust.partial.sigma.part2"),"branching.Rdata")
}
 
