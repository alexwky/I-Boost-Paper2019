library(R.utils)
library(Hmisc)
library(survival)
library(data.table)
library(survIDINRI)

rawdata <- data.frame(fread("../Data/TCGA_8cancer_rmmis.csv",header=TRUE,sep=","),row.names=1)
x <- list(outcome=t(rawdata[1:2,]),xd=rawdata[-c(1:2),],fnames=rownames(rawdata)[-c(1,2)],snames=colnames(rawdata))

Clinical <- grep("Clinical_",x$fnames)
Clinical <- setdiff(Clinical, which(x$fnames=="Clinical_BRCA"))
CNV <- grep("CNV_",x$fnames)
miRNA <- grep("miRNA_",x$fnames)
Mutation <- grep("Mutation_",x$fnames)
Module <- grep("Module_",x$fnames)
Protein <- grep("Protein_",x$fnames)
GeneExp <- grep("Gene_",x$fnames)
data.types <- list(Clinical=Clinical,GeneExp=GeneExp,Module=Module,Protein=Protein,miRNA=miRNA,CNV=CNV,Mutation=Mutation)
data.types.names <- list("Clinical","Gene","Module","Protein","miRNA","CNV","Mutation")

features.rows <- list(Clinical,CNV,Mutation,miRNA,Protein)
features <- c("Clinical","CNV","Mutation","miRNA","Protein")
m <- length(features)

FEATURES <- list()
NAMES <- rep(NA, 2^m)
for (i in 0:(2^m-1)) {
	a <- unlist(strsplit(intToBin(i), split=""))
	a <- c(rep(0,m-length(a)), a)
	a <- rev(a)
	index <- c("")
	name <- c("")
	feature.tmp <- c()
	for (j in 1:m) {
		if (a[j] == "1") {
			index <- paste(index,j+1,sep="")
			name <- paste(name,features[j],sep=".")
			feature.tmp <- c(feature.tmp, features.rows[[j]])
		}
	}
	NAMES[i+1] <- paste(index,name,sep="")
	FEATURES <- c(FEATURES, list(feature.tmp))
}
NAMES[1] <- c(".")

insert.string <- function(name, index, feature) {
	a <- unlist(strsplit(name, split="\\."))
	if (length(a) == 0) a <- c("")
	out <- paste(c(paste(index,a[1],sep=""),feature,a[-1]), collapse=".")
	c(out)
}
combine <- function(x, object) c(object, x)

NAMES.95 <- as.vector(c(sapply(NAMES, insert.string, index=0, feature="GeneExp"), sapply(NAMES, insert.string, index=1, feature="Module"), NAMES[-1]))
FEATURES.95 <- c(lapply(FEATURES, combine, object=GeneExp), lapply(FEATURES, combine, object=Module), FEATURES[-1])
n.model <- length(NAMES.95)
FEATURES.95 <- FEATURES.95[-grep("GeneExp",NAMES.95)]
NAMES.95 <- setdiff(NAMES.95, grep("GeneExp",NAMES.95,value=TRUE))

get.size <- function(x) {
	a <- unlist(strsplit(x,split="\\."))[1]
	c(nchar(a))
}

split <- read.table("../Data/DataSplit.csv",sep=",",header=FALSE,row.names=1)
split <- split[match(x$snames,rownames(split)),]
mean.rmmis <- function(x) mean(x[!is.na(x)])

measure <- "NRI"
thres <- 36
all.methods.res <- list()
methods <- c("en","IBcv","IBperm")
for (method in methods) {
	indir <- paste("../Results/Prediction_",method,sep="")
	res <- list()
	for (type in c(1,2,3)) {
		type.name <- ifelse(type == 1, "LUAD", ifelse(type == 2, "KIRC", "all"))
		out <- c()
		out <- c(out, name)
		if (type == 1) filter.t <- (x$xd["Clinical_LUAD",] == 1)
		if (type == 2) filter.t <- (x$xd["Clinical_KIRC",] == 1)
		if (type == 3) filter.t <- rep(TRUE,ncol(x$xd))
		rm.features <- c()
		for (k in Clinical) if (length(unique(as.numeric(x$xd[k,filter.t]))) == 1) rm.features <- c(rm.features, k)

		accum <- c()
		select.nested <- c()
		type.size <- sapply(NAMES.95, get.size)
		best.measure <- c(0)
		best.upper <- c(0)
		best.lower <- c(0)
		model.size <- c()
		abs.C <- c()
		for (k in 1:6) {
			new.preds <- list()
			subset <- NAMES.95[type.size==k]
			for (grepname in accum) subset <- grep(grepname, subset, value=TRUE)
			if (k == 1) {
				newmodel <- c("2.Clinical")
				new.preds <- list()
				features <- FEATURES.95[[match("2.Clinical",NAMES.95)]]
				msize.s <- c()
				C.s <- c()
				for (s in 1:30) {
					features.t <- setdiff(features, rm.features)
					test.modules <- t(x$xd[features.t, filter.t & split[,s] == 1])
					test.outcome <- Surv(x$outcome[filter.t&split[,s]==1,2],x$outcome[filter.t&split[,s]==1,1])
	
					beta <- try(read.table(paste("../Results/Prediction_Cox/Model_split",s,"-",type.name,"-2.Clinical.csv",sep=""), header=FALSE, sep=",", row.names=1))
					if (class(beta) == "try-error") {
						new.preds.i <- c(new.preds.i, list(rep(NA, nrow(test.modules))))
						next
					}
					if (nrow(beta) > 0) preds <- as.matrix(test.modules[,match(rownames(beta),x$fnames[features.t])]) %*% beta[,1] else preds <- rep(0, nrow(test.modules))
					new.preds <- c(new.preds, list(preds))
					msize.s <- c(msize.s, nrow(beta))
					C.s <- c(C.s, 1-rcorr.cens(preds,test.outcome)[[1]])
				}
				abs.C <- c(abs.C, mean(C.s))
				newmodelsize <- mean(msize.s)
			} else {
				subset.NRI <- c()
				subset.lower.NRI <- c()
				subset.upper.NRI <- c()
				subset.Cindex <- c()
				subset.abs.C <- c()
				subset.msize <- c()
				for (name in subset) {
					features <- FEATURES.95[[match(name,NAMES.95)]]
					new.preds.i <- list()
					NRI <- c()
					upper.NRI <- c()
					lower.NRI <- c()
					Cindex <- c()
					msize.s <- c()
					C.s <- c()
					for (s in 1:30) {
						features.t <- setdiff(features, rm.features)
						test.modules <- t(x$xd[features.t, filter.t & split[,s] == 1])
						test.outcome <- Surv(x$outcome[filter.t&split[,s]==1,2],x$outcome[filter.t&split[,s]==1,1])
	
						nl <- try(readLines(paste(indir, "/Model_split",s,"-",type.name,"-",name,".csv",sep="")))
						if (class(nl) == "try-error") {
							new.preds.i <- c(new.preds.i, list(rep(NA, nrow(test.modules))))
							next
						} else nl <- length(nl)
						if (nl == 0) {
							new.preds.i <- c(new.preds.i, list(rep(NA, nrow(test.modules))))
							next
						}
						beta <- try(read.table(paste(indir, "/Model_split",s,"-",type.name,"-",name,".csv",sep=""), header=FALSE, sep=",", row.names=1))
						if (nl > 0) preds <- as.matrix(test.modules[,match(rownames(beta),x$fnames[features.t])]) %*% beta[,1] else preds <- rep(0, nrow(test.modules))
						new.preds.i <- c(new.preds.i, list(preds))
						if (any(is.na(old.preds[[s]]))) next
	
						IDINRI <- IDI.INF(test.outcome,old.preds[[s]],preds,thres,npert=300)
						NRI <- c(NRI, IDINRI$m2[1])
						lower.NRI <- c(lower.NRI, IDINRI$m2[2])
						upper.NRI <- c(upper.NRI, IDINRI$m2[3])
						Cindex <- c(Cindex, rcorr.cens(old.preds[[s]],test.outcome)[[1]] - rcorr.cens(preds,test.outcome)[[1]])
						msize.s <- c(msize.s, nrow(beta))
						C.s <- c(C.s, 1-rcorr.cens(preds,test.outcome)[[1]])
					}
					subset.NRI <- c(subset.NRI, mean.rmmis(NRI))
					subset.upper.NRI <- c(subset.upper.NRI, mean.rmmis(upper.NRI))
					subset.lower.NRI <- c(subset.lower.NRI, mean.rmmis(lower.NRI))
					subset.Cindex <- c(subset.Cindex, mean.rmmis(Cindex))
					new.preds <- c(new.preds, list(new.preds.i))
					subset.msize <- c(subset.msize, mean(msize.s))
					subset.abs.C <- c(subset.abs.C, mean(C.s))
				}
				newmodel <- subset[which.max(subset.NRI)]
				best.measure <- c(best.measure, max(subset.NRI))
				best.upper <- c(best.upper, subset.upper.NRI[which.max(subset.NRI)])
				best.lower <- c(best.lower, subset.lower.NRI[which.max(subset.NRI)])
				new.preds <- new.preds[[which.max(subset.NRI)]]
				newmodelsize <- subset.msize[which.max(subset.NRI)]
				abs.C <- c(abs.C, subset.abs.C[which.max(subset.NRI)])
			}
			old.preds <- new.preds
			select.nested <- c(select.nested,newmodel)
			model.size <- c(model.size, newmodelsize)
			added <- setdiff(unlist(strsplit(newmodel,split="\\."))[-1],accum)
			accum <- c(accum,added)
		}
		res <- cbind(res, cbind(accum,model.size,best.measure,best.lower,best.upper,abs.C))
	}
	all.methods.res <- c(all.methods.res, list(res))
}

res.LUAD <- lapply(all.methods.res,function(r){x<-r[,2:6];x<-apply(x,2,as.numeric);rownames(x)<-r[,1];x})
res.KIRC <- lapply(all.methods.res,function(r){x<-r[,8:12];x<-apply(x,2,as.numeric);rownames(x)<-r[,7];x})
res.all <- lapply(all.methods.res,function(r){x<-r[,14:18];x<-apply(x,2,as.numeric);rownames(x)<-r[,13];x})

col1 <- "magenta"
col2 <- "blue"
col3 <- "yellow"

plot.nested <- function(res,col=col,ylab=measure,pos=12,up=1,newplot=TRUE,outtype="measure",yrange=c(),text.pos=NULL,lty=1,...) {
	if (outtype == "measure") {
		range <- range(res[,2:4])
		if (length(yrange) == 0) yrange <- range(res[,2:4])
		if (newplot) plot(1:6,res[,2],col=col,ylab=measure,xlab="Number of Data Types",ylim=yrange,pch=16,xlim=c(0.5,6.5)) else
		points(1:6,res[,2],col=col,pch=16)
		lines(1:6,res[,2],col=col,lwd=2)
		arrows(x0=1:6,y0=res[,3],x1=1:6,y1=res[,4],length=0.05,angle=90,code=3,col=col,lty=lty)
		abline(h=0,lwd=2)
		if (length(text.pos) == 0) text.pos <- res[,2]+up*(range[2]-range[1])/pos
		text(1:6+0.1,text.pos,rownames(res),col=col)
	} else if (outtype == "size") {
		if (length(yrange) == 0) yrange <- c(0,max(range(res[,1]))+20)
		if (newplot) plot(1:6,res[,1],col=col,ylab="Model Size",xlab="Number of Data Types",ylim=yrange,pch=16,xlim=c(0.5,6.5)) else
		points(1:6,res[,1],col=col,pch=16)
		lines(1:6,res[,1],col=col,lwd=2)
		range <- range(res[,1])
		if (length(text.pos) == 0) text.pos <- res[,1]+up*(range[2]-range[1])/pos
		text(1:6+0.1,text.pos,rownames(res),col=col)
	}
}

plot.C.size <- function(res,col="red",ylab="C-Index",newplot=TRUE,yrange=range(res[,5]),xrange=range(res[,1]),lty=1,...) {
	if (newplot) plot(res[,1],res[,5],xlab="Model size",ylab="C-Index",pch=16,ylim=yrange,xlim=xrange,col=col) else
	points(res[,1],res[,5],pch=16,col=col)
	x0 <- res[-6,1]
	y0 <- res[-6,5]
	x1 <- res[-1,1]
	y1 <- res[-1,5]
	arrows(x0,y0,x1,y1,lty=lty,length=0.05,lwd=1,col=col)
}

reset <- function(j=0) {
        par(mfrow=c(1,1), oma=rep(0,4), mar=rep(j,4), new=TRUE)
        plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}

add_legend <- function(...) {
        opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 5, 0), mar=c(1, 0, 3, 0), new=TRUE)
        on.exit(par(opar))
        plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
        legend(...)
}

pdf("../Plots/Figure6.pdf")
par(mfrow=c(3,2),mar=c(7,4,1,1),oma=c(1,0,2,0))
cols <- c(col1,col2,col3)
for (j in 1:length(methods)) {
	yrange <- range(unlist(lapply(res.LUAD,function(r)r[,2:4])))
	yrange <- c(min(yrange),max(yrange)+0.45)
	text.pos <- seq(0.9,length=3,by=-0.15)[j]
	plot.nested(res.LUAD[[j]],col=cols[j],ylab=measure,newplot=(j==1),yrange=yrange,text.pos=text.pos,lty=1)
}
for (j in 1:length(methods)) {
	yrange <- range(unlist(lapply(res.LUAD,function(r)r[,5])))
	yrange <- c(min(yrange),max(yrange)+0.02)
	xrange <- range(unlist(lapply(res.LUAD,function(r)r[,1])))
	xrange <- c(min(xrange),max(xrange)+5)
	plot.C.size(res.LUAD[[j]],col=cols[j],newplot=(j==1),yrange=yrange,xrange=xrange)
}
for (j in 1:length(methods)) {
	yrange <- range(unlist(lapply(res.KIRC,function(r)r[,2:4])))
	yrange <- c(min(yrange),max(yrange)+0.55)
	text.pos <- seq(1,length=3,by=-0.16)[j]
	plot.nested(res.KIRC[[j]],col=cols[j],ylab=measure,newplot=(j==1),yrange=yrange,text.pos=text.pos,lty=1)
}
for (j in 1:length(methods)) {
	yrange <- range(unlist(lapply(res.KIRC,function(r)r[,5])))
	yrange <- c(min(yrange),max(yrange)+0.02)
	xrange <- range(unlist(lapply(res.KIRC,function(r)r[,1])))
	xrange <- c(min(xrange),max(xrange)+5)
	plot.C.size(res.KIRC[[j]],col=cols[j],newplot=(j==1),yrange=yrange,xrange=xrange)
}
for (j in 1:length(methods)) {
	yrange <- range(unlist(lapply(res.all,function(r)r[,2:4])))
	yrange <- c(min(yrange),max(yrange)+0.5)
	text.pos <- seq(0.8,length=3,by=-0.12)[j]
	plot.nested(res.all[[j]],col=cols[j],ylab=measure,newplot=(j==1),yrange=yrange,text.pos=text.pos,lty=1)
}
for (j in 1:length(methods)) {
	yrange <- range(unlist(lapply(res.all,function(r)r[,5])))
	yrange <- c(min(yrange),max(yrange)+0.02)
	xrange <- range(unlist(lapply(res.all,function(r)r[,1])))
	xrange <- c(min(xrange),max(xrange)+5)
	plot.C.size(res.all[[j]],col=cols[j],newplot=(j==1),yrange=yrange,xrange=xrange)
}
reset()
text(0.52,1.03,"LUAD",cex=1.1)
text(0.52,0.68,"KIRC",cex=1.1)
text(0.52,0.33,"Pan-Cancer",cex=1.1)

shape <- 16
add_legend("bottom", c("Elastic Net","I-Boost-CV","I-Boost-Permutation"),bty="n",
		pch=c(shape,shape,shape), col=cols, horiz=TRUE, pt.cex=1, cex=1)
dev.off()
