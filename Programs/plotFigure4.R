library(R.utils)
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

thres <- 36
all.methods.NRI <- list()
for (method in c("IBcv","IBperm")) {
	all.NRI <- lapply(1:3,function(j)matrix(NA,95,30))
	for (s in 1:30) {
		split <- read.table("../Data/DataSplit.csv",sep=",",header=FALSE,row.names=1)
		split <- split[match(x$snames,rownames(split)),]
		for (i in seq(1,95)) {
			name <- NAMES.95[i]
			features <- FEATURES.95[[i]]
			if (!grepl("Clinical",name)) next
			for (type in c(1,2,3)) {
			        type.name <- ifelse(type == 1, "LUAD", ifelse(type == 2, "KIRC", "all"))
			        indir <- paste("../Results/Prediction_",method,sep="")
				if (type == 1) filter.t <- (x$xd["Clinical_LUAD",] == 1)
				if (type == 2) filter.t <- (x$xd["Clinical_KIRC",] == 1)
				if (type == 3) filter.t <- rep(TRUE,ncol(x$xd))
				rm.features <- c()
				for (k in Clinical) if (length(unique(as.numeric(x$xd[k,filter.t]))) == 1) rm.features <- c(rm.features, k)
				features.t <- setdiff(features, rm.features)
				test.modules <- t(x$xd[features.t, filter.t & split[,s] == 1])
				test.outcome <- Surv(x$outcome[filter.t&split[,s]==1,2],x$outcome[filter.t&split[,s]==1,1])
			
				nl <- try(readLines(paste(indir, "/Model_split",s,"-",type.name,"-",name,".csv",sep="")))
				if (class(nl) == "try-error") next else nl <- length(nl)
				if (nl == 0) next
				beta <- try(read.table(paste(indir, "/Model_split",s,"-",type.name,"-",name,".csv",sep=""), header=FALSE, sep=",", row.names=1))
				if (nl > 0) preds <- as.matrix(test.modules[,match(rownames(beta),x$fnames[features.t])]) %*% beta[,1] else preds <- rep(0, nrow(test.modules))
			
				clinical.beta <- read.table(paste("../Results/Prediction_Cox/Model_split",s,"-",type.name,"-2.Clinical.csv",sep=""), header=FALSE, sep=",", row.names=1)
				if (nrow(clinical.beta) > 0) clinical.preds <- as.matrix(test.modules[,match(rownames(clinical.beta),x$fnames[features.t])]) %*% clinical.beta[,1] else clinical.preds <- rep(0, nrow(test.modules))
			
				IDINRI <- IDI.INF(test.outcome,clinical.preds,preds,thres,npert=1)
				all.NRI[[type]][i,s] <- IDINRI$m2[1]
			}
		}
	}
	mean.rmmis <- function(x) mean(x[!is.na(x)])
	all.mean.NRI <- lapply(all.NRI,function(X) apply(X,1,mean.rmmis))
	result.NRI <- Reduce(cbind,all.mean.NRI)
	colnames(result.NRI) <- c("LUAD","KIRC","all")
	rownames(result.NRI) <- NAMES.95

	select <- grep("Clinical",NAMES.95)
	all.methods.NRI <- c(all.methods.NRI, list(result.NRI[match(NAMES.95[select],rownames(result.NRI)),]))
}

plot.box.strip <- function(data, var1, var2, var3, formula, label, col, ...) {
	data[,var2] <- factor(data[,var2],unique(data[,var2]))
	boxplot(formula, data=data, ...)
	stripchart(formula, method="jitter",
			jitter=0.1, vertical=T, pch=1,
			col=col, bg="bisque", add=TRUE, data=data)
	if (label) {
		k <- 0
		for (name in unique(data[,var2])) {
			Cindex.subgroup <- data[data[,var2]==name,var1]
			text(k+0.7,Cindex.subgroup, data[data[,var2]==name,var3],cex=0.7)
			k <- k+1
		}
	}
}

pdf("../Plots/Figure4.pdf", width=7, height=14) 
par(mfrow=c(2,1))

y.range <- c(-0.3,0.3)
for (j in 1:2) {
	nri <- all.methods.NRI[[j]]
	nri <- nri[-grep("2.Clinical",rownames(nri)),]
	nri.long <- as.vector(as.matrix(nri[,1:3]))
	nri.long <- as.data.frame(nri.long)
	nri.long$group <- rep(c("LUAD","KIRC","Pan-Cancer"),each=nrow(nri))
	colnames(nri.long) <- c("NRI","group")
	tmp <- c(nri.long$group)

	plot.box.strip(data=nri.long, var1="NRI", var2="group", var3="short.names",
		formula=NRI ~ group, label=F, col=2:4, ylab="NRI", outline=F, cex.axis=0.8, ylim=y.range,
		main=c("I-Boost-CV","I-Boost-Permutation")[j])
	abline(h=0)
}
dev.off()
