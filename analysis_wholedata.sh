mkdir Results > /dev/null 2> /dev/null
mkdir Results/WholeData > /dev/null 2> /dev/null
mkdir Programs > /dev/null 2> /dev/null
for method in "lasso" "en" "IBcv" "IBperm"; do
	outdir=Prediction_$method
        mkdir Results/WholeData/$outdir > /dev/null 2> /dev/null
	cat > Programs/WholeDataAnalysis-$method.R << EOF
library(Hmisc)
library(survC1)
library(glmnet)
library(R.utils)
library(survival)
library(data.table)
library(IBoost)

rawdata <- data.frame(fread("../Data/TCGA_8cancer_rmmis.csv",header=TRUE,sep=","),row.names=1)
x <- list(outcome=t(rawdata[1:2,]),xd=rawdata[-c(1:2),],fnames=rownames(rawdata)[-c(1,2)],snames=colnames(rawdata))

Clinical <- grep("Clinical_",x\$fnames)
Clinical <- setdiff(Clinical, which(x\$fnames=="Clinical_BRCA"))
CNV <- grep("CNV_",x\$fnames)
miRNA <- grep("miRNA_",x\$fnames)
Mutation <- grep("Mutation_",x\$fnames)
Module <- grep("Module_",x\$fnames)
Protein <- grep("Protein_",x\$fnames)
GeneExp <- grep("Gene_",x\$fnames)
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
	a <- unlist(strsplit(name, split="\\\\."))
	if (length(a) == 0) a <- c("")
	out <- paste(c(paste(index,a[1],sep=""),feature,a[-1]), collapse=".")
	c(out)
}
combine <- function(x, object) c(object, x)

NAMES.95 <- as.vector(c(sapply(NAMES, insert.string, index=0, feature="GeneExp"), sapply(NAMES, insert.string, index=1, feature="Module"), NAMES[-1]))
FEATURES.95 <- c(lapply(FEATURES, combine, object=GeneExp), lapply(FEATURES, combine, object=Module), FEATURES[-1])
n.model <- length(NAMES.95)

name<-NAMES.95[64]
features<-FEATURES.95[[64]]

for (type in c(1,2,3)) {
	type.name <- ifelse(type == 1, "LUAD", ifelse(type == 2, "KIRC", "all"))
	dir <- c("Results/WholeData/$outdir")
	out <- c()
	out <- c(out, name)

	if (type == 1) filter.t <- (x\$xd["Clinical_LUAD",] == 1)
	if (type == 2) filter.t <- (x\$xd["Clinical_KIRC",] == 1)
	if (type == 3) filter.t <- rep(TRUE,ncol(x\$xd))
	rm.features <- c()
	for (k in Clinical) if (length(unique(as.numeric(x\$xd[k,filter.t]))) == 1) rm.features <- c(rm.features, k)
	features.t <- setdiff(features, rm.features)
	all.modules <- t(x\$xd[features.t, filter.t])
	all.outcome <- Surv(x\$outcome[filter.t,2],x\$outcome[filter.t,1])
	train.modules <- t(x\$xd[features.t, filter.t])
	train.outcome <- Surv(x\$outcome[filter.t,2],x\$outcome[filter.t,1])

	train.modules <- t(t(train.modules) - apply(train.modules,2,mean))
	train.modules.sd <- apply(train.modules,2,function(x)ifelse(sd(x)>0,sd(x),1))
	train.modules <- t(t(train.modules) / train.modules.sd)

	if ("$method" == "lasso" | "$method" == "en") {
		n.repeat.cv <- 5
	        if ("$method" == "lasso") alpha <- c(1) else alpha <- c(0.05, seq(0.1,1,by=0.1))
	        foldids <- list()
	        for (j in 1:n.repeat.cv) {
	                set.seed(i * 31415 + 1938274 * j)
	                foldid <- rep(cumsum(rep(1,5)), length=nrow(train.modules))
	                foldid <- sample(foldid, replace=F)
	                foldids <- c(foldids, list(foldid))
	        }
	        min.cvm <- 1e8
	        out.cvm <- c()
	        for (a in alpha) {
	                cvm <- rep(0, 100)
	                for (j in 1:n.repeat.cv) {
	                        cv.fit.en <- cv_glmnet(train.modules, train.outcome, family="cox", alpha=a, foldid=foldids[[j]], error="cox", naive=TRUE)
	                        short <- min(length(cvm),length(cv.fit.en\$lambda))
	                        cvm <- cvm[1:short] + cv.fit.en\$cvm[1:short] / n.repeat.cv
	                }
	                out.cvm <- c(out.cvm, min(cvm))
	                if (min(cvm) < min.cvm) {
	                        best.alpha = a
	                        best.lambda = cv.fit.en\$lambda[1:which.min(cvm)]
	                        min.cvm <- min(cvm)
	                }
	        }
	        fit.en <- glmnet(train.modules, train.outcome, family="cox", alpha=best.alpha, lambda=best.lambda)
		beta <- fit.en\$beta[,length(best.lambda)]
		tuning <- c(best.alpha,tail(best.lambda,1))
	} else if ("$method" == "IBcv" | "$method" == "IBperm") {
		sub.data.types <- lapply(data.types,function(x)match(intersect(x,features.t),features.t))
		sub.data.types <- sub.data.types[as.numeric(lapply(sub.data.types,function(x)ifelse(length(x)==0,0,1)))==1]

		alpha.series <- c(0.05,seq(0.1,1,by=0.1))
		if ("$method" == "IBperm") {
			v <- 0.1
			method <- "permute"
		} else {
			v <- 0.1
			method <- "CV"
		}
		ib <- IBoost(train.modules,train.outcome,sub.data.types,method=method,iter.max=2000,v=v,alpha.series=alpha.series,seed=123987*i,permN=250)
		beta <- ib\$beta
		it <- ib\$iter.no
		tuning <- c(it,NA)
	}
	out <- c(out, tuning)
	selected <- x\$fnames[features.t[which(beta != 0)]]
	selected.pos <- which(beta != 0)
	beta.nonzero <- beta[which(beta != 0)]
	out <- c(out, length(beta.nonzero))

	if (length(grep("GeneExp", name)) > 0) GeneExpsize <- length(grep("Gene_", selected)) else GeneExpsize <- c(NA)
	if (length(grep("Module", name)) > 0) Modulesize <- length(grep("Module_", selected)) else Modulesize <- c(NA)
	if (length(grep("Clinical", name)) > 0) Clinicalsize <- length(grep("Clinical_", selected)) else Clinicalsize <- c(NA)
	if (length(grep("CNV", name)) > 0) DNAsize <- length(grep("CNV_", selected)) else DNAsize <- c(NA)
	if (length(grep("Mutation", name)) > 0) Mutationsize <- length(grep("Mutation_", selected)) else Mutationsize <- c(NA)
	if (length(grep("miRNA", name)) > 0) miRNAsize <- length(grep("miRNA_", selected)) else miRNAsize <- c(NA)
	if (length(grep("Protein", name)) > 0) RPPAsize <- length(grep("Protein_", selected)) else RPPAsize <- c(NA)
	out <- c(out, GeneExpsize, Modulesize, Clinicalsize, DNAsize, Mutationsize, miRNAsize, RPPAsize)

	outmodel <- cbind(selected[order(beta.nonzero)], beta.nonzero[order(beta.nonzero)])
	write.table(outmodel, file=paste("../",dir, "/Model-",type.name,"-",name,".csv",sep=""), col.names=F, row.names=F, sep=",", quote=F)

	gc()
}
EOF

done
