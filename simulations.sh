mkdir Programs > /dev/null 2> /dev/null
mkdir SimulationResults > /dev/null 2> /dev/null
mkdir SimulationSettings > /dev/null 2> /dev/null
mkdir SimulationData > /dev/null 2> /dev/null
cd Programs
	cat > inputdata.txt << EOF
library(data.table)
data <- data.frame(fread("../Data/TCGA_8cancer_rmmis.csv",header=TRUE,sep=","),row.names=1)
fnames <- rownames(data)

Clinical <- grep("Clinical_",fnames)
Clinical <- setdiff(Clinical, which(fnames=="Clinical_BRCA"))
CNV <- grep("CNV_",fnames)
miRNA <- grep("miRNA_",fnames)
Mutation <- grep("Mutation_",fnames)
Module <- grep("Module_",fnames)
Protein <- grep("Protein_",fnames)
GeneExp <- grep("Gene_",fnames)
data.types <- list(Clinical=Clinical,GeneExp=GeneExp,Module=Module,Protein=Protein,miRNA=miRNA,CNV=CNV,Mutation=Mutation)
data.types.names <- list("Clinical","Gene","Module","Protein","miRNA","CNV","Mutation")

features.name <- c("Clinical","GeneExp","Module","Protein","miRNA","Mutation","CNV")
features.pos <- list(Clinical,GeneExp,Module,Protein,miRNA,Mutation,CNV)

select.feature <- c(1,3,4,5,6,7)
pos <- unlist(features.pos[select.feature])
pos.type <- lapply(features.pos[select.feature],function(k) match(k,pos))

all.data <- data[pos,]
all.data <- apply(t(all.data),2,as.numeric)
all.data <- t(t(all.data) - apply(all.data,2,mean))
all.data <- t(t(all.data) / apply(all.data,2,function(x)ifelse(sd(x)>0,sd(x),1)))

EOF
	cat > output-summary.txt << EOF
all.corr <- cor(as.matrix(all.modules))
n.select <- sum(beta!=0)
n.true.select <- length(intersect(which(beta!=0),param[,1]))

preds <- as.numeric(all.modules[,param[,1]] %*% param[,2])
risks <- as.numeric(all.modules %*% beta)
corr.xb <- cor(preds,-risks)

true.beta <- rep(0,nrow(all.data))
true.beta[param[,1]] <- param[,2]
size.type <- unlist(lapply(pos.type,function(post) length(intersect(which(beta!=0),post))))
mse.type <- unlist(lapply(pos.type,function(post) sum((beta[post] + true.beta[post])^2) ))

cat(c(N,n.select,n.true.select,corr.xb,size.type,mse.type),sep=",")
cat("\\n")

}
sink()
EOF
	cat > sim_setting.R << EOF
library(data.table)

data <- data.frame(fread("../Data/TCGA_8cancer_rmmis.csv",sep=",",header=TRUE,data.table=FALSE),row.names=1)
fnames <- rownames(data)

Clinical <- grep("Clinical_",fnames)
Clinical <- setdiff(Clinical, which(fnames=="Clinical_BRCA"))
CNV <- grep("CNV_",fnames)
miRNA <- grep("miRNA_",fnames)
Mutation <- grep("Mutation_",fnames)
Module <- grep("Module_",fnames)
Protein <- grep("Protein_",fnames)
GeneExp <- grep("Gene_",fnames)

features.name <- c("Clinical","GeneExp","Module","Protein","miRNA","Mutation","CNV")
features.pos <- list(Clinical,GeneExp,Module,Protein,miRNA,Mutation,CNV)

select.feature <- c(1,3,4,5,6,7)
pos <- unlist(features.pos[select.feature])

all.modules <- data[pos,]
all.modules <- apply(t(all.modules),2,as.numeric)
all.modules <- t(t(all.modules) - apply(all.modules,2,mean))
all.modules <- t(t(all.modules) / apply(all.modules,2,function(x)ifelse(sd(x)>0,sd(x),1)))

library(stats)

signal <- 1.2
var.no <- c(5,20,20,10,10,10)
total.var <- sum(var.no)
dist <- list(c(0.7,0.2,0.05,0.05,0,0),c(0.5,0.5,0,0,0,0),c(0.5,0.1,0.1,0.1,0.1,0.1))

for (tv in total.var) {
for (setting in 1:3) {

res <- matrix(0,0,2)
count <- 0
for (j in select.feature) {
	count <- count + 1
	feature.signal <- signal * dist[[setting]][count]
	feature.novar <- var.no[count]
	if (feature.signal == 0) next

	matched.pos <- match(features.pos[[j]],pos)
	dissimilarity <- 1 - abs(cor(all.modules[,matched.pos]))
	distance <- as.dist(dissimilarity)

	hc <- hclust(distance, method="ave")
	memb <- cutree(hc, k = feature.novar)
	table(memb)
	cent <- NULL
	for(k in 1:feature.novar){
		cent <- cbind(cent, rowMeans(all.modules[,matched.pos][,memb == k,drop = FALSE]))
	}
	corr.centroid <- matrix(0,length(matched.pos),feature.novar)
	for (i in 1:feature.novar) {
		for (k in 1:length(matched.pos)) {
			corr.centroid[k,i] <- abs(cor(cent[,i],all.modules[,matched.pos[k]]))
		}
	}
	var.index <- apply(corr.centroid,2,which.max)
	beta <- sqrt(feature.signal / sum(var(all.modules[,matched.pos[var.index]])))
	res <- rbind(res,cbind(matched.pos[var.index],rep(beta,length(var.index))))

}
write.table(file=paste("../SimulationSettings/param-setting",setting,"tv",tv,".csv",sep=""),res,row.names=F,col.names=F,quote=F,sep=",")

}
}
EOF
	cat > gendata.R << EOF
library(data.table)
library(survival)
data <- data.frame(fread("../Data/TCGA_8cancer_rmmis.csv",header=TRUE,sep=","),row.names=1)
fnames <- rownames(data)

Clinical <- grep("Clinical_",fnames)
Clinical <- setdiff(Clinical, which(fnames=="Clinical_BRCA"))
CNV <- grep("CNV_",fnames)
miRNA <- grep("miRNA_",fnames)
Mutation <- grep("Mutation_",fnames)
Module <- grep("Module_",fnames)
Protein <- grep("Protein_",fnames)
GeneExp <- grep("Gene_",fnames)
data.types <- list(Clinical=Clinical,GeneExp=GeneExp,Module=Module,Protein=Protein,miRNA=miRNA,CNV=CNV,Mutation=Mutation)
data.types.names <- list("Clinical","Gene","Module","Protein","miRNA","CNV","Mutation")

features.name <- c("Clinical","GeneExp","Module","Protein","miRNA","Mutation","CNV")
features.pos <- list(Clinical,GeneExp,Module,Protein,miRNA,Mutation,CNV)

select.feature <- c(1,3,4,5,6,7)
pos <- unlist(features.pos[select.feature])
pos.type <- lapply(features.pos[select.feature],function(k) match(k,pos))

all.data <- data[pos,]
all.data <- apply(t(all.data),2,as.numeric)
all.data <- t(t(all.data) - apply(all.data,2,mean))
all.data <- t(t(all.data) / apply(all.data,2,function(x)ifelse(sd(x)>0,sd(x),1)))

tv <- 75
for (setting in 1:3) {
	param <- read.table(paste("../SimulationSettings/param-setting",setting,"tv",tv,".csv",sep=""),sep=",",header=F)
	for (N in 1:1000) {
		set.seed(N*269+1035)
		samp <- sample(1:nrow(all.data),500,replace=FALSE)
		all.modules <- all.data[samp,]
		Time <- sqrt(-2*log(runif(nrow(all.modules)))*exp(all.modules[,param[,1]] %*% param[,2]))
		C <- -log(runif(nrow(all.modules)))*1.7
		Y <- pmin(Time,C)
		CI <- as.numeric(Y==Time)
		all.outcome <- Surv(c(Y),CI)
		write.table(file=paste("../SimulationData/Data-setting",setting,"-",N,".csv",sep=""),cbind(Y,CI,samp),sep=",",row.names=FALSE,col.names=FALSE)
	}
}
EOF

for method in "Lasso" "ElasticNet" ; do
	for setting in 1 2 3; do
		programname=sim-$method-s$setting
		cat inputdata.txt > $programname.R
		cat >> $programname.R << EOF
library(glmnet)
library(IBoost)
library(survival)
setting <- $setting
tv <- 75

methodname <- c("$method")

dir <- "SimulationResults"
sink(paste("../",dir,"/setting",setting,"tv",tv,methodname,".csv",sep=""))
param <- read.table(paste("../SimulationSettings/param-setting",setting,"tv",tv,".csv",sep=""),sep=",",header=F)

for (N in 1:1000) {

set.seed(N*269+1035)
samp <- sample(1:nrow(all.data),500,replace=FALSE)
all.modules <- all.data[samp,]
Time <- sqrt(-2*log(runif(nrow(all.modules)))*exp(all.modules[,param[,1]] %*% param[,2]))
C <- -log(runif(nrow(all.modules)))*1.7
Y <- pmin(Time,C)
CI <- as.numeric(Y==Time)
all.outcome <- Surv(c(Y),CI)
model.family <- c("cox")

n.repeat.cv <- 5
if (methodname == "Lasso") alpha <- c(1) else alpha <- c(0.05, seq(0.1,1,by=0.1))
foldids <- list()
for (j in 1:n.repeat.cv) {
        set.seed(N * 31415 + 1938274 * j)
        foldid <- rep(cumsum(rep(1,5)), length=nrow(all.modules))
        foldid <- sample(foldid, replace=F)
        foldids <- c(foldids, list(foldid))
}
min.cvm <- 1e8
out.cvm <- c()
for (a in alpha) {
	cvm <- rep(0, 100)
	for (j in 1:n.repeat.cv) {
		cv.fit.en <- cv_glmnet(all.modules, all.outcome, family="cox", alpha=a, foldid=foldids[[j]], error="cox", naive=TRUE)
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
fit.en <- glmnet(all.modules, all.outcome, family=model.family, alpha=best.alpha, lambda=best.lambda)
beta <- fit.en\$beta[,length(best.lambda)]

EOF
		cat output-summary.txt >> $programname.R
	done
done

for method in "cv" "perm"; do
	for setting in 1 2 3 ; do
		programname=sim-Iboost-$method-s$setting
		cat inputdata.txt > $programname.R
		cat >> $programname.R << EOF
library(IBoost)
library(survival)
setting <- $setting
tv <- 75

methodname <- c("iboost-$method")

dir <- "SimulationResults"
sink(paste("../",dir,"/setting",setting,"tv",tv,methodname,".csv",sep=""), append=TRUE)
param <- read.table(paste("../SimulationSettings/param-setting",setting,"tv",tv,".csv",sep=""),sep=",",header=F)

for (N in 1:1000) {

set.seed(N*269+1035)
samp <- sample(1:nrow(all.data),500,replace=FALSE)
all.modules <- all.data[samp,]
Time <- sqrt(-2*log(runif(nrow(all.modules)))*exp(all.modules[,param[,1]] %*% param[,2]))
C <- -log(runif(nrow(all.modules)))*1.7
Y <- pmin(Time,C)
CI <- as.numeric(Y==Time)
all.outcome <- Surv(c(Y),CI)
model.family <- c("cox")

sub.data.types  <- list()
for (s in select.feature) sub.data.types <- c(sub.data.types,list(match(fnames[features.pos[[s]]],colnames(all.modules))))

alpha.series <- c(0.05,seq(0.1,1,by=0.1))
if (methodname == "iboost-perm") {
	v <- 0.1
	method <- "permute"
} else {
	v <- 0.2
	method <- "CV"
}
ib <- IBoost(all.modules,all.outcome,sub.data.types,alpha.series=alpha.series,v=v,method=method,seed=123987*N,iter.max=2000)
beta <- ib\$beta

EOF
		cat output-summary.txt >> $programname.R
	done
done
