all.methods.res <- list()
methods <- c("Lasso","ElasticNet","iboost-cv","iboost-perm")

for (method in methods) {

all.res <- list()
for (setting in c(1,2,3)) {

res <- matrix(nrow=0,ncol=16)
dir <- "../SimulationResults"
filename <- paste("setting",setting,"tv75",method,sep="")
files <- list.files(dir)
files <- grep(filename, files, value=TRUE)
for (f in files) {
	if (length(readLines(paste(dir,"/",f,sep=""))) == 0) next
	res <- rbind(res, read.table(paste(dir,"/",f,sep=""),sep=","))
}
res <- res[!duplicated(res[,1]),]
all.res <- c(all.res, list(res))

}
all.methods.res <- c(all.methods.res, list(all.res))

}

names(all.methods.res) <- methods

all.methods.summary <- list()
for (method in methods) {
	all.res <- all.methods.res[[method]]
	all.methods.summary <- c(all.methods.summary, list(Reduce(rbind, lapply(all.res,function(x) apply(x,2,mean)))))
}
names(all.methods.summary) <- methods

set.col <- c(2,4:16)
types <- c("clinical","modules","protein","miRNA","mutation","cnv")
measure.names <- c("ns","corr",paste("size.",types,sep=""),paste("mse.",types,sep=""))
all.methods.summary <- lapply(all.methods.summary, function(X) {
				X <- X[,set.col];
				colnames(X) <- measure.names;
				rownames(X) <- paste("Setting",1:3,sep="");
				X})

require(plotrix)
library(ggplot2)
library(gridExtra)
library(grid)
library(plotrix)

reset <- function() {
        par(mfrow=c(1,1), oma=rep(0,4), mar=rep(0,4), new=TRUE)
        plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}

all.plots <- list()

setting <- 1:3
methods.label <- c("Lasso","Elastic Net","I-Boost-CV","I-Boost-Permutation")

measure <- "corr"
res <- matrix(0,length(methods),length(setting))
i <- 0
for (method in methods) {
	j <- 0
	i <- i + 1
	for (s in setting) {
		j <- j + 1
		res[i,j] <- all.methods.summary[[method]][paste("Setting",s,sep=""),measure]
	}
}
colnames(res) <- c(paste("Setting",setting,sep=""))
rownames(res) <- methods.label

long.res <- as.data.frame(c(res))
colnames(long.res) <- c("measure")
long.res$Setting <- rep(setting,each=length(methods))
long.res$method <- factor(rep(rownames(res),length(setting)),levels=rownames(res))

cols.met <- c("Elastic Net"="magenta","Lasso"="green","I-Boost-CV"="blue","I-Boost-Permutation"="yellow")
p1 <- ggplot(data=long.res, aes(x=Setting, y=measure, fill=method)) +
		geom_bar(stat="identity",position=position_dodge(),colour="black") +
		scale_fill_manual(values = cols.met, guide = guide_legend(title = NULL)) +
		theme(axis.title.x = element_text(vjust = 0,face="plain",size=12)) +
		theme(axis.text.x = element_text(vjust = 0,face="plain",size=12)) +
		theme(axis.title.y = element_text(vjust = 1,face="plain",size=12)) +
		theme(axis.text.y = element_text(vjust = 0.3,face="plain",size=12)) +
		theme(plot.margin = unit(c(0.2,5,0.2,5),"lines")) +
		coord_cartesian(ylim=c(0,1.0)) +
		scale_y_continuous("Risk Correlation") +
		theme(legend.position="right")

types <- c("clinical","modules","protein","miRNA","mutation","cnv")
longtypes <- c("Clinical","Modules","Protein","miRNA","Mutation","CNV")

all.plot.size <- list()
for (s in 1:3) {
	res <- matrix(0,length(methods),length(types))
	i <- 0
	for (method in methods) {
		j <- 0
		i <- i + 1
		for (t in paste("size.",types,sep="")) {
			j <- j + 1
			res[i,j] <- all.methods.summary[[method]][paste("Setting",s,sep=""),t]
		}
	}
	colnames(res) <- paste("size.",types,sep="")
	rownames(res) <- methods.label

	long.res <- as.data.frame(c(res))
	colnames(long.res) <- c("measure")
	long.res$type <- rep(longtypes,each=length(methods))
	long.res$method <- factor(rep(rownames(res),length(types)),levels=rownames(res))

	cols.type <- c("Clinical"="purple","Modules"="green", 
	            "Protein"="orange", "miRNA"="tomato", "Mutation"="yellow","CNV"="royalblue")
	p2 <- ggplot(data=long.res, aes(x=method, y=measure, fill=type)) +
			geom_bar(stat="identity",colour="black") +
			scale_fill_manual(values = cols.type, guide = guide_legend(title = NULL)) +
			theme(axis.title.x = element_text(vjust = 0,face="plain",size=12)) +
			theme(axis.text.x = element_text(angle=45, vjust = 1.0,hjust=1,size=8)) +
			theme(axis.title.y = element_text(vjust = 1,face="plain",size=12)) +
			theme(axis.text.y = element_text(vjust = 0.3,face="plain",size=12)) +
			theme(plot.margin = unit(c(0.2,0.2,0.2,0.2),"lines")) +
			coord_cartesian(ylim=c(0,240)) +
			scale_y_continuous("Model Size") +
			labs(title=paste("Setting",s)) +
			theme(plot.title = element_text(hjust = 0.5)) +
			theme(legend.position="right")
	all.plot.size <- c(all.plot.size, list(p2))
}

all.plot.mse <- list()
for (s in 1:3) {
	res <- matrix(0,length(methods),length(types))
	i <- 0
	for (method in methods) {
		j <- 0
		i <- i + 1
		for (t in paste("mse.",types,sep="")) {
			j <- j + 1
			res[i,j] <- all.methods.summary[[method]][paste("Setting",s,sep=""),t]
		}
	}
	colnames(res) <- paste("mse.",types,sep="")
	rownames(res) <- methods.label

	long.res <- as.data.frame(c(res))
	colnames(long.res) <- c("measure")
	long.res$type <- rep(longtypes,each=length(methods))
	long.res$method <- factor(rep(rownames(res),length(types)),levels=rownames(res))

	p3 <- ggplot(data=long.res, aes(x=method, y=measure, fill=type)) +
			geom_bar(stat="identity",colour="black") +
			scale_fill_manual(values = cols.type, guide = guide_legend(title = NULL)) +
			theme(axis.title.x = element_text(vjust = 0,face="plain",size=12)) +
			theme(axis.text.x = element_text(angle = 45, vjust=1.0 ,size=8,hjust=1)) +
			theme(axis.title.y = element_text(vjust = 1,face="plain",size=12)) +
			theme(axis.text.y = element_text(vjust = 0.3,face="plain",size=12)) +
			theme(plot.margin = unit(c(0.2,0.2,0.2,0.2),"lines")) +
			coord_cartesian(ylim=c(0,0.85)) +
			scale_y_continuous("Mean-squared Error") +
			labs(title=paste("Setting",s)) +
			theme(plot.title = element_text(hjust = 0.5)) +
			theme(legend.position="right")
	all.plot.mse <- c(all.plot.mse, list(p3))
}

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 


signal <- rbind(c(5,20,20,10,10,10),
				c("70%","20%","5%","5%","0%","0%"),
				c("50%","50%","0%","0%","0%","0%"),
				c("50%","10%","10%","10%","10%","10%"))
rownames(signal) <- c("No. of Variables",paste("Setting",c(1:3)))
colnames(signal) <- c("Clinical","Modules","Protein","miRNA","Mutation","CNV")
t1 <- ttheme_default(colhead=list(fg_params=list(fontface="plain")),
				rowhead=list(fg_params=list(fontface="plain")))

signal.table <- tableGrob(signal,theme=t1)

legend1 <- g_legend(p1)
legend2 <- g_legend(all.plot.size[[1]])
legend3 <- g_legend(all.plot.mse[[1]])

blank <- rectGrob(gp=gpar(col="white"))

ggsave(grid.arrange(main=textGrob("a",x=0.1,gp = gpar(fontsize = 12,fontface ="bold")),blank,
		p1+theme(legend.position="none"),legend1,
		main=textGrob("b",x=0.1,gp = gpar(fontsize = 12,fontface ="bold")),blank,
		all.plot.size[[1]]+theme(legend.position="none"),all.plot.size[[2]]+theme(legend.position="none"),all.plot.size[[3]]+theme(legend.position="none"),legend2,
		main=textGrob("c",x=0.1,gp = gpar(fontsize = 12,fontface ="bold")),blank,
		all.plot.mse[[1]]+theme(legend.position="none"),all.plot.mse[[2]]+theme(legend.position="none"),all.plot.mse[[3]]+theme(legend.position="none"),legend2,
		main=textGrob("d",x=0.1,gp = gpar(fontsize = 12,fontface ="bold")),blank,
		signal.table,
		nrow=8,
		ncol=4,widths = unit(c(0.3,0.3,0.3,0.2),rep("null",4)),heights=unit(c(0.03,0.5,0.03,0.5,0.03,0.5,0.03,0.17),"null"),
		layout_matrix = matrix(c(1,2,2,2,
						3,3,3,4,
						5,6,6,6,
						7,8,9,10,
						11,12,12,12,
						13,14,15,16,
						17,18,18,18,
						19,19,19,19),8,4,byrow=TRUE)),
		file="../Plots/Figure1.pdf",width=10,height=15,dpi=300)
