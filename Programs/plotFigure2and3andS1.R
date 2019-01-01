# code has to be run with R 3.1.3

library(R.utils)
library(MASS)
library(ggplot2)
library(gridExtra)
library(reshape)
library(grid)
library(gtable)

features <- c("Clinical", "CNV", "Mutation", "miRNA", "Protein")
m <- length(features)

NAMES <- rep(NA, 2^m)
for (i in 0:(2^m-1)) {
        a <- unlist(strsplit(intToBin(i), split=""))
        a <- c(rep(0,m-length(a)), a)
        a <- rev(a)
        index <- c("")
        name <- c("")
        for (j in 1:m) {
                if (a[j] == "1") {
                        index <- paste(index,j+1,sep="")
                        name <- paste(name,features[j],sep=".")
                }
        }
        NAMES[i+1] <- paste(index,name,sep="")
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
ord <- order(sapply(NAMES.95,FUN=function(x)paste(unlist(strsplit(x,split="\\.")),collapse="7")))
NAMES.95 <- NAMES.95[ord]

n.model <- length(NAMES.95)

rm.na.mean <- function(x) {
	x <- x[!is.na(x)]
	c(mean(x))
}

index.plot <- function(group.names) {
	m <- length(group.names)
	composit.mat <- as.data.frame(matrix(0,0,2))
	colnames(composit.mat) <- c("Model","Type")
	all.types <- c("GeneExp","Module","Clinical","CNV",
						"Mutation","miRNA","Protein")
	k <- 1
	for (j in 1:m) {
		for (t in unlist(strsplit(group.names[j],split="\\."))) {
			composit.mat[k,1] <- j
			composit.mat[k,2] <- t
			k <- k + 1
		}
	}
	composit.mat <- composit.mat[length(composit.mat[,1]):1,]
	composit.mat$shape <- 4 ###if no GeneExp and Module in the model, the index panel somehow has a problematic format
	if (all(!grepl("GeneExp",composit.mat[,2])) & all(!grepl("Module",composit.mat[,2]))) {
		nr <- nrow(composit.mat)
		composit.mat[nr+1,1] <- 1
		composit.mat[nr+1,2] <- "Module"
		composit.mat[nr+1,3] <- NA
	}
	composit.mat[,1] <- m - composit.mat[,1] + 2
	composit.mat <- rbind(c(1,"",NA),composit.mat)
	composit.mat <- rbind(composit.mat,c(m+2,"",NA))
	class(composit.mat[,1]) <- "numeric"
	index.table <- ggplot(composit.mat, aes(y=Type,x=Model))+
	geom_point(shape=as.numeric(composit.mat$shape),size=1)+
	theme_bw()+
	scale_x_discrete("Model", limits = 1:(m+2))+
	scale_y_discrete("Type", limits = c("Clinical","GeneExp","Module","CNV","Mutation","miRNA","Protein")) +
	coord_flip() +
	theme(axis.text.y=element_blank()) +
	theme(axis.title.y=element_text(size=12,face="plain")) +
	theme(axis.text.x = element_text(angle=90,hjust=0.8,vjust=0),
        axis.title.x = element_blank())+
	theme(axis.ticks.y = element_blank())
	index.table
}

theme(axis.text.x = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))

cindex.plot <- function(cindex, pch, col, ylim) {
	m <- length(cindex)
	c <- as.data.frame(cbind(1:m,rev(cindex)))
	colnames(c) <- c("Model","Cindex")
	p <- ggplot(c, aes(y=Cindex,x=Model)) +
			geom_point(shape=pch,size=2.5,col=col) +
			scale_x_discrete("Model", breaks = 1:m) +
			scale_y_continuous("C-index Value", limits = ylim) +
			coord_flip() +
			theme_bw() +
			theme(axis.text.y=element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())+
			theme(plot.margin = unit(c(1,1,2.35,0),"lines")) +
			theme(axis.title.x = element_text(vjust = -2)) +
			theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())
	p
}

cindex.add.plot <- function(p, cindex, pch, col) {
	m <- length(cindex)
	c <- as.data.frame(cbind(1:m,rev(cindex)))
	colnames(c) <- c("Model","Cindex")
	p <- p + geom_point(data=c,aes(y=Cindex,x=Model),shape=pch,size=2.5,col=col)+coord_flip()
	p
}


add_legend <- function(...) {
	opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 2.2, 0), mar=c(0, 0, 0, 0), new=F)
	on.exit(par(opar))
	plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
	legend(...)
}

plot.cindex <- function(dir,label,name,outfile,n.split=rep(30,length(dir)),inc.clinical=TRUE) {

met.color <- c("green","magenta","blue","yellow")[match(label,c("LASSO","Elastic Net","I-Boost-CV","I-Boost-Permutation"))]

C.index.combined <- list()
for (jj in 1:length(name)) {
	C.index <- list()
	for (type in 1:3) {
		typename <- ifelse(type == 1, "LUAD", ifelse(type == 2, "KIRC", "all"))
		C.index.type <- matrix(NA, nrow=length(NAMES.95), ncol=n.split[jj])
		rownames(C.index.type) <- NAMES.95
		for (split in 1:n.split[jj]) {
			filename <- paste("../Results/",dir[jj],"summary_split",split,"-",typename,".csv",sep="")
			if (!file.exists(filename)) next
			res <- try(read.table(filename,sep=",",header=T))
			if (class(res) != "try-error" & nrow(res) > 0) {
				res[is.na(res[,"C.Index"]),"C.Index"] <- 0.5
				C.index.type[,split] <- res[match(NAMES.95,res[,1]),"C.Index"]
			}
		}
		C.index <- c(C.index, list(C.index.type))
	}
	C.index.combined <- c(C.index.combined,list(C.index))
}

check.NA <- function(x) {
	c(all(!is.na(x)))
}

cut.head <- function(x) {
	c(sapply(x, FUN=function(a) paste(unlist(strsplit(a,split="\\."))[-1],collapse=".")))
}

select <- 1:95
if (inc.clinical) select <- select[grep("Clinical",NAMES.95[select])] else select <- select[-grep("Clinical",NAMES.95[select])]

par(mar=c(2,12,4,1))

index <- index.plot(cut.head(NAMES.95[select])) 
index <- index + theme(plot.margin = unit(c(2.43,0.05,0.0,0.0), "line"))

all.avail <- list()
for (jj in 1:length(name)) {
	avail <- lapply(C.index.combined[[jj]],function(x)apply(x,2,function(y)!is.na(y)))
	all.avail <- c(all.avail,list(avail))
}
common.avail <- list()
for (j in 1:length(avail)) {
	j.all.avail <- lapply(all.avail,function(x)x[[j]])
	com <- Reduce('+',j.all.avail)
	com[com < length(name)] <- FALSE
	com[com == length(name)] <- TRUE
	com <- apply(com,2,as.logical)
	common.avail <- c(common.avail,list(com))
}

all.C <- list()
for (i in c(1,2,3)) {
	tmp2 <- as.data.frame(matrix(0,0,4))
	for (jj in 1:length(name)) {
		tmp <- as.data.frame(matrix(0,length(select),4))
		colnames(tmp) <- c("Cindex","Model","Method","Type")
		for (j in 1:length(select)) tmp$Cindex[j] <- mean(C.index.combined[[jj]][[i]][select[j],common.avail[[i]][select[j],]])
		tmp$Model <- length(select):1
		tmp$Method <- rep(paste("m",jj,sep=""),length(select))
		tmp$Type <- rep(as.character(i),length(select))
		tmp2 <- rbind(tmp2,tmp)
	}
	all.C <- c(all.C,list(tmp2))
}


k <- 1
all.C.append <- all.C[[k]]

min.nona <- function(x) min(x[!is.na(x)])
max.nona <- function(x) max(x[!is.na(x)])
min.all.C <- min.nona(c(all.C[[1]][,"Cindex"],all.C[[2]][,"Cindex"],all.C[[3]][,"Cindex"]))
max.all.C <- max.nona(c(all.C[[1]][,"Cindex"],all.C[[2]][,"Cindex"],all.C[[3]][,"Cindex"]))
all.C.append$Model <- all.C[[k]]$Model + 1
all.C.append <- rbind(c(NA,max(all.C[[k]]$Model)+1,NA,NA),all.C[[k]],c(NA,0,NA,NA))
points.LUAD <- ggplot(all.C.append, aes(y=Cindex,x=Model,color=Method,shape=Type)) +
		geom_point(size=2.5) +
		scale_colour_manual(values = met.color, breaks=paste("m",seq(1,length(name)),sep=""), labels=label) +
		scale_shape_manual(values = 16, breaks=paste("m",seq(1,length(name)),sep=""), labels=label) +
		scale_x_discrete("Model", limits=1:max(all.C[[k]]$Model)) +
		scale_y_continuous("C-index Value",limits=c(min.all.C,max.all.C)) +
		coord_flip() +
		theme_bw() +
		theme(axis.text.y=element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())+
		theme(plot.margin = unit(c(1.35,0.1,1.28,0.1),"lines")) +
		theme(axis.title.x = element_text(vjust = -1,face="plain",size=12)) +
		theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank()) +
		theme(legend.position="none") + 
		ggtitle("LUAD") + theme(plot.title = element_text(size=12, face="plain",vjust=2))
k <- 2
all.C.append <- all.C[[k]]
all.C.append$Model <- all.C[[k]]$Model + 1
all.C.append <- rbind(c(NA,max(all.C[[k]]$Model)+1,NA,NA),all.C[[k]],c(NA,0,NA,NA))
points.KIRC <- ggplot(all.C.append, aes(y=Cindex,x=Model,color=Method,shape=Type)) +
		geom_point(size=2.5) +
		scale_colour_manual(values = met.color, breaks=paste("m",seq(1,length(name)),sep=""), labels=label) +
		scale_shape_manual(values = 16, breaks=paste("m",seq(1,length(name)),sep=""), labels=label) +
		scale_x_discrete("Model", limits=1:max(all.C[[k]]$Model)) +
		scale_y_continuous("C-index Value",limits=c(min.all.C,max.all.C)) +
		coord_flip() +
		theme_bw() +
		theme(axis.text.y=element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())+
		theme(plot.margin = unit(c(1.35,0.1,1.28,0.1),"lines")) +
		theme(axis.title.x = element_text(vjust = -1,face="plain",size=12)) +
		theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank()) +
		theme(legend.position="none") + 
		ggtitle("KIRC") + theme(plot.title = element_text(size=12, face="plain",vjust=2))
k <- 3
all.C.append <- all.C[[k]]
all.C.append$Model <- all.C[[k]]$Model + 1
all.C.append <- rbind(c(NA,max(all.C[[k]]$Model)+1,NA,NA),all.C[[k]],c(NA,0,NA,NA))
points.all <- ggplot(all.C.append, aes(y=Cindex,x=Model,color=Method,shape=Type)) +
		geom_point(size=2.5) +
		scale_colour_manual(values = met.color, breaks=paste("m",seq(1,length(name)),sep=""), labels=label) +
		scale_shape_manual(values = 16, breaks=paste("m",seq(1,length(name)),sep=""), labels=label) +
		scale_x_discrete("Model", limits=1:max(all.C[[k]]$Model)) +
		scale_y_continuous("C-index Value",limits=c(min.all.C,max.all.C)) +
		coord_flip() +
		theme_bw() +
		theme(axis.text.y=element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())+
		theme(axis.title.x = element_text(vjust = -1,face="plain",size=12)) +
		theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank()) +
		theme(legend.box = "vertical", legend.position = unit(c(15.6,30.8),"lines"), legend.direction = "vertical") +
		guides(color = guide_legend(override.aes = list(shape = rep(16,length(name))), keyheight=0.7)) +
		theme(legend.title=element_blank()) +
		theme(legend.key = element_blank()) +
		theme(plot.margin = unit(c(1.35,0.1,1.28,.1), "line")) + 
		ggtitle("Pan-Cancer") + theme(plot.title = element_text(size=12, face="plain",vjust=2))

blank <- blank<-rectGrob(gp=gpar(col="white"))

ggsave(grid.arrange(index,points.LUAD,points.KIRC,points.all,blank,nrow=1,
		ncol=5,widths = unit(c(0.4,0.5,0.5,0.5,0.02),rep("null",5)),
		layout_matrix = matrix(c(1,2,3,4,5),1,5)),file=outfile,,height=7, width=14.8, units='in',dpi=300)
}

dir <- c("Prediction_lasso/","Prediction_en/")
label <- c("LASSO","Elastic Net")
name <- c("lasso","en")
plot.cindex(dir,label,name,"../Plots/Figure2.pdf")

dir <- c("Prediction_en/","Prediction_IBcv/","Prediction_IBperm/")
label <- c("Elastic Net","I-Boost-CV","I-Boost-Permutation")
name <- c("en","ibcv","ibperm")
plot.cindex(dir,label,name,"../Plots/Figure3.pdf")

dir <- c("Prediction_IBcv/")
label <- c("I-Boost-CV")
name <- c("ibcv")
plot.cindex(dir,label,name,"../Plots/FigureS1.pdf",inc.clinical=FALSE)
