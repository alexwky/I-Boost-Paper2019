setwd("../Results")
all.dir <- list.dirs(recursive=FALSE)

for (dir in all.dir) {
	setwd(dir)
	all.files <- list.files()
	for (split in 1:30) {
		for (type in c("LUAD","KIRC","all")) {
			filename <- paste("summary_split",split,"-",type,sep="")
			sink(paste("summary_split",split,"-",type,".csv",sep=""))
			cat("Features,alpha,beta,Total Size,GeneExp Size,Module Size,Clinical Size,DNA Size,Mutation Size,miRNA Size,RPPA Size,Test N,C-Index")
			files <- all.files[grep(filename,all.files)]
			files <- setdiff(files,paste(filename,".csv",sep=""))
			for (f in files) {
				cat("\n")
				cat(readLines(f)[1])
			}
			sink()
		}
	}
	setwd("../")
}
