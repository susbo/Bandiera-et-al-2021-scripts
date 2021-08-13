# Import cage data from step1
cage = read.csv("tmp/hg19.tss+fantom.cnt",sep=" ",header=FALSE,stringsAsFactors = FALSE)
colnames(cage) = c("EnsemblID","Reads")

# Import pre-defined gene sets to be compared to all genes
groups=system(paste("ls genesets/*.txt",sep=""),intern=TRUE)

all = list()
for (g in 1:length(groups)) {
	group = groups[g]
	genes = read.table(group,header=FALSE,stringsAsFactors = FALSE)$V1
	include = grepl(paste(genes, collapse="|"), cage$EnsemblID)
	
	all[[g]] = log2(cage$Reads[include]+1)
}

# Make figure
pdf("Figure3K.pdf",width=4,height=3)
par(mfrow=c(1,1))
par(oma = c(0,0,0,0) + 0.5)
par(mar = c(2,2,0,0)+0.5, mgp=c(1.5,0.6,0))
x=seq(0,100,100/(length(cage$Reads)-1));
plot(x,sort(log2(cage$Reads+1),decr=T),type="l",ylab="log2(CAGE)",xlab="% of genes",xaxs="i",lwd=2,col="grey",bty="n",xpd=NA,xlim=c(0,20))
x=seq(0,100,100/(length(all[[1]])-1));
lines(x,sort(all[[1]],decr=T),type="l",lwd=2,col="blue",xpd=NA)
x=seq(0,100,100/(length(all[[2]])-1));
lines(x,sort(all[[2]],decr=T),type="l",col="red",lwd=2,xpd=NA)
legend("topright",c("All","Down","Up"),col=c("grey","blue","red"),lwd=2,seg.len=1,bty="n")
dev.off()

# Calculate fraction of genes with CAGE evidence >0 (column 1)
# and use Fisher's exact test to determine significance (column 2)
for (i in 1:2) {
	cat(sum(all[[i]]>0)/length(all[[i]]),"\t")
	cat(fisher.test(matrix(c(sum(all[[i]]>0),sum(all[[i]]==0),sum(cage$Reads>0),sum(cage$Reads==0)),ncol=2))$p.value,"\n")
}
# 0.1546961 	3.090544e-05 
# 0.03606557 	0.0008888161 