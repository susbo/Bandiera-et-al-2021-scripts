# Import list of 16867 detected genes
expressed = read.table("expressed.txt",stringsAsFactors = FALSE)$V1

# Read matrix file produced by deepTools createMatrix function
dat = read.table(gzfile("matrix.allTSS.gz"),sep="\t",skip=1,stringsAsFactors = FALSE)
names = sapply(strsplit(dat$V4,"\\."),function(x) head(x,n=1)) # Remove version number from gene names

# Set the boundary between the up, all and down groups; this is currently done manually based on the header information: "group_boundaries":[0,305,20650,21012]
names = names[306:20650] # Subset "names" to 'all' group to avoid double-counting any genes
counts.up = dat[1:305,7:dim(dat)[2]]
counts.all = dat[306:20650,7:dim(dat)[2]]
counts.down = dat[20651:21012,7:dim(dat)[2]]

# Calculate average of the two replicates
counts.up = (counts.up[,1:600]+counts.up[,601:1200])/2
counts.all = (counts.all[,1:600]+counts.all[,601:1200])/2
counts.down = (counts.down[,1:600]+counts.down[,601:1200])/2

# Exclude genes with NaN values
names = names[!is.nan(rowMeans(counts.all))]
counts.up = counts.up[!is.nan(rowMeans(counts.up)),]
counts.all = counts.all[!is.nan(rowMeans(counts.all)),]
counts.down = counts.down[!is.nan(rowMeans(counts.down)),]

# Remove genes on chrM
chrM = sapply(strsplit(dat$V4[dat$V1=="chrM"],"\\."),function(x) head(x,n=1))
remove = which(names %in% chrM)
counts.all = counts.all[-remove,]
names = names[-remove]

# Subset the counts matrix to only expressed genes
counts.expr = counts.all[names %in% expressed,]

N = 10000 # Set the number of bootstrap replicates
limits = c(dim(counts.up)[1],dim(counts.down)[1]) # 305 and 358; the size of up and down-regulated genes
numbers = sample(limits[1]:limits[2],N,replace=TRUE) # Decide how many genes to select in each bootstrap replicate (limited by up and down group)
all = c()

# Calculate mean profile for each bootstrap replicate
# This will take a long time if N is large
for (i in 1:N) {

	selection = sample(1:dim(counts.all)[1],numbers[i])
	profile = colMeans(counts.all[selection,])
	all = rbind(all,profile)
	
}

# Calculate lower and upper bound of CI for each position
perc = 0.001 # Significance level for CI
max = c(); min = c() # Upper and lower limit of CI
for (i in 1:600) {
	sorted = sort(all[,i])
	min[i] = sorted[length(sorted)*perc/2]
	max[i] = sorted[length(sorted)*(1-perc/2)]
}

# Make final figure
pdf("Figure3L.pdf",width=4,height=3)
par(oma = c(0,0,0,0) + 0.5)
par(mar = c(2.5,2.5,0.5,0.5), mgp=c(1.7,0.6,0))
plot(colMeans(counts.all),type="l",lwd=2,ylim=c(0,550),bty="n",yaxs="i",xaxs="i",xaxt="n",ylab="Norm. ATAC-seq",xlab="Position",col="#FFFFFF")
polygon(c(1:600, 600:1), c(min, rev(max)), density=NULL,col="#CCCCCC",border="#CCCCCC")
lines(colMeans(counts.up),col="#FF0000CC",lwd=2)
lines(colMeans(counts.down),col="#0000FFCC",lwd=2)
axis(1,at=c(1,101,201,300,400,500,600),label=c("-3kb","","","TSS","","","+3kb"),xpd=NA)
legend("topright",legend=c("Up","Down"),col=c("red","blue"),lwd=2,bty="n",seg.len=1)
dev.off()

