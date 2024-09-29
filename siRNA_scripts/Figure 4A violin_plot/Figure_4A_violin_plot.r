
library(ggplot2)
library(reshape2)

################################# input data
#################################
args=c("siRNA_rpkm.txt", "sample_sheet.txt")
dat=read.table(args[1], header=T)
sample_sheet = read.table(args[2], sep="\t", header=T)

################################# color setting
#################################
col =  c("gray","#96C7CA","#E74A8F","#6DBF52","#F39A26","#F06724","gray","#96C7CA","#696EA9","#E6E755")
dat.melt = melt(dat, id=c("Cluster"))
names(dat.melt) = c("Cluster", "Sample", "rpkm")
sample_sheet2 = sample_sheet[match(as.character(dat.melt$Sample), as.character(sample_sheet$Sample)),]
dat.melt$genotype = sample_sheet2$Condition
dat.melt$log2rpkm = log(dat.melt$rpkm +1, 2)
dat.melt$Sample = factor(dat.melt$Sample, levels=unique(as.character(sample_sheet$Sample)))
dat.melt$genotype = factor(dat.melt$genotype, levels=unique(as.character(sample_sheet$Condition)))

################################# boxplot
#################################
p1 = ggplot(dat.melt) +
	geom_violin(aes(x=Sample, y=log2rpkm, fill=genotype), scale="width") +
	geom_boxplot(aes(x=Sample, y=log2rpkm, fill=genotype), width = 0.3, outlier.shape=NA, lwd=0.2) +
	scale_fill_manual(values=col)+
	labs(y = 'log2(rpkm+1)') +
	labs(x ="") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	theme(legend.position="none") +
	theme(plot.title=element_text(hjust=0.5)) +
	theme(plot.margin=margin(1,1,1,1,"cm"))	

pdf("Figure_4A.pdf", width=4+(length(unique(dat.melt$Sample)))/5, height=4)
	plot(p1)
dev.off()
