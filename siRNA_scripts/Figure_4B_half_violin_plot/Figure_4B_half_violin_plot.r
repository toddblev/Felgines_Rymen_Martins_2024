
library(ggplot2)
library(gghalves)

################################# input files
#################################
args = c("DESeq_FoldChange.txt", "subgroup_bed_files.txt")

fileName=read.table(args[1],sep="\t",header=T)
genotype_id=unique(as.character(fileName$genotype_id))
subset_sheet = read.table(args[2], sep="\t", header=T)

################################# combine different input files
#################################
data.rbind=NULL
for (i in 1:length(genotype_id)){
 print(i)
 genotype.one=genotype_id[i]
 # input data
 data.one=read.csv(as.character(fileName$csv.file.name)[i],sep="\t", header=T)
 data.one$geno=genotype_id[i]
 data.rbind=rbind(data.rbind,data.one)
 }

data.rbind = data.rbind[!is.na(data.rbind$log2FoldChange),]

################################# plot in the same figure with facet_wrap
#################################
dat.rbind.subAll = NULL
for (j in 1:nrow(subset_sheet)){
  print(j)
  subset.one.name = as.character(subset_sheet$genotypeID)[j]
  subset.one=read.table(as.character(subset_sheet$fileName)[j], sep="\t", header=F)
  subset.one.id=as.character(subset.one$V4)
  subset.one.number=nrow(subset.one)
  print(subset.one.number)
  data.sub = subset(data.rbind, Cluster %in% subset.one.id)
  #data.sub$sub=paste(subset.one.name,"-dep 24 nt siRNA clusters (n=", subset.one.number, ")", sep="")
  data.sub$sub=subset.one.name
  dat.rbind.subAll = rbind(dat.rbind.subAll, data.sub)
}

################################# FC distribution setting
#################################
FC.dat = dat.rbind.subAll
FC.dat$geno=factor(FC.dat$geno,levels=as.character(fileName$genotype_id))
FC.dat$sub=factor(FC.dat$sub, levels=as.character(subset_sheet$genotypeID))
col =  c("#F16623","#FBB040","#E7E755","#696EA9","#E94A91","#B3DDC5","#B3DDC5")

################################# FC distribution plot
#################################
facet_lab = c("All 24nt siRNA clusters (n=12,939)", "clsy12-dep 24nt siRNA clusters (n=6,617)", "clsy34-dep 24nt siRNA clusters (n=1,634)", "synergistic 24nt siRNA clusters (n=4,577)", "clsy1-dep 24nt siRNA clusters (n=1,947)", "clsy2-dep 24nt siRNA clusters (n=58)", "clsy3-dep 24nt siRNA clusters (n=802)", "clsy4-dep 24nt siRNA clusters (n=726)")
names(facet_lab) = c("all", "clsy12", "clsy34", "c12_c34_syner", "clsy1", "clsy2", "clsy3", "clsy4")
col =  rep(c("#B3DDC5","#B3DDC5","#E94A91","#696EA9","#E7E755","#FBB040","#F16623"),10)

p1=ggplot(FC.dat,aes(geno, log2FoldChange, fill=geno)) +
    geom_half_violin(data = FC.dat, aes(fill=geno, color=geno), position = position_nudge(x = 0), side="r", trim=TRUE, scale="area", draw_quantiles=c(0.25, 0.5, 0.75), show.legend=TRUE, width=1.5, alpha=0.8) +
    scale_fill_manual(values=col) +
    scale_color_manual(values=col) +
    xlab("")+
    ylab("log2(FC)")+
    coord_flip() +
    #scale_y_continuous(breaks = seq(floor(min(FC.dat$log2FoldChange)),0,2)) +
    scale_y_continuous(breaks = seq(-15,5,5)) +
    theme_bw() +
    facet_wrap(~sub, ncol=1, dir="v", labeller=as_labeller(facet_lab)) +
    theme(plot.margin = margin(1, 1, .5, .5, "cm"))+
    #labs(title=paste(as.character(subset_sheet$genotypeID)[j], "-dep 24 nt siRNA clusters (n=", subset.one.number, ")", sep="")) + 
    #theme(plot.title = element_text(hjust = 0.5, size=10)) +
    ylim(-15,5) +
    theme(strip.background = element_rect(colour="transparent", fill="transparent"))

pdf("Figure_4B.pdf",width=7,height=16)
  plot(p1)
dev.off()


