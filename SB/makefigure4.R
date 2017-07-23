## makefigure4.R
library(data.table) 
df.chimera.in <- fread("../chimera_fusions_detected_inBodyMap_samples_common_toMACH.corrected.by.additionalcalcsfile.tab")

df.chimera.in[,gene1:=sapply(lapply(lapply(QueryName,strsplit,split="-"),"[[",1),"[",1)]
df.chimera.in[,gene2:=sapply(lapply(lapply(QueryName,strsplit,split="-"),"[[",1),"[",2)]
df.chimera.in[,nicenames:= paste(gene1,gene2,sep="-")]
df.chimera <- data.frame(nicenames=df.chimera.in$nicenames,counts=df.chimera.in$counts)



color.indices <- c(399,642,473,142,561,181,137,657,524,139,254,263,490,567,552,653,637,470,498,115,116,562,258)
## color.indices <- sample.int(657, n.types)
color.choices <- colors()[color.indices]

barplot.heights <- df.chimera$counts
barplot.names <- df.chimera$nicenames
barplot.colors <- c(rep(color.choices[1], length(df.chimera$counts)))
barplot.spaces <- rep(0.2,length(barplot.heights))
barplot.title <- "A. Number of Body Map Samples in Which\nChimerSeq Fusions Found"


## barplot for smachete:
sm.heights <- counts.for.genenames.in.called.cancers.and.also.in.bodymap.sbt
sm.names <- genenames.in.called.cancers.and.also.in.bodymap.sbt
sm.colors <- c(rep(color.choices[1], length(sm.names)))
sm.spaces <- rep(0.2,length(sm.heights))
sm.title <- "B. Number of Body Map\nSamples in Which\nsMACHETE Fusions Found"




## plot four plots

pdf(file="../figure4.four.plots.pdf", onefile=T, paper='USr')

layout(matrix(c(1,1,2,3,3,4,4,4,4,4), nrow=2, ncol=5, byrow = TRUE))

## barplot for ChimerSeq:

op <- par(mar= c(5.5,4,4,2) + 0.1)
## par(mgp=c(3,0.75,0))
## par(oma=c(15,0,0,0))
plot1 <- barplot(barplot.heights, names.arg=barplot.names, horiz=F, las=2, cex.names= 0.55, cex.axis = 0.7, space=barplot.spaces, col=barplot.colors, main=barplot.title, cex.main = 0.85, ylim=c(0,17))

## doing mtext location manually
mtext.location <- 5.8
mtext.padj <- 0.5


mtext(text=paste("FUSIONS"), side=1, line=5, outer=FALSE, cex=0.7, at=mtext.location, padj=mtext.padj)


## barplot for smachete:

op <- par(mar= c(5.5,4,4,2) + 0.1)
## par(mgp=c(3,0.75,0))
## par(oma=c(15,0,0,0))
plot2 <- barplot(sm.heights, names.arg=sm.names, horiz=F, las=2, cex.names= 0.55, cex.axis = 0.7, space=sm.spaces, col=sm.colors, main=sm.title, cex.main = 0.8, ylim=c(0,17))

## doing mtext location manually
mtext.location <- 2
mtext.padj <- 0.5


mtext(text=paste("FUSIONS"), side=1, line=5, outer=FALSE, cex=0.7, at=mtext.location, padj=mtext.padj)







## par(mar= c(5,4,4,2) + 0.1)
## par(mgp=c(3,1,0))
## par(oma=c(0,0,0,0))

thisplot.cex.row.and.colnames <- 0.65

plot.new()
## par(mar=c(1,1,1,1) + 0.1)
plotpars <- par("usr")
xvals <- plotpars[1] + (plotpars[2]-plotpars[1])*c(0,.35,.75,1)
yvals <- plotpars[3] + (plotpars[4]-plotpars[3])*c(0,.6,.78,1)
similar.gene.offset <- .07*(yvals[2]-yvals[1])
abline(h=yvals[2])
abline(v=xvals[3])
lines(x=c(xvals[2],xvals[4]), y=rep(yvals[3],2))
lines(x=rep(xvals[2],2), y=c(yvals[1],yvals[3]))
text(x=mean(xvals[2:3]),y=mean(yvals[3:4]),labels=chart.colnames[1], font=2, cex=thisplot.cex.row.and.colnames)
text(x=mean(xvals[3:4]),y=mean(yvals[3:4]),labels=chart.colnames[2], font=2, cex=thisplot.cex.row.and.colnames)
text(x=mean(xvals[1:2]),y=mean(yvals[2:3]),labels=chart.rownames[1], font=2, cex=thisplot.cex.row.and.colnames)
text(x=mean(xvals[1:2]),y=mean(yvals[1:2]),labels=chart.rownames[2], font=2, cex=thisplot.cex.row.and.colnames)
text(x=mean(xvals[2:3]),y=mean(yvals[2:3]),labels=known.gold.standard.amls.in.chimera, cex=thisplot.cex.row.and.colnames)
text(x=mean(xvals[3:4]),y=mean(yvals[2:3]),labels=known.gold.standard.amls.in.smachete, cex=thisplot.cex.row.and.colnames)
title(main="C. Performance of sMACHETE compared to ChimerSeq\nin LAML on gold standard LAML fusions", cex.main= 0.85)

x.similar.genes <- rep(mean(xvals[2:3]), length(aml.similar.gene.names))
y.similar.genes <- rev(yvals[2] - seq(from=similar.gene.offset, by= similar.gene.offset, length.out = length(aml.similar.gene.names)))
text(x=x.similar.genes,y=y.similar.genes,labels=aml.similar.gene.names, cex=0.6)


## op <- par(mar = c(5,4,4,12) + 0.1)
par(mar = c(5,4,4,12) + 0.1)

## data (da) is calculated in  ChimeraVsMachete.r

# par(mfrow=c(2,1))
# par(new=FALSE)	
# colnames(da)[3:4]=c("# detected ChimerSeq","# detected sMACHETE")
colnames(da)[3:4]=c("Ratio for ChimerSeq","Ratio for sMACHETE")
## od=order(da[,3]/da[,2])
#barplot(t(as.matrix(da[od,3:4]/da[od,2])),beside=T, names.arg=da[od,1], cex.names=.5, legend=T, main="TCGA types ordered by chimerSeq")
od=order(da[,4]/da[,2])
#tp53=fread("TP53_mut")
#dat=merge(da,tp53,by="tcgaprefix",all=T)
#dat=dat[which(!is.na(dat$mutFreq)),]


## do these manually:
x.upper <- 28
y.upper <- 4.5
## this next one should be an integer:
scaling.for.mut.freq <- 4


matrix.to.plot <- t(as.matrix(da[od,3:4]/da[od,2]))


bar.out <- barplot(matrix.to.plot,beside=T, names.arg=da[od,1], cex.names=.9, legend=FALSE, main="D. Fusion and TP53 Prevalence, by TCGA types ordered by ratio for sMACHETE", col=c("red","cyan"),ylab="# Fusions / # Samples", las=2, xlim=c(0,x.upper), ylim=c(0,y.upper), xlab=NULL)


x.vals.tp53 <- colMeans(bar.out)

tp53rates.in = fread("../TP53mutationrates.tab")
tp53rates <- tp53rates.in[tcgaprefix %in% da$tcgaprefix,][od]
tp53rates[,mutFreq:=mutFreq/100]

## par(new=TRUE)	

y.vals.tp53 <- scaling.for.mut.freq*tp53rates$mutFreq

points(x.vals.tp53,y.vals.tp53, xlim=c(0,x.upper), ylim=c(0,y.upper), pch=22, bg="black", xlab="", ylab="")



tick.marks.heights.rhs <- c(0:scaling.for.mut.freq)
tick.marks.labels.rhs <- (100/scaling.for.mut.freq)*tick.marks.heights.rhs
axis(4,at=tick.marks.heights.rhs, labels=tick.marks.labels.rhs, las=2, cex.axis=0.8)
mtext("Percent of Samples\nwith TP53 Mutation", side = 4, line = 3, cex = 0.8)


par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
these.usr.coords <- par("usr")
legend(x=.8, y= .3, c(rownames(matrix.to.plot), "Percent w/ TP53"), xpd = TRUE, inset = c(0, 0), bty = "o", fill = c("red","cyan","black"), cex = 0.6)


mtext(text="Fig 4", at=c(0), adj=0, line =-3, side=3, outer=TRUE)


dev.off()

index.prad.in.da <- which(da$tcgaprefix=="PRAD")
index.prad.in.od <- which(od==index.prad.in.da)

print("\nWithout PRAD, which had samples picked to include fusions, PEARSON:")

matrix.without.prad.for.plotting <- matrix.to.plot[,-index.prad.in.od]
tp53mutfreq.without.prad <- tp53rates[tcgaprefix!="PRAD",]$mutFreq
print("sMACHETE, PEARSON:")
print('cor.test(matrix.without.prad.for.plotting[2,],tp53mutfreq.without.prad, method="pearson"')
print(cor.test(matrix.without.prad.for.plotting[2,],tp53mutfreq.without.prad, method="pearson"))
## chimerseq
print("ChimerSeq, PEARSON:")
print('cor.test(matrix.without.prad.for.plotting[1,],tp53mutfreq.without.prad, method="pearson")')
print(cor.test(matrix.without.prad.for.plotting[1,],tp53mutfreq.without.prad, method="pearson"))

##
print("\n\n\n\n\nWith PRAD, PEARSON:")
## smachete:
print("sMACHETE, PEARSON:")
print('cor.test(matrix.to.plot[2,],tp53rates$mutFreq, method="pearson"')
print(cor.test(matrix.to.plot[2,],tp53rates$mutFreq, method="pearson"))
## chimerseq
print("ChimerSeq, PEARSON:")
print('cor.test(matrix.to.plot[1,],tp53rates$mutFreq, method="pearson")')
print(cor.test(matrix.to.plot[1,],tp53rates$mutFreq, method="pearson"))









print("\nWithout PRAD, which had samples picked to include fusions, SPEARMAN:")

print("sMACHETE, SPEARMAN:")
print('cor.test(matrix.without.prad.for.plotting[2,],tp53mutfreq.without.prad, method="spearman"')
print(cor.test(matrix.without.prad.for.plotting[2,],tp53mutfreq.without.prad, method="spearman"))
## chimerseq
print("ChimerSeq, SPEARMAN:")
print('cor.test(matrix.without.prad.for.plotting[1,],tp53mutfreq.without.prad, method="spearman")')
print(cor.test(matrix.without.prad.for.plotting[1,],tp53mutfreq.without.prad, method="spearman"))

##
print("\n\n\n\n\nWith PRAD, SPEARMAN:")
## smachete:
print("sMACHETE, SPEARMAN:")
print('cor.test(matrix.to.plot[2,],tp53rates$mutFreq, method="spearman"')
print(cor.test(matrix.to.plot[2,],tp53rates$mutFreq, method="spearman"))
## chimerseq
print("ChimerSeq, SPEARMAN:")
print('cor.test(matrix.to.plot[1,],tp53rates$mutFreq, method="spearman")')
print(cor.test(matrix.to.plot[1,],tp53rates$mutFreq, method="spearman"))



par(op)


## OLD:
## pdf(file="../bodymap.fusions.found.chimera.pdf", width =5, height =5)
## par(mar= c(6.5,4,4,2) + 0.1)
## par(mgp=c(3,0.75,0))
## par(oma=c(1,0,0,0))
## plot1 <- barplot(barplot.heights, names.arg=barplot.names, horiz=F, las=2, cex.names= 0.55, cex.axis = 0.7, space=barplot.spaces, col=barplot.colors, main=barplot.title, cex.main = 0.85)
## mtext(text=paste("CHIMERA"), side=1, line=-1, outer=TRUE, cex=0.7, at=mtext.location, padj=mtext.padj)
## dev.off()
## par(mar= c(5,4,4,2) + 0.1)
## par(mgp=c(3,1,0))
## par(oma=c(0,0,0,0))




