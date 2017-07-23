## goals of this script: process scalpel output
## define which fusions are cosmic and or kinase annotated
## compare it to chimeradb, after subsetting on samples scalpel has been run on
## chimeraKB (the curated fusions) are also processed
## plots are mde
## ##########################################################

## This is designed to be sourced by ChimeraVSMachete.r
## counts.file, name.to.sample.ids.dir, and report.paths.with.meta.file
## were defined in those:

library(data.table)
total.scalpel.counts=fread(counts.file)
setnames(total.scalpel.counts,"V1","tcga")
setnames(total.scalpel.counts,"V3","TotalMachSamplesRun")
setnames(total.scalpel.counts,"V4","TotalBloomSamples")
total.scalpel.counts[,tcganame:=paste("TCGA-", toupper(tcga), sep=""), by=tcga]
bloomlist=total.scalpel.counts[!is.na(TotalBloomSamples)]$tcga
per.tumor.dir=name.to.sample.ids.dir
counter=0

duplicate.per.sample=""

#initialize

total.scalpel.counts[,samples.w.fusion:=0]
total.scalpel.counts[,mean.fusions:=0]
total.scalpel.counts[,upper.CI.fusions:=0]
total.scalpel.counts[,lower.CI.fusions:=NA]
total.scalpel.counts[,kinase.fusions:=NA]
total.scalpel.counts[,cosmic.fusions:=NA]

total.scalpel.counts=data.frame(total.scalpel.counts)
mach.per.tumor=list.files(per.tumor.dir)[which(list.files(per.tumor.dir) %like% "scalpel.samples.only")]

#for (myf in (mach.per.tumor)){

for ( tumor.name in bloomlist[bloomlist %like% ""]){

    if (tumor.name == "body"){tumor.name="bodymap"}
    myf=paste("matrix.name.to.sample.ids.machete.samples.only.",tumor.name,".",querydate,".csv",sep="")
    dfh=data.table(read.csv(paste(per.tumor.dir,myf,sep="")))

    n.tumors= length(unique(dfh$sample.id))

    ## reformat
    lower.pan.cancer=0
    passed.list=paste(full.final.list[pan_cancer>lower.pan.cancer]$genepair,"_full",sep="")

    passed.list.reformat =gsub("([ ])", "-", paste(passed.list))
    passed.fus.pos=dfh[which(is.sequence.present>0 & !is.na(match(dfh$sequence.name,passed.list.reformat))),]

    passed.fus.pos[,gene1 := as.character(lapply(strsplit(paste(passed.fus.pos$sequence.name), split="-"), "[", 1))]
    passed.fus.pos[,gene2 := as.character(lapply(strsplit(paste(passed.fus.pos$sequence.name), split="-"), "[", 2))]
    passed.fus.pos[,genepair:=paste(gene1,gene2,sep="-")]

    ## get counts of which fusions have gene1 paired with multiple gene2s in the same sample.
    passed.fus.pos[,g1samesample:=length(unique(gene2)),by=list(sample.id,gene1)]
    passed.fus.pos[,g2samesample:=length(unique(gene1)),by=list(sample.id,gene2)] 

    duplicate.per.sample=c(duplicate.per.sample,passed.fus.pos[g1samesample+g2samesample>2]$genepair)


    ## get meta info 
    meta.info=unique(allfdr[,list(g1prod,g2prod,gene1,gene2,is.cosmic)])
    meta.info[,genepair:=paste(gene1,gene2,sep="-")]

    passed.fus.pos=merge(passed.fus.pos,meta.info,by="genepair")

    passed.fus.pos[,is.kinase:=(paste(g1prod,g2prod) %like% "kinase")]
    passed.fus.pos[,is.onco:=(is.cosmic|is.kinase)]

    passed.fus.pos$is.kinase <- as.integer(passed.fus.pos$is.kinase)
    passed.fus.pos$is.cosmic <- as.integer(passed.fus.pos$is.cosmic)
    passed.fus.pos$is.onco <- as.integer(passed.fus.pos$is.onco)

    ## passed.fus.pos[is.kinase==TRUE,is.kinase:=1]
    ## passed.fus.pos[is.kinase!=TRUE,is.kinase:=0]
    ## passed.fus.pos[is.cosmic==TRUE,is.cosmic:=1]
    ## passed.fus.pos[is.cosmic!=TRUE,is.cosmic:=0]
    ## passed.fus.pos[is.onco==TRUE,is.onco:=1]
    ## passed.fus.pos[is.onco!=TRUE,is.onco:=0]

    passed.fus.pos = unique(passed.fus.pos[,list(sample.id,genepair,is.kinase,is.cosmic,is.onco)])

    passed.fus.pos[,nFus:=length(unique(genepair)),by=sample.id]
    passed.fus.pos[,n.is.kinase:=sum(is.kinase),by=sample.id]
    passed.fus.pos[,n.is.cosmic:=sum(is.cosmic),by=sample.id]
    passed.fus.pos[,n.is.onco:=sum(is.onco),by=sample.id]
    ## p.onco is the probability of either cosmic or kinase
    p.onco=length(unique(passed.fus.pos[is.onco==TRUE]$genepair))/length(unique(passed.fus.pos[]$genepair))
    passed.fus.pos[,prob.zero.onco:=(1-p.onco)^nFus]


    ## counter and concatenate all tumor/sample ids
    if (counter>0){allg=rbind(allg,passed.fus.pos[,scalpelname:=tumor.name])}
    if (counter==0){allg=passed.fus.pos[,scalpelname:=tumor.name]}
    counter=counter+1

    ## print (passed.fus.pos)

    ## Note that all of these variables are defined by sample.id, so
    ## this is the same as getting a list of unique sample.id's and
    ## the corresponding variables:
    per.tumor.fusion.freq=unique(passed.fus.pos[,list(sample.id,nFus, n.is.kinase, n.is.cosmic,n.is.onco,prob.zero.onco)])

    observed.zero= sum(per.tumor.fusion.freq$n.is.onco==0)
    expected.zero=sum(per.tumor.fusion.freq$prob.zero.onco)
    print (paste("obs zero ",observed.zero,"exp. zero",expected.zero))
    ## if similar, implies that conditional on fusions being cosmic or onco, no bias


    ## total fusions, total tumors
    ## number of zeros= total.scalpel.counts$TotalMachSamplesRun

    print (myf)
    n.tumors.w.fusions=length(unique(passed.fus.pos$sample.id))

    ## check:
    stopifnot(dim(per.tumor.fusion.freq)[1]==n.tumors.w.fusions)

    num.samples.with.onco.fusions <- sum((per.tumor.fusion.freq$n.is.kinase+per.tumor.fusion.freq$n.is.cosmic)>0)
    print(paste0("Number of samples with onco-fusions: ", num.samples.with.onco.fusions, " samples out of ", n.tumors , " total samples, which is ", round((100*(num.samples.with.onco.fusions/n.tumors)),2), "%; Number of samples with any fusions: ",n.tumors.w.fusions,"; Fraction w/ any fusion: ", n.tumors.w.fusions/n.tumors,"; Fraction of samples with any fusion that have an onco-fusion: ", round((100*(num.samples.with.onco.fusions/n.tumors.w.fusions)),2),"%\n"))

    ## onco.num=sum((per.tumor.fusion.freq$n.is.kinase+per.tumor.fusion.freq$n.is.cosmic)>0)

    nrow = which(total.scalpel.counts$tcga==tumor.name)
    ## fixing slight issue, because it did not match before:
    if (tumor.name == "bodymap"){
        nrow = which(total.scalpel.counts$tcga=="body")
    }

    total.kinase.fusions.per.tumortype=sum(passed.fus.pos$is.kinase)
    total.cosmic.fusions.per.tumortype=sum(passed.fus.pos$is.cosmic)

    n.tumors.w.kinase.fusions=length(unique(passed.fus.pos[is.kinase>0,]$sample.id))
    n.tumors.w.cosmic.fusions=length(unique(passed.fus.pos[is.cosmic>0,]$sample.id))

    
    total.scalpel.counts[nrow,6] = n.tumors.w.fusions
    total.scalpel.counts[nrow,7] = num.samples.with.onco.fusions
    total.scalpel.counts[nrow,8] = sum(per.tumor.fusion.freq$nFus)/n.tumors.w.fusions
    total.scalpel.counts[nrow,9] = n.tumors.w.kinase.fusions
    total.scalpel.counts[nrow,10] = n.tumors.w.cosmic.fusions
    colnames(total.scalpel.counts)[6:10]=c("n.tumors.w.fusions","num.samples.with.onco.fusions","mean.per.tumors.among.wfus","n.tumors.w.kinase.fusions","n.tumors.w.cosmic.fusions")

    ## mean fusion among those w/ fusion

    tumor.fusions=unique(passed.fus.pos[,list(genepair,is.onco)])
    total.uniq.fusions=length(unique(tumor.fusions$genepair))
    total.onco.fusions=length(unique(tumor.fusions[is.onco==TRUE]$genepair))

    # Note that these are numbers of fusions (not numbers of samples):
    print (paste("total onco", total.onco.fusions, "total uniq", total.uniq.fusions, "frac",total.onco.fusions/total.uniq.fusions))

    ##
    print ("onco by chance\n\n")

}

##

## plotting

## save total.scalpel.counts as is for plotting before changes below
total.scalpel.counts.for.plotting <- total.scalpel.counts

plot.total.scalpel.counts.version.a=total.scalpel.counts.for.plotting[which(!is.na(total.scalpel.counts[,4])),]

pdf(file="../Fig_5_barplots_scatterplots.pdf")

## par(mfrow=c(2,2))
layout(matrix(c(1,1,3,3,4,2,2,5,5,6), 2, 5, byrow = TRUE))
## requires running script analysis.r before plot:
library(calibrate)
title.cex.val <- 1
plot(full.final.list$maxMAFreq, full.final.list$maxCompFreq, pch=22, bg="black", main="A. Pan-cancer Discovery vs. Test freq", xlab="A. Discovery set Freq", ylab="SBT frequency", cex.main=title.cex.val, xlim=c(0,0.9), ylim=c(0,0.9))
abline(0,1)

tp=full.final.list[tcganame %like% "OV"]
plot(tp$maxMAFreq, tp$maxCompFreq, pch=22, bg="black", main="B. Discovery vs. Test Freq, OV only",xlab="Discovery set Freq", ylab="SBT frequency", xlim=c(0,.2),ylim=c(0,.2), cex.main=title.cex.val)
high.req=tp[maxMAFreq+maxCompFreq>.05]
nlength=dim(high.req)[1]
## myadd=sort(runif(nlength,0,.01))
myadd.x = c(.028,.018,.028,.028)
myadd.y = c(0,-.005,0,0)
## textxy(high.req$maxMAFreq+myadd.x, high.req$maxCompFreq+myadd.y, paste(high.req$gene1,high.req$gene2))
text(x=high.req$maxMAFreq+myadd.x, y=high.req$maxCompFreq+myadd.y, labels=paste(high.req$gene1,high.req$gene2), cex=0.5)
abline(0,1)


plot.total.scalpel.counts.for.plot1.version.b =plot.total.scalpel.counts.version.a[order(plot.total.scalpel.counts.version.a[,1]),]

## colors:

color.indices <- c(399,642,473,142,561,181,137,657,524,139,254,263,490,567,552,653,637,470,498)
## myinds=c(3,6,7)
mycols.plot1=c("TotalMachSamplesRun", "n.tumors.w.fusions", "num.samples.with.onco.fusions")

#
## barplot(t(as.matrix(plot.total.scalpel.counts.for.plot1[,myinds])) , names=toupper(paste(plot.total.scalpel.counts.for.plot1[,1])), beside=T, las=2,ylab="Number of individuals", legend.text=toupper(colnames(plot.total.scalpel.counts.for.plot1[,myinds])), col=color.indices[1:length(myinds)] ,main=paste("C. Chimeras detected in discovery set"), cex.main=title.cex.val)


plot.total.scalpel.counts.for.plot1 = plot.total.scalpel.counts.for.plot1.version.b[mycols.plot1]
barnames.plot1 <- toupper(plot.total.scalpel.counts.for.plot1.version.b$tcga)

## par(mar = c(5,4,4,12) + 0.1)


barplot(t(as.matrix(plot.total.scalpel.counts.for.plot1)) , names=barnames.plot1, beside=T, las=2,ylab="Number of samples", legend.text=FALSE, col=color.indices[1:length(mycols.plot1)] ,main=paste("C. Fusions detected in discovery set"), cex.main=title.cex.val)

## plot.new()
## plot.window(0:1, 0:1) 
## these.usr.coords <- par("usr")

labels.barplot.1 <- c("Samples Analyzed", "Samples w/ Fusions", "Samples w/ Kinase or COSMIC Fusions")
legend(x=(.1*par("usr")[1] + .9*par("usr")[2]), y= mean(par("usr")[3:4]), legend=labels.barplot.1, xpd = NA, inset = c(0, 0), bty = "o", fill = color.indices[1:length(mycols.plot1)], cex = 0.7)
## legend(x=0, y= .5, legend=labels.barplot.1, bty="n", fill=color.indices[1:length(mycols.plot1)], cex = 0.6)
## legend(x=.5, y= .5, legend=labels.barplot.1, xpd = TRUE, inset = c(0, 0), bty="n", fill=color.indices[1:length(mycols.plot1)], cex = 0.6)

plot.new()


# 8 is mean per tumor given 9 and 10 are cosmic and kinsase

## plot high counts 
## myinds=c(3,7,8,9)
mycols.plot2=c("TotalMachSamplesRun", "num.samples.with.onco.fusions", "n.tumors.w.kinase.fusions", "n.tumors.w.cosmic.fusions")
#

plot.total.scalpel.counts.for.plot2.version.a=total.scalpel.counts.for.plotting[which(total.scalpel.counts[,3]>45),]
plot.total.scalpel.counts.for.plot2 <- plot.total.scalpel.counts.for.plot2.version.a[mycols.plot2]
barnames.plot2 <- toupper(plot.total.scalpel.counts.for.plot2.version.a$tcga)

## barplot(t(as.matrix(plot.total.scalpel.counts.for.plot2[,myinds])) , names=toupper(paste(plot.total.scalpel.counts.for.plot2[,1])), beside=T, las=2, ylab="Number of individuals", legend.text=toupper(colnames(plot.total.scalpel.counts.for.plot2[,myinds])), col=color.indices[1:length(myinds)],main=paste("D. Chimera features in highly sampled tumors"), cex.main=title.cex.val)
barplot(t(as.matrix(plot.total.scalpel.counts.for.plot2)) , names=barnames.plot2, beside=T, las=2, ylab="Number of samples", legend.text=FALSE, col=color.indices[1:length(mycols.plot2)],main=paste("D. Fusion features in highly sampled tumors"), cex.main=title.cex.val)


## plot.new()
## plot.window(0:1, 0:1) 
## these.usr.coords <- par("usr")

labels.barplot.2 <- c("Samples Analyzed",  "Samples w/ Kinase or COSMIC Fusions", "Samples w/ Kinase Fusions", "Samples w/ COSMIC Fusions")

legend(x=(.05*par("usr")[1] + .95*par("usr")[2]), y= mean(par("usr")[3:4]), legend=labels.barplot.2, xpd = NA, inset = c(0, 0), bty = "o", fill = color.indices[1:length(mycols.plot2)], cex = 0.7)

## slegend(x=0, y= .6, legend=labels.barplot.2, bty="n", fill=color.indices[1:length(mycols.plot2)], cex = 0.6)

## legend(x="right", legend=labels.barplot.2, bty="n", fill=color.indices[1:length(mycols.plot2)])

plot.new()



par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
these.usr.coords <- par("usr")

mtext(text="Fig 5", at=c(0), adj=0, line =-3, side=3, outer=TRUE)




dev.off()

########### process unique onco files -- for ONCO vs Machete

## make list consistent:

onco=data.table(read.csv("../onc2014406x2.csv",skip=1))
onco[,genenames:=paste(Gene_A,Gene_B,sep="-")]

#merged.onco.scalpel=merge(full.final.list,onco, by="genenames",all.x=T)

meta.ids=fread(report.paths.with.meta.file)
## contains onco file name/ aliquot id

onco=onco[!is.na(match(onco$Filename,meta.ids$aliquot_id)),]

same.tumors.as.onco=allg#[!(is.na(match(sample.id,onco.tcgas)))]
setnames(same.tumors.as.onco,"genepair","genenames")


 

total.scalpel.counts=data.table(total.scalpel.counts)
setnames(total.scalpel.counts,"n.tumors.w.fusions","n.tumors.w.chimeras")
setnames(total.scalpel.counts,"num.samples.with.onco.fusions","num.w.kinase.or.COSMIC.ann")
setnames(total.scalpel.counts,"n.tumors.w.kinase.fusions","total.samples.w.kinase.chimeras")
setnames(total.scalpel.counts,"n.tumors.w.cosmic.fusions","total.samples.w.cosmic.chimeras")
total.scalpel.counts[,mean.per.tumors.among.wfus:=NULL]
total.scalpel.counts[,cosmic.fusions:=NULL]
total.scalpel.counts[,mean.per.tumors.among.wfus:=NULL]
total.scalpel.counts=data.frame(total.scalpel.counts)
	
write.table(file="../Table_Y_total_SCALPEL_Counts.tab", total.scalpel.counts, sep="\t", quote=F)

## analysis of enrichment in chimeras:
singles=full.final.list[panBloomCount==1]
## exclude multiples:
singles=singles[is.na(match(singles$genenames, duplicate.per.sample))]

singles[,nGene1partner:=length(unique(gene2)),by=gene1]
singles[,nGene2partner:=length(unique(gene1)),by=gene2]

unique.gene1s= unique(singles[nGene1partner>1 ,list(tcganame,genenames,nGene1partner,panBloomCount,gene1,gene2)][order(nGene1partner)]$gene1)
unique.gene2s= unique(singles[nGene2partner>1 ,list(tcganame,genenames,nGene1partner,panBloomCount,gene1,gene2)][order(nGene1partner)]$gene2)
  
print (paste("number all private fusions where if in a single tumor, gene1 is paired w/ gene a and b, gene1-geneA and gene1-geneB are removed-- length of gene1 and gene2 that are unique as partners" ,length(unique.gene2s),"gene2partner>1 for chimeras","gene1partner>1 for chimeras",length(unique.gene1s)))

## RALA
Tnm=dim(unique(singles[,list(gene1,gene2)]))[1]
n=20000 # number of genes
c=2
nu2=c(1:16)
cat("RALA: n=", n,";Tnm=, ", Tnm, ";c=", c, "\n")
t=Tnm/sqrt(n)
print(1-exp(-(1/c)*(t^c))*sum((((1/c)*t^c)^(nu2-1))/(factorial(nu2-1))))

