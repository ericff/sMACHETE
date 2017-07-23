# install.packages("calibrate", repos = "http://cran.cnr.berkeley.edu/")

##
require(data.table)
source("analysis.r") ## required
## machete ids used
mach.ids=fread(report.paths.with.meta.file)

setnames(mach.ids,"sample_id","BarcodeID")
## only using information from cancers with full "cMACHETE" ie those with bloom filters built and therefore only cancers in "full.final.list" so also restricting ids 
##
full.final.list[,tcgaprefix:=as.character(lapply(strsplit(paste(full.final.list$tcganame), split="-"), "[", 2))]

#chimeraseq fusions
chimera.info=fread("../ChimerDB3.0_ChimerSeq.txt")

## these are the samples that are in chimera and also were run through machete
## to subset tumors for figure

mach.ids.no.bodymap.no.normals <- mach.ids[(!is.na(mach.ids$sample_type) & mach.ids$sample_type!= "Solid Tissue Normal"),]

dim(mach.ids.no.bodymap.no.normals)

## these are chimerseq fusions which are found in one of the samples
## in the discovery set of cancers used for sMACHETE
chimera.and.mach = chimera.info[which(!(is.na(match(chimera.info$BarcodeID,mach.ids.no.bodymap.no.normals$BarcodeID))) & (!(is.na(match(chimera.info$Cancertype_or_disease,full.final.list$tcgaprefix)))))]

##  data used for table
#chimera.and.mach = chimera.info[which(!(is.na(match(chimera.info$BarcodeID,mach.ids$BarcodeID))) ]

chimeraBloom = fread(sbt.results.bodymap.chimera)

## use only fusions that are present in the samples screened in machete
## setting up names

setnames(chimeraBloom,"V1","QueryName")
chimera.and.mach[,QueryName:= paste(H_gene,"-",T_gene,"-",  H_position,"-",T_position,"_full_ChimeraScan_",sep="")]
chimera.and.mach[,genepair:= paste(H_gene,T_gene, H_position,T_position,sep=" ")]
chimera.and.mach[,genenames:= paste(H_gene,T_gene,sep="-")]

merged.chimera.mach=merge(full.final.list,chimera.and.mach,by="genenames",all=T)

chimeraBloom=chimeraBloom[which(!is.na(match(chimeraBloom$QueryName,chimera.and.mach$QueryName)))]

setnames(chimeraBloom,"V2","counts")

write.table(file="../chimera_fusions_detected_inBodyMap_samples_common_toMACH.tab",chimeraBloom[,list(QueryName,counts)], quote=F,sep="\t",row.names=FALSE)

write.table(file="../Table8_chimera_fusions_detected_inBodyMap_samples_common_toMACH.tab",chimeraBloom[,list(QueryName,counts)], quote=F,sep="\t",row.names=FALSE)

write.table(file="../chimera_fusions_detected_in_samples_common_toMACH.tab",chimera.and.mach, quote=F,sep="\t",row.names=FALSE)

write.table(file="../Table7_chimera_fusions_detected_in_samples_common_toMACH.tab",chimera.and.mach, quote=F,sep="\t",row.names=FALSE)


## output: chimera


#> length(unique(chimera.info$Fusion_pair))
#[1] 30001
#> length(unique(chimera.and.mach$Fusion_pair))
#[1] 1014
# includes splice variants: > length(unique(chimera.and.mach$QueryName))
#[1] 1289

#chimeraBloom[counts>0]
## figures used:

## Chimera per tumor vs. scalpel per tumor statistics

## Chimera per tumorvs. scalpel per tumor disease statistics
## run R script for machete results
source("analysis.r") 

## Do this again:
full.final.list[,tcgaprefix:=as.character(lapply(strsplit(paste(full.final.list$tcganame), split="-"), "[", 2))]


## generate figures
## output summaries for bodymap data by machete
## used for discovery via machete

## fusions_with_good_mq

## total fusions across all samples
cat("total fusions across all samples:\n")
cat(sum(unique(full.final.list[,list(genepair,maxMACount)])$maxMACount),"\n")
#>900
## total fusions across all bodymap
cat("total fusions across all bodymap:\n")
cat(sum(unique(full.final.list[tcganame %like% "BOD",list(genepair,maxMACount)])$maxMACount),"\n")
#[1] 9 each sample has 1 fusion, all except 1 are from KNIFE

chimeraBloom[counts>0]
## fusions detected in other tumors also in bodymap-- many are obvious FP
#QueryName counts     V3 V4

# no fusions detected by SCALPEL in tumors are found in bodymap.


## OFF BY 1
# only machete> xm[genepair %like% "SEPT"]

xm=merge(final.list,chimera.and.mach,by="genepair", all.x=T,all.y=T)

chimpub=fread("../ChimerDB3.0_ChimerPub.txt")
chimpub[,genenames:=paste(H_gene,T_gene,sep="-")]

## 
## problem bc LANCL-SEPT is not curated, and same for many MLL/ELL fusions that are actually described

merged.chimera.mach=merge(full.final.list,chimera.and.mach,by="genenames",all=T)
## merged.chimera.mach[,is.curated := !(is.na(match(chimera.and.mach$genenames,chimpub$genenames)))]

for (tcgatype in unique(merged.chimera.mach$Cancertype_or_disease)){
    write.table(file=paste("../Supp_table_1_chimera_seq_vs_scalpel_",tcgatype,".tab",sep=""),unique(merged.chimera.mach[(paste(Cancertype_or_disease,tcganame) %like% tcgatype),list(genenames,tcganame,Cancertype_or_disease)])[order(tcganame,Cancertype_or_disease)], quote=F,sep="\t", row.names=F)
}


chimera.and.mach.one.row.for.each.gene.pair.disease.combination <- unique(chimera.and.mach[,.(genenames,Cancertype_or_disease)])

chimera.and.mach.one.row.for.each.gene.pair.disease.combination[,NCperType:=length(unique(genenames)),by=Cancertype_or_disease]

nine.diseases.in.common <- unique(chimera.and.mach$Cancertype_or_disease)

full.final.list.in.diseases.common.to.chimseq.one.row.for.each.gene.pair.disease.combination <- unique(full.final.list[tcgaprefix %in% nine.diseases.in.common,.(genenames,tcgaprefix)])

allChimeraMach = unique(merged.chimera.mach[,list(genenames,tcganame,Cancertype_or_disease)])
allChimeraMach [,tcgaprefix:=as.character(lapply(strsplit(paste(allChimeraMach$tcganame), split="-"), "[", 2))]
	    
allChimeraMach[,NMperType:=length(unique(genenames)),by=tcganame]
allChimeraMach[,NCperType:=length(unique(genenames)),by=Cancertype_or_disease]

full.final.list.in.diseases.common.to.chimseq.one.row.for.each.gene.pair.disease.combination[,NMperType:= length(unique(genenames)),by=tcgaprefix]

## OLD: MTypeBar=unique(allChimeraMach[!is.na(tcganame),list(tcgaprefix,NMperType)])
## !is.na part seems not to be necessary, but keep in just in case:
MTypeBar=unique(full.final.list.in.diseases.common.to.chimseq.one.row.for.each.gene.pair.disease.combination[!is.na(tcgaprefix),.(tcgaprefix,NMperType)])

## OLD: CTypeBar=unique(allChimeraMach[!is.na(Cancertype_or_disease),list(Cancertype_or_disease,NCperType)])
## !is.na part seems not to be necessary, but keep in just in case:
CTypeBar=unique(chimera.and.mach.one.row.for.each.gene.pair.disease.combination[!is.na(Cancertype_or_disease),list(Cancertype_or_disease,NCperType)])

setnames(CTypeBar,"Cancertype_or_disease","tcgaprefix")
da = merge(CTypeBar,MTypeBar,by="tcgaprefix")
total.machete.counts[,tcgaprefix:=as.character(lapply(strsplit(paste(total.machete.counts$tcganame), split="-"), "[", 2))]
da=data.frame(merge(total.machete.counts[,list(TotalMachSamplesRun,tcgaprefix)],da,by="tcgaprefix"))


## Calculate how many samples chimerseq is run on; reasonable estimate
## but worth noting that there could be some samples in which
## no tumors are found
## so estimates of the ratio of fusions/sample are overestimates for
## chimerseq



## plot total chimera counts- NORMALIZED by sample numbers
## MOVING THE MAKING OF THE GRAPH WE ACTUALLY USE to makefigure4.R
## MOVING THE MAKING OF THE GRAPH WE ACTUALLY USE to makefigure4.R
## MOVING THE MAKING OF THE GRAPH WE ACTUALLY USE to makefigure4.R
## MOVING THE MAKING OF THE GRAPH WE ACTUALLY USE to makefigure4.R
# jpeg(file="total_mach_chimera_counts.jpeg")
pdf(file="../total_mach_chimera_counts.pdf", onefile=T, paper='USr')

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))

plot.new()

plot.new()

op <- par(mar = c(5,4,4,12) + 0.1)

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


bar.out <- barplot(matrix.to.plot,beside=T, names.arg=da[od,1], cex.names=.9, legend=FALSE, main="TCGA types ordered by ratio for sMACHETE", col=c("red","cyan"),ylab="#chimeras/#samples", las=2, xlim=c(0,x.upper), ylim=c(0,y.upper), xlab=NULL)


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
legend(x=.7, y= .3, c(rownames(matrix.to.plot), "Percent w/ TP53"), xpd = TRUE, inset = c(0, 0), bty = "o", fill = c("red","cyan","black"), cex = 0.6)


dev.off()
par(op)

## plot for AML

aml=data.frame(unique(merged.chimera.mach[(paste(Cancertype_or_disease,tcganame) %like% "AML"),list(genenames,tcganame,Cancertype_or_disease)]))

## rownames(aml)=aml[,1]
## aml=aml[,-1]
for (i in 1:2){
    aml[is.na(aml[,i]),i]=0
    aml[which(aml[,i]!=0 & !is.na(aml[,i])),i]=1
    aml[,i]=as.numeric(as.vector(aml[,i]))
}

nejm=fread("../NEJM_LAML_comparison.csv")
aml=cbind(aml,nejm[match(rownames(aml),genenames),]$is.Known.fusionEvent)
colnames(aml)[3]="is.annotated.aml.fusion"

## two panels
jpeg(file="../fig4_AML_NEJM.jpeg")
par(mfrow=c(1,2))
amlONE=aml[which(aml[,2]+aml[,1]==1),]
amlONE=as.matrix(amlONE[order(amlONE[,1]) , ])

## image(t(as.matrix(amlONE)), axes=F)

naml=dim(amlONE)[1]
plot(c(0,.01),c(0,naml-2),cex=.0000001, axes=F, ylab="", xlab="")
library(calibrate)			
textxy(X=c(rep(0,naml)), Y=1:naml,labs=rownames(amlONE), cex=1.1)
write.table(file="../FIG4_AML_names.txt",(amlONE),quote=F) ## used in ppt
dev.off()
