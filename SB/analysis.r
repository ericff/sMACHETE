require(data.table)
lower.mq=45 # could be more stringent by increasing bound; is a parameter of course

## auxiliary generating tables
allfdr=data.table(read.delim("../all_fdr_post_filter_all.tab",stringsAsFactors=FALSE))

## admittedly, genepair is not the clearest variable name,
##  but don't want to change it now.
##  genepair INCLUDES the positions.
allfdr[,genepair:=paste(gene1,gene2,pos1,pos2)]
## num.allfdr is the number of unique samples w/ this pair

## If investigation is one of the bodymaps, e.g.
## ERR030872, then change to TCGA-BODY
allfdr[grep(pattern="^ERR",investigation),investigation:= "TCGA-BODY"]




## numMacheteDiscoveryByTCGA is the number of different samples
##  in which the specific combination of disease with genepair 
##    is found
##  in the particular investigation by the MACHETE in the discovery set
##  for that particular disease
## So note that we could have, e.g., some fusion found by MACHETE in PAAD
##   and then found only in OV by the SBT queries
##   In the latter case (OV), we would get numMacheteDiscoveryByTCGA
##   equal to NA, and we later assign NA values to be 0.
allfdr[,numMacheteDiscoveryByTCGA:=length(unique(sType)),by=list(investigation, genepair)]
## OLD: note: this is a TCGA specific number, bodymap excluded
## OLD: print ("numMacheteDiscoveryByTCGA is pan TCGA specific number, bodymap excluded")

num.in=data.table(read.delim("../num.poly_leq7_nGenesum_geq1.tab"))
setnames(num.in,"fuspair","genepair")
num.in[,num.unfiltered:=length(unique(sType)),by=list(investigation,genepair)]

# create unique ids for merging
num.in[,pair.sample:=paste(investigation,genepair)]
allfdr[,pair.sample:=paste(investigation,genepair)]

## test set-- machete inputs

total.machete.counts=fread(counts.file)
setnames(total.machete.counts,"V1","tcga")
setnames(total.machete.counts,"V3","TotalMachSamplesRun")
setnames(total.machete.counts,"V4","TotalBloomSamples")

total.machete.counts[,tcganame:=paste("TCGA-", toupper(tcga), sep=""), by=tcga]

## test for whether expression of 3' 5' predicts fusion and decoy
## test for indels 

nr=100
## read in all files 
dir=querydir
myfilelist=list.files(dir)

myfile=myfilelist[c(myfilelist %like% "summary.bt.results" )]

i=0
for (filein in myfile){
    type=strsplit(filein,"_")[[1]][1]
    ## for summary.results type=substr(filein,17,19)
    type=substr(filein,20,23)
    ##type=substr(filein,20,23)
    g=data.table(fread(paste(dir,filein,sep=""), header=F))
    setnames(g,"V1","junction")
    setnames(g,"V2","count")
    setnames(g,"V3","freq")

    setnames(g,"V4","MAcount")
    setnames(g,"V5","MAfreq")
    setnames(g,"V6","COMPcount")
    setnames(g,"V7","COMPfreq")
    if ((type) %like% "ov"){type="ov"} ## QUIRK OF TCGA NAMES

    g[,cType:=type]
    g=g[count>0]


    if (i==0){allg=g}
    if (i>0){

        allg=rbind(g,allg)}
    i=i+1
    print (head(allg))
}


## split out fields
allg[,junction:=gsub("([_])", "-", paste(junction))]
allg[,gene1 := as.character(lapply(strsplit(paste(allg$junction), split="-"), "[", 1))]
allg[,gene2 := as.character(lapply(strsplit(paste(allg$junction), split="-"), "[", 2))]
allg[,start1 := as.character(lapply(strsplit(paste(allg$junction), split="-"), "[", 3))]
allg[,start2 := as.character(lapply(strsplit(paste(allg$junction), split="-"), "[", 4))]

allg[,is.bad:= (junction %like% "mut")|(junction %like% "decoy")|(junction %like% "del")|(junction %like% "one"), by=junction]

allg[,typeinfo:=paste(gene1,gene2,start1,start2,cType,sep="_"),by=list(cType,junction)]

allg[,fivetypeinfo:=paste(gene1,start1,cType,sep="_"),by=list(cType,gene1,start1,cType)]
allg[,maxfive:=max(freq,na.rm=T),by=fivetypeinfo]
allg[,threetypeinfo:=paste(gene2,start2,cType,sep="_"),by=list(cType,gene1,start1,cType)]
allg[,maxthree:=max(freq,na.rm=T),by=threetypeinfo]

f=allg[junction %like% "five",list(freq,gene1,start1,cType,fivetypeinfo)]

allg[,fiveExp:= round(nr*f[match(allg$fivetypeinfo,f$fivetypeinfo)]$freq)/nr]

f=allg[junction %like% "three",list(freq,threetypeinfo)]
nr=100
allg[,threeExp:= round(nr*f[match(allg$threetypeinfo,f$threetypeinfo)]$freq)/nr]

#all decoys -- these are 'decoy' fusions
nr=1000
allg[,freq:= round(nr*freq)/nr]

tcga = data.table(read.delim("../all_fdr_post_filter_all.tab")[,])# was 1,2,9
tcga[,genepair:=paste(gene1,gene2,pos1,pos2,sep=" ")]
allg[,genepair:=paste(gene1,gene2,start1,start2, sep=" ")]

## use only gene pairs id in TCGA

mallg=unique(allg[(!is.na(match(allg$genepair,tcga$genepair))),])

## test for multiple maps
threep=data.table(read.delim("../test.bowtie.out.3prime", header=F)[,c(1,5)])
setnames(threep,"V1","junction")
setnames(threep,"V5","mq3")
fivep=data.table(read.delim("../test.bowtie.out.5prime", header=F)[,c(1,5)])

setnames(fivep,"V1","junction")
setnames(fivep,"V5","mq5")
mq=cbind(fivep,threep[match(fivep$junction,threep$junction),mq3])
setnames(mq,"V2","mq3")
mq[,junction:=gsub("([_])", "-", paste(junction))]

fallg=mallg[  is.bad==F & junction %like% "full" ]

## merge hist(mq$mq3+mq$mq5)
mqallg= merge(mallg,mq[junction %like% "full",],by="junction") ## merge hist(mq$mq3+mq$mq5)

fusions_with_good_mq=unique(mallg[(!is.na(match( paste(mallg$junction),paste(mq[mq3+mq5>lower.mq]$junction))))])
## or genepair is like high fusion ev

fusions_with_good_mq=rbind(fusions_with_good_mq, unique(mallg))
fusions_with_good_mq=fusions_with_good_mq[junction %like% "full"] ## many decoy fusions are included in queries but not used for this paper

cosmic_fusions=data.table(read.csv("../COSMIC_fusions.csv"))
setnames(cosmic_fusions,"Gene.Symbol","gene")
unique(fallg[which(!is.na(match(fallg$gene1,cosmic_fusions$gene)))]$gene1)
cat("length(unique(fallg[which(!is.na(match(fallg$gene1,cosmic_fusions$gene)))]$gene1)):\n")
print(length(unique(fallg[which(!is.na(match(fallg$gene1,cosmic_fusions$gene)))]$gene1)))
cat("\n")

    
## take max freq in case different variants have differentfrequencies
## OLD COMMENT: max frequencies are computed by gene pair, not by 'chimera' which is defined by both gene pair and splice site pair
fusions_with_good_mq[,maxFreq:=max(freq),by=list(genepair,cType)]
fusions_with_good_mq[,maxCount:=max(count),by=list(genepair,cType)]
fusions_with_good_mq[,maxCompFreq:=max(COMPfreq),by=list(genepair,cType)]
fusions_with_good_mq[,maxCompCount:=max(COMPcount),by=list(genepair,cType)]
fusions_with_good_mq[,maxMAFreq:=max(MAfreq),by=list(genepair,cType)]

fusions_with_good_mq[,maxMACount:=max(MAcount),by=list(genepair,cType)]

final= unique(fusions_with_good_mq  [order(-freq),list(genepair,gene1,gene2,cType,maxFreq,maxCount,maxCompFreq,maxCompCount,maxMAFreq,maxMACount)][order(gene1)])

final[,g1sum:=length(unique(gene2)),by=gene1]
final[,g2sum:=length(unique(gene1)),by=gene2]

## RENAME INPUT FILES SO THEY HAVE NAMES FOR TCGA; COMPUTE PER-INVESTIGATION RATES ON ALLFDR
## FIND PREDICTED INTERVALS TO IDENTIFY FALSE POSITIVES



final[,pan_cancer:=sum(maxCount), by=genepair]

## OLD: final[,pan_cancer:=sum(maxCount), by=genepair]
print (paste("dim final", dim(final)))

####### naming conventions for these cancers are different- not 4 letters
final[,tcganame:=paste("TCGA-",toupper(cType),sep="")]
final[cType=="ov.n", tcganame:="TCGA-OV",]
final[cType=="gbm.", tcganame:="TCGA-GBM",]
final[,pair.sample:=paste(tcganame,genepair)]

# merge  with allfdr table and raw inputs from machete
## merge by genepair-disease combination (not by genepair.no.site-disease combination)
## WAS GENEPAIR
final=merge(final,unique(allfdr[,list(pos1,pos2,chr1,chr2,is.cosmic,numMacheteDiscoveryByTCGA, strand1,strand2,pair.sample)]), by="pair.sample", all.x=T)

strandinfo= data.table ( allfdr[,list(gene1,strand1)] )
final[ , strand1 := strandinfo [ match(final$gene1,strandinfo$gene1 )]$strand1]

strandinfo= data.table ( allfdr[,list(gene2,strand2)] )
final[ , strand2 := strandinfo [ match(final$gene2,strandinfo$gene2 )]$strand2]


final=merge(final,unique(num.in[,list( pair.sample)]), by="pair.sample", all.x=T)

allfdr[,minSpliceDist:=min(abs(pos1-pos2)),by=genepair]		       
uniq.ann = unique(allfdr[,list(genepair,chr1,chr2,strand1,strand2,minSpliceDist,is.cosmic )])

final=merge(final, total.machete.counts[,list(tcganame,TotalMachSamplesRun)], by="tcganame")

final[,MAfreq:=numMacheteDiscoveryByTCGA/TotalMachSamplesRun]

final[is.na(MAfreq),MAfreq:=0]
## print ("final here")


n.machete.counts.excluding.bodymap <- sum(total.machete.counts$TotalMachSamplesRun[total.machete.counts$tcganame!="TCGA-BODY"])
 
## panMAfreq is the ratio of the number of different samples
##  in which the specific genepair is found
##  in ANY investigation by the MACHETE in the (cancer) discovery set
##  to the total number of samples in the (cancer) discovery set;
##  it does NOT count any times the genepair is found in a
##  bodymap data set by MACHETE.
 

final[tcganame!="TCGA-BODY",panMAfreq:=sum(numMacheteDiscoveryByTCGA, na.rm=T)/sum(TotalMachSamplesRun), by=genepair]


## OLD: final[,panMAfreq:=sum(numMacheteDiscoveryByTCGA,na.rm=T)/sum(total.machete.counts$TotalMachSamplesRun), by=genepair]

## This does not appear to be used later in the script.
final[,panDetectedTumorType:=length(unique(tcganame)), by=genepair]

# require that either the machete detects at least once or the fusion is in only one tumor, 

## if this were required, would require at least one tcga sample w/chimera 
#final=final[panMAfreq>0]

##
library(binom)
## lower CI for machete prediction
ep=.0001
final[is.na(numMacheteDiscoveryByTCGA),numMacheteDiscoveryByTCGA:=0]

final[,lower.MApred:=binom.confint(numMacheteDiscoveryByTCGA,TotalMachSamplesRun, method="exact", conf.level=1-ep)$lower]

## remove fusions considered statistically inconsistent

## remove positions, rejoin
final[,chr1:=NULL]
final[,chr2:=NULL]

final=merge(final,unique(allfdr[,list(genepair)]), by="genepair", all.x=T)

x=merge(final,uniq.ann[,list(genepair, minSpliceDist,is.cosmic,chr1,chr2)], by="genepair", all.x=T)
## upper CI for machete prediction

x[,upper.MA:=binom.confint(numMacheteDiscoveryByTCGA,TotalMachSamplesRun, method="exact", conf.level=(1-ep))$upper]
#x[,sumR:=sum(num.unfiltered)/sum(numMacheteDiscoveryByTCGA) ,by=genepair]

## numMacheteDiscoveryByTCGA is across all samples,

## Define False Negative rate in case in which "the number of samples with a junction detected by the SBT is greater than the number of samples with a junction detected by the MACHETE or KNIFE":


## if fewer discoveries by MACHETE (or same), define fn rate this
## way; also 
## define it in general, just so it is a variable in the data table
## in next few lines, we will define it if fewer discoveries by SBT
x[,FNrate:=( maxMACount - numMacheteDiscoveryByTCGA )/ maxMACount]

## OLD: x[,FNrate:=( maxMACount - numMacheteDiscoveryByTCGA )/ maxMACount, by=pair.sample] #and take maximum at the very end

## if fewer discoveries by SBT, define fn rate as NaN
x[maxMACount<numMacheteDiscoveryByTCGA,FNrate:= NaN]

## OLD: x[maxMACount<numMacheteDiscoveryByTCGA,FNrate:= NaN, by=pair.sample] #and take maximum at the very end


## positions of gene pair are reaassigned because of original merging and because some fusions are detected in only tumor 1 by machete but tumor 1 and 2 by bloom filter

x[,pos1 := as.numeric(as.character(lapply(strsplit(paste(x$genepair), split=" "), "[", 3)))]
x[,pos2 := as.numeric(as.character(lapply(strsplit(paste(x$genepair), split=" "), "[", 4)))]

x[,pan_cancer:=sum(maxCount), by=genepair]

## Already did this:
## x[,upper.MA:=binom.confint(numMacheteDiscoveryByTCGA,TotalMachSamplesRun, method="exact", conf.level=1-ep)$upper]

bad.fusions=c()
bad.fusions= unique(x[, list(maxFreq,maxCount, upper.MA,genepair, upper.MA-maxFreq)][upper.MA<maxFreq]$genepair)
##
if (min(x$numMacheteDiscoveryByTCGA)==0){
    no.mach= x[numMacheteDiscoveryByTCGA==0]
    no.mach[,total:=sum(TotalMachSamplesRun),by=genepair] ## 0 below by definition

    no.mach[,upper.MA:=binom.confint(0,total, method="exact", conf.level=(1-ep))$upper]

    ## note:

    bad.fusions=c(bad.fusions,unique(no.mach[upper.MA<maxFreq]$genepair))

}

######### add pan-cancer bad fusions:
x[,panMachDisc:=sum(numMacheteDiscoveryByTCGA),by=genepair]
x[,panBloomFreq:=sum(maxCount)/ sum(total.machete.counts$TotalBloomSamples,na.rm=T),by=genepair] ## total.machete.counts is misnomer-- this includes all samples 
x[,panBloomCount:=sum(maxCount),by=genepair]
N.tot=sum(total.machete.counts$TotalMachSamplesRun[total.machete.counts$tcganame!="TCGA-BODY"])
## print("ESTIMATED TOTAL TCGA SCREENED")

ep=.0001
x[,upper.panMA:=binom.confint(panMachDisc,N.tot, method="exact", conf.level=(1-ep))$upper]
bad.pan.fusions=unique(x[panBloomFreq>upper.panMA]$genepair)

bad.fusions=c(bad.pan.fusions,bad.fusions)

final.list= data.table(unique(data.frame(x[ pan_cancer>0&is.na(match(genepair,bad.fusions)), list(genepair,pan_cancer,gene1,gene2,strand1,strand2,chr1,chr2,pos1,pos2, abs(pos1-pos2),FNrate,maxMACount,genepair,upper.MA, maxFreq, maxCount,FNrate,numMacheteDiscoveryByTCGA, panBloomFreq,panBloomCount,upper.panMA)][order(-pan_cancer)])))

final.list[,nGene1partner:=length(unique(gene2)),by=gene1]
final.list[,nGene2partner:=length(unique(gene1)),by=gene2]

write.table(file="../pan_cancer_list.tab",final.list,quote=F,sep="\t", row.names=F)


## Note that pan_cancer>0 ensures that the SBT query finds the fusion
##  at least once.
## Note that this has a separate row for each splice variant.
full.final.list= x[pan_cancer>0&is.na(match(genepair,bad.fusions)), list(genepair,gene1,gene2,pan_cancer,strand1,strand2,chr1,chr2,pos1,pos2,abs(pos1-pos2),FNrate,maxMACount,upper.MA, maxFreq, maxCount,numMacheteDiscoveryByTCGA, tcganame, panBloomFreq,panBloomCount,upper.panMA, maxCompFreq, maxMAFreq)][order(-pan_cancer)]

full.final.list[,nGene1partner:=length(unique(gene2)),by=gene1]
full.final.list[,nGene2partner:=length(unique(gene1)),by=gene2]

names(full.final.list)[11]="AbsPos1Pos2Diff"
write.table(file="../pan_cancer_list_w_info.tab", full.final.list[],quote=F,sep="\t", row.names=F)

write.table(file="../Table_X_pan_cancer_list_w_info.tab", full.final.list[,list(gene1,gene2,pan_cancer,strand1,strand2,chr1,chr2,AbsPos1Pos2Diff,FNrate,tcganame, maxFreq,maxCompFreq, maxMAFreq)],quote=F,sep="\t", row.names=F)

## OLD, has positions:
## write.table(file="../Table_X_pan_cancer_list_w_info.tab", full.final.list[,list(genepair,gene1,gene2,pan_cancer,strand1,strand2,chr1,chr2,pos1,pos2,FNrate,tcganame, maxFreq,maxCompFreq, maxMAFreq)],quote=F,sep="\t", row.names=F)

## compute the false negative rate only among by isoform

x[,FNrate:=( maxMACount - numMacheteDiscoveryByTCGA )/maxMACount, by=pair.sample] #and take maximum at the very end

## if no discoveries or fewer by SBT (maxMACcount=0), define fn rate as 0
x[maxMACount<numMacheteDiscoveryByTCGA,FNrate:=NaN, by=pair.sample] #and take maximum at the very end


### per-tumor rates of kinase and fusion incidence using bloom queries

dfh=read.csv(sequences.present.laml.file)

rownames(dfh)=dfh$sequence.names
dfh=dfh[,-1]
## example use
dfh[which(row.names(dfh) %like% "BCR"),]
## reformat
passed.list=paste(full.final.list$genepair,"_full",sep="")

passed.list.reformat =gsub("([ ])", "-", paste(passed.list))

passed=dfh[which(!is.na(match(row.names(dfh),passed.list.reformat))),]

#cbind(rownames(passed[which(rowSums(passed)>0),]),rowSums(passed)>0

## chimeraseq

chimseq=data.table(read.csv("../ChimerDB3.0_ChimerSeq.csv"))
chimseq[Cancertype_or_disease %like% "AML"]$Fusion_pair

sort(unique(chimseq[Cancertype_or_disease %like% "" & Fusion_pair %like% "ACT"]$Fusion_pair))

## figure N:
 par(mfrow=c(2,2))

plot(final$maxMAFreq, final$maxCompFreq, pch=19, main="machete discovery vs test freq", xlab="machete freq", ylab="test frequency")

tp=final[tcganame %like% "AML"]
tp=final[tcganame %like% "OV"]
plot(tp$maxMAFreq, tp$maxCompFreq, pch=19, main="machete discovery vs test freq", xlab="machete freq", ylab="test frequency")

abline(0,1)

#unique.genepairs=unique(full.final.list[(gene1!=gene2)& (!(tcganame %like% "LAML" | tcganame %like% "GBM"|gene1== "ERG")),list(gene1,gene2)])
unique.genepairs=unique(full.final.list[pan_cancer<2,list(gene1,gene2)])
n.fusion.genes=2*dim(unique.genepairs)[1]
pair.nums=tapply(rep(1,n.fusion.genes), c(unique.genepairs$gene1,unique.genepairs$gene2),sum)

# an is coupon number
#wn is the time when the an-th distinct coupon has been sampled
## an is number of distinct fusions-- does this occur too early or too late?
## if there is a pref. to collecting same 'coupon'/'fusion', wn's observed value will be larger than its expected value
## approximation that 2 fusion pairs cannot exist in same gene
## for pair of genes w/o specifying  to 5' or 3' partner

wn=n.fusion.genes ## time-- what set of an's (unique coupons) are compatible
an=length(pair.nums) ## an= number of distinct fusions

pair.nums.fiveprime=tapply(rep(1,n.fusion.genes/2), c(unique.genepairs$gene1),sum)
pair.nums.threeprime=tapply(rep(1,n.fusion.genes/2), c(unique.genepairs$gene2),sum)
wn=n.fusion.genes/2 ## time-- what set of an's (unique coupons) are compatible

#an=length(pair.nums.fiveprime) ## an= number of distinct fusions
an=length(pair.nums.threeprime) ## an= number of distinct fusions

n=15000 # number genes-- large underestimate
# an is the number of coupons to collect
bn=n-an
nvec=bn:n
mu.n=n*(sum(1/nvec))
sig.nsq=n* sum( ( n-nvec)/(nvec^2))
## formulas 2 and 3
#http://www.jstor.org/stable/pdf/2239126.pdf
z=(wn-mu.n)/sqrt(sig.nsq)
print("(wn-mu.n)/sqrt(sig.nsq)")
print (z)

hist(pair.nums,breaks=20, ylim=c(0,100),main="recurrent fusion pairs")
# wn is larger than expected: [1] 23.99578
## repeat w/ 5' pair and 3' pair

onco=data.table(read.csv("../onc2014406x2.csv",skip=1))

onco[,genenames:=paste(Gene_A,Gene_B,sep="-")]
full.final.list[,genenames:=paste(gene1,gene2,sep="-")]
merged.onco.scalpel=merge(full.final.list,onco, by="genenames",all.x=T)
