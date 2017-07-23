lower.mq=45
## convention for highly expressed intragenic 'scrambled junction' because many scrambled junctions within a gene will be circles
min.intra.gene=100

require(data.table)
cosmic=data.table(read.csv("cancer_gene_census.csv")[,c(1,8)])
setnames(cosmic,"Gene.Symbol","gene1")

#par(mfrow=c(2,2))
require(data.table)
radius=10000
lower.num.reads=1

#####################
## static defition of KNIFEbad: junctions that are circular but agree with linear junctions at the boundary
## define KNIFEbad: bowtie2 -p 8 --trim5 138 --trim3 138  -f  -x hg19_junctions_reg -U hg19_junctions_scrambled.fa > KNIFEtest
## ie junctions that are "circular" but have reasonable homology at the boundary to linear junctions
## read in junctions called "bad" and process them

knifebad=data.table(read.delim("knife_badFJ", skip=1, header=F))
setnames(knifebad, "V1","ijunction")
knifebad [ ,junction:=gsub("([|])", ":", paste(ijunction))]
knifebad[,gene1:= as.character(lapply(strsplit(paste(knifebad$junction), split=":"), "[", 2))]
knifebad[,pos1:= as.character(lapply(strsplit(paste(knifebad$junction), split=":"), "[", 3))]

knifebad[,gene2:= as.character(lapply(strsplit(paste(knifebad$junction), split=":"), "[", 4))]
knifebad[,pos2:= as.character(lapply(strsplit(paste(knifebad$junction), split=":"), "[", 5))]

## similarly, define the mapping quality for each junction
## script reads in filtered data -- from "secondary_process.r" which relies on tcga_process.r
x=data.table(read.delim("num.poly_leq7_nGenesum_geq1.tab"))

sam=unique(x$sType)
x[,sampleid:=match(sType,sam)]

print ("generating query strings, mapping with bowtie, summary criteria has been done by wrapper script, files below")

## read in 


## test for multiple maps
threep=data.table(read.delim("test.bowtie.out.3prime", header=F)[,c(1,5)])
setnames(threep,"V1","junctionF")
setnames(threep,"V5","mq3")
fivep=data.table(read.delim("test.bowtie.out.5prime", header=F)[,c(1,5)])

setnames(fivep,"V1","junctionF")
setnames(fivep,"V5","mq5")

#mq=merge(threep,fivep,by="junctionF")
mq=cbind(fivep,threep[match(fivep$junctionF,threep$junctionF),mq3])
setnames(mq,"V2","mq3")
mq[,junctionF:=gsub("([_])", "-", paste(junctionF))]

x[,junctionF:=paste(gene1,gene2,pos1,pos2, "full",sep="-")]

x=cbind(x,mq[match(x$junctionF,mq$junctionF),list(mq5,mq3)])

x=x[((pos1-pos2)>radius & (numReads>min.intra.gene |paste(gene1)!=paste(gene2)))| (strand1!=strand2) | (chr1!=chr2)]

x [ ,anomRatio:=anomSum/numReads]
#normals=x[sample_type %like% "ormal"]
allfdr=x[!(sample_type %like% "ormal")]
## uses benjamini-hochberg corr_emp_p-- corr_emp_p is a corrected p value

i=0

## not used information on case ids
x1=unique(allfdr[sample_type %like% "ormal" & nGene1+nGene2>1,list(case_id,substr(sample_type,1,20),investigation,fuspair)][order(fuspair)]$case_id)
x2=unique(x[sample_type %like% "ormal" & nGene1+nGene2>1,list(case_id,substr(sample_type,1,20),investigation,fuspair)][order(fuspair)]$case_id)
y1=unique(allfdr[!(sample_type %like% "ormal") & nGene1+nGene2>1,list(case_id,substr(sample_type,1,20),investigation,fuspair)][order(fuspair)]$case_id)
y2=unique(x[!(sample_type %like% "ormal") & nGene1+nGene2>1,list(case_id,substr(sample_type,1,20),investigation,fuspair)][order(fuspair)]$case_id)

print ("bypassed use.corr and tests of different allfdrs") 

allfdr=x[(!(tolower(sample_type) %like% "normal")) & seq1!="" & seq2!="" ]

## remove gene pairs with evidence of gene homology defined as: knife far junctions align to knife badfj

allfdr[,g1g2:=paste(gene1,gene2,pos1,pos2)]
knifebad[,g1g2:=paste(gene1,gene2,pos1,pos2)]
allfdr=allfdr[is.na(match(allfdr$g1g2,knifebad$g1g2)),]

## define an approximate confidence interval (upper confidence interval) on the fraction of anomaly reads compared to reads; for numerical approximation, set the p-hat for anomaly reads to the point estimate plus epsilon=.01

epsilon=.01 #
allfdr[,AnomPhat:=(epsilon+sum(anomSum)/(sum(numReads)+sum(anomSum))), by=list(fuspair,pos1,pos2)]
## make sure that the AnomPhat (in this case) is the same across the listed variables
allfdr[,AnomSD:=sqrt((1/(sum(numReads+anomSum)))*AnomPhat*(1-AnomPhat)), by=list(fuspair,pos1,pos2)]
allfdr[,upperCI:=AnomPhat + 2*AnomSD]

## text parsing
 allfdr[,NONE.1:= NULL]
 allfdr[,NONE.2:= NULL]
 allfdr[,NONE.3:= NULL]
 allfdr[,NONE.4:= NULL]

## impose lower limit on the number of anomaly reads as a fraction of 'real' reads as well as number of reads observed

## statistical filter used below
allfdr=allfdr[(mq3+mq5)>lower.mq]
th.allfdr=allfdr[ (numReads>lower.num.reads)& (upperCI<.25 ) ]

## start of text processing, though not used for generating lists, only for screening out "normals" as definied by TCGA, which are not molecularly normal per their documentation, so can and do contain cancer cells

#normals=th.allfdr[sample_type %like% "ormal"]
## cancers will be either ## NEED TO ASSIGN CASEID AND INVESTIATION

cancers=th.allfdr[!(sample_type %like% "ormal")]
bodymap=th.allfdr[sinfo %like% "1_ERR0"]
bodymap[,case_id:=sType]
bodymap[,investigation:=sType]
cancers=rbind(cancers,bodymap)

## cancers[,nFus:=length(unique(case_id)),by=paste(gene1,gene2)]
## cancers[,nFusDonor:=length(unique(case_id)),by=paste(gene1)]
## cancers[,nFusAcceptor:=length(unique(case_id)),by=paste(gene2)]
cancers[,nFus:=length(unique(case_id)),by=list(gene1,gene2)]
cancers[,nFusDonor:=length(unique(case_id)),by=gene1]
cancers[,nFusAcceptor:=length(unique(case_id)),by=gene2]

cancers=merge(cancers,cosmic, by="gene1", all.x=T)
setnames(cosmic,"gene1","gene2")
cancers=merge(cancers,cosmic, by="gene2", all.x=T)

cancers[,is.cosmic:=(!is.na(Tumour.Types.Somatic..x)|!is.na(Tumour.Types.Somatic..y))]

cancers[,nCan:=length(unique(sType)), by=fuspair]

x[,info:=tolower(paste(g1prod,g2prod))]
cancers[,info:=tolower(paste(g1prod,g2prod))]

print (names(cancers))

## print and cull information about each cancer/sample type, not used for paper

for (cancer_type in unique(x$investigation)){
    cancer_sub  = cancers[ investigation %like% cancer_type] 
    onco=cancers[( info %like% "oncogene"| fuspair %like% "EGF"| info %like% "histone"| info %like% "cyclin" | info %like% "growth" | info %like% "kinase"  | info %like% "ras-"|info %like% "homeobox"  |info %like% "serine" ) & investigation %like% cancer_type] ## | info %like% "transcription") 
    onco_ids=unique(onco$case_id)
    unonco_ids=cancer_sub[ which(is.na(match(paste(cancer_sub$case_id) , paste(onco_ids))))]
    print(cancer_type)
    print (length(unique(unonco_ids$case_id))/ (length(unique(unonco_ids$case_id)) +length(onco_ids)))

    ## output fusions by sample time (not used)


    write.table(file=paste("fdr_post_filter_",cancer_type,".tab",sep=""),cancers[investigation %like% cancer_type & (chr1!=chr2 | strand1!=strand2| abs(pos1-pos2)>radius) ,list(gene1,gene2,abs(pos1-pos2),numReads,maxGeneExp.x, maxGeneExp.y,case_id,investigation,fuspair,info,chr1,chr2,pos1,pos2,strand1,strand2,anomSum,seq1,seq2,Tumour.Types.Somatic..x,Tumour.Types.Somatic..y)],quote=F,sep="\t", row.names=F)

    ## output all junctions by sample type (not used)

    write.table(file=paste("sbprimer_in_post_filter_",cancer_type,".tab",sep=""),unique(cancers[investigation %like% cancer_type & (chr1!=chr2 | strand1!=strand2| abs(pos1-pos2)>radius & (nFusDonor+ nFusAcceptor)>0) ,list(gene1,gene2,chr1,chr2,pos1,pos2)]),quote=F,sep="\t", row.names=F)

    ##summary statistics (not used)
    print (names(cancers))

    print (paste("cosmic ratio",sum(cancers[investigation %like% cancer_type]$is.cosmic)/dim(cancers[investigation %like% cancer_type])[1]))

    ## is the number of cosmic gene fusions 'equally distributed' among case ids
    ##> dim(cancers[is.cosmic==T,])
    ##[1] 364  74
    ##> dim(cancers)
    ##[1] 2455   74
    ##unique(cancers[is.cosmic==T,]$sType)
    ## 234
    ## for each sType, probablity of having no cosmic fusions is p^k where p is the probability that a fusion is a cosmic one, which is 595/L and k is the number of fusions per sample

    per.cancers=cancers[investigation %like% cancer_type,]
    per.cancers[, nFusPerSample:=length(unique(fuspair)),by=sType]
    L = length(unique(c(per.cancers$gene1,per.cancers$gene2)))
    n=8000
    ## because of 16K expressed genes and approximating as 2* probability of hitting a cosmic gene
    ## so expected number of cancers w/o
    p=595/n

    per.cancers[,p.notCosmic := (1-p)^nFusPerSample,by=sType]
    expected.nocosmic=sum(unique(per.cancers[,list(sType,p.notCosmic)])$p.notCosmic)

    total.cancers= length(unique(per.cancers$sType))
    cancers.w.cosmic=length(unique(per.cancers[is.cosmic==T,]$sType))

    print (paste("cancersw.cosmic",cancers.w.cosmic, "fractotal cancers", cancers.w.cosmic/total.cancers, "versus expected", total.cancers-expected.nocosmic))
}


write.table(file=paste("sbprimer_in_post_filter_all.tab",sep=""),unique(data.frame(cancers[  (chr1!=chr2 | strand1!=strand2| abs(pos1-pos2)>radius) ,list(gene1,gene2,chr1,chr2,pos1,pos2)])),quote=F,sep="\t", row.names=F)
#removes cosmic for (i in 43:57){cancers[,paste(names(cancers)[i]):=NULL]}

write.table(file=paste("all_fdr_post_filter_all.tab",sep=""), unique(data.frame(cancers[  (chr1!=chr2 | strand1!=strand2| abs(pos1-pos2)>radius) & !(is.na(case_id)) ,list(gene1,gene2,abs(pos1-pos2),numReads,maxGeneExp.x, maxGeneExp.y,case_id,investigation,fuspair,info,chr1,chr2,pos1,pos2,strand1,strand2,anomSum,seq1,seq2,corr_emp_p,sType,Tumour.Types.Somatic..x,Tumour.Types.Somatic..y,sType,upperCI,is.cosmic,g1prod,g2prod) ])),quote=F,sep="\t", row.names=F)
write.table(file=paste("Table_4_all_fdr_post_filter_all.tab",sep=""), unique(data.frame(cancers[  (chr1!=chr2 | strand1!=strand2| abs(pos1-pos2)>radius) & !(is.na(case_id)) ,list(gene1,gene2,abs(pos1-pos2),numReads,maxGeneExp.x, maxGeneExp.y,case_id,investigation,fuspair,info,chr1,chr2,pos1,pos2,strand1,strand2,anomSum,seq1,seq2,corr_emp_p,sType,Tumour.Types.Somatic..x,Tumour.Types.Somatic..y,sType,upperCI,is.cosmic,g1prod,g2prod) ])),quote=F,sep="\t", row.names=F)

## write table WITHOUT sequences and WITHOUT positions, for publication:
write.table(file=paste("Table_4_NO_SEQUENCES_all_fdr_post_filter_all.tab",sep=""), unique(data.frame(cancers[  (chr1!=chr2 | strand1!=strand2| abs(pos1-pos2)>radius) & !(is.na(case_id)) ,list(gene1,gene2,abs(pos1-pos2),numReads,maxGeneExp.x, maxGeneExp.y,case_id,investigation,fuspair,info,chr1,chr2,strand1,strand2,anomSum,corr_emp_p,sType,Tumour.Types.Somatic..x,Tumour.Types.Somatic..y,sType,upperCI,is.cosmic,g1prod,g2prod) ])),quote=F,sep="\t", row.names=F)


# summary statistic plots

## par(mfrow=c(2,2))
## for (j in c(2,4,5,6)){
##  plot(ecdf(x$emp_p))
## # plot(ecdf(normals[numReads==j]$emp_p), col="red",add=TRUE)
##  plot(ecdf(x[numReads==j]$emp_p),col="blue",add=TRUE)

## }

## for (j in c(2,4,5,6)){
##  plot(ecdf(x$corr_emp_p))
## # plot(ecdf(normals[numReads==j]$corr_emp_p), col="red",add=TRUE)
##  plot(ecdf(x[numReads==j]$corr_emp_p),col="blue",add=TRUE)

## }

