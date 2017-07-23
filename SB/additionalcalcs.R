## additional calcs.R
## to be used after running everything in R_scripts_run.r

## assumes you are in SB directory

cat("Number of samples WITH bodymap:", sum(total.scalpel.counts$TotalMachSamplesRun),"\n")

cat("Number of samples WITHOUT bodymap:", sum(total.scalpel.counts$TotalMachSamplesRun[total.scalpel.counts$tcga!="body"]),"\n")

cat("dim(meta.ids):",dim(meta.ids),"\n")
print(table(meta.ids$sample_type, exclude=NULL))





## A large fraction (642) of the 997 gene fusions (847 unique gene fusions, as some occur multiple times) identified by sMACHETE in the TCGA tumor set are observed only once in our set of profiled tumors (i.e., they are private).
## 1006 and 849 and 760:

library(data.table)
## tablex = data.table(read.table("../pan_cancer_list_w_info.tab", sep ="\t", colClasses="character", stringsAsFactors=FALSE, header=TRUE))
tablex = fread("../pan_cancer_list_w_info.tab")

print("Number of UNIQUE gene fusions identified by sMACHETE in the TCGA tumor set, counting fusions as different if they have different positions:\n")
print(length(unique(tablex$genepair[tablex$tcganame!="TCGA-BODY"])))

print("Number of (nonunique) gene fusions identified by sMACHETE in the TCGA tumor set:\n")
print(dim(tablex[tablex$tcganame!="TCGA-BODY",])[1])


tablex[,genenames:=paste(gene1,gene2,sep="-")]

print("Number of UNIQUE gene fusions identified by sMACHETE in the TCGA tumor set, NOT counting fusions as different if they have different positions:\n")
print(length(unique(tablex$genenames[tablex$tcganame!="TCGA-BODY"])))










print(dim(tablex[tablex$tcganame=="TCGA-BODY",])[1])


## A large fraction (642) of the 997 gene fusions identified by sMACHETE in the TCGA tumor set are observed only once in our set of profiled tumors (i.e., they are private).

## myfile=myfilelist[c(myfilelist %like% "summary.bt.results" )]
##     setnames(g,"V2","count")
## So "count" is the number of times it's found in the SBT

## x[,pan_cancer:=sum(maxCount), by=genepair]

tablex$genenames <- paste(tablex$gene1, tablex$gene2, sep="-")

tablex[,pancancer.by.genenames:=sum(maxCount),by=genenames]

private.fusions <- tablex[pancancer.by.genenames==1,]
## OLD: print("length(unique(private.fusions$genenames))")
## OLD: print(length(unique(private.fusions$genenames)))
## ## OLD: 616
## OLD: Tnm= length(unique(private.fusions$genenames))

duplicate.genenames <- tablex$genenames[duplicated(tablex$genenames)]


single.genenames <- setdiff(unique(tablex$genenames),duplicate.genenames)
length(single.genenames)
## 639
## But these could have pan_cancer > 1

indices.single.genenames <- which(tablex$genenames %in% single.genenames)
table(tablex$pan_cancer[indices.single.genenames], exclude=NULL)
table(full.final.list$pan_cancer, exclude=NULL)


dim(unique(singles[,list(gene1,gene2)]))
dim(singles[,list(gene1,gene2)])

singles[,list(gene1,gene2)][1:10,]

## tablex[tablex$pan_cancer[indices.single.genenames]>1,1:10]

setdiff(singles$genenames, single.genenames)


## tablex[tablex$gene1=="ARFGEF2",]
## singles[singles$gene1=="ARFGEF2",]

## mallg[gene1=="ARFGEF2" & gene2=="NCOA3",]
## allg[genenames=="ARFGEF2-NCOA3",]


## Redoing this to account for possibility of this:
##   If two splice variants occur in the same sample, and only
##   in that same sample, we would still call that a private fusion.

## A large fraction (642) of the 997 gene fusions identified by sMACHETE in the TCGA tumor set are observed only once in our set of profiled tumors (i.e., they are private).

## Create data table with one row for every
## (splice variant)-sample combination

## First redo stuff from analysis.r:

## .add as suffix to ensure variables don't get mixed up with
## previous calculations

total.scalpel.counts.add=fread(counts.file)
setnames(total.scalpel.counts.add,"V1","tcga")
setnames(total.scalpel.counts.add,"V3","TotalMachSamplesRun")
setnames(total.scalpel.counts.add,"V4","TotalBloomSamples")
total.scalpel.counts.add[,tcganame:=paste("TCGA-", toupper(tcga), sep=""), by=tcga]
bloomlist.add=total.scalpel.counts.add[!is.na(TotalBloomSamples)]$tcga
per.tumor.dir.add=name.to.sample.ids.dir

matrix.files.add <- vector("list", length=length(bloomlist.add))
ii <- 0
for (tumor.name in bloomlist.add){
    ii <- ii + 1
    cat("\rWorking on ",ii)
    if (tumor.name == "body"){tumor.name="bodymap"}
    myf=paste("matrix.name.to.sample.ids.all.",tumor.name,".",querydate,".csv",sep="")
    matrix.files.add[[ii]]=data.table(read.csv(paste(per.tumor.dir.add,myf,sep=""), stringsAsFactors=FALSE))
}

## Add in investigation while doing next part; will need this elsewhere.
bloomlist.investigations <- paste("TCGA-", toupper(bloomlist.add),sep="")


mf.add.only.present.sequences <- list()
for (ii in 1:length(matrix.files.add)){
    temp.matrix <- matrix.files.add[[ii]][is.sequence.present==1,]
    mf.add.only.present.sequences[[ii]] <- data.frame(temp.matrix, investigation=rep(bloomlist.investigations[ii], times=dim(temp.matrix)[1]))
}



## Combine into data frame:

df.sv.sample <- mf.add.only.present.sequences[[1]]
for (ii in 2:length(mf.add.only.present.sequences)){
    df.sv.sample <- rbind(df.sv.sample, mf.add.only.present.sequences[[ii]])
}

## change to data table:
df.sv.sample <- data.table(df.sv.sample)

## as check
df.sv.sample[,is.bad:= (sequence.name %like% "mut")|(sequence.name %like% "decoy")|(sequence.name %like% "del")|(sequence.name %like% "one"), by=sequence.name]

stopifnot(sum(df.sv.sample[,is.bad])==0)

df.sv.sample$is.bad <- NULL

## collapse that to a data frame with one row for every
##  gene1-gene2-sample combination
## Note that this may not be precisely the correct way to
##  get gene1 and gene2 in df.sv.sample; see below re this.

df.sv.sample[,gene1:= sapply(lapply(lapply(sequence.name, FUN=strsplit, split="-"),"[[",1), "[",1)]
df.sv.sample[,gene2:= sapply(lapply(lapply(sequence.name, FUN=strsplit, split="-"),"[[",1), "[",2)]
df.sv.sample[,gene1.gene2.sample.combo:= paste(gene1,gene2,sample.id,sep=" ")]
df.sv.sample[,gene1.gene2:= paste(gene1,gene2,sep="-")]

## In df.sv.sample, there are some pairs of genes with a dash
##  in them, e.g. this sequence has two such genes:
##  LOC100132832-DTX2P1-UPK3BP1-PMS2P11


## Restrict to sequences that have a good mq:
good.mq.junctions <- gsub(pattern="-full", replacement="_full", x=mq[mq3+mq5>lower.mq]$junction)

## Check that full is in all of these, as a precaution:
stopifnot(identical(grep("full", good.mq.junctions),c(1:length(good.mq.junctions))))

df.sv.sample.with.good.mq = df.sv.sample[sequence.name %in% good.mq.junctions]

## Only keep fusions that pass this other filter
## the ".add" in tcga.add just means it's for this file
tcga.add = data.table(read.delim("../all_fdr_post_filter_all.tab")[,])
tcga.add[,sequence.name:=paste(paste(gene1,gene2,pos1,pos2,sep="-"),"_full",sep="")]

df.sv.sample.with.good.mq.and.pass.filter = df.sv.sample.with.good.mq[sequence.name %in% tcga.add$sequence.name]


## restrict to sv's that are in tablex
df.sv.sample.also.in.tablex <- df.sv.sample.with.good.mq.and.pass.filter[gene1.gene2 %in% tablex$genenames]

## Also check that there are no cases in which gene1.gene2
## matches part of the sequence in df.sv.sample, but it's not
## actually the full name of the pair of genes. E.g. something like
## AA-BB-CC-DD-11-22_full is the sequence name in df.sv.sample
## and AA-BB is the pair of genes in tablex. It would be
## surprising, but check to be sure.
## Check that when splitting these with dashes, also get 4 items:
lengths.strsplit.sequence <- sapply(lapply(lapply(df.sv.sample.also.in.tablex$sequence.name, FUN=strsplit, split="-"),"[[",1), FUN=length)
stopifnot(sum(lengths.strsplit.sequence!=4)==0)

## Extra check that there are no genes with dashes in them in
## tablex. Not really necessary, but good to double-check.
gene1.with.dashes.in.tablex <- tablex[grep(pattern="-",gene1),gene1]
gene2.with.dashes.in.tablex <- tablex[grep(pattern="-",gene2),gene2]
stopifnot((length(gene1.with.dashes.in.tablex) + length(gene2.with.dashes.in.tablex))==0)


write.table(file="../tablex.with.sequence.names.tab", df.sv.sample.also.in.tablex, sep="\t", quote=F, row.names=FALSE, append=FALSE)

## If we have two splice variants in one sample, we still want
## to count that as a private fusion
## So look at list of unique gene1-gene2-sample combinations
##  Then see how many of those have only 1 gene1-gene2 pair

df.sv.sample.also.in.tablex.no.bodymap <- df.sv.sample.also.in.tablex[investigation!="TCGA-BODY",]

unique.gene1.gene2.sample.combos <- unique(df.sv.sample.also.in.tablex.no.bodymap$gene1.gene2.sample.combo)


gene1.within.unique.gene1.gene2.sample.combos<- sapply(lapply(lapply(unique.gene1.gene2.sample.combos, FUN=strsplit,split=" "),"[[",1),"[",1)
gene2.within.unique.gene1.gene2.sample.combos<- sapply(lapply(lapply(unique.gene1.gene2.sample.combos, FUN=strsplit,split=" "),"[[",1),"[",2)
gene1.gene2.combos.within.unique.gene1.gene2.sample.combos <- paste(gene1.within.unique.gene1.gene2.sample.combos,gene2.within.unique.gene1.gene2.sample.combos,sep="-")
unique.gene1.gene2.combos.within.unique.gene1.gene2.sample.combos <- unique(gene1.gene2.combos.within.unique.gene1.gene2.sample.combos)
print("length(unique.gene1.gene2.combos.within.unique.gene1.gene2.sample.combos)")
print(length(unique.gene1.gene2.combos.within.unique.gene1.gene2.sample.combos))

df.temp <- data.table(unique.gene1.gene2.sample.combos,gene1.gene2.combos.within.unique.gene1.gene2.sample.combos)
n.df.temp <- df.temp[,.N,by=gene1.gene2.combos.within.unique.gene1.gene2.sample.combos][order(-N)]
sum(n.df.temp$N==1)
private.fusions.v1 <- n.df.temp[N==1,gene1.gene2.combos.within.unique.gene1.gene2.sample.combos]
print("Number of private fusions, without restrictions for geneA-geneB/geneA-geneC occurring in same sample issue.:\n")
print("length(private.fusions.v1):")
n.private.fusions.v1 <- length(private.fusions.v1)
print(length(private.fusions.v1))
## 660


## Now get data table of df.sv that are in the private fusions list:
## might have multiple splice variants per fusion
df.sv.sample.in.private.fusions.v1 <- df.sv.sample.also.in.tablex[gene1.gene2 %in% private.fusions.v1,]

## Now take this data table and pick one splice variant for each
## fusion; just do first one
## df.sv.sample.in.private.fusions.one.sv.per.fusion
list.for.constructing.df <- vector("list", length=n.private.fusions.v1)

for (ii in 1:n.private.fusions.v1){
    rows.for.this.fusion <- df.sv.sample.in.private.fusions.v1[gene1.gene2==private.fusions.v1[ii],]
    stopifnot(dim(rows.for.this.fusion)[1]>=1)
    list.for.constructing.df[[ii]] <- rows.for.this.fusion[1,]
}

df.sv.sample.in.private.fusions.one.sv.per.fusion <- list.for.constructing.df[[1]]
for (ii in 2:n.private.fusions.v1){
    df.sv.sample.in.private.fusions.one.sv.per.fusion <- rbind(df.sv.sample.in.private.fusions.one.sv.per.fusion, list.for.constructing.df[[ii]])
}





## Now consider the private fusions we want to look at
##  for the sake of finding genes recurrently present in
##  gene fusions.
## 
## First: are there any cases in which
## (g1,g2) and (g1,g3) occur in same sample:
counts.g1.sample <- df.sv.sample.in.private.fusions.one.sv.per.fusion[,.N,by=.(gene1,sample.id)]
setorder(counts.g1.sample,-N)

counts.g2.sample <- df.sv.sample.in.private.fusions.one.sv.per.fusion[,.N,by=.(gene2,sample.id)]
setorder(counts.g2.sample,-N)

## Find cases in which (g1,g2) and (g1,g3) occur in same sample
## AND (g1,gx) does NOT occur in another sample (where gx could
##  equal g2,g3 or some other gene)
## We will remove these.
indices.g1.counts.geq.2 <- which(counts.g1.sample$N>=2)
n.g1.counts.geq.2 <- length(indices.g1.counts.geq.2)
gene1.values.with.more.than.one.per.sample <- counts.g1.sample[N>=2,gene1]
gene1.values.with.more.than.one.per.sample.and.not.in.other.samples <- character(0)

for (ii in indices.g1.counts.geq.2){
    ## How many fusions are there in total that have gene1 equal g1?
    gene1.matches <- df.sv.sample.in.private.fusions.one.sv.per.fusion[gene1==counts.g1.sample[indices.g1.counts.geq.2[ii],gene1],]
    if (dim(gene1.matches)[1] == counts.g1.sample$N[indices.g1.counts.geq.2[ii]]){
        gene1.values.with.more.than.one.per.sample.and.not.in.other.samples <- c(gene1.values.with.more.than.one.per.sample.and.not.in.other.samples,counts.g1.sample[indices.g1.counts.geq.2[ii],gene1])
    }
    print(counts.g1.sample[indices.g1.counts.geq.2[ii],gene1])
    print(dim(gene1.matches))
    print("\n")
}


## Do same thing for gene 2 now

## Find cases in which (g1,g2) and (g1,g3) occur in same sample
## AND (g1,gx) does NOT occur in another sample (where gx could
##  equal g2,g3 or some other gene)
## We will remove these.
indices.g2.counts.geq.2 <- which(counts.g2.sample$N>=2)
n.g2.counts.geq.2 <- length(indices.g2.counts.geq.2)
gene2.values.with.more.than.one.per.sample <- counts.g2.sample[N>=2,gene2]
gene2.values.with.more.than.one.per.sample.and.not.in.other.samples <- character(0)

df.for.recurrent.genes.first.chop <- df.sv.sample.in.private.fusions.one.sv.per.fusion[!(gene1 %in% gene1.values.with.more.than.one.per.sample.and.not.in.other.samples),]


for (ii in indices.g2.counts.geq.2){
    ## How many fusions are there in total that have gene2 equal g2?
    gene2.matches <- df.sv.sample.in.private.fusions.one.sv.per.fusion[gene2==counts.g2.sample[indices.g2.counts.geq.2[ii],gene2],]
    if (dim(gene2.matches)[1] == counts.g2.sample$N[indices.g2.counts.geq.2[ii]]){
        gene2.values.with.more.than.one.per.sample.and.not.in.other.samples <- c(gene2.values.with.more.than.one.per.sample.and.not.in.other.samples,counts.g2.sample[indices.g2.counts.geq.2[ii],gene2])
    }
    gene2.matches.after.first.chop <- df.for.recurrent.genes.first.chop[gene2==counts.g2.sample[indices.g2.counts.geq.2[ii],gene2],]
    print(counts.g2.sample[indices.g2.counts.geq.2[ii],gene2])
    print(dim(gene2.matches))
    print("\n")
    ##
    print(dim(gene2.matches.after.first.chop))
    print("\n")
}




df.for.recurrent.genes <- df.for.recurrent.genes.first.chop[!(gene2 %in% gene2.values.with.more.than.one.per.sample.and.not.in.other.samples),]

print("dim(df.for.recurrent.genes.first.chop)[1]:\n")
print(dim(df.for.recurrent.genes.first.chop)[1])

print("dim(df.for.recurrent.genes)[1]:\n")
print(dim(df.for.recurrent.genes)[1])

Tnm <- dim(df.for.recurrent.genes)[1]

n.recurrent.gene1 <- df.for.recurrent.genes[,.N,by=gene1][order(-N)]
n.recurrent.gene2 <- df.for.recurrent.genes[,.N,by=gene2][order(-N)]

print("n.recurrent.gene1[1:10]")
print(n.recurrent.gene1[1:10])
print("n.recurrent.gene2[1:10]")
print(n.recurrent.gene2[1:10])

num.five.prime <- sum(n.recurrent.gene1$N>=2)
num.three.prime <- sum(n.recurrent.gene2$N>=2)
print("num.five.prime:")
print(num.five.prime)
print("num.three.prime:")
print(num.three.prime)



## dim(full.final.list)
## table(full.final.list$panBloomCount,exclude=NULL)



## The probability that no chimera is recurrent is computed using the poisson approximation with P(X=0)~=  e^-lambda where lambda= (n(n-1)/2) / (g(g-1)).  For n=644 and g=2*10^5, 1-exp(-lambda)  ~= 0.00051
nch = length(unique(tablex$genenames[tablex$tcganame!="TCGA-BODY"]))
print("nch:")
print(nch)
gch = 20000
lambdach = ((nch*(nch-1))/2)/(gch*(gch-1))
print("The probability that no chimera is recurrent is computed using the poisson approximation; calculation:\n")
print(1 - exp(-lambdach))
cat("\n")





## uses Corollary 2.3 of N. Henze
## Assumes m is   fixed and c = 2





## sMACHETE reports 31 5’ partners (p-value of 1.5x10-7) and 37 3’ partners, (p-value of 1.0x10-10), which are highly statistically significant findings.

## OLD: Tnm=dim(unique(singles[,list(gene1,gene2)]))[1]
n=20000 # number of genes
c=2
m=num.five.prime
nu2=c(1:m)
cat("RALA: n=", n,";Tnm=, ", Tnm, ";c=", c, ";m=", m, "\n")
t=Tnm/sqrt(n)
print(1-exp(-(1/c)*(t^c))*sum((((1/c)*t^c)^(nu2-1))/(factorial(nu2-1))))

cat("\n\n")
m=num.three.prime
nu2=c(1:m)
cat("RALA: n=", n,";Tnm=, ", Tnm, ";c=", c, ";m=", m, "\n")
t=Tnm/sqrt(n)
print(1-exp(-(1/c)*(t^c))*sum((((1/c)*t^c)^(nu2-1))/(factorial(nu2-1))))



cat("\n\n")
m=17
nu2=c(1:m)
cat("RALA: n=", n,";Tnm=, ", Tnm, ";c=", c, ";m=", m, "\n")
t=Tnm/sqrt(n)
print(1-exp(-(1/c)*(t^c))*sum((((1/c)*t^c)^(nu2-1))/(factorial(nu2-1))))


print("Tnm/sqrt(n):")
print(t)


## This does the same thing in a different way:

summand.nu <- function(t,nu){
    (((t^2)/2)^nu )/(factorial(nu))
}

prob.T.over.root.n.leq.t <- function(t, m) {
    summands <- vector("numeric", length = m)
    summands[] <- NA
    for (i.t in 0:(m-1)){
        summands[(i.t+1)] <- summand.nu(t, i.t)
    }
    1 - ((exp((-t^2)/2))* sum(summands))
}

## m = 31
## n = 20000
## 641 fusions (i.e. "balls")

## So set t = T/sqrt(n) = 641/sqrt(20000)
n.fusions <- Tnm

t0 = n.fusions/sqrt(20000)
m1 = num.five.prime

print("prob.T.over.root.n.leq.t(t=t0, m=m1)")
print(prob.T.over.root.n.leq.t(t=t0, m=m1))



t0 = n.fusions/sqrt(20000)
m2 = num.three.prime
print("prob.T.over.root.n.leq.t(t=t0, m=m2)")
print(prob.T.over.root.n.leq.t(t=t0, m=m2))



## RALA, a Ras-family G-protein and known oncogene (Lim et al., 2005), has three distinct partners found in OV and GBM
## Check both the df that takes out the weird g1-g2,g1-g3 cases and
## the one that doesn't:
print("df.for.recurrent.genes[gene1==\"RALA\"|gene2==\"RALA\",]")
print(df.for.recurrent.genes[gene1=="RALA"|gene2=="RALA",])
print("df.sv.sample.in.private.fusions.one.sv.per.fusion[gene1==\"RALA\"|gene2==\"RALA\",]")
print(df.sv.sample.in.private.fusions.one.sv.per.fusion[gene1=="RALA"|gene2=="RALA",])
## equal for RALA

## also see this to see that they are indeed OV and GBM
print("tablex[gene1==\"RALA\"|gene2==\"RALA\",]")
print(tablex[gene1=="RALA"|gene2=="RALA",])




##  A fourth fusion involving RALA, RALA-YAE1D1, was identified by sMACHETE as a recurrent gene fusion in OV
## Note that pancancer.by.genenames is 3 for YAE1D1
print("tablex[gene1==\"RALA\"|gene2==\"RALA\",]")
print(tablex[gene1=="RALA"|gene2=="RALA",])



## ZBTB20, another known oncogene (Lim et al., 2005; Zhao, Ren, and Tang, 2014) that has also been identified as a gene fusion partner in tumors, was also recovered purely on the basis of participating in private fusions.
print("tablex[gene1==\"ZBTB20\"|gene2==\"ZBTB20\",]")
print(tablex[gene1=="ZBTB20"|gene2=="ZBTB20",])

df.for.recurrent.genes[gene1=="ZBTB20"|gene2=="ZBTB20",]
df.sv.sample.in.private.fusions.one.sv.per.fusion[gene1=="ZBTB20"|gene2=="ZBTB20",]


## UVRAG, a tumor supressor with activating oncogene activity (He and Liang, 2015) had the highest diversity of 5’ and 3’ partners and has previously not been reported as participating in fusions.
print("tablex[gene1==\"UVRAG\"|gene2==\"UVRAG\",]")
print(tablex[gene1=="UVRAG"|gene2=="UVRAG",])

print("df.for.recurrent.genes[gene1==\"UVRAG\"|gene2==\"UVRAG\",]")
print(df.for.recurrent.genes[gene1=="UVRAG"|gene2=="UVRAG",])
print("df.sv.sample.in.private.fusions.one.sv.per.fusion[gene1==\"UVRAG\"|gene2==\"UVRAG\",]")
print(df.sv.sample.in.private.fusions.one.sv.per.fusion[gene1=="UVRAG"|gene2=="UVRAG",])

gene1.among.private.fusions <- sapply(lapply(lapply(private.fusions.v1,strsplit,split="-"),"[[",1),"[",1)
gene2.among.private.fusions <- sapply(lapply(lapply(private.fusions.v1,strsplit,split="-"),"[[",1),"[",2)
unique.genes.among.private.fusions <- unique(union(gene1.among.private.fusions,gene2.among.private.fusions))    

## What has highest diversity among private fusions?
matrix.diversity.private.fusions <- matrix(NA,nrow=length(unique.genes.among.private.fusions),ncol=2)

for (ii in 1:length(unique.genes.among.private.fusions)){
    this.gene <- unique.genes.among.private.fusions[ii]
    matrix.diversity.private.fusions[ii,1]<- dim(df.for.recurrent.genes[gene1==this.gene|gene2==this.gene,])[1]
    matrix.diversity.private.fusions[ii,2]<- dim(df.sv.sample.in.private.fusions.one.sv.per.fusion[gene1==this.gene|gene2==this.gene,])[1]
}


## Are columns 1 and 2 identical? presumably not:
unequal.entries.12 <- which(matrix.diversity.private.fusions[,1]!=matrix.diversity.private.fusions[,2])

df.diversity.private.fusions <- data.table(matrix.diversity.private.fusions,unique.genes.among.private.fusions)


## biggest 5 by first column
print("biggest 5 by first column")
print(df.diversity.private.fusions[order(-V1)][1:5])

## biggest 5 by second column
print("biggest 5 by second column")
print(df.diversity.private.fusions[order(-V2)][1:5])



## SORL1

print("tablex[gene1==\"SORL1\"|gene2==\"SORL1\",]")
print(tablex[gene1=="SORL1"|gene2=="SORL1",])





## sMACHETE detected 7 splice variants of TMPRSS2-ERG
print("dim(full.final.list[gene1==\"TMPRSS2\" & gene2==\"ERG\" & tcganame==\"TCGA-PRAD\",.(genepair)])[1]\n")
print(dim(full.final.list[gene1=="TMPRSS2" & gene2=="ERG" & tcganame=="TCGA-PRAD",.(genepair)])[1])


## the fusion in GBM with the highest incidence is FGFR3-TACC3
print("full.final.list[tcganame==\"TCGA-GBM\",.(genepair,maxCount)][order(maxCount, decreasing=TRUE)][1:3]\n")
print(full.final.list[tcganame=="TCGA-GBM",.(genepair,maxCount)][order(maxCount, decreasing=TRUE)][1:3])





## _ unique gene pairs expressed as fusions are called by sMACHETE,
y <- copy(full.final.list)
y[,genenames:= paste0(gene1,"-",gene2)]
print("y<- full.final.list; y[,genenames:= paste0(gene1,\"-\",gene2)]\n")
print("length(unique(y$genenames):")
print(length(unique(y$genenames)))
rm(y)


## 768 unique fusions are called by ChimerSeq, ..., and only 215 are common to both algorithms

common.samples.chimera.mach <- unique(chimera.and.mach$BarcodeID)

## Number of samples in common:
print("Number of files in common: length(common.samples.chimera.mach):\n")
print(length(common.samples.chimera.mach))

## (approx) number of unique fusions from ChimerSeq on common samples:
print("length(unique(chimera.and.mach$Fusion_pair)):\n")
print(length(unique(chimera.and.mach$Fusion_pair)))


## How many does sMACHETE find in these samples?
machete.output.table = fread("../all_fdr_post_filter_all.tab")

tablex[,genes.positions.investigation:= paste(gene1,gene2,pos1,pos2,tcganame)]
machete.output.table[,genes.positions.investigation:= paste(gene1,gene2,pos1,pos2,investigation)]


## Check that these are all unique, i.e. that tablex does not have
## a repeat within it.
stopifnot(sum(duplicated(tablex[,genes.positions.investigation]))==0)

print("length(unique(mach.ids.no.bodymap.no.normals$BarcodeID)):\n")
print(length(unique(mach.ids.no.bodymap.no.normals$BarcodeID)))

## endings.common.samples <- sapply(lapply(lapply(lapply(common.samples.chimera.mach, FUN=strsplit, split=""),FUN="[[",1), tail, n=4),paste,collapse="")

## double-check that chopping off the ending of common.samples
##  won't lead to problems; of course we've already ensured that
##  there shouldn't be multiple samples in the discovery set with the
##  same case id
stopifnot(sum(duplicated(mach.ids.no.bodymap.no.normals$case_id))==0)

## Make case_id's for common.samples
common.samples.case.ids <- gsub(pattern="-0[123][ABC]$", x=common.samples.chimera.mach, replacement="")

## There should not be any duplicates (especially as we already checked
##  above for the superset) but check again.
stopifnot(sum(duplicated(common.samples.case.ids))==0)

## Keep those things in machete.output.table with samples in
## common.samples.case.ids

machete.output.table.restricted.to.common.samples <- machete.output.table[case_id %in% common.samples.case.ids,]

## Now see which gene1-gene2-pos1-pos2-disease combinations in
## tablex (output of all sMACHETE samples) came from a sample in the
## common samples set, i.e. which of these matches with a
## gene1-gene2-pos1-pos2-disease combination in
## machete.output.table.restricted.to.common.samples

genes.positions.investigation.stemming.from.common.samples <- unique(machete.output.table.restricted.to.common.samples$genes.positions.investigation)

tablex.part.stemming.from.common.samples <- tablex[genes.positions.investigation %in% genes.positions.investigation.stemming.from.common.samples,]

print("dim(tablex.part.stemming.from.common.samples):\n")
print(dim(tablex.part.stemming.from.common.samples))

## 525 unique gene pairs expressed as fusions are called on this set by sMACHETE
print("length(unique(tablex.part.stemming.from.common.samples$genenames)):\n")
print(length(unique(tablex.part.stemming.from.common.samples$genenames)))

print("Number of fusions on common sample, called by both algorithms:\n")
print(length(intersect(gsub(pattern="-", replacement="_", unique(tablex.part.stemming.from.common.samples$genenames)),unique(chimera.and.mach$Fusion_pair))))

## Of note, among this list, 12 distinct gene fusions involving HLA or ribosomal protein subunit genes, proxies for likely false positives due to their high expression, are called by ChimerSeq while none are called by sMACHETE.
## ChimerSeq
print("Chimerseq: length(unique(chimera.and.mach[Fusion_pair %like% \"HLA\",Fusion_pair]))")
print(length(unique(chimera.and.mach[Fusion_pair %like% "HLA",Fusion_pair])))

print("sMACHETE: length(tablex.part.stemming.from.common.samples[genenames %like% \"HLA\",genenames])")
print(length(tablex.part.stemming.from.common.samples[genenames %like% "HLA",genenames]))

## TFG-GPR128 (Chase et al., 2010) and two predicted ovarian-specific recurrent fusions, METTL3-TM4SF1 and RCC1-UBE2D2. Consistent with the range of previous reports of the prevalence of TFG-GPR128 in the population (3/120 as reported in Chase et al., 2010, 0.5%-7.1%), sMACHETE estimates its frequency in TCGA data to be <1% in sarcoma (SARC), 2.8% in pancreatic adenocarcinoma (PAAD), and 1.5% in ovarian carcinoma (OV) (see Table 5).
## can also read these off of maxFreq from table 5

print("full.final.list[gene1==\"TFG\" & gene2==\"GPR128\" & tcganame==\"TCGA-SARC\",maxFreq]:\n")
print(full.final.list[gene1=="TFG" & gene2=="GPR128" & tcganame=="TCGA-SARC",maxFreq])

print("full.final.list[gene1==\"TFG\" & gene2==\"GPR128\" & tcganame==\"TCGA-PAAD\",maxFreq]:\n")
print(full.final.list[gene1=="TFG" & gene2=="GPR128" & tcganame=="TCGA-PAAD",maxFreq])
# .022d

print("full.final.list[gene1==\"TFG\" & gene2==\"GPR128\" & tcganame==\"TCGA-OV\",maxFreq]:\n")
print(full.final.list[gene1=="TFG" & gene2=="GPR128" & tcganame=="TCGA-OV",maxFreq])





print("full.final.list[gene1==\"METTL3\" & gene2==\"TM4SF1\" & tcganame==\"TCGA-SARC\",maxFreq]:\n")
print(full.final.list[gene1=="METTL3" & gene2=="TM4SF1" & tcganame=="TCGA-SARC",maxFreq])

print("full.final.list[gene1==\"METTL3\" & gene2==\"TM4SF1\" & tcganame==\"TCGA-PAAD\",maxFreq]:\n")
print(full.final.list[gene1=="METTL3" & gene2=="TM4SF1" & tcganame=="TCGA-PAAD",maxFreq])

print("full.final.list[gene1==\"METTL3\" & gene2==\"TM4SF1\" & tcganame==\"TCGA-OV\",maxFreq]:\n")
print(full.final.list[gene1=="METTL3" & gene2=="TM4SF1" & tcganame=="TCGA-OV",maxFreq])



print("full.final.list[gene1==\"RCC1\" & gene2==\"UBE2D2\" & tcganame==\"TCGA-OV\",maxFreq]:\n")
print(full.final.list[gene1=="RCC1" & gene2=="UBE2D2" & tcganame=="TCGA-OV",maxFreq])



## RCC1-UBE2D2 is predicted to be specific to ovarian tumors
print("tablex[genenames==\"RCC1-UBE2D2\",tcganame]:\n")
print(tablex[genenames=="RCC1-UBE2D2",tcganame])

## The chimera METTL3-TM4SF1 of METTL3, a methyltransferase-like protein involved in splicing, and TM4SF1, a transmembrane protein of unknown function, was seen in 5.9% of tumors and also specific to ovarian cancer.
print("tablex[genenames==\"METTL3-TM4SF1\",tcganame]:\n")
print(tablex[genenames=="METTL3-TM4SF1",tcganame])








## Of the fusions called for bodymap, are they also called
##  for a cancer?

tablex$genenames[tablex$tcganame=="TCGA-BODY"]

dim(tablex[genenames %in% tablex$genenames[tablex$tcganame=="TCGA-BODY"],])

tablex[genenames %in% tablex$genenames[tablex$tcganame=="TCGA-BODY"],.(genenames,tcganame)]


## sMACHETE predicted 100 recurrent gene fusions 
## look at df.sv.sample.also.in.tablex
## and see how many gene1.gene2 values have at least two values of sample.id
length.unique <- function(x){ length(unique(x)) }
recurrent.fusions.calc <- df.sv.sample.also.in.tablex[,lapply(.SD,length.unique),by=.(gene1.gene2),.SDcols=c("sample.id")][order(-sample.id)]

recurrent.fusions.add <- recurrent.fusions.calc[sample.id>1,]

print("How many recurrent gene fusions?\n")
print(dim(recurrent.fusions.add)[1])
print("Check on how many private; doing in different way from before:\n")
print(sum(recurrent.fusions.calc$sample.id==1))

recurrent.fusions.vec <- recurrent.fusions.add$gene1.gene2


## Strawberry Notched Homolog, SBNO2, in the putative fusion product SBNO2-SERINC2
print(recurrent.fusions.add[gene1.gene2=="SBNO2-SERINC2",])
print(df.sv.sample.also.in.tablex[gene1.gene2=="SBNO2-SERINC2",])
print(tablex[genenames=="SBNO2-SERINC2",])

## RPS6KB1-VMP1 fusion, previously identified as a recurrent fusion in breast cancers (BRCA) (Inaki, et al., 2011), was detected for the first time in other cancer types, such as LUAD and OV (Table 5).
print(recurrent.fusions.add[gene1.gene2=="RPS6KB1-VMP1",])
print(df.sv.sample.also.in.tablex[gene1.gene2=="RPS6KB1-VMP1",])
print(tablex[genenames=="RPS6KB1-VMP1",])

### PAAD, which had previously lacked reports of recurrent fusions, were found to harbor a handful of rare recurrent fusions when all cancer types were used to estimate recurrence.
print("How many PAAD recurrent fusions?\n")
print(dim(tablex[genenames %in% recurrent.fusions.vec & tcganame=="TCGA-PAAD",]))

print('tablex[genenames %in% recurrent.fusions.vec & tcganame=="TCGA-PAAD",.(gene1,gene2,pancancer.by.genenames)][order(-pancancer.by.genenames)]')
print(tablex[genenames %in% recurrent.fusions.vec & tcganame=="TCGA-PAAD",.(gene1,gene2,pancancer.by.genenames)][order(-pancancer.by.genenames)])


## Some of these rare recurrent gene fusions were present across tumor types in addition to PAAD; for example, ERBB2-PPP1R1B was detected in just two total tumors across TCGA and, if validated, could be targetable with current drugs (Table 5).
print(recurrent.fusions.add[gene1.gene2=="ERBB2-PPP1R1B",])
print(df.sv.sample.also.in.tablex[gene1.gene2=="ERBB2-PPP1R1B",])
print(tablex[genenames=="ERBB2-PPP1R1B",])




## The sMACHETE algorithm called 12 fusions that were nominated by the MACHETE from the Body Map samples; 7 of these were found by the SBT queries only in the Body Map query, and each of these was only counted once in the query. The fusion C10orf68-CCDC7 and HTATSF1-BRS3 were each counted once in the Body Map query and 93 and 6 times, respectively, in total when summing over all the counts in each of the 10 cancer SBT queries.


## Fusions NOMINATED by bodymap:
tcga.nominated.by.body <- tcga.add[investigation %like% "ERR",]
tcga.nominated.by.body[,genepair:=paste(gene1,gene2,pos1,pos2,sep=" ")]
## which of called fusions come from these:
tablex.nominated.by.body <- tablex[genepair %in% tcga.nominated.by.body$genepair,]
bodymap.fusions <- unique(tablex.nominated.by.body$genenames)
print("how many fusions called that were nominated by bodymap, i.e. found by MACHETE in it?")
print(bodymap.fusions)

## How many samples is each one in?
## OLD: bodymap.fusions <- tablex[tcganame %like% "BODY"]$genenames

## Information on each of the fusions called that were nominated by bodymap:
tablex.nominated.by.body[,.(genepair, gene1, gene2, pan_cancer, strand1, strand2, chr1, chr2, pos1, pos2, AbsPos1Pos2Diff, tcganame, pancancer.by.genenames)]

## In df.sv.sample.also.in.tablex.for.bodymap.fusions, the
## investigation variable tells you what SBT it comes from
df.sv.sample.also.in.tablex.for.bodymap.fusions <- df.sv.sample.also.in.tablex[gene1.gene2 %in% bodymap.fusions,]

## each identified chimera was found in exactly one sample
n.bodymap.fusions.in.bodymap <- df.sv.sample.also.in.tablex.for.bodymap.fusions[investigation=="TCGA-BODY",lapply(.SD,length.unique),by=.(gene1.gene2),.SDcols=c("gene1.gene2.sample.combo")]
print("For this output, the investigation variable tells you what SBT it comes from. So this is the number of times fusions nominated by bodymap were found in bodymap by the SBT.")
print("\n\n n.bodymap.fusions.in.bodymap:\n")
print(n.bodymap.fusions.in.bodymap)




n.bodymap.fusions.in.cancers <- df.sv.sample.also.in.tablex.for.bodymap.fusions[investigation!="TCGA-BODY",lapply(.SD,length.unique),by=.(gene1.gene2),.SDcols=c("gene1.gene2.sample.combo")]
print("For this output, the investigation variable tells you what SBT it comes from. So this is the number of times fusions nominated by bodymap were found in cancers by the SBT.")
print("n.bodymap.fusions.in.cancers <- df.sv.sample.also.in.tablex.for.bodymap.fusions[investigation!=\"TCGA-BODY\",lapply(.SD,length.unique),by=.(gene1.gene2),.SDcols=c(\"gene1.gene2.sample.combo\")]")
print("n.bodymap.fusions.in.cancers")
print(n.bodymap.fusions.in.cancers)


print("fusions nominated by bodymap and not found in bodymap by SBTs:")
print(setdiff(bodymap.fusions,n.bodymap.fusions.in.bodymap$gene1.gene2))

## all fusions were intrachromosomal within 1 MB.

df.sv.sample.also.in.tablex.for.bodymap.fusions[,pos1:= as.numeric(sapply(lapply(lapply(lapply(lapply(df.sv.sample.also.in.tablex.for.bodymap.fusions$sequence.name, FUN=strsplit, split="_"),FUN="[[",1),FUN=strsplit, split="-"),FUN="[[",1),"[",3))]

df.sv.sample.also.in.tablex.for.bodymap.fusions[,pos2:= as.numeric(sapply(lapply(lapply(lapply(lapply(df.sv.sample.also.in.tablex.for.bodymap.fusions$sequence.name, FUN=strsplit, split="_"),FUN="[[",1),FUN=strsplit, split="-"),FUN="[[",1),"[",4))]


df.sv.sample.also.in.tablex.for.bodymap.fusions[,abspos1pos2diff:= abs(pos1-pos2)]
print("Within 1 MB?")
print("table(df.sv.sample.also.in.tablex.for.bodymap.fusions$abspos1pos2diff)")
print(table(df.sv.sample.also.in.tablex.for.bodymap.fusions$abspos1pos2diff))

print("Are all of these intrachromosomal")
print("tablex.nominated.by.body[,.(chr1,chr2)]")
print(tablex.nominated.by.body[,.(chr1,chr2)])





## how many are called by a TCGA and called by body map, include counts of how often found in SBTs
called.bodymap.fusions <- unique(tablex[tcganame=="TCGA-BODY",]$genenames)
print("called.bodymap.fusions")
print(called.bodymap.fusions)

genenames.in.called.bodymap.fusions.and.also.in.a.cancer <- unique(tablex[tcganame!="TCGA-BODY" & genenames %in% called.bodymap.fusions,]$genenames)
print("genenames.in.called.bodymap.fusions.and.also.in.a.cancer")
print(genenames.in.called.bodymap.fusions.and.also.in.a.cancer)
print("tablex[genenames %in% genenames.in.called.bodymap.fusions.and.also.in.a.cancer,.(genepair, gene1, gene2, strand1, strand2, chr1, chr2, pos1, pos2, AbsPos1Pos2Diff, tcganame, maxCount, pan_cancer, pancancer.by.genenames)]")
print(tablex[genenames %in% genenames.in.called.bodymap.fusions.and.also.in.a.cancer,.(genepair, gene1, gene2, strand1, strand2, chr1, chr2, pos1, pos2, AbsPos1Pos2Diff, tcganame, maxCount, pan_cancer, pancancer.by.genenames)])



## Prevalent recurrent fusions were not detected, with the exception of a handful
## this step should be unnecessary:
## tablex[,.(genepair,genenames,maxCount,tcganame)]
tablex[,count.genenames:=sum(maxCount), by=.(genenames,tcganame)]
tablex[genenames=="C10orf68-CCDC7" & tcganame=="TCGA-OV",]

## print("unique(tablex[tcganame==\"TCGA-OV\",.(gene1,gene2,count.genenames)][order(-count.genenames)])[1:20]")
## print(unique(tablex[tcganame=="TCGA-OV",.(gene1,gene2,count.genenames)][order(-count.genenames)])[1:20])


print('unique(tablex[tcganame=="TCGA-OV",.(gene1,gene2,pancancer.by.genenames)][order(-pancancer.by.genenames)])[1:20]')
print(unique(tablex[tcganame=="TCGA-OV",.(gene1,gene2,pancancer.by.genenames)][order(-pancancer.by.genenames)])[1:20])




## sMACHETE’s false positive (FP) rate by this method was very small: only one of the fusions identified in primary tumors in the discovery set was found in the Illumina Body Map data set (see Table 5).
## C10orf68 CCDC7
## sMACHETE’s false positive (FP) rate by this method was zero: none of the fusions identified in primary tumors in the discovery set were found in the Illumina Body Map data set (see Table 5)
tcga.nominated.by.tumors <- tcga.add[!(investigation %like% "ERR"),]
tcga.nominated.by.tumors[,genepair:=paste(gene1,gene2,pos1,pos2,sep=" ")]
## which of called fusions from SBT for bodymap come from fusions nominated by tumors:
tablex.from.bodymap.sbt.nominated.by.tumors <- tablex[genepair %in% tcga.nominated.by.tumors$genepair & tcganame=="TCGA-BODY",]
print("Which of called fusions from SBT for bodymap come from fusions nominated by tumors:")
print(tablex.from.bodymap.sbt.nominated.by.tumors)



## In contrast to sMACHETE, nine distinct fusions found by ChimerSeq were also detected in the Body Map data, including several interchromosomal fusions and four that are present in all samples (see Fig. 4), strongly suggestive of being false-positives.
## In contrast to sMACHETE, nine distinct fusions found by ChimerSeq were also detected in the Body Map data, including several interchromosomal fusions and present in all samples (see Fig. 4)
chimeraBloom.add = fread(sbt.results.bodymap.chimera)
setnames(chimeraBloom.add,"V1","QueryName")
setnames(chimeraBloom.add,"V2","counts")

chimera.info.add=fread("../ChimerDB3.0_ChimerSeq.txt")

chimera.and.mach.add = chimera.info.add[which(!(is.na(match(chimera.info.add$BarcodeID,mach.ids.no.bodymap.no.normals$BarcodeID))))]

chimeraBloom.add = fread("/scratch/PI/horence/sb/bloom/chimera/summary.bt.results.bodymap.chimera.txt")
setnames(chimeraBloom.add,"V1","QueryName")
chimera.and.mach.add[,QueryName:= paste(H_gene,"-",T_gene,"-",  H_position,"-",T_position,"_full_ChimeraScan_",sep="")]
chimeraBloom.add=chimeraBloom.add[which(!is.na(match(chimeraBloom.add$QueryName,chimera.and.mach.add$QueryName)))]

setnames(chimeraBloom.add,"V2","counts")
chimerabloom.add.nonzero <- chimeraBloom.add[counts>0]
print("chimerabloom.add.nonzero")
print(chimerabloom.add.nonzero)

print("interchromosomal fusions??")
print(chimera.and.mach.add[QueryName %in% chimerabloom.add.nonzero$QueryName,.(H_chr,T_chr,QueryName)])

print("How many present in all samples? There are 16 bodymap samples.")
print(chimerabloom.add.nonzero$counts)

write.table(file="../chimera_fusions_detected_inBodyMap_samples_common_toMACH.corrected.by.additionalcalcsfile.tab",chimerabloom.add.nonzero[,list(QueryName,counts)], quote=F,sep="\t",row.names=FALSE)


print("Information about these fusions:")
print(chimera.and.mach.add[QueryName %in% chimerabloom.add.nonzero$QueryName,.(Fusion_pair, H_gene, H_chr, H_position, H_strand, T_gene, T_chr, T_position, T_strand, Breakpoint_Type, Cancertype_or_disease, BarcodeID, Seed_reads_num, Spanning_pairs_num, Frame, Chr_info, H_locus, T_locus, Kinase, Oncogene, QueryName)])




## if one restricts to the subset of the sMACHETE discovery set in the 10 cancers for which an SBT was built, then three distinct fusions found by ChimerSeq were also detected in the Body Map data.

chimera.and.mach.add.only.10.cancers.with.sbts = chimera.info.add[which(!(is.na(match(chimera.info.add$BarcodeID,mach.ids.no.bodymap.no.normals$BarcodeID))) & (!(is.na(match(chimera.info.add$Cancertype_or_disease,full.final.list$tcgaprefix)))))]

chimera.and.mach.add.only.10.cancers.with.sbts[,QueryName:= paste(H_gene,"-",T_gene,"-",  H_position,"-",T_position,"_full_ChimeraScan_",sep="")]


chimerabloom2 = fread("/scratch/PI/horence/sb/bloom/chimera/summary.bt.results.bodymap.chimera.txt")
setnames(chimerabloom2,"V1","QueryName")

chimerabloom.add.only.10.cancers.with.sbts=chimerabloom2[which(!is.na(match(chimerabloom2$QueryName,chimera.and.mach.add.only.10.cancers.with.sbts$QueryName)))]

setnames(chimerabloom.add.only.10.cancers.with.sbts,"V2","counts")

chimerabloom.nonzero.add.only.10.cancers.with.sbts <- chimerabloom.add.only.10.cancers.with.sbts[counts>0]

print("restricted to the subset of the sMACHETE discovery set in the 10 cancers for which an SBT was built:")
print(chimerabloom.nonzero.add.only.10.cancers.with.sbts)



## how many are called by Chimerseq and called by body map, include counts of how often found in SBTs; probably don't need this

chimera.and.mach.add[,genenames:= paste(H_gene,T_gene,sep="-")]

print("chimera.and.mach.add[genenames %in% called.bodymap.fusions,]")
print(chimera.and.mach.add[genenames %in% called.bodymap.fusions,])

genenames.in.called.bodymap.fusions.and.also.in.chimerseq <- chimera.and.mach.add[genenames %in% called.bodymap.fusions,]$genenames
print('tablex[genenames %in% genenames.in.called.bodymap.fusions.and.also.in.chimerseq & tcganame=="TCGA-BODY",.(genepair, gene1, gene2, strand1, strand2, chr1, chr2, pos1, pos2, AbsPos1Pos2Diff, tcganame, maxCount, pan_cancer, pancancer.by.genenames)]')
print(tablex[genenames %in% genenames.in.called.bodymap.fusions.and.also.in.chimerseq & tcganame=="TCGA-BODY",.(genepair, gene1, gene2, strand1, strand2, chr1, chr2, pos1, pos2, AbsPos1Pos2Diff, tcganame, maxCount, pan_cancer, pancancer.by.genenames)])



## ## and show that sMACHETE’s has sensitive detection of gold standard positive controls, such as TMPRSS2-ERG in prostate cancer 

## dim(df.sv.sample.also.in.tablex[gene1.gene2=="TMPRSS2-ERG",]))

print("number of PRAD samples in SBT in which some splice variant of TMPRSS2-ERG is found, and total samples in this SBT, and the fraction:")
print(length(unique(df.sv.sample.also.in.tablex[gene1.gene2=="TMPRSS2-ERG" & investigation=="TCGA-PRAD",]$sample.id)))
print(total.scalpel.counts$TotalBloomSamples[total.scalpel.counts$tcganame=="TCGA-PRAD"])
print(length(unique(df.sv.sample.also.in.tablex[gene1.gene2=="TMPRSS2-ERG" & investigation=="TCGA-PRAD",]$sample.id))/total.scalpel.counts$TotalBloomSamples[total.scalpel.counts$tcganame=="TCGA-PRAD"])

## and show that sMACHETE’s has sensitive detection of gold standard positive controls, such as TMPRSS2-ERG in prostate adenocarcinoma (PRAD) (Tomlins et al., 2008) (which we find in 192 out of 497 prostate cancer cases, or 39%),
## How many in ChimerSeq

## from above:
## chimera.info.add=fread("../ChimerDB3.0_ChimerSeq.txt")

## Need to first get the 497 ids in PRAD that we use, don't want to
## use normals.

this.prad.index <- which(bloomlist.add=="prad")
stopifnot(length(this.prad.index)==1)
sample.ids.prad.sbt <- unique(matrix.files.add[[this.prad.index]]$sample.id)

chimera.prad = chimera.info.add[chimera.info.add$BarcodeID %in% sample.ids.prad.sbt,]

print('length(unique(chimera.prad[H_gene=="TMPRSS2" & T_gene=="ERG",]$BarcodeID))')
print(length(unique(chimera.prad[H_gene=="TMPRSS2" & T_gene=="ERG",]$BarcodeID)))



## Does Chimerseq find FGFR3-TACC3 and if so, what is the prevalence it calls vs. Machete
print('length(unique(df.sv.sample.also.in.tablex[gene1.gene2=="FGFR3-TACC3" & investigation=="TCGA-GBM",]$sample.id))')
print(length(unique(df.sv.sample.also.in.tablex[gene1.gene2=="FGFR3-TACC3" & investigation=="TCGA-GBM",]$sample.id)))
print('total.scalpel.counts$TotalBloomSamples[total.scalpel.counts$tcganame=="TCGA-GBM"]')
print(total.scalpel.counts$TotalBloomSamples[total.scalpel.counts$tcganame=="TCGA-GBM"])
print("Ratio for GBM for sMACHETE:")
print((length(unique(df.sv.sample.also.in.tablex[gene1.gene2=="FGFR3-TACC3" & investigation=="TCGA-GBM",]$sample.id)))/(total.scalpel.counts$TotalBloomSamples[total.scalpel.counts$tcganame=="TCGA-GBM"]))

this.gbm.index <- which(bloomlist.add=="gbm")
stopifnot(length(this.gbm.index)==1)
sample.ids.gbm.sbt <- unique(matrix.files.add[[this.gbm.index]]$sample.id)

chimera.gbm = chimera.info.add[chimera.info.add$BarcodeID %in% sample.ids.gbm.sbt,]

print('length(unique(chimera.gbm[H_gene=="FGFR3" & T_gene=="TACC3",]$BarcodeID))')
print(length(unique(chimera.gbm[H_gene=="FGFR3" & T_gene=="TACC3",]$BarcodeID)))

print("length(sample.ids.gbm.sbt)")
print(length(sample.ids.gbm.sbt))

print("Ratio for GBM for ChimerSeq:")
print((length(unique(chimera.gbm[H_gene=="FGFR3" & T_gene=="TACC3",]$BarcodeID)))/(length(sample.ids.gbm.sbt)))





# figure 4B

chimera.and.mach.add = chimera.info.add[which(!(is.na(match(chimera.info.add$BarcodeID,mach.ids.no.bodymap.no.normals$BarcodeID))))]

chimera.and.mach.add[,genenames:= paste(H_gene,T_gene,sep="-")]
chimera.and.mach.add[,genepair:= paste(H_gene,T_gene, H_position,T_position,sep=" ")]

merged.chimera.mach.genepairs=merge(full.final.list,chimera.and.mach.add,by="genepair",all=T)


aml.from.chimera= data.table(unique(chimera.and.mach.add[Cancertype_or_disease %like% "AML",.(genenames,Cancertype_or_disease)]))

aml.from.smachete = data.table(unique(tablex[tcganame %like% "AML",.(genenames,tcganame)]))

aml.add=data.table(unique(merged.chimera.mach.genepairs[(paste(Cancertype_or_disease,tcganame) %like% "AML"),list(genepair,genenames.x,genenames.y,tcganame,Cancertype_or_disease)]))

## rownames(aml)=aml[,1]
## aml=aml[,-1]
for (i in 1:2){
    aml[is.na(aml[,i]),i]=0
    aml[which(aml[,i]!=0 & !is.na(aml[,i])),i]=1
    aml[,i]=as.numeric(as.vector(aml[,i]))
}

amlgenenames <- vector("character", length=dim(aml.add)[1])
for (ii in 1:length(amlgenenames)){
    if (is.na(aml.add$genenames.x[ii])){
        amlgenenames[ii] <- aml.add$genenames.y[ii]
    } else {
        amlgenenames[ii] <- aml.add$genenames.x[ii]
    }
}
stopifnot(sum(is.na(amlgenenames))==0)

aml.add$genenames <- amlgenenames
aml.add$genenames.x <- NULL
aml.add$genenames.y <- NULL

nejm.add=fread("../NEJM_LAML_comparison.csv")
aml.nejm=merge(x=aml.add,y=nejm.add,by="genenames")
## aml=cbind(aml,nejm[match(rownames(aml),genenames),]$is.Known.fusionEvent)

stopifnot(sum(duplicated(aml.from.chimera$genenames))==0)
stopifnot(sum(duplicated(aml.from.smachete$genenames))==0)
aml.chimera.nejm <- merge(x=aml.from.chimera,y=nejm.add, by="genenames")
aml.smachete.nejm <- merge(x=aml.from.smachete,y=nejm.add, by="genenames")

sum(aml.chimera.nejm$is.NEJM)

sum(aml.smachete.nejm$is.NEJM)

sum(aml.chimera.nejm$is.Known.fusionEvent)

sum(aml.smachete.nejm$is.Known.fusionEvent)

ratio.aml.chimera.nejm <- round(sum(aml.chimera.nejm$is.NEJM)/(dim(aml.from.chimera)[1]),2)
ratio.aml.smachete.nejm <- round(sum(aml.smachete.nejm$is.NEJM)/(dim(aml.from.smachete)[1]),2)

ratio.aml.chimera.known.nejm <- round(sum(aml.chimera.nejm$is.Known.fusionEvent)/(dim(aml.from.chimera)[1]),2)
ratio.aml.smachete.known.nejm <- round(sum(aml.smachete.nejm$is.Known.fusionEvent)/(dim(aml.from.smachete)[1]),2)


known.gold.standard.amls.in.chimera <- paste0(sum(aml.chimera.nejm$is.Known.fusionEvent), " (", 100*ratio.aml.chimera.known.nejm, "%)")
known.gold.standard.amls.in.smachete <- paste0(sum(aml.smachete.nejm$is.Known.fusionEvent), " (", 100*ratio.aml.smachete.known.nejm, "%)")

nejm.gold.standard.amls.in.chimera <- paste0(sum(aml.chimera.nejm$is.NEJM), " (", 100*ratio.aml.chimera.nejm, "%)")
nejm.gold.standard.amls.in.smachete <- paste0(sum(aml.smachete.nejm$is.NEJM), " (", 100*ratio.aml.smachete.nejm, "%)")


print("known.gold.standard.amls.in.chimera:")
print(known.gold.standard.amls.in.chimera)
print("known.gold.standard.amls.in.smachete:")
print(known.gold.standard.amls.in.smachete)

print("nejm.gold.standard.amls.in.chimera:")
print(nejm.gold.standard.amls.in.chimera)
print("nejm.gold.standard.amls.in.smachete:")
print(nejm.gold.standard.amls.in.smachete)


### Similar gene names, first three letters of gene name


aml.from.chimera.with.more.info= data.table(unique(chimera.and.mach.add[Cancertype_or_disease %like% "AML",.(genenames,Cancertype_or_disease,H_gene,T_gene)]))

aml.from.chimera.with.more.info[,gene1.first3:= substr(H_gene,1,3)]
aml.from.chimera.with.more.info[,gene2.first3:= substr(T_gene,1,3)]

first.pass.aml.chimera.similar.gene.names <- aml.from.chimera.with.more.info[gene1.first3==gene2.first3,]

## Remove RUNX1 and MLL cases:

aml.similar.gene.names <- sort(first.pass.aml.chimera.similar.gene.names[(H_gene != "RUNX1" & gene1.first3!="MLL"),]$genenames)

## OLD (same): aml.similar.gene.names <- c("BOD1L1-BOD1","EIF3I-EIF3IP1","HLA-A-HLA-E",  "HLA-C-HLA-B", "HLA-DPA1-HLA-DPB1", "HLA-E-HLA-B", "NBPF1-NBPF15", "SAFB2-SAFB")

chart.colnames <- c("Only detected\nby\nChimerSeq", "Only detected\nby\nsMACHETE")
chart.rownames <- c("Total Gold\nStandards (rate)", "Fusions \nbetween genes with\nsimilar gene names\n(presumed False\nPositives)")
## chart.rownames <- c("Total Gold\nStandards (rate)", "Non-AML fusions \nbetween genes with\nsimilar gene names\n(presumed False\nPositives)")


## make the pdf in makefigure4.R

## plot.new()
## plotpars <- par("usr")
## xvals <- plotpars[1] + (plotpars[2]-plotpars[1])*c(0,.35,.75,1)
## yvals <- plotpars[3] + (plotpars[4]-plotpars[3])*c(0,.6,.78,1)
## similar.gene.offset <- .07*(yvals[2]-yvals[1])
## abline(h=yvals[2])
## abline(v=xvals[3])
## lines(x=c(xvals[2],xvals[4]), y=rep(yvals[3],2))
## lines(x=rep(xvals[2],2), y=c(yvals[1],yvals[3]))
## text(x=mean(xvals[2:3]),y=mean(yvals[3:4]),labels=chart.colnames[1], font=2)
## text(x=mean(xvals[3:4]),y=mean(yvals[3:4]),labels=chart.colnames[2], font=2)
## text(x=mean(xvals[1:2]),y=mean(yvals[2:3]),labels=chart.rownames[1], font=2)
## text(x=mean(xvals[1:2]),y=mean(yvals[1:2]),labels=chart.rownames[2], font=2)
## text(x=mean(xvals[2:3]),y=mean(yvals[2:3]),labels=known.gold.standard.amls.in.chimera)
## text(x=mean(xvals[3:4]),y=mean(yvals[2:3]),labels=known.gold.standard.amls.in.smachete)

## x.similar.genes <- rep(mean(xvals[2:3]), length(aml.similar.gene.names))
## y.similar.genes <- rev(yvals[2] - seq(from=similar.gene.offset, by= similar.gene.offset, length.out = length(aml.similar.gene.names)))
## text(x=x.similar.genes,y=y.similar.genes,labels=aml.similar.gene.names)



colnames(aml)[3]="is.annotated.aml.fusion"


### Calculations re figure 4C

fusions.called.by.a.cancer <- unique(tablex[tcganame!="TCGA-BODY",]$genenames)

this.bodymap.index <- which(bloomlist.add=="body")
stopifnot(length(this.bodymap.index)==1)


## look for rows with the above fusions.called.by.a.cancer

bodymapdt <- data.table(mf.add.only.present.sequences[[this.bodymap.index]])
bodymapdt[,gene1:= sapply(lapply(lapply(sequence.name, FUN=strsplit, split="-"),"[[",1), "[",1)]
bodymapdt[,gene2:= sapply(lapply(lapply(sequence.name, FUN=strsplit, split="-"),"[[",1), "[",2)]
bodymapdt[,genenames:= paste(gene1,gene2,sep="-")]

bodydt.in.fusions.called.by.cancers <- bodymapdt[genenames %in% fusions.called.by.a.cancer,]

genenames.in.called.cancers.and.also.in.bodymap.sbt <- bodydt.in.fusions.called.by.cancers$genenames

print("genenames.in.called.cancers.and.also.in.bodymap.sbt:")
print(genenames.in.called.cancers.and.also.in.bodymap.sbt)


dt.counts.for.genenames.in.called.cancers.and.also.in.bodymap.sbt <- bodydt.in.fusions.called.by.cancers[,.N,by=genenames]

print("dt.counts.for.genenames.in.called.cancers.and.also.in.bodymap.sbt")
print(dt.counts.for.genenames.in.called.cancers.and.also.in.bodymap.sbt)


## check that order is the same, as a precaution:
stopifnot(identical(dt.counts.for.genenames.in.called.cancers.and.also.in.bodymap.sbt$genenames,genenames.in.called.cancers.and.also.in.bodymap.sbt))

counts.for.genenames.in.called.cancers.and.also.in.bodymap.sbt<- dt.counts.for.genenames.in.called.cancers.and.also.in.bodymap.sbt$N


## get info about these:
print("tablex[genenames %in% genenames.in.called.cancers.and.also.in.bodymap.sbt,.(genepair, gene1, gene2, strand1, strand2, chr1, chr2, pos1, pos2, AbsPos1Pos2Diff, tcganame, maxCount, pan_cancer, pancancer.by.genenames)]")
print(tablex[genenames %in% genenames.in.called.cancers.and.also.in.bodymap.sbt,.(genepair, gene1, gene2, strand1, strand2, chr1, chr2, pos1, pos2, AbsPos1Pos2Diff, tcganame, maxCount, pan_cancer, pancancer.by.genenames)])










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




