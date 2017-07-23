
radius=10 ### radius for gene fusion splice donor/acceptor


require(data.table)

args = commandArgs(trailingOnly = TRUE)
circ=args[1]
lin=args[2]
fusion=args[3]
wdir = args[4] ## 

filetypes=c("circ","lin","fusion")
file.list=c(paste(circ),paste(lin),paste(fusion))

i=0
##for (j in c(3,1,2)){
for (j in c(3,1,2)){
    print (paste("J",j))

    ## would like to add in decoy and anomaly reads

    type=filetypes[j]

    m=data.table(read.delim( file.list[j]))
    m[,dirname:= strsplit(paste(file.list[j]),"/")[[1]][7] ]
    m[,sampleType:=type ]
    niceName=file.list[4] ##strsplit(paste(file.list[j]),"/")[[1]][8]
    m[,niceName:=strsplit(paste(file.list[j]),"/")[[1]][8] ]
    m[,sample:=file.list[j]]

    ## setnames fields for KNIFE:
    ##junction	numReads	logsum	productPhat.x	junction_cdf.x

    if (j!=3){
        setnames(m, "productPhat.x","productPhat")
        setnames(m, "junction_cdf.x","junction_cdf")
    }
    
    if (j==3){ ## j=3 corresponds to MACHETE, this script uses only use some info: 
        setnames(m, "productPhat.y","productPhat")
        setnames(m, "junction_cdf_lower.y","junction_cdf_lower")
        setnames(m, "junction_cdf.y","junction_cdf")
        setnames(m, "numReads.y","numReads")
        setnames(m, names(m)[1],"junction")
        setnames(m, names(m)[19],"badfj1is1")
        setnames(m, names(m)[20],"badfj2is1")
        m[,anomSum:=junc.anom+reg.anomaly+genome.anomaly+ FarJunc.anom, by = list ( junction, niceName) ]
    }

    ## set data in correct format
    g=data.frame(m[  productPhat!="-",])
    g$junction_cdf=as.numeric(as.vector(g$junction_cdf))
    g$numReads = as.numeric(as.vector(g$numReads))
    g$productPhat=as.numeric(as.vector(g$productPhat))

    ## ############ assign names to datatable##
    ## now parse strings


    ## included sample
    used.names=c("sample","chr1","gene1","pos1","strand1","chr2","gene2","pos2","strand2","junction","numReads","productPhat","junction_cdf","junction_cdf_lower", "niceName", "dirname","sampleType","badfj2is1","badfj1is1", "anomSum","maxPosStrandOverlap","maxNegStrandOverlap","type")

    ## retain only the names desired
    g=g[,-which(is.na(match(names(g) , used.names)))]

    g=data.table(g)
    ## substitute special characters in junction name
    g[,junction:=gsub("([|])", ":", paste(junction))]
    if (j==3){
        ## resulting junction name
        ## example: chr19:DAZAP1:99999999:+:chr22:SEPT5:99999999:+:fusion      0
        g[,chr1:= as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 1))]
        g[,gene1:= as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 2))]
        g[,pos1:= as.numeric(as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 3)))]
        g[,strand1:= as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 4))]

        g[,chr2:= as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 5))]
        g[,gene2:= as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 6))]
        g[,pos2:= as.numeric(as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 7)))]
        g[,strand2:= as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 8))]

        g[,type:= as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 9))]

        ## add in names of g if doesn't exist:

    }

    if (j<3){
        g[,chr1:= as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 1))]
        g[,gene1:= as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 2))]
        g[,pos1:= as.numeric(as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 3)))]
        g[,strand1:= as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 4))]
        g[,chr2:= chr1]
        g[,gene2:= as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 4))]
        g[,pos2:= as.numeric(as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 5)))]
        g[,strand1:= as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 7))]
        g[,type:= as.character(lapply(strsplit(paste(g$junction), split=":"), "[", 6))]
        g[,strand2:=strand1]
        ## ideallyg[,anomSum:=decoy]
        g[,anomSum:=NA]

    }
###########################################
#####################
    print ("starting text processing for naming ")
    for (rn in used.names){
        if (is.na(match(rn,names(g)))){
            ## add this as NA
            g[,paste(rn):=NA]
        }
    } 

    print (names(g))
    if (i==0){
        allg=data.table(g)
        rm(g)
    }

    if (i>0){
        allg=rbind(g,allg)
        print (dim(allg))
        rm(g)
    }

    dirsuf=niceName

    print ("i=i+1")
    i=i+1

} ## end for (j in c(3,1,2)){

## done with loop
## make data frame consistent, and impose definitions of fusion

allg=data.table(allg)
allg= allg[ (chr1!=chr2 | abs(pos1-pos2)>radius | strand1!=strand2 )]

##############################

ii=0
for (my.sample in unique (allg$sample)){
    print(my.sample)

    if (my.sample %like% "Appended"){is.machete=1}
    if (!(my.sample %like% "Appended")){is.machete=0}

    allg[,mtype:=paste(is.machete,type), by="junction"]

###################################
    ## assign badfj1 and 2 to be 0 if they are na
    allg[is.na(badfj1is1),badfj1is1:=0]
    allg[is.na(badfj2is1),badfj2is1:=0]

    train.good=allg [ paste(sample)==my.sample & (badfj1is1+badfj2is1)==0]

    ## these are junctions that are bad whereas fj2is could just be linear transcripts

    train.bad=allg[paste(sample)==my.sample & badfj1is1>0 ]
    ## ################################
    machetesample=my.sample ## assumption
    if (is.machete==0){
        ##redefine train.good and bad

        train.bad=allg[sampleType == "fusion" & badfj1is1>0 ]
    }
    INC=100

    train.good[, discrete_p:=round(INC*junction_cdf)/INC]
    train.good[,emp_p:=1] ## will get reassigned later

    for (th in c(0:INC)/INC){
        ## assign a threshold
        tot.bad=length(train.bad$junction_cdf)
        tot.good=length(train.good$junction_cdf)
        prop=tot.bad/tot.good

        current.p=sum(train.bad$junction_cdf>th)/tot.bad

        ## assigning empirical p value as described in machete paper Hsieh et al
        train.good[discrete_p==th,emp_p:=current.p]
    }
    if (ii>0 & dim(train.good)[1]>0){
        print (ii)
        all.passed=rbind(all.passed,train.good)
    }
    if (ii==0){
        all.passed=train.good
    }
    ii=ii+1
} ## end for (my.sample in unique (allg$sample)){

print ("COMPLETED")
#

## falsely called over all called:
## find the good threshold, then report posterior and fdr
## if there are no junctions, set emp_p=0

all.passed[emp_p=="NaN", emp_p:=0]

all.passed[,sinfo:=paste(sample,sampleType)]
all.passed[,fusionInfo:=junction]

## establish cut-offs if used
#min.reads=1
#min.reads.high.cdf=10
#p.thresh=.1
#lower.jun.cdf=.2
test=all.passed

sam=unique(test$sample)
test[,sampleid:=match(sample,sam)]
test[,sinfo:=paste(sample,sampleType)]
test[,fusionInfo:=junction]


test.fj=test[sampleType %like% "fusion",]
test.cir=test[sampleType %like% "cir",]
testf=rbind(test.cir,test.fj)

## want to merge max gene expression in circle and fusion w/ linear from each partner

## now merge to all samples with max linear expression:

tlinears=allg[sample %like% "lin",]
tlinears[,maxGeneExp:=max(as.numeric(numReads)),by=list(niceName,gene1)]
tlinears[,merger:=paste(niceName,gene1,sep="_")]
linears=unique(tlinears[,list(merger,maxGeneExp)])
testf[,merger:=paste(niceName,gene1,sep="_")]
setkey(testf,"merger")
setkey(linears,"merger")

## merge data
quicktest=merge(testf,linears[,list( maxGeneExp,merger)],by="merger", all.x=TRUE)

quicktest[,merger:=paste(niceName,gene2,sep="_")]
setkey(testf,"merger")
fquicktest=merge(quicktest,linears[,list( maxGeneExp,merger)],by="merger",all.x=TRUE)
quicktest=unique(fquicktest[,merger:=NULL])
# 
sam=unique(quicktest$niceName)
quicktest[,sampleid:=match(niceName,sam)]


## allinst ofnicename replacing stem
tablepref=allg$niceName[1]
quicktest[,sample:=NULL]

output_files_dir = file.path(wdir, "output_files")

prefixdirsuf= file.path(output_files_dir, tablepref)
write.table(file=paste(prefixdirsuf,"quicktcga",sep="_"),unique(quicktest[order(gene2)]),sep="\t",quote=F, row.names=F)
write.table(file=paste(prefixdirsuf,"diff_quicktcga",sep="_"),unique(quicktest[gene1!=gene2][order(gene2)]),sep="\t",quote=F, row.names=F)
write.table(file=paste(prefixdirsuf,"same_genequicktcga",sep="_"),unique(quicktest[gene1==gene2][order(gene2)]),sep="\t",quote=F, row.names=F)

write.table(file=paste(prefixdirsuf,"all",sep="_"),unique(quicktest),sep="\t",quote=F,row.names=F)
