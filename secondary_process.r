# adds crude gene annotations
## run perl on each file to get the SEQUENCES and add them to allg
## The file hg19_crude_geneann should be in the working directory wdir.
args = commandArgs(trailingOnly = TRUE)
metadatafile=args[1]
wdir=args[2]

outputfilesdir = file.path(wdir, "output_files")

library(stats)
require(data.table)
remake=1

metainfo=data.table(read.csv(metadatafile))

setnames(metainfo,"nicename","sType")

setwd(wdir)

g=data.table(read.delim("hg19_crude_geneann"))
# refflat product https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=502825941_NQQWFDm7G51vKlIgkPhbm9a4N3N4&hgta_doSchemaDb=hg19&hgta_doSchemaTable=refLink)

setwd(wdir)
myfiles=list.files(wdir)
i=0


if (remake==1){

    screenlist=myfiles[myfiles %like% "all_merged_w_anomalies"]


    for (filename in screenlist){

        print ("filename passed")
        niceName=strsplit(filename, "__all")[[1]][1]

        print (niceName)
        if (exists("m")){
            rm(m)
        }
        ## read in file, process if it is larger than 1 row
        m= read.delim(filename) 
        print (dim(m))
        if (dim(m)[1]>1){

            print (paste("i is ",i))
            print (head(m,n=1))

            m=cbind(m,rep(niceName,dim(m)[1]))

            m=m[,-match("niceName",colnames(m))]
            colnames(m)[dim(m)[2]]="sType"


            if (i==0){oldlength=(length(unique(m$niceName)))}
            if (i>0) {oldlength=(length(unique(allg$niceName)))}
            if (i==0){allg=m}

            if (i>0) {allg=rbind(allg,m)}
            i=i+1
            newlength=(length(unique(allg$niceName)))

        }
    }


    ## join to gene
    totalcount=i
    allg=data.table(unique(allg))

    write.table(file="consolidated.tab",allg,sep="\t",quote=F, row.names=F)
    allg=merge(allg,metainfo[,list(investigation,case_id,sample_type,sType)],by="sType", all.x=TRUE,allow.cartesian=TRUE)

    ##setnames(allg,"investigation.x", "investigation")

    print("allg merged with meta information on sample")
    allg[,sampleCount:=length(unique(sType)),by=investigation]
    allg[,investigation.x:=NULL]
    ##setnames(allg,"investigation.y","investigation")

    ## assign meta information
    gene1info=g[match(tolower(allg$gene1),tolower(g$X.name)),]$product
    gene2info=g[match(tolower(allg$gene2),tolower(g$X.name)),]$product
    allg[,g1prod:=gene1info]
    allg[,g2prod:=gene2info]

    allg[,fuspair:=paste(gene1,gene2)]
    setkey(allg,"fuspair")

    ## counts the number of unique samples per fusion pair
    allg[,fuspaircount:=length(unique(sType)),by=fuspair]

    ## defines corr_emp_p which is the empirical p value times the maximum gene expression in each gene participating in fusion
    allg[,corr_emp_p:= emp_p*sum(maxGeneExp.y,maxGeneExp.x, na.rm=T) ,by=list(sType,chr1,chr2,pos1,pos2)]

    ## strand processing:
    ## if the last of strand2 and seq 2 are resp.
    ## calculates poly-a strings for the sequences at the 3' and 5' ends of seq1 and the 5' end of seq 2

    ## num.poly -- ie 5' end of seq 1
    allg[, num.poly:=sum("A"==strsplit(substr(seq2,0,8), "")[[1]]), by=list(gene1,gene2,chr1,chr2,pos1,pos2)]
    allg[, num.poly.s15:=sum("A"==strsplit(substr(seq1,0,8), "")[[1]]), by=list(gene1,gene2,chr1,chr2,pos1,pos2)]

    allg[, num.polyT:=sum("T"==strsplit(substr(seq2,0,8), "")[[1]]), by=list(gene1,gene2,chr1,chr2,pos1,pos2)]
    allg[, num.polyTs15:=sum("T"==strsplit(substr(seq1,0,8), "")[[1]]), by=list(gene1,gene2,chr1,chr2,pos1,pos2)]

    allg[, num.poly.s2:=sum("A"==strsplit(substr(seq2,33,41), "")[[1]]), by=list(gene1,gene2,chr1,chr2,pos1,pos2)]
    allg[, num.poly.s1:=sum("A"==strsplit(substr(seq1,33,41), "")[[1]]), by=list(gene1,gene2,chr1,chr2,pos1,pos2)]


    write.table(file="temp_allg.tab",allg,row.names=F,quote=F,sep="\t")
} ## end if (remake == 1)
maxPoly=7

if (remake==0){
    print ("reading in")
    allg=data.table(read.delim("temp_allg.tab"))
}
print ("done reading")

allg[,no.poly:= num.poly<maxPoly , by=list(gene1,gene2,chr1,chr2,pos1,pos2)]
type2candidates=allg[(no.poly) & (numReads>1) & (emp_p<.1)&(corr_emp_p<5) ][order(fuspair),]

print ("finished type2cands")

type2candidates[,nFus:=length(unique(sType)),by=list(fuspair)]
type2candidates[,nFusPerSample:=length(unique(fuspair)),by=list(sType)]
tempCountGene=data.table(unique(data.frame(type2candidates[(!(gene1 %like% gene2)), list(gene1,gene2)])))

## # extra information added to reports but not used
gene.freq=tapply(rep(1,2*dim(tempCountGene)[1]), c(paste(tempCountGene$gene1),paste(tempCountGene$gene2)),sum)

type2candidates[,nGene1:= gene.freq[match(type2candidates$gene1,names(gene.freq))]]
type2candidates[,nGene2:= gene.freq[match(type2candidates$gene2,names(gene.freq))]]

tot.num=1
write.table(file=paste("num.poly_leq7_nGenesum_geq",tot.num,".tab",sep=""),type2candidates[  (nGene2+nGene1)>tot.num ,],quote=F,sep="\t", row.names=F)
write.table(file=paste("TABLE_Z_RAW_MACHETE.tab",tot.num,".tab",sep=""),type2candidates[  (nGene2+nGene1)>tot.num ,],quote=F,sep="\t", row.names=F)



