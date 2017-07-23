library(data.table)
i=0
for (name in c("ov_tcga_pub","luad_tcga_pub","laml_tcga_pub","gbm_tcga_pub2013","prad_tcga_pub","brca_tcga_pub2015","paad.tcga")){
    i=i+1
    print (paste(i,name))
    mym=data.frame(fread(paste("../",name,".data_mutations_extended.txt",sep="")))[,c(1,16)]
    mym[,2]=substr(mym[,2],0,12)
    mym=data.table(mym)
    mym[,id:=paste(name),]
    print (head(mym))
    if (i==1){mut=mym}
    if (i>1){
        mut=rbind(mut,mym)
        print (name)
    }
}

names(mut)[2]="Tumor_Sample_Barcode"

y.mut <- fread("../tablex.with.sequence.names.tab")

length.seqname <- sapply(lapply(lapply(y.mut$sequence.name, FUN=strsplit, split="-"),FUN="[[",1), FUN=length)
## these should all have length 4, so no genes with dashes in names;
## this should have already been checked in previous script anyway.
stopifnot(sum(length.seqname!=4)==0)

y.mut[,pos1 := as.numeric(sapply(lapply(lapply(lapply(lapply(y.mut$sequence.name, FUN=strsplit, split="_"),FUN="[[",1),FUN=strsplit, split="-"),FUN="[[",1),"[",3))]

y.mut[,pos2 := as.numeric(sapply(lapply(lapply(lapply(lapply(y.mut$sequence.name, FUN=strsplit, split="_"),FUN="[[",1),FUN=strsplit, split="-"),FUN="[[",1),"[",4))]

y.mut = y.mut[gene1!=gene2]

y.mut[,case.id:= substr(sample.id,0,12)]

y.mut=unique(y.mut[,list(case.id,investigation,gene1,gene2)])


names(mut)[2]="case.id"
names(mut)[1]="gene1"
mut[,gene2:=gene1]

n=length(unique(mut$case.id))
mut[,idn:=length(unique(case.id)),by=id]
## mut[,prevalence:=length(unique(case.id))/n,by=gene1]
## prevalence is now computed by disease
mut[,prevalence:=length(unique(case.id))/idn,by=list(gene1,id)]

ymNoCase=unique(merge(y.mut,mut,by=c("gene1"))[,list(gene1,gene2.y)])

y.mut=cbind(y.mut, mut[match(y.mut$gene1,mut$gene1)]$prevalence)
setnames(y.mut,"V2","prevG1")
y.mut=cbind(y.mut, mut[match(y.mut$gene2,mut$gene2)]$prevalence)
setnames(y.mut,"V2","prevG2")


y.mut[is.na(y.mut)]=0

y.mut[,nFus:=length(unique(case.id)),by=paste(gene1,gene2, investigation)]
## y.mut[,expe:=(prevG1+prevG2-prevG1*prevG2)*nFus]
y.mut[,expe:=(prevG1+prevG2)*nFus]
## prevalence of major fusions (report/describe)
unique(y.mut[gene1=="TMPRSS2",list(investigation,case.id)])[!(investigation %like% "PRAD")]

ym1=merge(mut,y.mut,by=c("gene1","case.id"))
ym2=merge(mut,y.mut,by=c("gene2","case.id"))

## below for the 'pan cancer rate:
print ("expected") ## expected based on mutqation frequences" excluding COL is conservative
texp=(sum(unique(y.mut[!(gene1%like% "COL"),list(gene1,gene2,expe)])$expe))
print ("observed")
## tobs=(length(unique(c(ym1[!(gene1 %like% "COL")]$case.id,ym2[!(gene2 %like% "COL")]$case.id))))
tobs=length(c(ym1[!(gene1 %like% "COL")]$case.id,ym2[!(gene2 %like% "COL")]$case.id))
print(tobs)

pval= sum(dpois(0:tobs,texp))
print (paste("pan cancer pval", pval))



## for (ty in unique(y.mut[investigation %like% "PRAD"]$investigation)){
for (ty in c("TCGA-LAML","TCGA-BRCA","TCGA-GBM","TCGA-LUAD","TCGA-OV","TCGA-PAAD","TCGA-PRAD")){
    temp=y.mut[investigation==ty]
    ym1=merge(mut,temp,by=c("gene1","case.id"))
    ym2=merge(mut,temp,by=c("gene2","case.id"))
    print (ty)
    ## expected based on mutation frequences".. would be biased downwards if 
    lambda.est=sum(unique(temp[!(gene1%like% "COL" | gene2%like% "COL"),list(gene1,gene2,expe)])$expe)
    ## observed=length(unique(c(ym1[!(gene1 %like% "COL")]$case.id,ym2[!(gene2 %like% "COL")]$case.id)))
    observed=length(c(ym1[!(gene1 %like% "COL")]$case.id,ym2[!(gene2 %like% "COL")]$case.id))
    pval= sum(dpois(0:observed,lambda.est))
    print (paste("observed",observed,"est",lambda.est,"pval",pval))
    print ("WITHOUT COL EXCLUSION")
    lambda.est=sum(unique(temp[(gene1 %like% ""),list(gene1,gene2,expe)])$expe)
    observed=length(c(ym1[(gene1 %like% "")]$case.id,ym2[(gene2 %like% "")]$case.id))
    pval= sum(dpois(0:observed,lambda.est))
    print (paste("observed",observed,"est",lambda.est,"pval",pval))
}


ty="TCGA-PRAD"
temp=y.mut[investigation==ty & !(gene1 %like% "TMPRSS2" & gene2 %like% "ERG") & !(gene1 %like% "COL" | gene2 %like% "COL"),]
ym1=merge(mut,temp,by=c("gene1","case.id"))
ym2=merge(mut,temp,by=c("gene2","case.id"))
print(ty)
print("PRAD WITHOUT TMPRSS2-ERG and WITH COL EXCLUSION:\n")
## expected based on mutation frequences".. would be biased downwards if 
lambda.est=sum(unique(temp[,list(gene1,gene2,expe)])$expe)
observed=length(unique(c(ym1$case.id,ym2$case.id)))
pval= sum(dpois(0:observed,lambda.est))
print (paste("observed",observed,"est",lambda.est,"pval",pval))








