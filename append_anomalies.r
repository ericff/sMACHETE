## simulated data anlysis
rm(list=ls())

radius=100 ### for Hsieh et al paper,used 1000000
#par(mfrow=c(2,4))
require(data.table)

args = commandArgs(trailingOnly = TRUE)
circ=args[1]
lin=args[2]
fusion=args[3]
naive_knife=args[4]
tablepref=args[5]
outputfilesdir=args[6]

## read in circ and linear and fusion if
if (!is.na(naive_knife)){
    ## read in output files and output new files

    prefixdirsuf=file.path(outputfilesdir, tablepref)

    filein=paste(prefixdirsuf,"all",sep="_")
    fin=data.table(read.delim(paste(filein)))
    ## fin junction has strand twice in middle that isn't there 

    ## need to split junction for g
    g=data.table(read.delim(paste(naive_knife), header=T, sep="\t", skip=1))
    setnames(g,names(g)[1],"junction")
    g[,junction:=gsub("([|])", ":", paste(junction))]
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

    setnames(g,names(g)[1],"fusionInfo2")
    g[,anomSum:=anomaly+decoy]
    ## add in

    setkey(fin,fusionInfo)

    ## wasg[,fusionInfo:=paste(chr1,gene1,pos1,strand1,chr2,gene2,pos2, strand2,type, sep=":"),]
    g[,fusionInfo:=paste(chr1,gene1,pos1,gene2,pos2, type,strand2, sep=":"),]
    setkey(g,fusionInfo)

    g[,fusionInfo:=gsub("([:])", "_", paste(fusionInfo))]
    fin[,fusionInfo:=gsub("([:])", "_", paste(fusionInfo))]

    ## g is naive reports, foregiven if doesn't exist
    merged=merge(fin,g[,list(fusionInfo,anomSum)], by="fusionInfo", all.x=T)	
    merged[is.na(anomSum.y), anomSum.y:=0]
    merged[is.na(anomSum.x), anomSum.x:=0]
    merged[,anomSum:=anomSum.x+anomSum.y]
    merged[,anomSum.x:=NULL]
    merged[,anomSum.y:=NULL]
} ## end if (!is.na(naive_knife)){

write.table(file=paste(prefixdirsuf,"_all_merged_w_anomalies",sep="_"),merged,sep="\t",quote=F, row.names=F)
