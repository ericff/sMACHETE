library(data.table)
tablex.count.file=fread("Table_X_pan_cancer_list_w_info.tab")
tablex.count.file[,Fusion_pair:=paste(gene1,gene2,sep="-")]

m=fread("../chimera_fusions_detected_in_samples_common_toMACH.tab", header=FALSE, skip=1)
setnames(m,"V32","Fusion_pair")
all.f=merge(tablex.count.file,m,by="Fusion_pair",all.x=T,all.y=T)


## 8 for chimerseq
unique(all.f[is.na(gene1) & Fusion_pair %like% "HLA"]$Fusion_pair)
print("length(unique(all.f[is.na(gene1) & Fusion_pair %like% \"HLA\"]$Fusion_pair))")
print(length(unique(all.f[is.na(gene1) & Fusion_pair %like% "HLA"]$Fusion_pair)))

## none for smachete
unique(all.f[is.na(V1) & Fusion_pair %like% "HLA"]$Fusion_pair)
print("length(unique(all.f[is.na(V1) & Fusion_pair %like% \"HLA\"]$Fusion_pair))")
print(length(unique(all.f[is.na(V1) & Fusion_pair %like% "HLA"]$Fusion_pair)))


## only called by chimerseq
print("length(unique(all.f[is.na(gene1) & Fusion_pair %like% \"\"]$Fusion_pair))")
print(length(unique(all.f[is.na(gene1) & Fusion_pair %like% ""]$Fusion_pair)))
## [1] 466

## only by sMACHETE:
print("length(unique(all.f[is.na(V1) & Fusion_pair %like% \"\"]$Fusion_pair))")
print(length(unique(all.f[is.na(V1) & Fusion_pair %like% ""]$Fusion_pair)))
## 550

## common
print("length(unique(all.f[(!(is.na(gene1) | is.na(V1))) & Fusion_pair %like% \"\"]$Fusion_pair)")
print(length(unique(all.f[(!(is.na(gene1) | is.na(V1))) & Fusion_pair %like% ""]$Fusion_pair)))
## 214

## called by machete:
print("length(unique(m$Fusion_pair)): ")
print(length(unique(m$Fusion_pair)))
## 680
print("length(unique(tablex.count.file$Fusion_pair))")
print(length(unique(tablex.count.file$Fusion_pair)))
## 764

