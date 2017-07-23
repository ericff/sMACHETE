args = commandArgs(trailingOnly = TRUE)
counts.file = args[1]
querydir = args[2]
sequences.present.laml.file = args[3]
report.paths.with.meta.file = args[4]
sbt.results.bodymap.chimera = args[5]
name.to.sample.ids.dir = args[6]
querydate = args[7]


#
source("ChimeraVSMachete.r")
# Above will automatically run analysis.r;

source("FigX_Mach_fusion_per_cancer.r")
source("summary_fig_plots.r")

## Clean up Table 4 for publication; remove identifying info.

library(data.table)
table4.in <- fread("../Table_4_NO_SEQUENCES_all_fdr_post_filter_all.tab")
## table4.in <- fread("../Table_4_NO_SEQUENCES_all_fdr_post_filter_all.tab", stringsAsFactors=FALSE, sep="\t"))
names(table4.in)[3]="AbsPos1Pos2Diff"

write.table(file="../Table_4_NO_IDENTIFYING_INFO_all_fdr_post_filter_all.tab", table4.in[,.(gene1,gene2,AbsPos1Pos2Diff,numReads,maxGeneExp.x, maxGeneExp.y,investigation,fuspair,info,chr1,chr2,strand1,strand2,anomSum,corr_emp_p,Tumour.Types.Somatic..x,Tumour.Types.Somatic..y,upperCI,is.cosmic,g1prod,g2prod) ],quote=F,sep="\t", row.names=F)

source("count_ever_discovered.r")

source("additionalcalcs.R")

source("mutation_merge.r")

# make figure 4A
source("makefigure4.R")
