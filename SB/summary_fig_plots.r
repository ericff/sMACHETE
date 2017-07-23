require(data.table)
total.machete.counts=fread(counts.file)
names(total.machete.counts)=c("short","name","macheteCount","bloomCount")

pdf(file="../plots_from_summary_fig_plots.pdf")
par(mfrow=c(2,2))
barplot(total.machete.counts$macheteCount, names.arg=total.machete.counts$short,horiz=T)
barplot(total.machete.counts$bloomCount, names.arg=total.machete.counts$short,horiz=T)
dev.off()
