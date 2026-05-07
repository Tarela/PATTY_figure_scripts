setwd("/Users/shengenhu/Dropbox/CUTTag/Results/noX_rerun/otherHM/")

data <- read.table("otherActiveHM_K562corEXP.txt",header=T)
rownames(data) <- data[,"GSMID"]
###### F5
useCol <- c( "red", "darkred","#DD864E")  # adjust as needed

bar3plot <- function(userow, YL){
  plotdata <- as.numeric(data[userow,c("m2","m2IgG","PATTY")])
  barplot(plotdata,las=2,col=useCol,ylim=YL)
  abline(h=0)
  box(lwd=1)
}
pdf(file="H3K4me3_raw.pdf",width=2*6,height=3)
par(mfrow=c(1,6))
bar3plot("GSM3536518", c(0,0.75))
bar3plot("GSM3680226", c(0,0.75))
bar3plot("GSM4308159", c(0,0.75))
bar3plot("GSM4308168", c(0,0.75))
bar3plot("GSM4835703", c(0,0.75))
bar3plot("GSM4842199", c(0,0.75))
dev.off()

pdf(file="H3K4me1_raw.pdf",width=2*5,height=3)
par(mfrow=c(1,5))
bar3plot("GSM3536516", c(0,0.65))
bar3plot("GSM3680223", c(0,0.65))
bar3plot("GSM3680224", c(0,0.65))
bar3plot("GSM4835702", c(0,0.65))
bar3plot("GSM4842197", c(0,0.65))
dev.off()

pdf(file="H3K36me3_raw.pdf",width=2*1,height=3)
par(mfrow=c(1,1))
bar3plot("GSM4835701", c(0,0.25))
dev.off()





