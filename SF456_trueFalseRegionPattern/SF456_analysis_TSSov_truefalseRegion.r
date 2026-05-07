## K562 H3K27me3
ChIP <- read.table("TSS3kb_K562/TSS3kSig_GSM608166_ChIP_H3K27me3.bed",row.names=5)
ATAC <- read.table("TSS3kb_K562/TSS3kSig_ATAC.bed",row.names=5)
CnT1raw <- read.table("TSS3kb_K562/TSS3kSig_GSM4308151_CUTTag_H3K27me3_raw.bed",row.names=5)
CnT2raw <- read.table("TSS3kb_K562/TSS3kSig_GSM4308152_CUTTag_H3K27me3_raw.bed",row.names=5)
CnT1score <- read.table("TSS3kb_K562/TSS3kSig_GSM4308151_CUTTag_H3K27me3_score_new.bed",row.names=5)
CnT2score <- read.table("TSS3kb_K562/TSS3kSig_GSM4308152_CUTTag_H3K27me3_score_new.bed",row.names=5)
ATACwithsig <- rownames(ATAC)[which( apply(ATAC[,6:605],1,mean) > 0)]


#true1kb <- intersect(ATACwithsig,rownames(read.table("hg38_TSS1kb_ov_H3K27me3_true.bed",row.names=5)))#intersect(intersect(rownames(activeGene_tmp)[order(activeGene_tmp[,4],decreasing=T)[1:4000]],ATACwithsig), noXchrom)
#false1kb <- intersect(ATACwithsig,rownames(read.table("hg38_TSS1kb_ov_H3K27me3_false.bed",row.names=5)))
true3kb <- intersect(ATACwithsig,rownames(read.table("hg38_TSS3kb_ov_H3K27me3_true.bed",row.names=5)))#intersect(intersect(rownames(activeGene_tmp)[order(activeGene_tmp[,4],decreasing=T)[1:4000]],ATACwithsig), noXchrom)
false3kb <- intersect(ATACwithsig,rownames(read.table("hg38_TSS3kb_ov_H3K27me3_false.bed",row.names=5)))



sitepro <- function(indata, M, COL,YLIM){
  plot(apply(indata,2,mean),ylim=YLIM,col=COL,main=M,type='l',xlab="relative position to TSS (bp)", ylab="aggregate signal",axes=F)
  axis(side=1,at=c(1, 300.5, 600), labels=c("-3k","TSS","3k"))
  axis(side=2,las=2)
  box()
}

heatmap <- function(data0, M, usecolor){
  data<-data0
  data <- c(as.matrix(data))
  data <- sort(data)
  temp<-data[round(c(0.010000,0.5,0.99)*length(data))]
  p20<-temp[1]
  p50<-temp[2]
  p80<-temp[3]
  zmin=p20
  zmax=p80
  #zmin=MIN
  #zmax=MAX
  ColorRamp <- colorRampPalette(c("white",usecolor), bias=1)(10000)   #color list
  ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence
  data0[data0<zmin] <- zmin
  data0[data0>zmax] <- zmax
  ColorRamp_ex <- ColorRamp[round( (min(data0)-zmin)*10000/(zmax-zmin) ) : round( (max(data0)-zmin)*10000/(zmax-zmin) )]
  image(1:ncol(data0), 1:nrow(data0), t(data0), axes=F, col=ColorRamp_ex, xlab="", ylab=paste(zmin,zmax,sep=":"),main=M,useRaster=T,cex.main=1)
  axis(side=2,las=2)
  axis(side=1,at=c(1, 300.5, 600), labels=c("-3k","TSS","3k"))
  box()
}

heatmapFIX <- function(data0, M, usecolor,MAX){
  zmin=0#MIN
  zmax=MAX
  ColorRamp <- colorRampPalette(c("white",usecolor), bias=1)(10000)   #color list
  ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence
  data0[data0<zmin] <- zmin
  data0[data0>zmax] <- zmax
  ColorRamp_ex <- ColorRamp[round( (min(data0)-zmin)*10000/(zmax-zmin) ) : round( (max(data0)-zmin)*10000/(zmax-zmin) )]
  image(1:ncol(data0), 1:nrow(data0), t(data0), axes=F, col=ColorRamp_ex, xlab="", ylab=paste(zmin,zmax,sep=":"),main=M,useRaster=T,cex.main=1)
  axis(side=2,las=2)
  axis(side=1,at=c(1, 300.5, 600), labels=c("-3k","TSS","3k"))
  box()
}



pdf(file="TSS3kb_sitepro_K562_H3K27me3true3kb.pdf",width=3*6,height=3)
par(mfrow=c(1,6),mar=c(4,4,2,2))
sitepro(CnT1raw[true3kb,6:605], "cut1raw true3kb","red",c(0,0.6))
sitepro(CnT2raw[true3kb,6:605], "cut2raw true3kb","red",c(0,0.6))
sitepro(ChIP[true3kb,6:605], "ChIPraw true3kb","blue",c(0,0.4))
sitepro(ATAC[true3kb,6:605], "ATAC true3kb","black",c(0,7))
sitepro(CnT1score[true3kb,6:605], "cut1score true3kb","#DD864E",c(0,1))
sitepro(CnT2score[true3kb,6:605], "cut2score true3kb","#DD864E",c(0,1))
dev.off()

pdf(file="TSS3kb_heatmap_K562_H3K27me3true3kb.pdf",width=3*6,height=6)
par(mfrow=c(1,6),mar=c(4,4,2,2))
heatmapFIX(CnT1raw[true3kb,6:605], "cut1raw true3kb","red",0.6)
heatmapFIX(CnT2raw[true3kb,6:605], "cut2raw true3kb","red",0.6)
heatmapFIX(ChIP[true3kb,6:605], "ChIPraw true3kb","blue",0.3)
heatmapFIX(ATAC[true3kb,6:605], "ATAC true3kb","black",12)
heatmapFIX(CnT1score[true3kb,6:605], "cut1score true3kb","#DD864E",1)
heatmapFIX(CnT2score[true3kb,6:605], "cut2score true3kb","#DD864E",1)
dev.off()



pdf(file="TSS3kb_sitepro_K562_H3K27me3false3kb.pdf",width=3*6,height=3)
par(mfrow=c(1,6),mar=c(4,4,2,2))
sitepro(CnT1raw[false3kb,6:605], "cut1raw false3kb","red",c(0,1.2))
sitepro(CnT2raw[false3kb,6:605], "cut2raw false3kb","red",c(0,1.2))
sitepro(ChIP[false3kb,6:605], "ChIPraw false3kb","blue",c(0,0.42))
sitepro(ATAC[false3kb,6:605], "ATAC false3kb","black",c(0,8))
sitepro(CnT1score[false3kb,6:605], "cut1score false3kb","#DD864E",c(0,1))
sitepro(CnT2score[false3kb,6:605], "cut2score false3kb","#DD864E",c(0,1))
dev.off()

pdf(file="TSS3kb_heatmap_K562_H3K27me3false3kb.pdf",width=3*6,height=6)
par(mfrow=c(1,6),mar=c(4,4,2,2))
heatmapFIX(CnT1raw[false3kb,6:605], "cut1raw false3kb","red",0.6)
heatmapFIX(CnT2raw[false3kb,6:605], "cut2raw false3kb","red",0.6)
heatmapFIX(ChIP[false3kb,6:605], "ChIPraw false3kb","blue",0.3)
heatmapFIX(ATAC[false3kb,6:605], "ATAC false3kb","black",12)
heatmapFIX(CnT1score[false3kb,6:605], "cut1score false3kb","#DD864E",1)
heatmapFIX(CnT2score[false3kb,6:605], "cut2score false3kb","#DD864E",1)
dev.off()


## K562 H3K27ac
ChIP <- read.table("TSS3kb_K562/TSS3kSig_GSM1652918_ChIP_H3K27ac.bed",row.names=5)
ATAC <- read.table("TSS3kb_K562/TSS3kSig_ATAC.bed",row.names=5)
CnT1raw <- read.table("TSS3kb_K562/TSS3kSig_GSM3536514_CUTTag_H3K27ac_raw.bed",row.names=5)
CnT1score <- read.table("TSS3kb_K562/TSS3kSig_GSM3536514_CUTTag_H3K27ac_score.bed",row.names=5)
ATACwithsig <- rownames(ATAC)[which( apply(ATAC[,6:605],1,mean) > 0)]

#true1kb <- intersect(ATACwithsig,rownames(read.table("hg38_TSS1kb_ov_H3K27ac_true.bed",row.names=5)))#intersect(intersect(rownames(activeGene_tmp)[order(activeGene_tmp[,4],decreasing=T)[1:4000]],ATACwithsig), noXchrom)
#false1kb <- intersect(ATACwithsig,rownames(read.table("hg38_TSS1kb_ov_H3K27ac_false.bed",row.names=5)))
true3kb <- intersect(ATACwithsig,rownames(read.table("hg38_TSS3kb_ov_H3K27ac_true.bed",row.names=5)))#intersect(intersect(rownames(activeGene_tmp)[order(activeGene_tmp[,4],decreasing=T)[1:4000]],ATACwithsig), noXchrom)
false3kb <- intersect(ATACwithsig,rownames(read.table("hg38_TSS3kb_ov_H3K27ac_false.bed",row.names=5)))


pdf(file="TSS3kb_sitepro_K562_H3K27actrue3kb.pdf",width=3*4,height=3)
par(mfrow=c(1,4),mar=c(4,4,2,2))
sitepro(CnT1raw[true3kb,6:605], "cut1raw true3kb","red",c(0,3))
sitepro(ChIP[true3kb,6:605], "ChIPraw true3kb","blue",c(0,3))
sitepro(ATAC[true3kb,6:605], "ATAC true3kb","black",c(0,7))
sitepro(CnT1score[true3kb,6:605], "cut1score true3kb","#DD864E",c(0,1))
dev.off()

pdf(file="TSS3kb_heatmap_K562_H3K27actrue3kb.pdf",width=3*4,height=6)
par(mfrow=c(1,4),mar=c(4,4,2,2))
heatmapFIX(CnT1raw[true3kb,6:605], "cut1raw true3kb","red",0.6)
heatmapFIX(ChIP[true3kb,6:605], "ChIPraw true3kb","blue",0.3)
heatmapFIX(ATAC[true3kb,6:605], "ATAC true3kb","black",12)
heatmapFIX(CnT1score[true3kb,6:605], "cut1score true3kb","#DD864E",1)
dev.off()



pdf(file="TSS3kb_sitepro_K562_H3K27acfalse3kb.pdf",width=3*4,height=3)
par(mfrow=c(1,4),mar=c(4,4,2,2))
sitepro(CnT1raw[false3kb,6:605], "cut1raw false3kb","red",c(0,0.5))
sitepro(ChIP[false3kb,6:605], "ChIPraw false3kb","blue",c(0,0.5))
sitepro(ATAC[false3kb,6:605], "ATAC false3kb","black",c(0,7))
sitepro(CnT1score[false3kb,6:605], "cut1score false3kb","#DD864E",c(0,1))
dev.off()

pdf(file="TSS3kb_heatmap_K562_H3K27acfalse3kb.pdf",width=3*4,height=6)
par(mfrow=c(1,4),mar=c(4,4,2,2))
heatmapFIX(CnT1raw[false3kb,6:605], "cut1raw false3kb","red",0.6)
heatmapFIX(ChIP[false3kb,6:605], "ChIPraw false3kb","blue",0.3)
heatmapFIX(ATAC[false3kb,6:605], "ATAC false3kb","black",12)
heatmapFIX(CnT1score[false3kb,6:605], "cut1score false3kb","#DD864E",1)
dev.off()



### K562 H3K9me3

ChIP <- read.table("TSS3kb_K562/TSS3kSig_K562_H3K9me3_ChIP1_raw.bed",row.names=5)
ATAC <- read.table("TSS3kb_K562/TSS3kSig_ATAC.bed",row.names=5)
CnT1raw <- read.table("TSS3kb_K562/TSS3kSig_K562_H3K9me3_CUTTag1_raw.bed",row.names=5)
CnT1score <- read.table("TSS3kb_K562/TSS3kSig_CUTTag1_H3K9me3_score.bed",row.names=5)
CnT2score <- read.table("TSS3kb_K562/TSS3kSig_CUTTag2_H3K9me3_score.bed",row.names=5)
ATACwithsig <- rownames(ATAC)[which( apply(ATAC[,6:605],1,mean) > 0)]

#true1kb <- intersect(ATACwithsig,rownames(read.table("hg38_TSS1kb_ov_H3K27ac_true.bed",row.names=5)))#intersect(intersect(rownames(activeGene_tmp)[order(activeGene_tmp[,4],decreasing=T)[1:4000]],ATACwithsig), noXchrom)
#false1kb <- intersect(ATACwithsig,rownames(read.table("hg38_TSS1kb_ov_H3K27ac_false.bed",row.names=5)))
true3kb <- intersect(ATACwithsig,rownames(read.table("hg38_TSS3kb_ov_H3K9me3_true.bed",row.names=5)))#intersect(intersect(rownames(activeGene_tmp)[order(activeGene_tmp[,4],decreasing=T)[1:4000]],ATACwithsig), noXchrom)
false3kb <- intersect(ATACwithsig,rownames(read.table("hg38_TSS3kb_ov_H3K9me3_false.bed",row.names=5)))


pdf(file="TSS3kb_sitepro_K562_H3K9me3true3kb.pdf",width=3*6,height=3)
par(mfrow=c(1,6),mar=c(4,4,2,2))
sitepro(CnT1raw[true3kb,6:605], "cut1raw true3kb","red",c(0,2.5))
sitepro(ChIP[true3kb,6:605], "ChIPraw true3kb","blue",c(0,0.4))
sitepro(ATAC[true3kb,6:605], "ATAC true3kb","black",c(0,12))
sitepro(CnT1score[true3kb,6:605], "cut1score true3kb","#DD864E",c(0,1))
dev.off()
 
pdf(file="TSS3kb_heatmap_K562_H3K9me3true3kb.pdf",width=3*6,height=6)
par(mfrow=c(1,6),mar=c(4,4,2,2))
heatmapFIX(CnT1raw[true3kb,6:605], "cut1raw true3kb","red",2)
heatmapFIX(ChIP[true3kb,6:605], "ChIPraw true3kb","blue",0.3)
heatmapFIX(ATAC[true3kb,6:605], "ATAC true3kb","black",12)
heatmapFIX(CnT1score[true3kb,6:605], "cut1score true3kb","#DD864E",1)
dev.off()



pdf(file="TSS3kb_sitepro_K562_H3K9me3false3kb.pdf",width=3*6,height=3)
par(mfrow=c(1,6),mar=c(4,4,2,2))
sitepro(CnT1raw[false3kb,6:605], "cut1raw false3kb","red",c(0,2.5))
sitepro(ChIP[false3kb,6:605], "ChIPraw false3kb","blue",c(0,0.4))
sitepro(ATAC[false3kb,6:605], "ATAC false3kb","black",c(0,12))
sitepro(CnT1score[false3kb,6:605], "cut1score false3kb","#DD864E",c(0,1))
dev.off()

pdf(file="TSS3kb_heatmap_K562_H3K9me3false3kb.pdf",width=3*6,height=6)
par(mfrow=c(1,6),mar=c(4,4,2,2))
heatmapFIX(CnT1raw[false3kb,6:605], "cut1raw false3kb","red",2.5)
heatmapFIX(ChIP[false3kb,6:605], "ChIPraw false3kb","blue",0.3)
heatmapFIX(ATAC[false3kb,6:605], "ATAC false3kb","black",12)
heatmapFIX(CnT1score[false3kb,6:605], "cut1score false3kb","#DD864E",1)
dev.off()




#### H3K9me3 K562
H3K9me3_CnT1raw <- read.table("TSS3kSig_K562_H3K9me3_CUTTag1_raw.bed",row.names=5)
H3K9me3_CnT2raw <- read.table("TSS3kSig_K562_H3K9me3_CUTTag2_raw.bed",row.names=5)
H3K9me3_CnT3raw <- read.table("TSS3kSig_K562_H3K9me3_CUTTag3_raw.bed",row.names=5)
H3K9me3_CnT4raw <- read.table("TSS3kSig_K562_H3K9me3_CUTTag4_raw.bed",row.names=5)
H3K9me3_ChIP1 <- read.table("TSS3kSig_K562_H3K9me3_ChIP1_raw.bed",row.names=5)



pdf(file="TSS3kb_H3K9me3_sitepro_active.pdf",width=3*5,height=3)
par(mfrow=c(1,5),mar=c(4,4,2,2))
sitepro(H3K9me3_CnT1raw[activeIDX1,6:605], "cut1raw active","red",c(0,0.42))
sitepro(H3K9me3_CnT2raw[activeIDX1,6:605], "cut2raw active","red",c(0,0.42))
sitepro(H3K9me3_CnT3raw[activeIDX1,6:605], "cut3raw active","red",c(0,0.42))
sitepro(H3K9me3_CnT4raw[activeIDX1,6:605], "cut4raw active","red",c(0,0.42))
sitepro(H3K9me3_ChIP1[activeIDX1,6:605], "ChIP1raw active","red",c(0,0.42))
dev.off()

pdf(file="TSS3kb_H3K9me3_heatmap_active.pdf",width=3*5,height=6)
par(mfrow=c(1,5),mar=c(4,4,2,2))
heatmapFIX(H3K9me3_CnT1raw[activeIDX1,6:605], "cut1raw active","red",0.6)
heatmapFIX(H3K9me3_CnT2raw[activeIDX1,6:605], "cut2raw active","red",0.6)
heatmapFIX(H3K9me3_CnT3raw[activeIDX1,6:605], "cut3raw active","red",0.6)
heatmapFIX(H3K9me3_CnT4raw[activeIDX1,6:605], "cut4raw active","red",0.6)
heatmapFIX(H3K9me3_ChIP1[activeIDX1,6:605], "ChIP1raw active","red",0.6)
dev.off()


pdf(file="TSS3kb_H3K9me3_sitepro_repress.pdf",width=3*5,height=3)
par(mfrow=c(1,5),mar=c(4,4,2,2))
sitepro(H3K9me3_CnT1raw[repressIDX2,6:605], "cut1raw repress","red",c(0,0.42))
sitepro(H3K9me3_CnT2raw[repressIDX2,6:605], "cut2raw repress","red",c(0,0.42))
sitepro(H3K9me3_CnT3raw[repressIDX2,6:605], "cut3raw repress","red",c(0,0.42))
sitepro(H3K9me3_CnT4raw[repressIDX2,6:605], "cut4raw repress","red",c(0,0.42))
sitepro(H3K9me3_ChIP1[repressIDX2,6:605], "ChIP1raw repress","red",c(0,0.42))
dev.off()

pdf(file="TSS3kb_H3K9me3_heatmap_repress.pdf",width=3*5,height=6)
par(mfrow=c(1,5),mar=c(4,4,2,2))
heatmapFIX(H3K9me3_CnT1raw[repressIDX2,6:605], "cut1raw repress","red",0.6)
heatmapFIX(H3K9me3_CnT2raw[repressIDX2,6:605], "cut2raw repress","red",0.6)
heatmapFIX(H3K9me3_CnT3raw[repressIDX2,6:605], "cut3raw repress","red",0.6)
heatmapFIX(H3K9me3_CnT4raw[repressIDX2,6:605], "cut4raw repress","red",0.6)
heatmapFIX(H3K9me3_ChIP1[repressIDX2,6:605], "ChIP1raw repress","red",0.6)
dev.off()

