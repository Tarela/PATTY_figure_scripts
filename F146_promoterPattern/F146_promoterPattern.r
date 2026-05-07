ChIP <- read.table("TSS3kSig_GSM608166_ChIP_H3K27me3.bed",row.names=5)
ATAC <- read.table("TSS3kSig_ATAC.bed",row.names=5)
CnT1raw <- read.table("TSS3kSig_GSM4308151_CUTTag_H3K27me3_raw.bed",row.names=5)
CnT2raw <- read.table("TSS3kSig_GSM4308152_CUTTag_H3K27me3_raw.bed",row.names=5)
CnT1score <- read.table("TSS3kSig_GSM4308151_CUTTag_H3K27me3_score_new.bed",row.names=5)
CnT2score <- read.table("TSS3kSig_GSM4308152_CUTTag_H3K27me3_score_new.bed",row.names=5)
GRO <- read.table("TSS3kbSig_GRO_GSM1480325_SEuniq.bed",row.names=5)
RNAPIIs5p <- read.table("TSS3kbSig_RNAPIIs5p_GSM803443_SEuniq.bed",row.names=5)
RNAPII1 <- read.table("TSS3kbSig_RNAPII_GSM2423412_rdSE.bed",row.names=5)
RNAPII2 <- read.table("TSS3kbSig_RNAPII_GSM2423413_rdSE.bed",row.names=5)

hs1raw <- read.table("TSS3kSig_GSM4308145_CUTTag_H3K27me3_raw.bed",row.names=5)
hs2raw <- read.table("TSS3kSig_GSM4308146_CUTTag_H3K27me3_raw.bed",row.names=5)
hs3raw <- read.table("TSS3kSig_GSM4308147_CUTTag_H3K27me3_raw.bed",row.names=5)
hs4raw <- read.table("TSS3kSig_GSM4308148_CUTTag_H3K27me3_raw.bed",row.names=5)

activeGene_tmp <- read.table("signal_on_activeGene/active_repress_gene/K562_RNA_8K_promoter1kb.bed",row.names=4)
repressGene_tmp <- read.table("signal_on_activeGene/active_repress_gene/K562_RNA_none_promoter1kb.bed",row.names=4)
TSSovPeak <- read.table("signal_on_activeGene/active_repress_gene/hg38geneTSS_ovH3K27me3Peak.bed",row.names=5)
TSS3kbovPeak <- read.table("signal_on_activeGene/active_repress_gene/hg38geneTSS3kb_ovH3K27me3Peak.bed",row.names=5)
TSS3kbovChIPATACPeak <- read.table("signal_on_activeGene/active_repress_gene/hg38geneTSS3kb_ovChIPATACpeak.bed",row.names=5)

TSS3kbovReads <- read.table("signal_on_activeGene/active_repress_gene/hg38geneTSS3kb_ovReads.bed",row.names=5)


ATACwithsig <- rownames(ATAC)[which( apply(ATAC[,6:605],1,mean) > 0)]

TSS3kb_noChIP_withATAC <- rownames(TSS3kbovChIPATACPeak)[which(TSS3kbovChIPATACPeak[,6] == 0 & TSS3kbovChIPATACPeak[,7]>0)]
TSS3kb_withH3K27me3 <- rownames(TSS3kbovPeak)[which(apply(TSS3kbovPeak[,6:8],1,sum) > 0)]

#TSS3kb_withChIPreads5 <- rownames(TSS3kbovReads)[which(TSS3kbovReads[,6] > 5)]
#TSS3kb_withChIPreads10 <- rownames(TSS3kbovReads)[which(TSS3kbovReads[,6] > 10)]
TSS3kb_withChIPreads20 <- rownames(TSS3kbovReads)[which(TSS3kbovReads[,6] > 20)]
TSS3kb_withReads5 <- rownames(TSS3kbovReads)[which(apply(TSS3kbovReads[,6:8],1,min) > 5)]
noXchrom <- rownames(TSS3kbovReads)[which(TSS3kbovReads[,1] != "chrX")]

#TSS_withChIP <- rownames(TSS3kbovChIPATACPeak)[which(TSS3kbovChIPATACPeak[,6] > 0)]

activeIDX1 <- intersect(intersect(rownames(activeGene_tmp)[order(activeGene_tmp[,4],decreasing=T)[1:4000]],ATACwithsig), noXchrom)
#activeIDX2 <- intersect(rownames(activeGene_tmp)[order(activeGene_tmp[,4],decreasing=T)[1:4200]],TSS3kb_noChIP_withATAC)
#repressIDX1 <- intersect(rownames(repressGene_tmp),TSS3kb_withChIPreads20)
repressIDX2 <- intersect(intersect(rownames(repressGene_tmp),TSS3kb_withReads5),noXchrom)


hg38_TSS1kb <- read.table("hg38_gene_annotation_geneID_LenOrder_TSS1kb.bed",row.names=5)
hg38_TSS3kb <- read.table("hg38_gene_annotation_geneID_LenOrder_TSS3kb.bed",row.names=5)
hg38_TSS5kb <- read.table("hg38_gene_annotation_geneID_LenOrder_TSS5kb.bed",row.names=5)


idx <- intersect(intersect(rownames(hg38_TSS5kb), rownames(repressGene_tmp)), noXchrom)

printGeneBed(hg38_TSS3kb[idx,], "K562_repressALL_TSS3kb.bed")
printGeneBed(hg38_TSS5kb[idx,], "K562_repressALL_TSS5kb.bed")


printGeneBed <- function(indata,outf){
  outdat <- cbind(indata[,1:4], rownames(indata), indata[,5])
  write.table(outdat, file=outf,sep="\t",quote=F,row.names=F,col.names=F)
}
printGeneBed(hg38_TSS1kb[activeIDX1,], "K562_active_TSS1kb.bed")
printGeneBed(hg38_TSS3kb[activeIDX1,], "K562_active_TSS3kb.bed")
printGeneBed(hg38_TSS5kb[activeIDX1,], "K562_active_TSS5kb.bed")
printGeneBed(hg38_TSS1kb[repressIDX2,], "K562_repress_TSS1kb.bed")
printGeneBed(hg38_TSS3kb[repressIDX2,], "K562_repress_TSS3kb.bed")
printGeneBed(hg38_TSS5kb[repressIDX2,], "K562_repress_TSS5kb.bed")



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


pdf(file="TSS3kb_sitepro_active_newPATTY.pdf",width=3*6,height=3)
par(mfrow=c(1,6),mar=c(4,4,2,2))
sitepro(CnT1raw[activeIDX1,6:605], "cut1raw active","red",c(0,0.42))
sitepro(CnT2raw[activeIDX1,6:605], "cut2raw active","red",c(0,0.42))
sitepro(ChIP[activeIDX1,6:605], "ChIPraw active","blue",c(0,0.42))
sitepro(ATAC[activeIDX1,6:605], "ATAC active","black",c(0,7))
sitepro(CnT1score[activeIDX1,6:605], "cut1score active","#DD864E",c(0,1))
sitepro(CnT2score[activeIDX1,6:605], "cut2score active","#DD864E",c(0,1))
dev.off()

pdf(file="TSS3kb_heatmap_active_newPATTY.pdf",width=3*6,height=6)
par(mfrow=c(1,6),mar=c(4,4,2,2))
heatmapFIX(CnT1raw[activeIDX1,6:605], "cut1raw active","red",0.6)
heatmapFIX(CnT2raw[activeIDX1,6:605], "cut2raw active","red",0.6)
heatmapFIX(ChIP[activeIDX1,6:605], "ChIPraw active","blue",0.3)
heatmapFIX(ATAC[activeIDX1,6:605], "ATAC active","black",12)
heatmapFIX(CnT1score[activeIDX1,6:605], "cut1score active","#DD864E",1)
heatmapFIX(CnT2score[activeIDX1,6:605], "cut2score active","#DD864E",1)
dev.off()


pdf(file="TSS3kb_sitepro_repress2_newPATTY.pdf",width=3*6,height=3)
par(mfrow=c(1,6),mar=c(4,4,2,2))
sitepro(CnT1raw[repressIDX2,6:605], "cut1raw repress","red",c(0,0.42))
sitepro(CnT2raw[repressIDX2,6:605], "cut2raw repress","red",c(0,0.42))
sitepro(ChIP[repressIDX2,6:605], "ChIPraw repress","blue",c(0,0.42))
sitepro(ATAC[repressIDX2,6:605], "ATAC repress","black",c(0,7))
sitepro(CnT1score[repressIDX2,6:605], "cut1score repress","#DD864E",c(0,1))
sitepro(CnT2score[repressIDX2,6:605], "cut2score repress","#DD864E",c(0,1))
dev.off()

pdf(file="TSS3kb_heatmap_repress2_newPATTY.pdf",width=3*6,height=6)
par(mfrow=c(1,6),mar=c(4,4,2,2))
heatmapFIX(CnT1raw[repressIDX2,6:605], "cut1raw repress","red",0.6)
heatmapFIX(CnT2raw[repressIDX2,6:605], "cut2raw repress","red",0.6)
heatmapFIX(ChIP[repressIDX2,6:605], "ChIPraw repress","blue",0.3)
heatmapFIX(ATAC[repressIDX2,6:605], "ATAC repress","black",12)
heatmapFIX(CnT1score[repressIDX2,6:605], "cut1score repress","#DD864E",1)
heatmapFIX(CnT2score[repressIDX2,6:605], "cut2score repress","#DD864E",1)
dev.off()



outdata1 <- rbind(apply(CnT1raw[activeIDX1,6:605],2,mean),
                 apply(CnT2raw[activeIDX1,6:605],2,mean),
                 apply(ChIP[activeIDX1,6:605],2,mean),
                 apply(ATAC[activeIDX1,6:605],2,mean),
                 apply(CnT1score[activeIDX1,6:605],2,mean),
                 apply(CnT2score[activeIDX1,6:605],2,mean))
rownames(outdata1) <- c("CnT1raw","CnT2raw","ChIP","ATAC","CnT1score","CnT2score")
write.table(outdata1, file="F4_A-F_sourceData.txt",row.names=T,col.names=F,sep="\t",quote=F)

outdata2 <- rbind(apply(CnT1raw[repressIDX2,6:605],2,mean),
                 apply(CnT2raw[repressIDX2,6:605],2,mean),
                 apply(ChIP[repressIDX2,6:605],2,mean),
                  apply(ATAC[repressIDX2,6:605],2,mean),
                 apply(CnT1score[repressIDX2,6:605],2,mean),
                 apply(CnT2score[repressIDX2,6:605],2,mean))
rownames(outdata2) <- c("CnT1raw","CnT2raw","ChIP","ATAC","CnT1score","CnT2score")
write.table(outdata2, file="SF3_A-F_sourceData.txt",row.names=T,col.names=F,,sep="\t",quote=F)



### for high salt

pdf(file="TSS3kb_sitepro_active_H3K27me3_highsalt.pdf",width=3*6,height=3)
par(mfrow=c(1,6),mar=c(4,4,2,2))
sitepro(hs1raw[activeIDX1,6:605], "hs1 active","red",c(0,0.1))
sitepro(hs2raw[activeIDX1,6:605], "hs2 active","red",c(0,0.1))
sitepro(hs3raw[activeIDX1,6:605], "hs3 active","red",c(0,0.1))
sitepro(hs4raw[activeIDX1,6:605], "hs4 active","red",c(0,0.1))
sitepro(ChIP[activeIDX1,6:605], "ChIPraw active","blue",c(0,0.1))
sitepro(ATAC[activeIDX1,6:605], "ATAC active","black",c(0,7))
dev.off()

pdf(file="TSS3kb_heatmap_active_H3K27me3_highsalt.pdf",width=3*6,height=6)
par(mfrow=c(1,6),mar=c(4,4,2,2))
heatmapFIX(hs1raw[activeIDX1,6:605], "hs1 active","red",0.03)
heatmapFIX(hs2raw[activeIDX1,6:605], "hs2 active","red",0.03)
heatmapFIX(hs3raw[activeIDX1,6:605], "hs3 active","red",0.03)
heatmapFIX(hs4raw[activeIDX1,6:605], "hs4 active","red",0.03)
heatmapFIX(ChIP[activeIDX1,6:605], "ChIPraw active","blue",0.3)
heatmapFIX(ATAC[activeIDX1,6:605], "ATAC active","black",12)
dev.off()










### for RNAPII







pdf(file="TSS3kb_sitepro_active_newPATTY_RNAPII.pdf",width=3*4,height=3)
par(mfrow=c(1,4),mar=c(4,4,2,2))
sitepro(GRO[activeIDX1,6:605], "GRO active","red",c(0,2*1.2))
sitepro(RNAPIIs5p[activeIDX1,6:605], "RNAPIIs5p active","darkblue",c(0,3*1.2))
sitepro(RNAPII1[activeIDX1,6:605], "RNAPII1 active","blue",c(0,4*1.2))
sitepro(RNAPII2[activeIDX1,6:605], "RNAPII2 active","blue",c(0,4*1.2))
dev.off()

pdf(file="TSS3kb_heatmap_active_newPATTY_RNAPII.pdf",width=3*4,height=6)
par(mfrow=c(1,4),mar=c(4,4,2,2))
heatmapFIX(GRO[activeIDX1,6:605], "GRO active","red",2)
heatmapFIX(RNAPIIs5p[activeIDX1,6:605], "RNAPIIs5p active","darkblue",3)
heatmapFIX(RNAPII1[activeIDX1,6:605], "RNAPII1 active","blue",4)
heatmapFIX(RNAPII2[activeIDX1,6:605], "RNAPII2 active","blue",4)
dev.off()



pdf(file="TSS3kb_sitepro_repress2_newPATTY_RNAPII.pdf",width=3*4,height=3)
par(mfrow=c(1,4),mar=c(4,4,2,2))
sitepro(GRO[repressIDX2,6:605], "GRO repress","red",c(0,2*1.2))
sitepro(RNAPIIs5p[repressIDX2,6:605], "RNAPIIs5p repress","darkblue",c(0,3*1.2))
sitepro(RNAPII1[repressIDX2,6:605], "RNAPII1 repress","blue",c(0,4*1.2))
sitepro(RNAPII2[repressIDX2,6:605], "RNAPII2 repress","blue",c(0,4*1.2))
dev.off()

pdf(file="TSS3kb_heatmap_repress2_newPATTY_RNAPII.pdf",width=3*4,height=6)
par(mfrow=c(1,4),mar=c(4,4,2,2))
heatmapFIX(GRO[repressIDX2,6:605], "GRO repress","red",2)
heatmapFIX(RNAPIIs5p[repressIDX2,6:605], "RNAPIIs5p repress","darkblue",3)
heatmapFIX(RNAPII1[repressIDX2,6:605], "RNAPII1 repress","blue",4)
heatmapFIX(RNAPII2[repressIDX2,6:605], "RNAPII2 repress","blue",4)
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







#### for H3K36me3, K562

H3K36me3_ChIP1 <- read.table("TSS3kSig_GSM646439_ChIP_H3K36me3_raw.bed",row.names=5)

H3K36me3_CnT1score <- read.table("TSS3kSig_GSM4835701_CUTTag_H3K36me3_score.bed",row.names=5)
H3K36me3_CnT2score <- read.table("TSS3kSig_GSM4308162_CUTTag_H3K36me3_score.bed",row.names=5)
H3K36me3_CnT3score <- read.table("TSS3kSig_GSM4308163_CUTTag_H3K36me3_score.bed",row.names=5)
H3K36me3_CnT4score <- read.table("TSS3kSig_GSM4308161_CUTTag_H3K36me3_score.bed",row.names=5)
H3K36me3_CnT5score <- read.table("TSS3kSig_GSM4308160_CUTTag_H3K36me3_score.bed",row.names=5)
H3K36me3_CnT6score <- read.table("TSS3kSig_GSM4842195_CUTTag_H3K36me3_score.bed",row.names=5)
H3K36me3_CnT7score <- read.table("TSS3kSig_GSM4797829_CUTTag_H3K36me3_score.bed",row.names=5)
H3K36me3_CnT8score <- read.table("TSS3kSig_GSM4797830_CUTTag_H3K36me3_score.bed",row.names=5)

H3K36me3_CnT1raw <- read.table("TSS3kSig_GSM4835701_CUTTag_H3K36me3_raw.bed",row.names=5)
H3K36me3_CnT2raw <- read.table("TSS3kSig_GSM4308162_CUTTag_H3K36me3_raw.bed",row.names=5)
H3K36me3_CnT3raw <- read.table("TSS3kSig_GSM4308163_CUTTag_H3K36me3_raw.bed",row.names=5)
H3K36me3_CnT4raw <- read.table("TSS3kSig_GSM4308161_CUTTag_H3K36me3_raw.bed",row.names=5)
H3K36me3_CnT5raw <- read.table("TSS3kSig_GSM4308160_CUTTag_H3K36me3_raw.bed",row.names=5)
H3K36me3_CnT6raw <- read.table("TSS3kSig_GSM4842195_CUTTag_H3K36me3_raw.bed",row.names=5)
H3K36me3_CnT7raw <- read.table("TSS3kSig_GSM4797829_CUTTag_H3K36me3_raw.bed",row.names=5)
H3K36me3_CnT8raw <- read.table("TSS3kSig_GSM4797830_CUTTag_H3K36me3_raw.bed",row.names=5)
#H3K36me3_ChIP2 <- read.table("TSS3kSig_GSM646438_ChIP_H3K36me3_raw.bed",row.names=5)
#H3K36me3_ChIP3 <- read.table("TSS3kSig_GSM733714_ChIP_H3K36me3_raw.bed",row.names=5)
#H3K36me3_ChIP4 <- read.table("TSS3kSig_GSM945302_ChIP_H3K36me3_raw.bed",row.names=5)
#H3K36me3_ChIP5 <- read.table("TSS3kSig_GSM1782705_ChIP_H3K36me3_raw.bed",row.names=5)


pdf(file="TSS3kb_sitepro_active_H3K36me3.pdf",width=3*8,height=3)
par(mfrow=c(1,8),mar=c(4,4,2,2))
sitepro(H3K36me3_CnTraw[activeIDX1,6:605], "cutraw active","red",c(0,0.42))
sitepro(ATAC[activeIDX1,6:605], "ATAC active","black",c(0,8))
sitepro(H3K36me3_CnTscore[activeIDX1,6:605], "cutscore active","#DD864E",c(0,1))
sitepro(H3K36me3_ChIP1[activeIDX1,6:605], "ChIP1raw active","blue",c(0,0.42))
sitepro(H3K36me3_ChIP2[activeIDX1,6:605], "ChIP2raw active","blue",c(0,0.42))
sitepro(H3K36me3_ChIP3[activeIDX1,6:605], "ChIP3raw active","blue",c(0,0.42))
sitepro(H3K36me3_ChIP4[activeIDX1,6:605], "ChIP4raw active","blue",c(0,0.42))
sitepro(H3K36me3_ChIP5[activeIDX1,6:605], "ChIP5raw active","blue",c(0,0.42))
dev.off()

pdf(file="TSS3kb_heatmap_active_H3K36me3.pdf",width=3*8,height=6)
par(mfrow=c(1,8),mar=c(4,4,2,2))
heatmapFIX(H3K36me3_CnTraw[activeIDX1,6:605], "cutraw active","red",0.6)
heatmapFIX(ATAC[activeIDX1,6:605], "ATAC active","black",12)
heatmapFIX(H3K36me3_CnTscore[activeIDX1,6:605], "cutscore active","#DD864E",1)
heatmapFIX(H3K36me3_ChIP1[activeIDX1,6:605], "ChIP1raw active","blue",0.3)
heatmapFIX(H3K36me3_ChIP2[activeIDX1,6:605], "ChIP2raw active","blue",0.3)
heatmapFIX(H3K36me3_ChIP3[activeIDX1,6:605], "ChIP3raw active","blue",0.3)
heatmapFIX(H3K36me3_ChIP4[activeIDX1,6:605], "ChIP4raw active","blue",0.3)
heatmapFIX(H3K36me3_ChIP5[activeIDX1,6:605], "ChIP5raw active","blue",0.3)
dev.off()



outdata1 <- rbind(apply(H3K36me3_CnT1raw[activeIDX1,6:605],2,mean),
                 apply(H3K36me3_ChIP1[activeIDX1,6:605],2,mean),
                 apply(ATAC[activeIDX1,6:605],2,mean),
                 apply(H3K36me3_CnT1score[activeIDX1,6:605],2,mean))
rownames(outdata1) <- c("CnT1raw","ChIP","ATAC","CnT1score")
write.table(outdata1, file="SF11_d-g_sourceData.txt",row.names=T,col.names=F,sep="\t",quote=F)



pdf(file="TSS3kb_sitepro_repress2_H3K36me3.pdf",width=3*8,height=3)
par(mfrow=c(1,8),mar=c(4,4,2,2))
sitepro(H3K36me3_CnTraw[repressIDX2,6:605], "cutraw repress","red",c(0,0.42))
sitepro(ATAC[repressIDX2,6:605], "ATAC repress","black",c(0,8))
sitepro(H3K36me3_CnTscore[repressIDX2,6:605], "cutscore repress","#DD864E",c(0,1))
sitepro(H3K36me3_ChIP1[repressIDX2,6:605], "ChIP1raw repress","blue",c(0,0.42))
sitepro(H3K36me3_ChIP2[repressIDX2,6:605], "ChIP2raw repress","blue",c(0,0.42))
sitepro(H3K36me3_ChIP3[repressIDX2,6:605], "ChIP3raw repress","blue",c(0,0.42))
sitepro(H3K36me3_ChIP4[repressIDX2,6:605], "ChIP4raw repress","blue",c(0,0.42))
sitepro(H3K36me3_ChIP5[repressIDX2,6:605], "ChIP5raw repress","blue",c(0,0.42))
dev.off()

pdf(file="TSS3kb_heatmap_repress2_H3K36me3.pdf",width=3*8,height=6)
par(mfrow=c(1,8),mar=c(4,4,2,2))
heatmapFIX(H3K36me3_CnTraw[repressIDX2,6:605], "cutraw repress","red",0.6)
heatmapFIX(ATAC[repressIDX2,6:605], "ATAC repress","black",12)
heatmapFIX(H3K36me3_CnTscore[repressIDX2,6:605], "cutscore repress","#DD864E",1)
heatmapFIX(H3K36me3_ChIP1[repressIDX2,6:605], "ChIP1raw repress","blue",0.3)
heatmapFIX(H3K36me3_ChIP2[repressIDX2,6:605], "ChIP2raw repress","blue",0.3)
heatmapFIX(H3K36me3_ChIP3[repressIDX2,6:605], "ChIP3raw repress","blue",0.3)
heatmapFIX(H3K36me3_ChIP4[repressIDX2,6:605], "ChIP4raw repress","blue",0.3)
heatmapFIX(H3K36me3_ChIP5[repressIDX2,6:605], "ChIP5raw repress","blue",0.3)
dev.off()







pdf(file="TSS3kb_sitepro_active_H3K36me3_ex.pdf",width=3*18,height=3)
par(mfrow=c(1,18),mar=c(4,4,2,2))
sitepro(ATAC[activeIDX1,6:605], "ATAC active","black",c(0,8))
sitepro(H3K36me3_ChIP1[activeIDX1,6:605], "ChIP1raw active","blue",c(0,0.42))
sitepro(H3K36me3_CnT1raw[activeIDX1,6:605], "cutraw GSM4835701","red",c(0,0.42))
sitepro(H3K36me3_CnT2raw[activeIDX1,6:605], "cutraw GSM4308162","red",c(0,0.42))
sitepro(H3K36me3_CnT3raw[activeIDX1,6:605], "cutraw GSM4308163","red",c(0,0.42))
sitepro(H3K36me3_CnT4raw[activeIDX1,6:605], "cutraw GSM4308161","red",c(0,0.42))
sitepro(H3K36me3_CnT5raw[activeIDX1,6:605], "cutraw GSM4308160","red",c(0,0.42))
sitepro(H3K36me3_CnT6raw[activeIDX1,6:605], "cutraw GSM4842195","red",c(0,0.42))
sitepro(H3K36me3_CnT7raw[activeIDX1,6:605], "cutraw GSM4797829","red",c(0,0.42))
sitepro(H3K36me3_CnT8raw[activeIDX1,6:605], "cutraw GSM4797830","red",c(0,0.42))
sitepro(H3K36me3_CnT1score[activeIDX1,6:605], "patty GSM4835701","#DD864E",c(0,1))
sitepro(H3K36me3_CnT2score[activeIDX1,6:605], "patty GSM4308162","#DD864E",c(0,1))
sitepro(H3K36me3_CnT3score[activeIDX1,6:605], "patty GSM4308163","#DD864E",c(0,1))
sitepro(H3K36me3_CnT4score[activeIDX1,6:605], "patty GSM4308161","#DD864E",c(0,1))
sitepro(H3K36me3_CnT5score[activeIDX1,6:605], "patty GSM4308160","#DD864E",c(0,1))
sitepro(H3K36me3_CnT6score[activeIDX1,6:605], "patty GSM4842195","#DD864E",c(0,1))
sitepro(H3K36me3_CnT7score[activeIDX1,6:605], "patty GSM4797829","#DD864E",c(0,1))
sitepro(H3K36me3_CnT8score[activeIDX1,6:605], "patty GSM4797830","#DD864E",c(0,1))
dev.off()

pdf(file="TSS3kb_heatmap_active_H3K36me3_ex.pdf",width=3*18,height=6)
par(mfrow=c(1,18),mar=c(4,4,2,2))
heatmapFIX(ATAC[activeIDX1,6:605], "ATAC active","black",12)
heatmapFIX(H3K36me3_ChIP1[activeIDX1,6:605], "ChIP1raw active","blue",0.3)
heatmapFIX(H3K36me3_CnT1raw[activeIDX1,6:605], "cutraw GSM4835701","red",0.6)
heatmapFIX(H3K36me3_CnT2raw[activeIDX1,6:605], "cutraw GSM4308162","red",0.6)
heatmapFIX(H3K36me3_CnT3raw[activeIDX1,6:605], "cutraw GSM4308163","red",0.6)
heatmapFIX(H3K36me3_CnT4raw[activeIDX1,6:605], "cutraw GSM4308161","red",0.6)
heatmapFIX(H3K36me3_CnT5raw[activeIDX1,6:605], "cutraw GSM4308160","red",0.6)
heatmapFIX(H3K36me3_CnT6raw[activeIDX1,6:605], "cutraw GSM4842195","red",0.6)
heatmapFIX(H3K36me3_CnT7raw[activeIDX1,6:605], "cutraw GSM4797829","red",0.6)
heatmapFIX(H3K36me3_CnT8raw[activeIDX1,6:605], "cutraw GSM4797830","red",0.6)
heatmapFIX(H3K36me3_CnT1score[activeIDX1,6:605], "patty GSM4835701","#DD864E",1)
heatmapFIX(H3K36me3_CnT2score[activeIDX1,6:605], "patty GSM4308162","#DD864E",1)
heatmapFIX(H3K36me3_CnT3score[activeIDX1,6:605], "patty GSM4308163","#DD864E",1)
heatmapFIX(H3K36me3_CnT4score[activeIDX1,6:605], "patty GSM4308161","#DD864E",1)
heatmapFIX(H3K36me3_CnT5score[activeIDX1,6:605], "patty GSM4308160","#DD864E",1)
heatmapFIX(H3K36me3_CnT6score[activeIDX1,6:605], "patty GSM4842195","#DD864E",1)
heatmapFIX(H3K36me3_CnT7score[activeIDX1,6:605], "patty GSM4797829","#DD864E",1)
heatmapFIX(H3K36me3_CnT8score[activeIDX1,6:605], "patty GSM4797830","#DD864E",1)
dev.off()





pdf(file="TSS3kb_sitepro_active_H3K36me3_3kb.pdf",width=3*4,height=3)
par(mfrow=c(1,4),mar=c(4,4,2,2))
sitepro(H3K36me3_CnT1raw[activeIDX1,6:605], "cutraw GSM4835701","red",c(0,0.42))
sitepro(H3K36me3_ChIP1[activeIDX1,6:605], "ChIP1raw active","blue",c(0,0.42))
sitepro(ATAC[activeIDX1,6:605], "ATAC active","black",c(0,8))
sitepro(H3K36me3_CnT1score[activeIDX1,6:605], "patty GSM4835701","#DD864E",c(0,1))
dev.off()

pdf(file="TSS3kb_heatmap_active_H3K36me3_3kb.pdf",width=3*4,height=6)
par(mfrow=c(1,4),mar=c(4,4,2,2))
heatmapFIX(H3K36me3_CnT1raw[activeIDX1,6:605], "cutraw GSM4835701","red",0.6)
heatmapFIX(H3K36me3_ChIP1[activeIDX1,6:605], "ChIP1raw active","blue",0.3)
heatmapFIX(ATAC[activeIDX1,6:605], "ATAC active","black",12)
heatmapFIX(H3K36me3_CnT1score[activeIDX1,6:605], "patty GSM4835701","#DD864E",1)
dev.off()


pdf(file="TSS3kb_sitepro_active_H3K36me3_2kb.pdf",width=3*4,height=3)
par(mfrow=c(1,4),mar=c(4,4,2,2))
sitepro(H3K36me3_CnT1raw[activeIDX1,106:505], "cutraw GSM4835701","red",c(0,0.42))
sitepro(H3K36me3_ChIP1[activeIDX1,106:505], "ChIP1raw active","blue",c(0,0.42))
sitepro(ATAC[activeIDX1,106:505], "ATAC active","black",c(0,8))
sitepro(H3K36me3_CnT1score[activeIDX1,106:505], "patty GSM4835701","#DD864E",c(0,1))
dev.off()

pdf(file="TSS3kb_heatmap_active_H3K36me3_2kb.pdf",width=3*4,height=6)
par(mfrow=c(1,4),mar=c(4,4,2,2))
heatmapFIX(H3K36me3_CnT1raw[activeIDX1,106:505], "cutraw GSM4835701","red",0.6)
heatmapFIX(H3K36me3_ChIP1[activeIDX1,106:505], "ChIP1raw active","blue",0.3)
heatmapFIX(ATAC[activeIDX1,106:505], "ATAC active","black",12)
heatmapFIX(H3K36me3_CnT1score[activeIDX1,106:505], "patty GSM4835701","#DD864E",1)
dev.off()


pdf(file="TSS3kb_sitepro_active_H3K36me3_1kb.pdf",width=3*4,height=3)
par(mfrow=c(1,4),mar=c(4,4,2,2))
sitepro(H3K36me3_CnT1raw[activeIDX1,206:405], "cutraw GSM4835701","red",c(0,0.42))
sitepro(H3K36me3_ChIP1[activeIDX1,206:405], "ChIP1raw active","blue",c(0,0.42))
sitepro(ATAC[activeIDX1,206:405], "ATAC active","black",c(0,8))
sitepro(H3K36me3_CnT1score[activeIDX1,206:405], "patty GSM4835701","#DD864E",c(0,1))
dev.off()

pdf(file="TSS3kb_heatmap_active_H3K36me3_1kb.pdf",width=3*4,height=6)
par(mfrow=c(1,4),mar=c(4,4,2,2))
heatmapFIX(H3K36me3_CnT1raw[activeIDX1,206:405], "cutraw GSM4835701","red",0.6)
heatmapFIX(H3K36me3_ChIP1[activeIDX1,206:405], "ChIP1raw active","blue",0.3)
heatmapFIX(ATAC[activeIDX1,206:405], "ATAC active","black",12)
heatmapFIX(H3K36me3_CnT1score[activeIDX1,206:405], "patty GSM4835701","#DD864E",1)
dev.off()



