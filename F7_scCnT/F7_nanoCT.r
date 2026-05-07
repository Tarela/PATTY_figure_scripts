library(mclust)
library(Seurat)
#library(Signac)
library(gridExtra)
#library(infotheo)
library(cluster)


ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=8)


#### nanoCT_H3K27me3
nanoCT_H3K27me3_read <- read.table("nanoCT_H3K27me3_allcells_highVar_read.txt",row.names=1,header=T,sep="\t")
nanoCT_H3K27me3_patty <- read.table("nanoCT_H3K27me3_allcells_highVar_patty.txt",row.names=1,header=T,sep="\t")
nanoCT_H3K27me3_obj <- readRDS("nanoCT/GSE198467_H3K27me3_Seurat_object_clustered_renamed.Rds")
nanoCT_H3K27me3_idx <- rownames(nanoCT_H3K27me3_read)
nanoCT_H3K27me3_UMAP <- nanoCT_H3K27me3_obj@reductions$umap@cell.embeddings[nanoCT_H3K27me3_idx,]
nanoCT_H3K27me3_labels <- cbind(nanoCT_H3K27me3_read[nanoCT_H3K27me3_idx, c("CTl2","clusterKM_CT2num_hv10k")],nanoCT_H3K27me3_patty[nanoCT_H3K27me3_idx,"clusterKM_CT2num_hv10k"] )
colnames(nanoCT_H3K27me3_labels) <- c("cellname","read","patty")
color_list <- ggplotColours(n=length(unique(nanoCT_H3K27me3_labels[,1])))
nanoCT_H3K27me3_cellname_list <- sort(unique(nanoCT_H3K27me3_labels[,1]))
nanoCT_H3K27me3_cluster_list <- seq(1, length(nanoCT_H3K27me3_cellname_list))#sort(unique(nanoCT_H3K27me3_labels[,1]))
nanoCT_H3K27me3_colors <- matrix(rep("black",nrow(nanoCT_H3K27me3_labels)*3),ncol=3)
colnames(nanoCT_H3K27me3_colors) <-  c("cellname","read","patty")
rownames(nanoCT_H3K27me3_colors) <- rownames(nanoCT_H3K27me3_labels)



confmat <- function(inmat){
  feature1 <- unique(inmat[,1])
  feature2 <- unique(inmat[,2])
  outmat <- matrix(rep(0, length(feature1)*length(feature2)),nrow=length(feature1),ncol=length(feature2))
  rownames(outmat) <- feature1
  colnames(outmat) <- feature2
  for(i in 1:nrow(inmat)){
    outmat[inmat[i,1], inmat[i,2]] <- outmat[inmat[i,1], inmat[i,2]] + 1
  }
  return(outmat)
}

confmat_CT_patty <- confmat(cbind(as.character(nanoCT_H3K27me3_labels[,1]),nanoCT_H3K27me3_labels[,3]))
write.table(confmat_CT_patty , file="nanoCT_H3K27me3_confmat_CT_patty.txt",row.names=T,col.names=T,sep="\t",quote=F)
confmat_CT_read <- confmat(cbind(as.character(nanoCT_H3K27me3_labels[,1]),nanoCT_H3K27me3_labels[,2]))
write.table(confmat_CT_read , file="nanoCT_H3K27me3_confmat_CT_read.txt",row.names=T,col.names=T,sep="\t",quote=F)
confmat_read_patty <- confmat(cbind(as.character(nanoCT_H3K27me3_labels[,2]),nanoCT_H3K27me3_labels[,3]))
write.table(confmat_read_patty , file="nanoCT_H3K27me3_confmat_read_patty.txt",row.names=T,col.names=T,sep="\t",quote=F)


outdata <- cbind(nanoCT_H3K27me3_UMAP[nanoCT_H3K27me3_idx,],
                 nanoCT_H3K27me3_labels[nanoCT_H3K27me3_idx,],
                 nanoCT_H3K27me3_colors[nanoCT_H3K27me3_idx,])
rownames(outdata) <- nanoCT_H3K27me3_idx
write.table(outdata, file="F7_a-c_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)




png(file="nanoCT_H3K27me3_UMAP_cellnameColor.png")
par(mar=c(0,0,0,0))
plot(nanoCT_H3K27me3_UMAP[nanoCT_H3K27me3_idx,],cex=3,pch=".",col=nanoCT_H3K27me3_colors[nanoCT_H3K27me3_idx,"cellname"], axes=F,xlab="",ylab="")
box()
dev.off()
png(file="nanoCT_H3K27me3_UMAP_readColor.png")
par(mar=c(0,0,0,0))
plot(nanoCT_H3K27me3_UMAP[nanoCT_H3K27me3_idx,],cex=3,pch=".",col=nanoCT_H3K27me3_colors[nanoCT_H3K27me3_idx,"read"], axes=F,xlab="",ylab="")
box()
dev.off()
png(file="nanoCT_H3K27me3_UMAP_pattyColor.png")
par(mar=c(0,0,0,0))
plot(nanoCT_H3K27me3_UMAP[nanoCT_H3K27me3_idx,],cex=3,pch=".",col=nanoCT_H3K27me3_colors[nanoCT_H3K27me3_idx,"patty"], axes=F,xlab="",ylab="")
box()
dev.off()
pdf(file="nanoCT_H3K27me3_colorLegend.pdf",width=16,height=16)
plot(1,1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=nanoCT_H3K27me3_cellname_list,col=color_list,pch=16,bty="n")
legend("right",legend=paste0("c",nanoCT_H3K27me3_cluster_list),col=color_list[c(14,6,5,15,8,3,13,4,7,9,11,10,1,2,12)],pch=16,bty="n")
dev.off()

pdf(file="nanoCT_H3K27me3_UMAP_cellnameColor.pdf")
par(mar=c(0,0,0,0))
plot(nanoCT_H3K27me3_UMAP[nanoCT_H3K27me3_idx,],cex=3,pch=".",col=nanoCT_H3K27me3_colors[nanoCT_H3K27me3_idx,"cellname"], axes=F,xlab="",ylab="")
box()
dev.off()
pdf(file="nanoCT_H3K27me3_UMAP_readColor.pdf")
par(mar=c(0,0,0,0))
plot(nanoCT_H3K27me3_UMAP[nanoCT_H3K27me3_idx,],cex=3,pch=".",col=nanoCT_H3K27me3_colors[nanoCT_H3K27me3_idx,"read"], axes=F,xlab="",ylab="")
box()
dev.off()
pdf(file="nanoCT_H3K27me3_UMAP_pattyColor.pdf")
par(mar=c(0,0,0,0))
plot(nanoCT_H3K27me3_UMAP[nanoCT_H3K27me3_idx,],cex=3,pch=".",col=nanoCT_H3K27me3_colors[nanoCT_H3K27me3_idx,"patty"], axes=F,xlab="",ylab="")
box()
dev.off()

#### nanoCT_H3K27ac
nanoCT_H3K27ac_read <- read.table("nanoCT_H3K27ac_allcells_highVar_read.txt",row.names=1,header=T,sep="\t")
nanoCT_H3K27ac_patty <- read.table("nanoCT_H3K27ac_allcells_highVar_patty.txt",row.names=1,header=T,sep="\t")
nanoCT_H3K27ac_obj <- readRDS("nanoCT/GSE198467_H3K27ac_Seurat_object_clustered_renamed.Rds")
nanoCT_H3K27ac_idx <- rownames(nanoCT_H3K27ac_read)
nanoCT_H3K27ac_UMAP <- nanoCT_H3K27ac_obj@reductions$umap@cell.embeddings[nanoCT_H3K27ac_idx,]
nanoCT_H3K27ac_labels <- cbind(nanoCT_H3K27ac_read[nanoCT_H3K27ac_idx, c("CTl2","clusterKM_CT2num_hv10k")],nanoCT_H3K27ac_patty[nanoCT_H3K27ac_idx,"clusterKM_CT2num_hv10k"] )
colnames(nanoCT_H3K27ac_labels) <- c("cellname","read","patty")
color_list <- ggplotColours(n=length(unique(nanoCT_H3K27ac_labels[,1])))
nanoCT_H3K27ac_cellname_list <- sort(unique(nanoCT_H3K27ac_labels[,1]))
nanoCT_H3K27ac_cluster_list <- seq(1, length(nanoCT_H3K27ac_cellname_list))#sort(unique(nanoCT_H3K27ac_labels[,1]))
nanoCT_H3K27ac_colors <- matrix(rep("black",nrow(nanoCT_H3K27ac_labels)*3),ncol=3)
colnames(nanoCT_H3K27ac_colors) <-  c("cellname","read","patty")
rownames(nanoCT_H3K27ac_colors) <- rownames(nanoCT_H3K27ac_labels)
for(i in nanoCT_H3K27ac_cluster_list){
  this_CT <- nanoCT_H3K27ac_cellname_list[i]
  this_color <- color_list[i]
  nanoCT_H3K27ac_colors[which(nanoCT_H3K27ac_labels[,1] == this_CT),1] <- this_color
  nanoCT_H3K27ac_colors[which(nanoCT_H3K27ac_labels[,2] == i),2] <- this_color
  #nanoCT_H3K27ac_colors[which(nanoCT_H3K27ac_labels[,3] == i),3] <- this_color
}

confmat_CT_patty <- confmat(cbind(as.character(nanoCT_H3K27ac_labels[,1]),nanoCT_H3K27ac_labels[,3]))
write.table(confmat_CT_patty , file="nanoCT_H3K27ac_confmat_CT_patty.txt",row.names=T,col.names=T,sep="\t",quote=F)
confmat_CT_read <- confmat(cbind(as.character(nanoCT_H3K27ac_labels[,1]),nanoCT_H3K27ac_labels[,2]))
write.table(confmat_CT_read , file="nanoCT_H3K27ac_confmat_CT_read.txt",row.names=T,col.names=T,sep="\t",quote=F)
confmat_read_patty <- confmat(cbind(as.character(nanoCT_H3K27ac_labels[,2]),nanoCT_H3K27ac_labels[,3]))
write.table(confmat_read_patty , file="nanoCT_H3K27ac_confmat_read_patty.txt",row.names=T,col.names=T,sep="\t",quote=F)




outdata <- cbind(nanoCT_H3K27ac_UMAP[nanoCT_H3K27ac_idx,],
                 nanoCT_H3K27ac_labels[nanoCT_H3K27ac_idx,],
                 nanoCT_H3K27ac_colors[nanoCT_H3K27ac_idx,])
rownames(outdata) <- nanoCT_H3K27ac_idx
write.table(outdata, file="F7_d-f_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)

png(file="nanoCT_H3K27ac_UMAP_cellnameColor.png")
par(mar=c(0,0,0,0))
plot(nanoCT_H3K27ac_UMAP[nanoCT_H3K27ac_idx,],cex=3,pch=".",col=nanoCT_H3K27ac_colors[nanoCT_H3K27ac_idx,"cellname"], axes=F,xlab="",ylab="")
box()
dev.off()
png(file="nanoCT_H3K27ac_UMAP_readColor.png")
par(mar=c(0,0,0,0))
plot(nanoCT_H3K27ac_UMAP[nanoCT_H3K27ac_idx,],cex=3,pch=".",col=nanoCT_H3K27ac_colors[nanoCT_H3K27ac_idx,"read"], axes=F,xlab="",ylab="")
box()
dev.off()
png(file="nanoCT_H3K27ac_UMAP_pattyColor.png")
par(mar=c(0,0,0,0))
plot(nanoCT_H3K27ac_UMAP[nanoCT_H3K27ac_idx,],cex=3,pch=".",col=nanoCT_H3K27ac_colors[nanoCT_H3K27ac_idx,"patty"], axes=F,xlab="",ylab="")
box()
dev.off()
pdf(file="nanoCT_H3K27ac_colorLegend.pdf",width=16,height=16)
plot(1,1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=nanoCT_H3K27ac_cellname_list,col=color_list,pch=16,bty="n")
legend("right",legend=paste0("c",nanoCT_H3K27ac_cluster_list),col=color_list[c(8,1,3,11,10,2,16,4,7,15,18,9,14,13,6,12,17,5)],pch=16,bty="n")
dev.off()


pdf(file="nanoCT_H3K27ac_UMAP_cellnameColor.pdf")
par(mar=c(0,0,0,0))
plot(nanoCT_H3K27ac_UMAP[nanoCT_H3K27ac_idx,],cex=3,pch=".",col=nanoCT_H3K27ac_colors[nanoCT_H3K27ac_idx,"cellname"], axes=F,xlab="",ylab="")
box()
dev.off()
pdf(file="nanoCT_H3K27ac_UMAP_readColor.pdf")
par(mar=c(0,0,0,0))
plot(nanoCT_H3K27ac_UMAP[nanoCT_H3K27ac_idx,],cex=3,pch=".",col=nanoCT_H3K27ac_colors[nanoCT_H3K27ac_idx,"read"], axes=F,xlab="",ylab="")
box()
dev.off()
pdf(file="nanoCT_H3K27ac_UMAP_pattyColor.pdf")
par(mar=c(0,0,0,0))
plot(nanoCT_H3K27ac_UMAP[nanoCT_H3K27ac_idx,],cex=3,pch=".",col=nanoCT_H3K27ac_colors[nanoCT_H3K27ac_idx,"patty"], axes=F,xlab="",ylab="")
box()
dev.off()





pdf(file="F7DH_raw.pdf",width=2*2,height=3)
par(mfrow=c(1,2))
plotdata <- c(adjustedRandIndex(nanoCT_H3K27me3_read[,"CTl2"],nanoCT_H3K27me3_read[,"clusterKM_CT2num_hv10k"]),
              adjustedRandIndex(nanoCT_H3K27me3_patty[,"CTl2"],nanoCT_H3K27me3_patty[,"clusterKM_CT2num_hv10k"]))
barplot(plotdata,las=2,col=c("red","#DD864E"))
abline(h=0)
box(lwd=1)

plotdata <- c(adjustedRandIndex(nanoCT_H3K27ac_read[,"CTl2"],nanoCT_H3K27ac_read[,"clusterKM_CT2num_hv10k"]),
              adjustedRandIndex(nanoCT_H3K27ac_patty[,"CTl2"],nanoCT_H3K27ac_patty[,"clusterKM_CT2num_hv10k"]))
barplot(plotdata,las=2,col=c("red","#DD864E"))
box(lwd=1)
dev.off()

outdata <- c(adjustedRandIndex(nanoCT_H3K27me3_read[,"CTl2"],nanoCT_H3K27me3_read[,"clusterKM_CT2num_hv10k"]),
  adjustedRandIndex(nanoCT_H3K27me3_patty[,"CTl2"],nanoCT_H3K27me3_patty[,"clusterKM_CT2num_hv10k"]),
  adjustedRandIndex(nanoCT_H3K27ac_read[,"CTl2"],nanoCT_H3K27ac_read[,"clusterKM_CT2num_hv10k"]),
  adjustedRandIndex(nanoCT_H3K27ac_patty[,"CTl2"],nanoCT_H3K27ac_patty[,"clusterKM_CT2num_hv10k"]))
names(outdata) <- c("uncorrected_H3K27me3","PATTY_H3K27me3","uncorrected_H3K27ac","PATTY_H3K27ac")
write.table(outdata, file="F7_g-h_sourceData.txt",row.names=T,col.names=F,sep="\t",quote=F)



### WNNconnect kmeans
nanoCT_WNNconnect_patty <- read.table("nanoCT_kmeansWNNconnect_patty.txt",row.names=1)
nanoCT_WNNconnect_read <- read.table("nanoCT_kmeansWNNconnect_read.txt",row.names=1)
idx <- rownames(nanoCT_WNNconnect_read)
plotdata_H3K27me3 <- c(adjustedRandIndex(nanoCT_H3K27me3_read[,"CTl2"],nanoCT_H3K27me3_read[,"clusterKM_CT2num_hv10k"]),
              adjustedRandIndex(nanoCT_H3K27me3_patty[,"CTl2"],nanoCT_H3K27me3_patty[,"clusterKM_CT2num_hv10k"]))
plotdata_H3K27ac <- c(adjustedRandIndex(nanoCT_H3K27ac_read[,"CTl2"],nanoCT_H3K27ac_read[,"clusterKM_CT2num_hv10k"]),
              adjustedRandIndex(nanoCT_H3K27ac_patty[,"CTl2"],nanoCT_H3K27ac_patty[,"clusterKM_CT2num_hv10k"]))
plotdata_WNN <- c(adjustedRandIndex(nanoCT_H3K27me3_read[idx,"CTl2"],nanoCT_WNNconnect_read[idx,1]),
              adjustedRandIndex(nanoCT_H3K27me3_read[idx,"CTl2"],nanoCT_WNNconnect_patty[idx,1]))

pdf(file="F7DH_kmeansWNNconnect.pdf",width=2,height=3)
par(mfrow=c(1,1))
barplot(plotdata_yesWNN,las=2,col=c("red","#DD864E"),main="wnn cluster")
abline(h=0)
box(lwd=1)
dev.off()


pdf(file="SF10_WNN_ARIcmp.pdf",width=4,height=3)
par(mfrow=c(1,1))
barplot(c(plotdata_H3K27me3,plotdata_WNN),las=2,col=c("red","#DD864E","red","#DD864E"),main="")
abline(h=0)
box(lwd=1)
dev.off()

outdata <- c(plotdata_H3K27me3,plotdata_WNN)
names(outdata) <- c("uncorrected_withoutWNN","PATTY_withoutWNN","uncorrected_withWNN","PATTY_withWNN")
write.table(outdata, file="SF10_g_sourceData.txt",row.names=T,col.names=F,sep="\t",quote=F)



### WNN kmeans
nanoCT_WNN_patty <- read.table("nanoCT_kmeansWNN_patty.txt",row.names=1)
nanoCT_WNN_read <- read.table("nanoCT_kmeansWNN_read.txt",row.names=1)
idx <- rownames(nanoCT_WNN_read)

pdf(file="F7DH_kmeansWNN.pdf",width=2,height=3)
par(mfrow=c(1,1))
plotdata <- c(adjustedRandIndex(nanoCT_H3K27me3_read[idx,"CTl2"],nanoCT_WNN_read[idx,1]),
              adjustedRandIndex(nanoCT_H3K27me3_read[idx,"CTl2"],nanoCT_WNN_patty[idx,1]))
barplot(plotdata,las=2,col=c("red","#DD864E"),main="wnn cluster")
abline(h=0)
box(lwd=1)
dev.off()



pdf(file="F7DH_ALL.pdf",width=2*4,height=3)
par(mfrow=c(1,4))

plotdata <- c(adjustedRandIndex(nanoCT_H3K27me3_read[,"CTl2"],nanoCT_H3K27me3_read[,"clusterKM_CT2num_hv10k"]),
              adjustedRandIndex(nanoCT_H3K27me3_patty[,"CTl2"],nanoCT_H3K27me3_patty[,"clusterKM_CT2num_hv10k"]))
barplot(plotdata,las=2,col=c("red","#DD864E"),main="K27me3 cluster",ylim=c(0,0.35))
abline(h=0)
box(lwd=1)

plotdata <- c(adjustedRandIndex(nanoCT_H3K27ac_read[,"CTl2"],nanoCT_H3K27ac_read[,"clusterKM_CT2num_hv10k"]),
              adjustedRandIndex(nanoCT_H3K27ac_patty[,"CTl2"],nanoCT_H3K27ac_patty[,"clusterKM_CT2num_hv10k"]))
barplot(plotdata,las=2,col=c("red","#DD864E"),main="K27ac cluster",ylim=c(0,0.35))
box(lwd=1)

plotdata <- c(adjustedRandIndex(nanoCT_H3K27me3_read[idx,"CTl2"],nanoCT_WNN_read[idx,1]),
              adjustedRandIndex(nanoCT_H3K27me3_read[idx,"CTl2"],nanoCT_WNN_patty[idx,1]))
barplot(plotdata,las=2,col=c("red","#DD864E"),main="wnn me3+ac cluster",ylim=c(0,0.35))
abline(h=0)
box(lwd=1)

plotdata <- c(adjustedRandIndex(nanoCT_H3K27me3_read[idx,"CTl2"],nanoCT_WNNconnect_read[idx,1]),
              adjustedRandIndex(nanoCT_H3K27me3_read[idx,"CTl2"],nanoCT_WNNconnect_patty[idx,1]))
barplot(plotdata,las=2,col=c("red","#DD864E"),main="wnn me3 connect ac cluster",ylim=c(0,0.35))
abline(h=0)
box(lwd=1)

dev.off()

















