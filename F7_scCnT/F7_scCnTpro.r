library(mclust)
library(Seurat)
#library(Signac)
library(gridExtra)
library(cluster)


ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=8)


#### scCnTpro_H3K27me3
scCnTpro_H3K27me3_read <- read.table("scCnTpro_H3K27me3_allcells_highVar_read.txt",row.names=1,header=T,sep="\t")
scCnTpro_H3K27me3_patty <- read.table("scCnTpro_H3K27me3_allcells_highVar_patty.txt",row.names=1,header=T,sep="\t")
scCnTpro_H3K27me3_obj <- readRDS("scCnTpro/H3K27me3.rds")
scCnTpro_H3K27me3_idx <- rownames(scCnTpro_H3K27me3_read)
scCnTpro_H3K27me3_UMAPwnn <- scCnTpro_H3K27me3_obj@reductions$wnn.umap@cell.embeddings[scCnTpro_H3K27me3_idx,]
scCnTpro_H3K27me3_UMAPadt <- scCnTpro_H3K27me3_obj@reductions$adt.umap@cell.embeddings[scCnTpro_H3K27me3_idx,]
scCnTpro_H3K27me3_UMAPcnt <- scCnTpro_H3K27me3_obj@reductions$cut.tag.umap@cell.embeddings[scCnTpro_H3K27me3_idx,]
scCnTpro_H3K27me3_UMAPazi <- scCnTpro_H3K27me3_obj@reductions$azimuth.umap@cell.embeddings[scCnTpro_H3K27me3_idx,]
scCnTpro_H3K27me3_labels <- cbind(scCnTpro_H3K27me3_read[scCnTpro_H3K27me3_idx, c("CTl1","clusterKM_CT1num_hv10k")],scCnTpro_H3K27me3_patty[scCnTpro_H3K27me3_idx,"clusterKM_CT1num_hv10k"] )
colnames(scCnTpro_H3K27me3_labels) <- c("cellname","read","patty")
color_list <- ggplotColours(n=length(unique(scCnTpro_H3K27me3_labels[,1])))
scCnTpro_H3K27me3_cellname_list <- sort(unique(scCnTpro_H3K27me3_labels[,1]))
scCnTpro_H3K27me3_cluster_list <- seq(1, length(scCnTpro_H3K27me3_cellname_list))#sort(unique(scCnTpro_H3K27me3_labels[,1]))
scCnTpro_H3K27me3_colors <- matrix(rep("black",nrow(scCnTpro_H3K27me3_labels)*3),ncol=3)
colnames(scCnTpro_H3K27me3_colors) <-  c("cellname","read","patty")
rownames(scCnTpro_H3K27me3_colors) <- rownames(scCnTpro_H3K27me3_labels)
for(i in scCnTpro_H3K27me3_cluster_list){
  this_CT <- scCnTpro_H3K27me3_cellname_list[i]
  this_color <- color_list[i]
  scCnTpro_H3K27me3_colors[which(scCnTpro_H3K27me3_labels[,1] == this_CT),1] <- this_color
}


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

confmat_CT_patty <- confmat(cbind(as.character(scCnTpro_H3K27me3_labels[,1]),scCnTpro_H3K27me3_labels[,3]))
write.table(confmat_CT_patty , file="scCnTpro_H3K27me3_confmat_CT_patty.txt",row.names=T,col.names=T,sep="\t",quote=F)
confmat_CT_read <- confmat(cbind(as.character(scCnTpro_H3K27me3_labels[,1]),scCnTpro_H3K27me3_labels[,2]))
write.table(confmat_CT_read , file="scCnTpro_H3K27me3_confmat_CT_read.txt",row.names=T,col.names=T,sep="\t",quote=F)
confmat_read_patty <- confmat(cbind(as.character(scCnTpro_H3K27me3_labels[,2]),scCnTpro_H3K27me3_labels[,3]))
write.table(confmat_read_patty , file="scCnTpro_H3K27me3_confmat_read_patty.txt",row.names=T,col.names=T,sep="\t",quote=F)



outdata <- cbind(scCnTpro_H3K27me3_UMAP[scCnTpro_H3K27me3_idx,],
      scCnTpro_H3K27me3_labels[scCnTpro_H3K27me3_idx,],
      scCnTpro_H3K27me3_colors[scCnTpro_H3K27me3_idx,])
write.table(outdata, file="SF9_a-c_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)

scCnTpro_H3K27me3_UMAP <- scCnTpro_H3K27me3_UMAPazi
png(file="scCnTpro_H3K27me3_UMAP_cellnameColor.png")
par(mar=c(0,0,0,0))
plot(scCnTpro_H3K27me3_UMAP[scCnTpro_H3K27me3_idx,],cex=3,pch=".",col=scCnTpro_H3K27me3_colors[scCnTpro_H3K27me3_idx,"cellname"], axes=F,xlab="",ylab="")
box()
dev.off()
png(file="scCnTpro_H3K27me3_UMAP_readColor.png")
par(mar=c(0,0,0,0))
plot(scCnTpro_H3K27me3_UMAP[scCnTpro_H3K27me3_idx,],cex=3,pch=".",col=scCnTpro_H3K27me3_colors[scCnTpro_H3K27me3_idx,"read"], axes=F,xlab="",ylab="")
box()
dev.off()
png(file="scCnTpro_H3K27me3_UMAP_pattyColor.png")
par(mar=c(0,0,0,0))
plot(scCnTpro_H3K27me3_UMAP[scCnTpro_H3K27me3_idx,],cex=3,pch=".",col=scCnTpro_H3K27me3_colors[scCnTpro_H3K27me3_idx,"patty"], axes=F,xlab="",ylab="")
box()
dev.off()
pdf(file="scCnTpro_H3K27me3_colorLegend.pdf",width=16,height=16)
plot(1,1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=scCnTpro_H3K27me3_cellname_list,col=color_list,pch=16,bty="n")
legend("right",legend=paste0("c",scCnTpro_H3K27me3_cluster_list),col=color_list[c(6,3,7,8,4,5,2,1)],pch=16,bty="n")
dev.off()

pdf(file="scCnTpro_H3K27me3_UMAP_cellnameColor.pdf")
par(mar=c(0,0,0,0))
plot(scCnTpro_H3K27me3_UMAP[scCnTpro_H3K27me3_idx,],cex=3,pch=".",col=scCnTpro_H3K27me3_colors[scCnTpro_H3K27me3_idx,"cellname"], axes=F,xlab="",ylab="")
box()
dev.off()
pdf(file="scCnTpro_H3K27me3_UMAP_readColor.pdf")
par(mar=c(0,0,0,0))
plot(scCnTpro_H3K27me3_UMAP[scCnTpro_H3K27me3_idx,],cex=3,pch=".",col=scCnTpro_H3K27me3_colors[scCnTpro_H3K27me3_idx,"read"], axes=F,xlab="",ylab="")
box()
dev.off()
pdf(file="scCnTpro_H3K27me3_UMAP_pattyColor.pdf")
par(mar=c(0,0,0,0))
plot(scCnTpro_H3K27me3_UMAP[scCnTpro_H3K27me3_idx,],cex=3,pch=".",col=scCnTpro_H3K27me3_colors[scCnTpro_H3K27me3_idx,"patty"], axes=F,xlab="",ylab="")
box()
dev.off()




#### scCnTpro_H3K27ac
scCnTpro_H3K27ac_read <- read.table("scCnTpro_H3K27ac_allcells_highVar_read.txt",row.names=1,header=T,sep="\t")
scCnTpro_H3K27ac_patty <- read.table("scCnTpro_H3K27ac_allcells_highVar_patty.txt",row.names=1,header=T,sep="\t")
scCnTpro_H3K27ac_obj <- readRDS("scCnTpro/H3K27ac.rds")
scCnTpro_H3K27ac_idx <- rownames(scCnTpro_H3K27ac_read)
scCnTpro_H3K27ac_UMAPwnn <- scCnTpro_H3K27ac_obj@reductions$wnn.umap@cell.embeddings[scCnTpro_H3K27ac_idx,]
scCnTpro_H3K27ac_UMAPadt <- scCnTpro_H3K27ac_obj@reductions$adt.umap@cell.embeddings[scCnTpro_H3K27ac_idx,]
scCnTpro_H3K27ac_UMAPcnt <- scCnTpro_H3K27ac_obj@reductions$cut.tag.umap@cell.embeddings[scCnTpro_H3K27ac_idx,]
scCnTpro_H3K27ac_UMAPazi <- scCnTpro_H3K27ac_obj@reductions$azimuth.umap@cell.embeddings[scCnTpro_H3K27ac_idx,]
scCnTpro_H3K27ac_labels <- cbind(scCnTpro_H3K27ac_read[scCnTpro_H3K27ac_idx, c("CTl1","clusterKM_CT1num_hv10k")],scCnTpro_H3K27ac_patty[scCnTpro_H3K27ac_idx,"clusterKM_CT1num_hv10k"] )
colnames(scCnTpro_H3K27ac_labels) <- c("cellname","read","patty")
color_list <- ggplotColours(n=length(unique(scCnTpro_H3K27ac_labels[,1])))
scCnTpro_H3K27ac_cellname_list <- sort(unique(scCnTpro_H3K27ac_labels[,1]))
scCnTpro_H3K27ac_cluster_list <- seq(1, length(scCnTpro_H3K27ac_cellname_list))#sort(unique(scCnTpro_H3K27ac_labels[,1]))
scCnTpro_H3K27ac_colors <- matrix(rep("black",nrow(scCnTpro_H3K27ac_labels)*3),ncol=3)
colnames(scCnTpro_H3K27ac_colors) <-  c("cellname","read","patty")
rownames(scCnTpro_H3K27ac_colors) <- rownames(scCnTpro_H3K27ac_labels)


pdf(file="scCnTpro_ARIbar.pdf",width=2*2,height=4)
par(mfrow=c(1,2))
plotdata <- c(adjustedRandIndex(scCnTpro_H3K27me3_read[,c("CTl1")], scCnTpro_H3K27me3_read[,c("clusterKM_CT1num_hv10k")]),
              adjustedRandIndex(scCnTpro_H3K27me3_patty[,c("CTl1")], scCnTpro_H3K27me3_patty[,c("clusterKM_CT1num_hv10k")]))
barplot(plotdata,las=2,col=c("red","#DD864E"),ylim=c(-0.02,0.14),main="H3K27me3",axes=F)
axis(side=2,at=seq(0,0.12,0.03),las=2)
abline(h=0)
box(lwd=1)

plotdata <- c(adjustedRandIndex(scCnTpro_H3K27ac_read[,c("CTl1")], scCnTpro_H3K27ac_read[,c("clusterKM_CT1num_hv10k")]),
              adjustedRandIndex(scCnTpro_H3K27ac_patty[,c("CTl1")], scCnTpro_H3K27ac_patty[,c("clusterKM_CT1num_hv10k")]))
barplot(plotdata,las=2,col=c("red","#DD864E"),ylim=c(-0.02,0.14),main="H3K27ac",axes=F)
axis(side=2,at=seq(0,0.12,0.03),las=2)
abline(h=0)
box(lwd=1)
dev.off()


outdata <- c(adjustedRandIndex(scCnTpro_H3K27me3_read[,c("CTl1")], scCnTpro_H3K27me3_read[,c("clusterKM_CT1num_hv10k")]),
             adjustedRandIndex(scCnTpro_H3K27me3_patty[,c("CTl1")], scCnTpro_H3K27me3_patty[,c("clusterKM_CT1num_hv10k")]),
             adjustedRandIndex(scCnTpro_H3K27ac_read[,c("CTl1")], scCnTpro_H3K27ac_read[,c("clusterKM_CT1num_hv10k")]),
             adjustedRandIndex(scCnTpro_H3K27ac_patty[,c("CTl1")], scCnTpro_H3K27ac_patty[,c("clusterKM_CT1num_hv10k")]))
names(outdata) <- c("uncorrected_H3K27me3","PATTY_H3K27me3","uncorrected_H3K27ac","PATTY_H3K27ac")
write.table(outdata, file="F7_i-j_sourceData.txt",row.names=T,col.names=F,sep="\t",quote=F)


confmat_CT_patty <- confmat(cbind(as.character(scCnTpro_H3K27ac_labels[,1]),scCnTpro_H3K27ac_labels[,3]))
write.table(confmat_CT_patty , file="scCnTpro_H3K27ac_confmat_CT_patty.txt",row.names=T,col.names=T,sep="\t",quote=F)
confmat_CT_read <- confmat(cbind(as.character(scCnTpro_H3K27ac_labels[,1]),scCnTpro_H3K27ac_labels[,2]))
write.table(confmat_CT_read , file="scCnTpro_H3K27ac_confmat_CT_read.txt",row.names=T,col.names=T,sep="\t",quote=F)
confmat_read_patty <- confmat(cbind(as.character(scCnTpro_H3K27ac_labels[,2]),scCnTpro_H3K27ac_labels[,3]))
write.table(confmat_read_patty , file="scCnTpro_H3K27ac_confmat_read_patty.txt",row.names=T,col.names=T,sep="\t",quote=F)


outdata <- cbind(scCnTpro_H3K27ac_UMAP[scCnTpro_H3K27ac_idx,],
                 scCnTpro_H3K27ac_labels[scCnTpro_H3K27ac_idx,],
                 scCnTpro_H3K27ac_colors[scCnTpro_H3K27ac_idx,])
write.table(outdata, file="SF9_d-f_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)

scCnTpro_H3K27ac_UMAP <- scCnTpro_H3K27ac_UMAPazi
png(file="scCnTpro_H3K27ac_UMAP_cellnameColor.png")
par(mar=c(0,0,0,0))
plot(scCnTpro_H3K27ac_UMAP[scCnTpro_H3K27ac_idx,],cex=3,pch=".",col=scCnTpro_H3K27ac_colors[scCnTpro_H3K27ac_idx,"cellname"], axes=F,xlab="",ylab="")
box()
dev.off()
png(file="scCnTpro_H3K27ac_UMAP_readColor.png")
par(mar=c(0,0,0,0))
plot(scCnTpro_H3K27ac_UMAP[scCnTpro_H3K27ac_idx,],cex=3,pch=".",col=scCnTpro_H3K27ac_colors[scCnTpro_H3K27ac_idx,"read"], axes=F,xlab="",ylab="")
box()
dev.off()
png(file="scCnTpro_H3K27ac_UMAP_pattyColor.png")
par(mar=c(0,0,0,0))
plot(scCnTpro_H3K27ac_UMAP[scCnTpro_H3K27ac_idx,],cex=3,pch=".",col=scCnTpro_H3K27ac_colors[scCnTpro_H3K27ac_idx,"patty"], axes=F,xlab="",ylab="")
box()
dev.off()
pdf(file="scCnTpro_H3K27ac_colorLegend.pdf",width=16,height=16)
plot(1,1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=scCnTpro_H3K27ac_cellname_list,col=color_list,pch=16,bty="n")
legend("right",legend=paste0("c",scCnTpro_H3K27ac_cluster_list),col=color_list[c(7,3,5,8,1,4,2,6)],pch=16,bty="n")
dev.off()

pdf(file="scCnTpro_H3K27ac_UMAP_cellnameColor.pdf")
par(mar=c(0,0,0,0))
plot(scCnTpro_H3K27ac_UMAP[scCnTpro_H3K27ac_idx,],cex=3,pch=".",col=scCnTpro_H3K27ac_colors[scCnTpro_H3K27ac_idx,"cellname"], axes=F,xlab="",ylab="")
box()
dev.off()
pdf(file="scCnTpro_H3K27ac_UMAP_readColor.pdf")
par(mar=c(0,0,0,0))
plot(scCnTpro_H3K27ac_UMAP[scCnTpro_H3K27ac_idx,],cex=3,pch=".",col=scCnTpro_H3K27ac_colors[scCnTpro_H3K27ac_idx,"read"], axes=F,xlab="",ylab="")
box()
dev.off()
pdf(file="scCnTpro_H3K27ac_UMAP_pattyColor.pdf")
par(mar=c(0,0,0,0))
plot(scCnTpro_H3K27ac_UMAP[scCnTpro_H3K27ac_idx,],cex=3,pch=".",col=scCnTpro_H3K27ac_colors[scCnTpro_H3K27ac_idx,"patty"], axes=F,xlab="",ylab="")
box()
dev.off()


