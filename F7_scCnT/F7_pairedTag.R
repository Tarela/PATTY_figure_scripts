ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

### Kmeans
PC_ac_patty  <- read.table("repCB_RNAlabel_PC_ac_patty.txt" ,row.names=1,header=T,check.names = FALSE)
PC_me3_patty <- read.table("repCB_RNAlabel_PC_me3_patty.txt",row.names=1,header=T,check.names = FALSE)
PC_ac_reads  <- read.table("repCB_RNAlabel_PC_ac_reads.txt" ,row.names=1,header=T,check.names = FALSE)
PC_me3_reads <- read.table("repCB_RNAlabel_PC_me3_reads.txt",row.names=1,header=T,check.names = FALSE)

rownames(PC_ac_patty ) <- sub("^([^.]+)\\.(.+)\\.([^.]+)$", "\\1:\\2-\\3", rownames(PC_ac_patty ))
rownames(PC_me3_patty) <- sub("^([^.]+)\\.(.+)\\.([^.]+)$", "\\1:\\2-\\3", rownames(PC_me3_patty))
rownames(PC_ac_reads ) <- sub("^([^.]+)\\.(.+)\\.([^.]+)$", "\\1:\\2-\\3", rownames(PC_ac_reads ))
rownames(PC_me3_reads) <- sub("^([^.]+)\\.(.+)\\.([^.]+)$", "\\1:\\2-\\3", rownames(PC_me3_reads))

H3K27me3_rawmeta <- read.table(paste0("../H3K27me3_meta.txt"),row.names=1,header=T,comment.char=";",sep="\t")[rownames(PC_me3_patty),]
H3K27ac_rawmeta <- read.table(paste0("../H3K27ac_meta.txt"),row.names=1,header=T,comment.char=";",sep="\t")[rownames(PC_ac_patty),]

pdf(file="UMAP_CT_RNAlabel.pdf",width=15,height=5)
par(mfrow=c(1,3),mar=c(4,4,2,2))
df <- H3K27me3_rawmeta
df$rna_Anno0.8 <- as.factor(df$rna_Anno0.8)
# assign a color per cell type
ct <- levels(df$rna_Anno0.8)
cols <- setNames(ggplotColours(n=length(ct)), ct)  # simple distinct colors
plot(df$int_umap_1, df$int_umap_2,
     pch = 16, cex = 0.5,
     col = cols[df$rna_Anno0.8],
     xlab = "UMAP-1", ylab = "UMAP-2",
     main = "H3K27me3 RNACT")

df <- H3K27ac_rawmeta
df$rna_Anno0.8 <- as.factor(df$rna_Anno0.8)
# assign a color per cell type
ct <- levels(df$rna_Anno0.8)
cols <- setNames(ggplotColours(n=length(ct)), ct)  # simple distinct colors
plot(df$int_umap_1, df$int_umap_2,
     pch = 16, cex = 0.5,
     col = cols[df$rna_Anno0.8],
     xlab = "UMAP-1", ylab = "UMAP-2",
     main = "H3K27ac RNACT")

plot(1,1,type="n",xlab="",ylab="",main="",axes=F)
legend("center", legend = ct, col = cols[ct], pch = 16, cex = 1, bty = "n")
dev.off()


###### Kmeans
CTnum <- 20
set.seed(1)
ac_kmeans_patty <- kmeans(PC_ac_patty , centers = CTnum,iter.max=100)$cluster
me3_kmeans_patty <- kmeans(PC_me3_patty, centers = CTnum,iter.max=100)$cluster
ac_kmeans_reads <- kmeans(PC_ac_reads , centers = CTnum,iter.max=100)$cluster
me3_kmeans_reads <- kmeans(PC_me3_reads, centers = CTnum,iter.max=100)$cluster

H3K27me3_plotdata <- c(adjustedRandIndex(H3K27me3_rawmeta[names(me3_kmeans_reads),"rna_Anno0.8"],me3_kmeans_reads),
                       adjustedRandIndex(H3K27me3_rawmeta[names(me3_kmeans_patty),"rna_Anno0.8"],me3_kmeans_patty))
H3K27ac_plotdata <- c(adjustedRandIndex(H3K27ac_rawmeta[names(ac_kmeans_reads),"rna_Anno0.8"],ac_kmeans_reads),
                       adjustedRandIndex(H3K27ac_rawmeta[names(ac_kmeans_patty),"rna_Anno0.8"],ac_kmeans_patty))

pdf(file="PairedTag_ARI_RNAlabel_fixY.pdf",width=3*2,height=4)
par(mfrow=c(1,2),mar=c(4,4,2,2))
barplot(H3K27me3_plotdata,las=2,col=c("red","#DD864E"),main="H3K27me3",ylim=c(0,0.1))
abline(h=0)
barplot(H3K27ac_plotdata,las=2,col=c("red","#DD864E"),main="H3K27ac",ylim=c(0,0.1))
abline(h=0)
dev.off()


outdata <- c(H3K27me3_plotdata,H3K27ac_plotdata)
names(outdata) <- c("uncorrected_H3K27me3","PATTY_H3K27me3","uncorrected_H3K27ac","PATTY_H3K27ac")
write.table(outdata, file="F7_k-l_sourceData.txt",row.names=T,col.names=F,sep="\t",quote=F)



confmat <- function( inmat){
  feature1 <- sort(unique(inmat[,1]))
  feature2 <- sort(unique(inmat[,2]))
  outmat <- matrix(rep(0, length(feature1)*length(feature2)),nrow=length(feature1),ncol=length(feature2))
  rownames(outmat) <- feature1
  colnames(outmat) <- feature2
  for(i in 1:nrow(inmat)){
    outmat[inmat[i,1], inmat[i,2]] <- outmat[inmat[i,1], inmat[i,2]] + 1
  }
  return(outmat)
}
confmatCT_patty_H3K27me3 <- confmat(cbind(H3K27me3_rawmeta[names(me3_kmeans_patty),"rna_Anno0.8"],paste0("C",me3_kmeans_patty)))
confmatCT_patty_H3K27ac <- confmat(cbind(H3K27ac_rawmeta[names(ac_kmeans_patty),"rna_Anno0.8"],paste0("C",ac_kmeans_patty)))

write.table(confmatCT_patty_H3K27ac,file="confmat_RNACT_patty_H3K27ac.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(confmatCT_patty_H3K27me3,file="confmat_RNACT_patty_H3K27me3.txt",row.names=T,col.names=T,sep="\t",quote=F)




H3K27me3_idx <- names(me3_kmeans_patty)
H3K27me3_labels <- cbind(H3K27me3_rawmeta[H3K27me3_idx, "rna_Anno0.8"], me3_kmeans_reads[H3K27me3_idx],me3_kmeans_patty[H3K27me3_idx] )
colnames(H3K27me3_labels) <- c("cellname","read","patty")
rownames(H3K27me3_labels) <- H3K27me3_idx#, c("CTl2")
color_list <- cols#ggplotColours(n=length(unique(H3K27me3_labels[,1])))
H3K27me3_cellname_list <- sort(unique(H3K27me3_labels[,1]))
H3K27me3_cluster_list <- seq(1, length(H3K27me3_cellname_list))#sort(unique(H3K27me3_labels[,1]))
H3K27me3_colors <- matrix(rep("black",nrow(H3K27me3_labels)*3),ncol=3)
colnames(H3K27me3_colors) <-  c("cellname","read","patty")
rownames(H3K27me3_colors) <- H3K27me3_idx#rownames(H3K27me3_labels)
H3K27me3_UMAP <- H3K27me3_rawmeta[H3K27me3_idx,c("int_umap_1","int_umap_2")]


outdata <- cbind(H3K27me3_UMAP[H3K27me3_idx,],
                 H3K27me3_labels[H3K27me3_idx,],
                 H3K27me3_colors[H3K27me3_idx,])
write.table(outdata, file="SF9_g-i_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)


png(file="pairedTag_RNAlabel_UMAP_cellnameColor.png")
par(mar=c(0,0,0,0))
plot(H3K27me3_UMAP[H3K27me3_idx,],cex=3,pch=".",col=H3K27me3_colors[H3K27me3_idx,"cellname"], axes=F,xlab="",ylab="")
box()
dev.off()
png(file="pairedTag_RNAlabel_UMAP_readColor.png")
par(mar=c(0,0,0,0))
plot(H3K27me3_UMAP[H3K27me3_idx,],cex=3,pch=".",col=H3K27me3_colors[H3K27me3_idx,"read"], axes=F,xlab="",ylab="")
box()
dev.off()
png(file="pairedTag_RNAlabel_UMAP_pattyColor.png")
par(mar=c(0,0,0,0))
plot(H3K27me3_UMAP[H3K27me3_idx,],cex=3,pch=".",col=H3K27me3_colors[H3K27me3_idx,"patty"], axes=F,xlab="",ylab="")
box()
dev.off()








###ac


H3K27ac_idx <- names(ac_kmeans_patty)
H3K27ac_labels <- cbind(H3K27ac_rawmeta[H3K27ac_idx, "rna_Anno0.8"], ac_kmeans_reads[H3K27ac_idx],ac_kmeans_patty[H3K27ac_idx] )
colnames(H3K27ac_labels) <- c("cellname","read","patty")
rownames(H3K27ac_labels) <- H3K27ac_idx#, c("CTl2")
color_list <- cols#ggplotColours(n=length(unique(H3K27ac_labels[,1])))
H3K27ac_cellname_list <- sort(unique(H3K27ac_labels[,1]))
H3K27ac_cluster_list <- seq(1, length(H3K27ac_cellname_list))#sort(unique(H3K27ac_labels[,1]))
H3K27ac_colors <- matrix(rep("black",nrow(H3K27ac_labels)*3),ncol=3)
colnames(H3K27ac_colors) <-  c("cellname","read","patty")
rownames(H3K27ac_colors) <- H3K27ac_idx#rownames(H3K27ac_labels)
H3K27ac_UMAP <- H3K27ac_rawmeta[H3K27ac_idx,c("int_umap_1","int_umap_2")]



outdata <- cbind(H3K27ac_UMAP[H3K27ac_idx,],
                 H3K27ac_labels[H3K27ac_idx,],
                 H3K27ac_colors[H3K27ac_idx,])
write.table(outdata, file="SF9_j-l_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)


png(file="pairedTag_RNAlabel_UMAP_cellnameColor_H3K27ac.png")
par(mar=c(0,0,0,0))
plot(H3K27ac_UMAP[H3K27ac_idx,],cex=3,pch=".",col=H3K27ac_colors[H3K27ac_idx,"cellname"], axes=F,xlab="",ylab="")
box()
dev.off()
png(file="pairedTag_RNAlabel_UMAP_readColor_H3K27ac.png")
par(mar=c(0,0,0,0))
plot(H3K27ac_UMAP[H3K27ac_idx,],cex=3,pch=".",col=H3K27ac_colors[H3K27ac_idx,"read"], axes=F,xlab="",ylab="")
box()
dev.off()
png(file="pairedTag_RNAlabel_UMAP_pattyColor_H3K27ac.png")
par(mar=c(0,0,0,0))
plot(H3K27ac_UMAP[H3K27ac_idx,],cex=3,pch=".",col=H3K27ac_colors[H3K27ac_idx,"patty"], axes=F,xlab="",ylab="")
box()
dev.off()



pdf(file="pairedTag_RNAlabel_colorLegend.pdf",width=32,height=16)
par(mfrow=c(1,2))
plot(1,1,type="n",axes=F,xlab="",ylab="",main="H3K27me3")
legend("left",legend=H3K27me3_cellname_list,col=cols,pch=16,bty="n",y.intersp = 0.7)
legend("right",legend=paste0("c",H3K27me3_cluster_list),col=color_list[c(4,19,13,18,11,3,20,5,10,12,7,14,16,6,2,1,8,9,15,17)],pch=16,bty="n",y.intersp = 0.7)

plot(1,1,type="n",axes=F,xlab="",ylab="",main="H3K27ac")
legend("left",legend=H3K27ac_cellname_list,col=cols,pch=16,bty="n",y.intersp = 0.7)
legend("right",legend=paste0("c",H3K27ac_cluster_list),col=color_list[c(9,8,10,11,15,1,16,20,6,12,5,7,18,2,13,3,14,19,17,4)],pch=16,bty="n",y.intersp = 0.7)
dev.off()



###### adjusted rand index
data <- read.table("PairedTag_kmeansCMPdata.txt",row.names=1)
colnames(data) <- c("reads","patty")
pdf(file="F7_PairedTag_ARI_fixY.pdf",width=3*2,height=2*4)
par(mfrow=c(2,3),mar=c(4,4,2,2))
for(i in c(2,4,6,1,3,5)){
  plotdata <- as.numeric(data[i,])
  usename <- rownames(data)[i]
  barplot(plotdata,las=2,col=c("red","#DD864E"),main=usename,ylim=c(0,0.25),cex.main=0.7)
}
dev.off()
pdf(file="F7_PairedTag_ARI_flexibleY.pdf",width=3*2,height=2*4)
par(mfrow=c(2,3),mar=c(4,4,2,2))
for(i in c(2,4,6,1,3,5)){
  plotdata <- as.numeric(data[i,])
  usename <- rownames(data)[i]
  barplot(plotdata,las=2,col=c("red","#DD864E"),main=usename,cex.main=0.7)
}
dev.off()





