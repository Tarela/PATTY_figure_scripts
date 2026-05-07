library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpattern)
library(tidyverse)


### H3K27me3
dNN <- read.table("crossValidation/K562_H3K27me3_modelAcc_noX.txt")
H3K27me3_dNN_best <- c()
for(METHOD in c("CNN","GRU","MLP","RNN")){
  for(Feature in c("CnTATACIgG","CnTATACseq","CnTATAC",'CnTIgG',"CnTonly","onehot")){
    thisdata <- dNN[which(dNN[,1] == METHOD & dNN[,2] == Feature),  1:11]
    outdata <- thisdata[order(apply(thisdata[,7:11],1,mean),decreasing=T),][1,]
    H3K27me3_dNN_best <- rbind(H3K27me3_dNN_best, outdata)
  }
}

other_H3K27me3 <- read.table("crossValidation/K562_H3K27me3_otherModels_crossValidationAcc_noX.txt")
cbdata <- rbind(as.matrix(H3K27me3_dNN_best[,c(1,2,7:11)]),as.matrix(other_H3K27me3[c(7:12,13:18,1:6),]))
rownames(cbdata) <- paste0(cbdata[,1],"_",cbdata[,2])
plotdata_H3K27me3 <- cbdata[,3:7]
plotdata_H3K27me3 <- matrix(as.numeric(plotdata_H3K27me3), nrow = nrow(plotdata_H3K27me3), ncol = ncol(plotdata_H3K27me3))
rownames(plotdata_H3K27me3) <- rownames(cbdata)

plotdata_H3K27me3 <- rbind(
  plotdata_H3K27me3[grep("CNN_",rownames(plotdata_H3K27me3)),],
  plotdata_H3K27me3[grep("RNN_",rownames(plotdata_H3K27me3)),],
  plotdata_H3K27me3[grep("GRU_",rownames(plotdata_H3K27me3)),],
  plotdata_H3K27me3[grep("MLP_",rownames(plotdata_H3K27me3)),],
  plotdata_H3K27me3[grep("LR_",rownames(plotdata_H3K27me3)),] ,
  plotdata_H3K27me3[grep("RF_",rownames(plotdata_H3K27me3)),] ,
  plotdata_H3K27me3[grep("GBM_",rownames(plotdata_H3K27me3)),])

plotdata_H3K27me3_cb <- rbind(
  as.numeric(plotdata_H3K27me3[grep("CNN_",rownames(plotdata_H3K27me3)),]),
  as.numeric(plotdata_H3K27me3[grep("RNN_",rownames(plotdata_H3K27me3)),]),
  as.numeric(plotdata_H3K27me3[grep("GRU_",rownames(plotdata_H3K27me3)),]),
  as.numeric(plotdata_H3K27me3[grep("MLP_",rownames(plotdata_H3K27me3)),]),
  as.numeric(plotdata_H3K27me3[grep("LR_",rownames(plotdata_H3K27me3)),] ),
  as.numeric(plotdata_H3K27me3[grep("RF_",rownames(plotdata_H3K27me3)),] ),
  as.numeric(plotdata_H3K27me3[grep("GBM_",rownames(plotdata_H3K27me3)),]))
rownames(plotdata_H3K27me3_cb) <- c("CNN","RNN","GRU","MLP","LR","RF","GBM")

colors <- c("#DD864E","#969696", "#E6E6E6","#C6DBEF", "#6BAED6", "#2171B5", "#08306B")

write.table(plotdata_H3K27me3, file="SF2_d_sourceData.txt",row.names=T,col.names=F,sep="\t",quote=F)


#### SF2D
mat <- plotdata_H3K27me3
# Extract model and input from row names
parts <- strsplit(rownames(mat), "_")
model <- sapply(parts, `[`, 1)
input <- sapply(parts, `[`, 2)
# Desired order
model_order <- c("LR","RF","GBM","CNN","MLP","RNN","GRU")
input_order <- c("CnTonly","CnTATAC","CnTIgG","CnTATACIgG","CnTATACseq","onehot")
# Reorder rows by section (input) then model
ord <- order(factor(input, levels=input_order), factor(model, levels=model_order))
mat <- mat[ord, ]
model <- model[ord]
input <- input[ord]

# Create labels for x-axis combining section and model
labels <- paste(model, input, sep="_")

# Make boxplot
pdf(file="SF2D_raw.pdf",width=12,height=6)
boxplot(
  t(mat),               # transpose so each row is a set of values
  at = seq_along(labels),las=2,
  col = rep(colors, times=6),
  xaxt = "n",
  outline = FALSE,
  ylab = "Score"
)
axis(1, at = seq_along(labels), labels = model, las = 2, cex.axis = 0.7)
# Add vertical lines to separate sections
section_sizes <- table(factor(input, levels=input_order))
section_ends <- cumsum(section_sizes)
abline(v = c(0.5, section_ends + 0.5), col = "grey50", lty = 2)
# Add section names in the middle of each group
section_starts <- c(1, head(section_ends+1, -1))
section_mids <- (section_starts + section_ends) / 2
mtext(input_order, side = 3, at = section_mids, line = 1, cex = 0.8, font = 2)
dev.off()


##### genome-wide cor
LR_EXP <- read.table("K562_modelCMP/H3K27me3CnTscore_vsEXPcor_LR.txt",row.names=1,header=T)
RF_EXP <- read.table("K562_modelCMP/H3K27me3CnTscore_vsEXPcor_RF.txt",row.names=1,header=T)
GBM_EXP <- read.table("K562_modelCMP/H3K27me3CnTscore_vsEXPcor_GBM.txt",row.names=1,header=T)
CNN_EXP <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsEXPcor_CNN.txt",row.names=1,header=T)
MLP_EXP <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsEXPcor_MLP.txt",row.names=1,header=T)
RNN_EXP <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsEXPcor_RNN.txt",row.names=1,header=T)
GRU_EXP <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsEXPcor_GRU.txt",row.names=1,header=T)

LR_HM <- read.table("K562_modelCMP/H3K27me3CnTscore_vsHMcor_LR.txt",row.names=1,header=T)
RF_HM <- read.table("K562_modelCMP/H3K27me3CnTscore_vsHMcor_RF.txt",row.names=1,header=T)
GBM_HM <- read.table("K562_modelCMP/H3K27me3CnTscore_vsHMcor_GBM.txt",row.names=1,header=T)
CNN_HM <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsHMcor_CNN.txt",row.names=1,header=T)
MLP_HM <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsHMcor_MLP.txt",row.names=1,header=T)
RNN_HM <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsHMcor_RNN.txt",row.names=1,header=T)
GRU_HM <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsHMcor_GRU.txt",row.names=1,header=T)

write.table(LR_EXP, file="F3_b_LR_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(RF_EXP, file="F3_b_RF_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(GBM_EXP, file="F3_b_GBM_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(CNN_EXP, file="F3_b_CNN_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(MLP_EXP, file="F3_b_MLP_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(RNN_EXP, file="F3_b_RNN_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(GRU_EXP, file="F3_b_GRU_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)


write.table(LR_HM, file="F3_c_LR_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(RF_HM, file="F3_c_RF_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(GBM_HM, file="F3_c_GBM_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(CNN_HM, file="F3_c_CNN_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(MLP_HM, file="F3_c_MLP_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(RNN_HM, file="F3_c_RNN_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(GRU_HM, file="F3_c_GRU_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)


#### F3b
colors <- c("#DD864E","#969696", "#E6E6E6","#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
df <- data.frame(
  Method = rep(c("LR", "RF", "GBM", "CNN", "MLP", "RNN", "GRU"),
               times = c(length(as.matrix(LR_EXP)),
                         length(as.matrix(RF_EXP)),
                         length(as.matrix(GBM_EXP)),
                         length(as.matrix(CNN_EXP)),
                         length(as.matrix(MLP_EXP)),
                         length(as.matrix(RNN_EXP)),
                         length(as.matrix(GRU_EXP)))),
  Value = c(as.numeric(as.matrix(LR_EXP)),
            as.numeric(as.matrix(RF_EXP)),
            as.numeric(as.matrix(GBM_EXP)),
            as.numeric(as.matrix(CNN_EXP)),
            as.numeric(as.matrix(MLP_EXP)),
            as.numeric(as.matrix(RNN_EXP)),
            as.numeric(as.matrix(GRU_EXP)))
)
# Factor to ensure order
df$Method <- factor(df$Method, levels = c("LR", "RF", "GBM", "CNN", "MLP", "RNN", "GRU"))
# Plot
ggplot(df, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  labs(title = "corEXP", y = "Value", x = "Method") +
  scale_y_continuous(limits = c(min(df[,2]),0)) +  # y-axis starts at 0
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none"  # remove color legend
  )
ggsave("F3B_raw.pdf",width=6,height=6)

#### F3C
df <- data.frame(
  Method = rep(c("LR", "RF", "GBM", "CNN", "MLP", "RNN", "GRU"),
               times = c(length(as.matrix(LR_HM)),
                         length(as.matrix(RF_HM)),
                         length(as.matrix(GBM_HM)),
                         length(as.matrix(CNN_HM)),
                         length(as.matrix(MLP_HM)),
                         length(as.matrix(RNN_HM)),
                         length(as.matrix(GRU_HM)))),
  Value = c(as.numeric(as.matrix(LR_HM)),
            as.numeric(as.matrix(RF_HM)),
            as.numeric(as.matrix(GBM_HM)),
            as.numeric(as.matrix(CNN_HM)),
            as.numeric(as.matrix(MLP_HM)),
            as.numeric(as.matrix(RNN_HM)),
            as.numeric(as.matrix(GRU_HM)))
)
# Factor to ensure order
df$Method <- factor(df$Method, levels = c("LR", "RF", "GBM", "CNN", "MLP", "RNN", "GRU"))
# Plot
ggplot(df, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  labs(title = "corHM", y = "Value", x = "Method") +
  scale_y_continuous(limits = c(min(df[,2]),max(df[,2]))) +  # y-axis starts at 0
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none"  # remove color legend
  )
ggsave("F3C_raw.pdf",width=6,height=6)


#### F3D

input_order <- c("CnTonly","CnTATAC","CnTIgG","CnTATACIgG","CnTATACseq","onehot")
# Convert LR_EXP to data frame in long format
df <- data.frame(
  Input = rep(colnames(LR_EXP), each = nrow(LR_EXP)),
  Value = as.numeric(as.matrix(LR_EXP)),
  stringsAsFactors = FALSE
)
# Keep only finite values
df <- df[is.finite(df$Value), ]
# Apply factor order
df$Input <- factor(df$Input, levels = input_order)
# Plot
ggplot(df, aes(x = Input, y = Value, color = Input)) +
  geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
  scale_color_manual(values = rep("#DD864E", length(input_order))) +
  scale_y_continuous(limits = c(min(df[,2]),0)) +  # y-axis starts at 0
  labs(title = "LR_EXP", y = "Value", x = "Input") +
  geom_jitter(width = 0.15,size = 2,alpha = 0.8) +#
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
ggsave("F3D_withdots.pdf",width=6,height=6)

### F3E
input_order <- c("CnTonly","CnTATAC","CnTIgG","CnTATACIgG","CnTATACseq","onehot")
# Convert LR_HM to data frame in long format
df <- data.frame(
  Input = rep(colnames(LR_HM), each = nrow(LR_HM)),
  Value = as.numeric(as.matrix(LR_HM)),
  stringsAsFactors = FALSE
)
# Keep only finite values
df <- df[is.finite(df$Value), ]
# Apply factor order
df$Input <- factor(df$Input, levels = input_order)
# Plot
ggplot(df, aes(x = Input, y = Value, color = Input)) +
  geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
  scale_color_manual(values = rep("#DD864E", length(input_order))) +
  scale_y_continuous(limits = c(min(df[,2]),0)) +  # y-axis starts at 0
  labs(title = "LR_HM", y = "Value", x = "Input") +  
  geom_jitter(width = 0.15,size = 2,alpha = 0.8) +#
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
ggsave("F3E_withdots.pdf",width=6,height=6)





#### SF2E
colors <- c("#DD864E","#969696", "#E6E6E6","#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
##### genome-wide cor
LR_EXP <- read.table("K562_modelCMP/H3K27me3CnTscore_vsEXPcor_LR.txt",row.names=1,header=T)
RF_EXP <- read.table("K562_modelCMP/H3K27me3CnTscore_vsEXPcor_RF.txt",row.names=1,header=T)
GBM_EXP <- read.table("K562_modelCMP/H3K27me3CnTscore_vsEXPcor_GBM.txt",row.names=1,header=T)
CNN_EXP <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsEXPcor_CNN.txt",row.names=1,header=T)
MLP_EXP <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsEXPcor_MLP.txt",row.names=1,header=T)
RNN_EXP <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsEXPcor_RNN.txt",row.names=1,header=T)
GRU_EXP <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsEXPcor_GRU.txt",row.names=1,header=T)
colnames(LR_EXP) <- paste0("LR_",colnames(LR_EXP))
colnames(RF_EXP) <- paste0("RF_",colnames(RF_EXP))
colnames(GBM_EXP) <- paste0("GBM_",colnames(GBM_EXP))
colnames(CNN_EXP) <- paste0("CNN_",colnames(CNN_EXP))
colnames(MLP_EXP) <- paste0("MLP_",colnames(MLP_EXP))
colnames(RNN_EXP) <- paste0("RNN_",colnames(RNN_EXP))
colnames(GRU_EXP) <- paste0("GRU_",colnames(GRU_EXP))
mat <- t(cbind(LR_EXP,RF_EXP,GBM_EXP,CNN_EXP,MLP_EXP,RNN_EXP,GRU_EXP))
# Extract model and input from row names
parts <- strsplit(rownames(mat), "_")
model <- sapply(parts, `[`, 1)
input <- sapply(parts, `[`, 2)
# Desired order
model_order <- c("LR","RF","GBM","CNN","MLP","RNN","GRU")
input_order <- c("CnTonly","CnTATAC","CnTIgG","CnTATACIgG","CnTATACseq","onehot")
# Reorder rows by section (input) then model
ord <- order(factor(input, levels=input_order), factor(model, levels=model_order))
mat <- mat[ord, ]
model <- model[ord]
input <- input[ord]
# Create labels for x-axis combining section and model
labels <- paste(model, input, sep="_")
# Make boxplot
pdf(file="SF2E_raw.pdf",width=12,height=6)
boxplot(
  t(mat),               # transpose so each row is a set of values
  at = seq_along(labels),las=2,
  col = rep(colors, times=6),
  xaxt = "n",
  outline = FALSE,
  ylab = "Score"
)
axis(1, at = seq_along(labels), labels = model, las = 2, cex.axis = 0.7)
# Add vertical lines to separate sections
section_sizes <- table(factor(input, levels=input_order))
section_ends <- cumsum(section_sizes)
abline(v = c(0.5, section_ends + 0.5), col = "grey50", lty = 2)
# Add section names in the middle of each group
section_starts <- c(1, head(section_ends+1, -1))
section_mids <- (section_starts + section_ends) / 2
mtext(input_order, side = 3, at = section_mids, line = 1, cex = 0.8, font = 2)
dev.off()


#### SF2F
LR_HM <- read.table("K562_modelCMP/H3K27me3CnTscore_vsHMcor_LR.txt",row.names=1,header=T)
RF_HM <- read.table("K562_modelCMP/H3K27me3CnTscore_vsHMcor_RF.txt",row.names=1,header=T)
GBM_HM <- read.table("K562_modelCMP/H3K27me3CnTscore_vsHMcor_GBM.txt",row.names=1,header=T)
CNN_HM <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsHMcor_CNN.txt",row.names=1,header=T)
MLP_HM <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsHMcor_MLP.txt",row.names=1,header=T)
RNN_HM <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsHMcor_RNN.txt",row.names=1,header=T)
GRU_HM <- read.table("K562_dnnCMP/H3K27me3CnTscore_vsHMcor_GRU.txt",row.names=1,header=T)

colnames(LR_HM) <- paste0("LR_",colnames(LR_HM))
colnames(RF_HM) <- paste0("RF_",colnames(RF_HM))
colnames(GBM_HM) <- paste0("GBM_",colnames(GBM_HM))
colnames(CNN_HM) <- paste0("CNN_",colnames(CNN_HM))
colnames(MLP_HM) <- paste0("MLP_",colnames(MLP_HM))
colnames(RNN_HM) <- paste0("RNN_",colnames(RNN_HM))
colnames(GRU_HM) <- paste0("GRU_",colnames(GRU_HM))
mat <- t(cbind(LR_HM,RF_HM,GBM_HM,CNN_HM,MLP_HM,RNN_HM,GRU_HM))
# Extract model and input from row names
parts <- strsplit(rownames(mat), "_")
model <- sapply(parts, `[`, 1)
input <- sapply(parts, `[`, 2)
# Desired order
model_order <- c("LR","RF","GBM","CNN","MLP","RNN","GRU")
input_order <- c("CnTonly","CnTATAC","CnTIgG","CnTATACIgG","CnTATACseq","onehot")
# Reorder rows by section (input) then model
ord <- order(factor(input, levels=input_order), factor(model, levels=model_order))
mat <- mat[ord, ]
model <- model[ord]
input <- input[ord]
# Create labels for x-axis combining section and model
labels <- paste(model, input, sep="_")
# Make boxplot
pdf(file="SF2F_raw.pdf",width=12,height=6)
boxplot(
  t(mat),               # transpose so each row is a set of values
  at = seq_along(labels),las=2,
  col = rep(colors, times=6),
  xaxt = "n",
  outline = FALSE,
  ylab = "Score"
)
axis(1, at = seq_along(labels), labels = model, las = 2, cex.axis = 0.7)
# Add vertical lines to separate sections
section_sizes <- table(factor(input, levels=input_order))
section_ends <- cumsum(section_sizes)
abline(v = c(0.5, section_ends + 0.5), col = "grey50", lty = 2)
# Add section names in the middle of each group
section_starts <- c(1, head(section_ends+1, -1))
section_mids <- (section_starts + section_ends) / 2
mtext(input_order, side = 3, at = section_mids, line = 1, cex = 0.8, font = 2)
dev.off()


#### F4 ij
allGSM <- c("GSM3536507","GSM3536508","GSM3536509","GSM3536510","GSM3536511","GSM3536512","GSM3536513","GSM3536515","GSM3560261","GSM3680215","GSM3680216","GSM3680217","GSM3680218","GSM3680219","GSM3680220","GSM3680221","GSM3680222")
all_corEXP <- c()
for(GSM in allGSM){
  thisdata <- read.table(paste0("PATTY_vs_macs2/K562_H3K27me3_corFeatures_",GSM,".txt"),row.names=1,header=T)
  all_corEXP <- rbind(all_corEXP, thisdata[,"corEXP"])
}
rownames(all_corEXP) <- allGSM
colnames(all_corEXP) <- c("patty","m2","m2IgG")

write.table(all_corEXP, file="F4_i_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)

# Convert to long format
df <- data.frame(Sample = rownames(all_corEXP), all_corEXP, check.names = FALSE)
df_long <- reshape2::melt(df, id.vars = "Sample", variable.name = "Method", value.name = "Value")
# Keep only finite values
df_long <- df_long[is.finite(df_long$Value), ]
# Order of boxes
method_order <- c( "m2", "m2IgG","patty")
df_long$Method <- factor(df_long$Method, levels = method_order)
# Colors for each box
colors <- c( "red", "darkred","#DD864E")  # adjust as needed
names(colors) <- method_order
ggplot(df_long, aes(x = Method, y = Value, color = Method)) +
  geom_boxplot( na.rm = TRUE) +
  scale_color_manual(values = colors) +
  scale_y_continuous(limits = c(min(df_long$Value, na.rm = TRUE), NA)) +
  labs(title = "all_corEXP", y = "Value", x = "Method") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
ggsave("F4I_raw.pdf",width=6,height=6)



all_corHM <- c()
for(GSM in allGSM){
  thisdata <- read.table(paste0("PATTY_vs_macs2/K562_H3K27me3_corFeatures_",GSM,".txt"),row.names=1,header=T)
  all_corHM <- rbind(all_corHM, thisdata[,"corHM"])
}
rownames(all_corHM) <- allGSM
colnames(all_corHM) <- c("patty","m2","m2IgG")
write.table(all_corHM, file="F4_j_sourceData.txt",row.names=T,col.names=T,sep="\t",quote=F)

# Convert to long format
df <- data.frame(Sample = rownames(all_corHM), all_corHM, check.names = FALSE)
df_long <- reshape2::melt(df, id.vars = "Sample", variable.name = "Method", value.name = "Value")
# Keep only finite values
df_long <- df_long[is.finite(df_long$Value), ]
# Order of boxes
method_order <- c( "m2", "m2IgG","patty")
df_long$Method <- factor(df_long$Method, levels = method_order)
# Colors for each box
colors <- c( "red", "darkred","#DD864E")  # adjust as needed
names(colors) <- method_order
ggplot(df_long, aes(x = Method, y = Value, color = Method)) +
  geom_boxplot( na.rm = TRUE) +
  scale_color_manual(values = colors) +
  scale_y_continuous(limits = c(min(df_long$Value, na.rm = TRUE), NA)) +
  labs(title = "all_corHM", y = "Value", x = "Method") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
ggsave("F4J_raw.pdf",width=6,height=6)



###### F5
colors <- c( "red", "darkred","#DD864E")  # adjust as needed

bar3plot <- function(dataname, Feature, YL){
  thisdata <- read.table(paste0("PATTY_vs_macs2/",dataname),row.names=1,header=T)
  plotdata <- thisdata[c("m2","m2IgG","PATTY"),Feature]
  barplot(plotdata,las=2,col=colors,ylim=YL)
  abline(h=0)
  box(lwd=1)
}
pdf(file="F5A_raw.pdf",width=2*4,height=3)
par(mfrow=c(1,4))
bar3plot("K562_H3K27me3_corFeatures_GSM4308145.txt","corEXP", c(-0.55,0))
bar3plot("K562_H3K27me3_corFeatures_GSM4308146.txt","corEXP", c(-0.55,0))
bar3plot("K562_H3K27me3_corFeatures_GSM4308147.txt","corEXP", c(-0.55,0))
bar3plot("K562_H3K27me3_corFeatures_GSM4308148.txt","corEXP", c(-0.55,0))
dev.off()

pdf(file="F5B_raw.pdf",width=2*4,height=3)
par(mfrow=c(1,4))
bar3plot("K562_H3K27me3_corFeatures_GSM4308145.txt","corHM",c(-0.42, 0.05))
bar3plot("K562_H3K27me3_corFeatures_GSM4308146.txt","corHM",c(-0.42, 0.05))
bar3plot("K562_H3K27me3_corFeatures_GSM4308147.txt","corHM",c(-0.42, 0.05))
bar3plot("K562_H3K27me3_corFeatures_GSM4308148.txt","corHM",c(-0.42, 0.05))
dev.off()

pdf(file="F5C_raw.pdf",width=2*9,height=3)
par(mfrow=c(1,9))
bar3plot("K562_H3K27me3_corFeatures_GSM5549925.txt","corEXP",c(-0.55,0.2))
bar3plot("K562_H3K27me3_corFeatures_GSM5549926.txt","corEXP",c(-0.55,0.2))
bar3plot("K562_H3K27me3_corFeatures_GSM5549927.txt","corEXP",c(-0.55,0.2))
bar3plot("K562_H3K27me3_corFeatures_GSM5549928.txt","corEXP",c(-0.55,0.2))
bar3plot("K562_H3K27me3_corFeatures_GSM5549929.txt","corEXP",c(-0.55,0.2))
bar3plot("K562_H3K27me3_corFeatures_GSM5549930.txt","corEXP",c(-0.55,0.2))
bar3plot("K562_H3K27me3_corFeatures_GSM5549931.txt","corEXP",c(-0.55,0.2))
bar3plot("K562_H3K27me3_corFeatures_GSM5549932.txt","corEXP",c(-0.55,0.2))
bar3plot("K562_H3K27me3_corFeatures_GSM5549933.txt","corEXP",c(-0.55,0.2))
dev.off()

pdf(file="F5D_raw.pdf",width=2*9,height=3)
par(mfrow=c(1,9))
bar3plot("K562_H3K27me3_corFeatures_GSM5549925.txt","corHM",c(-0.45,0.25))
bar3plot("K562_H3K27me3_corFeatures_GSM5549926.txt","corHM",c(-0.45,0.25))
bar3plot("K562_H3K27me3_corFeatures_GSM5549927.txt","corHM",c(-0.45,0.25))
bar3plot("K562_H3K27me3_corFeatures_GSM5549928.txt","corHM",c(-0.45,0.25))
bar3plot("K562_H3K27me3_corFeatures_GSM5549929.txt","corHM",c(-0.45,0.25))
bar3plot("K562_H3K27me3_corFeatures_GSM5549930.txt","corHM",c(-0.45,0.25))
bar3plot("K562_H3K27me3_corFeatures_GSM5549931.txt","corHM",c(-0.45,0.25))
bar3plot("K562_H3K27me3_corFeatures_GSM5549932.txt","corHM",c(-0.45,0.25))
bar3plot("K562_H3K27me3_corFeatures_GSM5549933.txt","corHM",c(-0.45,0.25))
dev.off()


pdf(file="F5E_raw.pdf",width=2*1,height=3)
par(mfrow=c(1,1))
bar3plot("K562_H3K27ac_corFeatures_GSM3536514.txt","corEXP",c(0,0.7))
dev.off()
pdf(file="F5F_raw.pdf",width=2*1,height=3)
par(mfrow=c(1,1))
bar3plot("K562_H3K27ac_corFeatures_GSM3536514.txt","corHM",c(-0.13,0))
dev.off()


pdf(file="F5G_raw.pdf",width=2*6,height=3)
par(mfrow=c(1,6))
bar3plot("H1_H3K27me3_corFeatures_GSM5549907.txt","corEXP",c(-0.25,0.12))
bar3plot("H1_H3K27me3_corFeatures_GSM5549908.txt","corEXP",c(-0.25,0.12))
bar3plot("H1_H3K27me3_corFeatures_GSM5549909.txt","corEXP",c(-0.25,0.12))
bar3plot("H1_H3K27me3_corFeatures_GSM5549910.txt","corEXP",c(-0.25,0.12))
bar3plot("H1_H3K27me3_corFeatures_GSM5549911.txt","corEXP",c(-0.25,0.12))
bar3plot("H1_H3K27me3_corFeatures_GSM5549912.txt","corEXP",c(-0.25,0.12))
dev.off()

pdf(file="F5H_raw.pdf",width=2*6,height=3)
par(mfrow=c(1,6))
bar3plot("H1_H3K27me3_corFeatures_GSM5549907.txt","corHM",c(-0.38,0.22))
bar3plot("H1_H3K27me3_corFeatures_GSM5549908.txt","corHM",c(-0.38,0.22))
bar3plot("H1_H3K27me3_corFeatures_GSM5549909.txt","corHM",c(-0.38,0.22))
bar3plot("H1_H3K27me3_corFeatures_GSM5549910.txt","corHM",c(-0.38,0.22))
bar3plot("H1_H3K27me3_corFeatures_GSM5549911.txt","corHM",c(-0.38,0.22))
bar3plot("H1_H3K27me3_corFeatures_GSM5549912.txt","corHM",c(-0.38,0.22))
dev.off()

pdf(file="F5I_raw.pdf",width=2*1,height=3)
par(mfrow=c(1,1))
bar3plotEXP("H1_H3K27ac_corFeatures_GSM3536497.txt")
dev.off()
pdf(file="F5J_raw.pdf",width=2*1,height=3)
par(mfrow=c(1,1))
bar3plotHM("H1_H3K27ac_corFeatures_GSM3536497.txt")
dev.off()

### F6 HI


pdf(file="F6H_raw.pdf",width=2*1,height=3)
par(mfrow=c(1,1))
thisdata <- read.table(paste0("PATTY_vs_macs2/HCT116_H3K9me3_corFeatures_1.txt"),row.names=1,header=T)
plotdata <- thisdata[c("m2","PATTY"),"corEXP"]
barplot(plotdata,las=2,col=colors[c(1,3)])
abline(h=0)
box(lwd=1)
dev.off()

pdf(file="F6I_raw.pdf",width=2*1,height=3)
par(mfrow=c(1,1))
thisdata <- read.table(paste0("PATTY_vs_macs2/HCT116_H3K9me3_corFeatures_1.txt"),row.names=1,header=T)
plotdata <- thisdata[c("m2","PATTY"),"corHM"]
barplot(plotdata,las=2,col=colors[c(1,3)])
abline(h=0)
box(lwd=1)
dev.off()




### get pearsonCor matrix
getMAT <- function(dataname, Feature,note){
  thisdata <- read.table(paste0("PATTY_vs_macs2/",dataname),row.names=1,header=T)
  plotdata <- thisdata[c("m2","m2IgG","PATTY"),Feature]
  return(c(note, dataname, Feature,plotdata))
}


outdata <- rbind(
getMAT("K562_H3K27me3_corFeatures_GSM4308145.txt","corEXP", "Figure5A"),
getMAT("K562_H3K27me3_corFeatures_GSM4308146.txt","corEXP", "Figure5A"),
getMAT("K562_H3K27me3_corFeatures_GSM4308147.txt","corEXP", "Figure5A"),
getMAT("K562_H3K27me3_corFeatures_GSM4308148.txt","corEXP", "Figure5A"),
getMAT("K562_H3K27me3_corFeatures_GSM4308145.txt","corHM","Figure5B"),
getMAT("K562_H3K27me3_corFeatures_GSM4308146.txt","corHM","Figure5B"),
getMAT("K562_H3K27me3_corFeatures_GSM4308147.txt","corHM","Figure5B"),
getMAT("K562_H3K27me3_corFeatures_GSM4308148.txt","corHM","Figure5B"),
getMAT("K562_H3K27me3_corFeatures_GSM5549925.txt","corEXP","Figure5C"),
getMAT("K562_H3K27me3_corFeatures_GSM5549926.txt","corEXP","Figure5C"),
getMAT("K562_H3K27me3_corFeatures_GSM5549927.txt","corEXP","Figure5C"),
getMAT("K562_H3K27me3_corFeatures_GSM5549928.txt","corEXP","Figure5C"),
getMAT("K562_H3K27me3_corFeatures_GSM5549929.txt","corEXP","Figure5C"),
getMAT("K562_H3K27me3_corFeatures_GSM5549930.txt","corEXP","Figure5C"),
getMAT("K562_H3K27me3_corFeatures_GSM5549931.txt","corEXP","Figure5C"),
getMAT("K562_H3K27me3_corFeatures_GSM5549932.txt","corEXP","Figure5C"),
getMAT("K562_H3K27me3_corFeatures_GSM5549933.txt","corEXP","Figure5C"),
getMAT("K562_H3K27me3_corFeatures_GSM5549925.txt","corHM","Figure5D"),
getMAT("K562_H3K27me3_corFeatures_GSM5549926.txt","corHM","Figure5D"),
getMAT("K562_H3K27me3_corFeatures_GSM5549927.txt","corHM","Figure5D"),
getMAT("K562_H3K27me3_corFeatures_GSM5549928.txt","corHM","Figure5D"),
getMAT("K562_H3K27me3_corFeatures_GSM5549929.txt","corHM","Figure5D"),
getMAT("K562_H3K27me3_corFeatures_GSM5549930.txt","corHM","Figure5D"),
getMAT("K562_H3K27me3_corFeatures_GSM5549931.txt","corHM","Figure5D"),
getMAT("K562_H3K27me3_corFeatures_GSM5549932.txt","corHM","Figure5D"),
getMAT("K562_H3K27me3_corFeatures_GSM5549933.txt","corHM","Figure5D"),
getMAT("K562_H3K27ac_corFeatures_GSM3536514.txt","corEXP","Figure5E"),
getMAT("K562_H3K27ac_corFeatures_GSM3536514.txt","corHM","Figure5F"),
getMAT("H1_H3K27me3_corFeatures_GSM5549907.txt","corEXP","Figure5G"),
getMAT("H1_H3K27me3_corFeatures_GSM5549908.txt","corEXP","Figure5G"),
getMAT("H1_H3K27me3_corFeatures_GSM5549909.txt","corEXP","Figure5G"),
getMAT("H1_H3K27me3_corFeatures_GSM5549910.txt","corEXP","Figure5G"),
getMAT("H1_H3K27me3_corFeatures_GSM5549911.txt","corEXP","Figure5G"),
getMAT("H1_H3K27me3_corFeatures_GSM5549912.txt","corEXP","Figure5G"),
getMAT("H1_H3K27me3_corFeatures_GSM5549907.txt","corHM","Figure5H"),
getMAT("H1_H3K27me3_corFeatures_GSM5549908.txt","corHM","Figure5H"),
getMAT("H1_H3K27me3_corFeatures_GSM5549909.txt","corHM","Figure5H"),
getMAT("H1_H3K27me3_corFeatures_GSM5549910.txt","corHM","Figure5H"),
getMAT("H1_H3K27me3_corFeatures_GSM5549911.txt","corHM","Figure5H"),
getMAT("H1_H3K27me3_corFeatures_GSM5549912.txt","corHM","Figure5H"),
getMAT("H1_H3K27ac_corFeatures_GSM3536497.txt","corHM","Figure5I"),
getMAT("H1_H3K27ac_corFeatures_GSM3536497.txt","corHM","Figure5J"))
write.table(outdata, file="F5_cor_table.txt",row.names=F,col.names=F,sep="\t",quote=F)




methods <- rownames(plotdata_H3K27me3_cb)
#colors <- c("#5A76B4","#5A76B4","#5A76B4","#5A76B4","#DD864E","#A7A6A6","#A7A6A6")
colors <- c("#DD864E","#969696", "#E6E6E6","#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
names(colors) <- rownames(plotdata_H3K27me3_cb)
  
df <- as.data.frame(plotdata_H3K27me3_cb)
#colnames(df) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
df$Method <- methods
df_long <- gather(df, Fold, Accuracy, -Method)
desired_order <- names(colors)#c("CNNseq", "RF_out", "RFseq_out", "GBM_out", "GBMseq_out", "LR_out", "LRseq_out")
df_long$Method <- factor(df_long$Method, levels = desired_order)
ggplot(df_long, aes(x=Method, y=Accuracy, fill=Method)) +
  geom_boxplot() +
  scale_fill_manual(values=colors) +
  theme_minimal() +
  scale_y_continuous(limits=c(0.7,0.95)) + # Setting the y limits here
  labs(title="Cross Validation Results", y="Accuracy", x="Methods") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.border = element_rect(color="black", fill=NA, size=1.5))
ggsave("H3K27me3_crossVal_accuracy_boxplot_methodCB.pdf")




### H3K27ac
dNN <- read.table("K562_H3K27ac_modelAcc.txt")
H3K27ac_dNN_best <- c()
for(METHOD in c("CNN","GRU","MLP","RNN")){
  for(Feature in c("CnTATACIgG","CnTATACseq","CnTATAC",'CnTIgG',"CnTonly","onehot")){
    thisdata <- dNN[which(dNN[,1] == METHOD & dNN[,2] == Feature),  1:11]
    outdata <- thisdata[order(apply(thisdata[,7:11],1,mean),decreasing=T),][1,]
    H3K27ac_dNN_best <- rbind(H3K27ac_dNN_best, outdata)
  }
}
other_H3K27ac <- read.table("K562_H3K27ac_otherModels_crossValidationAcc.txt")
cbdata <- rbind(as.matrix(H3K27ac_dNN_best[,c(1,2,7:11)]),as.matrix(other_H3K27ac[c(7:12,13:18,1:6),]))
rownames(cbdata) <- paste0(cbdata[,1],"_",cbdata[,2])
plotdata_H3K27ac <- cbdata[,3:7]
plotdata_H3K27ac <- matrix(as.numeric(plotdata_H3K27ac), nrow = nrow(plotdata_H3K27ac), ncol = ncol(plotdata_H3K27ac))
rownames(plotdata_H3K27ac) <- rownames(cbdata)

plotdata_H3K27ac <- rbind(
  plotdata_H3K27ac[grep("CNN_",rownames(plotdata_H3K27ac)),],
  plotdata_H3K27ac[grep("RNN_",rownames(plotdata_H3K27ac)),],
  plotdata_H3K27ac[grep("GRU_",rownames(plotdata_H3K27ac)),],
  plotdata_H3K27ac[grep("MLP_",rownames(plotdata_H3K27ac)),],
  plotdata_H3K27ac[grep("LR_",rownames(plotdata_H3K27ac)),] ,
  plotdata_H3K27ac[grep("RF_",rownames(plotdata_H3K27ac)),] ,
  plotdata_H3K27ac[grep("GBM_",rownames(plotdata_H3K27ac)),])

plotdata_H3K27ac_cb <- rbind(
  as.numeric(plotdata_H3K27ac[grep("CNN_",rownames(plotdata_H3K27ac)),]),
  as.numeric(plotdata_H3K27ac[grep("RNN_",rownames(plotdata_H3K27ac)),]),
  as.numeric(plotdata_H3K27ac[grep("GRU_",rownames(plotdata_H3K27ac)),]),
  as.numeric(plotdata_H3K27ac[grep("MLP_",rownames(plotdata_H3K27ac)),]),
  as.numeric(plotdata_H3K27ac[grep("LR_",rownames(plotdata_H3K27ac)),] ),
  as.numeric(plotdata_H3K27ac[grep("RF_",rownames(plotdata_H3K27ac)),] ),
  as.numeric(plotdata_H3K27ac[grep("GBM_",rownames(plotdata_H3K27ac)),]))
rownames(plotdata_H3K27ac_cb) <- c("CNN","RNN","GRU","MLP","LR","RF","GBM")




methods <- rownames(plotdata_H3K27ac)
#colors <- c("#5A76B4","#5A76B4","#5A76B4","#5A76B4","#DD864E","#A7A6A6","#A7A6A6")
colors <- rep(c("#C6DBEF", "#6BAED6", "#2171B5", "#08306B","#DD864E","#969696", "#E6E6E6"),each=6)
names(colors) <- rownames(plotdata_H3K27ac)

df <- as.data.frame(plotdata_H3K27ac)
#colnames(df) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
df$Method <- methods
df_long <- gather(df, Fold, Accuracy, -Method)
desired_order <- names(colors)#c("CNNseq", "RF_out", "RFseq_out", "GBM_out", "GBMseq_out", "LR_out", "LRseq_out")
df_long$Method <- factor(df_long$Method, levels = desired_order)
ggplot(df_long, aes(x=Method, y=Accuracy, fill=Method)) +
  geom_boxplot() +
  scale_fill_manual(values=colors) +
  theme_minimal() +
  scale_y_continuous(limits=c(0.78,0.95)) + # Setting the y limits here
  labs(title="Cross Validation Results", y="Accuracy", x="Methods") + guides(fill = FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.border = element_rect(color="black", fill=NA, size=1.5))
ggsave("H3K27ac_crossVal_accuracy_boxplot_methodFeature.pdf",width=10)




methods <- rownames(plotdata_H3K27ac_cb)
#colors <- c("#5A76B4","#5A76B4","#5A76B4","#5A76B4","#DD864E","#A7A6A6","#A7A6A6")
colors <- c("#C6DBEF", "#6BAED6", "#2171B5", "#08306B","#DD864E","#969696", "#E6E6E6")
names(colors) <- rownames(plotdata_H3K27ac_cb)

df <- as.data.frame(plotdata_H3K27ac_cb)
#colnames(df) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
df$Method <- methods
df_long <- gather(df, Fold, Accuracy, -Method)
desired_order <- names(colors)#c("CNNseq", "RF_out", "RFseq_out", "GBM_out", "GBMseq_out", "LR_out", "LRseq_out")
df_long$Method <- factor(df_long$Method, levels = desired_order)
ggplot(df_long, aes(x=Method, y=Accuracy, fill=Method)) +
  geom_boxplot() +
  scale_fill_manual(values=colors) +
  theme_minimal() +
  scale_y_continuous(limits=c(0.78,0.95)) + # Setting the y limits here
  labs(title="Cross Validation Results", y="Accuracy", x="Methods") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.border = element_rect(color="black", fill=NA, size=1.5))
ggsave("H3K27ac_crossVal_accuracy_boxplot_methodCB.pdf")





#### 5-fold cross validation on strictly T/F peaks
LR_out <- c(0.8195488721804511, 0.8045112781954887, 0.8327067669172933, 0.8571428571428571, 0.8116760828625236)
LRseq_out<- c(0.7255639097744361, 0.7518796992481203, 0.7669172932330827, 0.768796992481203, 0.7288135593220338)
CNNseq <- c(0.9285714030265808,0.9248120188713074,0.9229323267936707,0.9210526347160339,0.9191729426383972)
RF_out <- c(0.8909774436090225, 0.9060150375939849, 0.9304511278195489, 0.9191729323308271, 0.9152542372881356)
RFseq_out<-c(0.8834586466165414, 0.9078947368421053, 0.9285714285714286, 0.9191729323308271, 0.9133709981167608)
GBM_out <- c(0.881578947368421, 0.8984962406015038, 0.9342105263157895, 0.924812030075188, 0.9152542372881356)
GBMseq_out <- c(0.8834586466165414,0.9078947368421053, 0.9360902255639098, 0.9285714285714286, 0.9133709981167608)

outdata <- rbind(CNNseq,LR_out, LRseq_out, RF_out,RFseq_out,GBM_out,GBMseq_out)
methods <- rownames(outdata)#c("CNNseq", "LR_out", "LRseq_out", "RF_out", "RFseq_out", "GBM_out", "GBMseq_out")




# Google-style color scheme
colors <- c(
  "LR_out" = "#BBDEFB", 
  "RF_out" = "#C8E6C9", 
  "GBM_out" = "#FFF9C4", 
  "LRseq_out" = "#2196F3",  # A slightly lighter shade of red
  "RFseq_out" = "#4CAF50",  # A slightly lighter shade of green
  "GBMseq_out" = "#FFEB3B",  # A slightly lighter shade of yellow
  "CNNseq" = "#F44336" 
)

  
df <- as.data.frame(outdata)
colnames(df) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
df$Method <- methods
df_long <- gather(df, Fold, Accuracy, -Method)
desired_order <- names(colors)#c("CNNseq", "RF_out", "RFseq_out", "GBM_out", "GBMseq_out", "LR_out", "LRseq_out")
df_long$Method <- factor(df_long$Method, levels = desired_order)
ggplot(df_long, aes(x=Method, y=Accuracy, fill=Method)) +
  geom_boxplot() +
  scale_fill_manual(values=colors) +
  theme_minimal() +
  scale_y_continuous(limits=c(0.7,0.95)) + # Setting the y limits here
  labs(title="5-Fold Cross Validation Results", y="Accuracy", x="Methods") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.border = element_rect(color="black", fill=NA, size=1.5))
ggsave("crossVal_accuracy_boxplot.pdf")


#### K562 vs Exp, vs K27ac sig
K562_ChIP_evaluation <- read.table("../K562_ChIP_evaluation/K562_H3K27me3_ChIP_macs2peak_evaluation.txt",row.names=1,header=T)
c("macs_chip",mean(K562_ChIP_evaluation[,"vsChIPsigK27ac"]),sd(K562_ChIP_evaluation[,"vsChIPsigK27ac"]))

K562_Pexp <- read.table("../CBmodelAVEsig_epochs10k/vsEXPpromoter_cor_summary.txt",row.names=1)
K562_K27ac <- read.table("../CBmodelAVEsig_epochs10k/vsH3K27acChIPsig_cor_summary.txt",row.names=1)

all_models <- c("CNN","GBM","RF","LR")
summary_colnames <- c("macs_enrich","macsQ",
                      paste0(rep(c("onehot_",
                                   "onehot2mer_",
                                   #"onehotRandom_",
                                   #"onehotPermute_",
                                   "simplex_",
                                   "simplex2mer_",
                                   "Tn5bias_",
                                   #"DNasebias_",
                                   "CnTATACIgG_",
                                   "CnTATAC_",
                                   "CnTIgG_",
                                   "CnTonly_",
                                   "CnTATACseq_"
                      ),each=length(all_models)),
                      rep(all_models, 9)))

colnames(K562_Pexp) <- summary_colnames
colnames(K562_K27ac) <- summary_colnames

colors <- c("macs_enrich"="#2962FF", 
            "CnTATACseq_CNN"="#FFD54F", 
            "CnTATACseq_RF"="#0F9D58", 
            "CnTATACseq_GBM"="#9C27B0", 
            "CnTATACseq_LR"="#FF5252")


feature_bar_CnT <- function(indata, outfile,YLIM,YLAB,usecolor ){
  data_long <- indata %>%
    gather(key="variable", value="value", names(usecolor) )
  # Calculate means and standard deviations for the selected columns
  summary_stats <- data_long %>% 
    group_by(variable) %>%
    summarise(mean=mean(value), sd=sd(value))
  summary_stats$variable <- factor(summary_stats$variable, 
                                   levels=names(usecolor))
  ggplot(summary_stats, aes(x=variable, y=mean, fill=variable)) +
    geom_bar(stat="identity", position="dodge", width=0.7) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.25, position=position_dodge(width=0.7)) +
    scale_fill_manual(values=usecolor) +
    scale_y_continuous(limits=YLIM) + # Setting the y limits here
    theme_minimal() +
    labs(y=YLAB, x="Method", fill="Method") +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          panel.border = element_rect(color="black", fill=NA, size=1.5))
  ggsave(outfile)
}
feature_bar_CnT(K562_Pexp,"K562_Pexp_barplots_section1_CnT.pdf",c(-0.3, 0) ,"pearson R vs expression",
            c("macs_enrich"="grey", 
              "CnTATACseq_LR"="#2196F3",
              "CnTATACseq_RF"="#4CAF50", 
              "CnTATACseq_GBM"="#FFEB3B", 
              "CnTATACseq_CNN"="#F44336" 
            )
)
feature_bar_CnT(K562_Pexp,"K562_Pexp_barplots_section2_CnT.pdf" ,c(-0.3, 0),"pearson R vs expression",
            c( "macs_enrich"="grey",
               "CnTonly_LR"="#E3F2FD",
               "CnTATAC_LR"="#BBDEFB",
               "CnTATACseq_LR"="#1976D2")
)

feature_bar_CnT(K562_K27ac,"K562_K27ac_barplots_section1_CnT.pdf",c(-0.35, 0.05) ,"pearson R vs H3K27ac sig",
            c("macs_enrich"="grey", 
              "CnTATACseq_LR"="#2196F3",
              "CnTATACseq_RF"="#4CAF50", 
              "CnTATACseq_GBM"="#FFEB3B", 
              "CnTATACseq_CNN"="#F44336" 
            )
)
feature_bar_CnT(K562_K27ac,"K562_K27ac_barplots_section2_CnT.pdf",c(-0.35, 0.05) ,"pearson R vs H3K27ac sig",
            c( "macs_enrich"="grey",
               "CnTonly_LR"="#E3F2FD",
               "CnTATAC_LR"="#BBDEFB",
               "CnTATACseq_LR"="#1976D2")
)


feature_bar <- function(indata, outfile,YLIM,YLAB,usecolor ,ChIPdata){
  data_long <- indata %>%
    gather(key="variable", value="value", names(usecolor)[-1] )
  # Calculate means and standard deviations for the selected columns
  summary_stats <- data_long %>% 
    group_by(variable) %>%
    summarise(mean=mean(value), sd=sd(value))
  summary_stats <- summary_stats %>%
    add_row(variable = ChIPdata[1], mean = as.numeric(ChIPdata[2]), sd = as.numeric(ChIPdata[3]))
  summary_stats$variable <- factor(summary_stats$variable, 
                                   levels=names(usecolor))
  #print(ChIPdata)
  print(summary_stats)
  # Plot using ggplot2
  ggplot(summary_stats, aes(x=variable, y=mean, fill=variable)) +
    geom_bar(stat="identity", position="dodge", width=0.7) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.25, position=position_dodge(width=0.7)) +
    scale_fill_manual(values=usecolor) +
    scale_y_continuous(limits=YLIM) + # Setting the y limits here
    theme_minimal() +
    labs(y=YLAB, x="Method", fill="Method") +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          panel.border = element_rect(color="black", fill=NA, size=1.5))
  ggsave(outfile)
}
feature_bar(K562_Pexp,"K562_Pexp_barplots_section1.pdf",c(-0.3, 0) ,"pearson R vs expression",
            c("macs2_chip"="black",
              "macs_enrich"="grey", 
              "CnTATACseq_LR"="#2196F3",
              "CnTATACseq_RF"="#4CAF50", 
              "CnTATACseq_GBM"="#FFEB3B", 
              "CnTATACseq_CNN"="#F44336" 
              ),
            c("macs2_chip",mean(K562_ChIP_evaluation[,"vsEXPpromoter"]),sd(K562_ChIP_evaluation[,"vsEXPpromoter"]))
)
feature_bar(K562_Pexp,"K562_Pexp_barplots_section2.pdf" ,c(-0.3, 0),"pearson R vs expression",
            c( "macs2_chip"="black",
               "macs_enrich"="grey",
               "CnTonly_LR"="#E3F2FD",
               "CnTATAC_LR"="#BBDEFB",
               "CnTATACseq_LR"="#1976D2"),
            c("macs2_chip",mean(K562_ChIP_evaluation[,"vsEXPpromoter"]),sd(K562_ChIP_evaluation[,"vsEXPpromoter"]))
)

feature_bar(K562_K27ac,"K562_K27ac_barplots_section1.pdf",c(-0.35, 0.05) ,"pearson R vs H3K27ac sig",
            c("macs2_chip"="black",
              "macs_enrich"="grey", 
              "CnTATACseq_LR"="#2196F3",
              "CnTATACseq_RF"="#4CAF50", 
              "CnTATACseq_GBM"="#FFEB3B", 
              "CnTATACseq_CNN"="#F44336" 
            ),
            c("macs2_chip",mean(K562_ChIP_evaluation[,"vsChIPsigK27ac"]),sd(K562_ChIP_evaluation[,"vsChIPsigK27ac"]))
)
feature_bar(K562_K27ac,"K562_K27ac_barplots_section2.pdf",c(-0.35, 0.05) ,"pearson R vs H3K27ac sig",
            c( "macs2_chip"="black",
               "macs_enrich"="grey",
               "CnTonly_LR"="#E3F2FD",
               "CnTATAC_LR"="#BBDEFB",
               "CnTATACseq_LR"="#1976D2"),
            c("macs2_chip",mean(K562_ChIP_evaluation[,"vsChIPsigK27ac"]),sd(K562_ChIP_evaluation[,"vsChIPsigK27ac"]))
)


                         
                         
H1_Pexp <- read.table("../CBmodelAVEsig_applyH1/vsEXPpromoter_cor_summary.txt",row.names=1)
H1_K27ac <- read.table("../CBmodelAVEsig_applyH1/vsH3K27acChIPsig_cor_summary.txt",row.names=1)

all_models <- c("CNN","CNNtune","LR")
summary_colnames <- c("macs_enrich","macsQ",
                      paste0(rep(c("CnTATAC_",
                                   "CnTonly_"
                      ),each=length(all_models)),
                      rep(all_models, 2)))

colnames(H1_Pexp) <- summary_colnames
colnames(H1_K27ac) <- summary_colnames



colors <- c("macs_enrich"="grey", 
            "CnTATAC_CNNtune"="#F44336", 
            "CnTATAC_LR"="#BBDEFB")

# Plot using ggplot2

data_long <- H1_Pexp %>%
  gather(key="variable", value="value", c("macs_enrich", "CnTATAC_CNNtune", "CnTATAC_LR"))
# Calculate means and standard deviations for the selected columns
summary_stats <- data_long %>% 
  group_by(variable) %>%
  summarise(mean=mean(value), sd=sd(value))

summary_stats$variable <- factor(summary_stats$variable, 
                                 levels=c("macs_enrich", "CnTATAC_CNNtune", "CnTATAC_LR"))

# Plot using ggplot2
ggplot(summary_stats, aes(x=variable, y=mean, fill=variable)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.25, position=position_dodge(width=0.7)) +
  scale_fill_manual(values=colors) +
  theme_minimal() +
  labs(y="Value", x="Variable", fill="Variable") +
  theme(axis.text.x=element_text(angle=45, hjust=1),panel.border = element_rect(color="black", fill=NA, size=1.5))
ggsave("H1_Pexp_barplots.pdf")
print(plot)





#### H1 K27ac signal
data_long <- H1_K27ac %>%
  gather(key="variable", value="value", c("macs_enrich", "CnTATAC_CNNtune", "CnTATAC_LR"))
# Calculate means and standard deviations for the selected columns
summary_stats <- data_long %>% 
  group_by(variable) %>%
  summarise(mean=mean(value), sd=sd(value))

summary_stats$variable <- factor(summary_stats$variable, 
                                 levels=c("macs_enrich", "CnTATAC_CNNtune", "CnTATAC_LR"))

# Plot using ggplot2
ggplot(summary_stats, aes(x=variable, y=mean, fill=variable)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.25, position=position_dodge(width=0.7)) +
  scale_fill_manual(values=colors) +
  theme_minimal() +
  labs(y="Value", x="Variable", fill="Variable") +
  theme(axis.text.x=element_text(angle=45, hjust=1),panel.border = element_rect(color="black", fill=NA, size=1.5))
ggsave("H1_K27acSig_barplots.pdf")


