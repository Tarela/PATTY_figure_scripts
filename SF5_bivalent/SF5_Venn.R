library(VennDiagram)
library(grid)
#txt <- 'ID\tpeak_total\tPATTYuniq\tm2uniq\tOV
#GSM3560257\t61061\t12607\t19584\t28870
#GSM3560258\t47387\t13477\t15830\t18080
#GSM3560259\t54404\t23294\t17777\t13333
#GSM3680214\t127304\t5170\t26040\t96094
#'
#df <- read.table(text = txt, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
#
col_patty <- "#D78551"
col_m2    <- "#E12A2D"


pdf("bivalent_venn.pdf")
grid.newpage()
vp <- draw.pairwise.venn(
  area1 = 152+838,
  area2 = 838+1013,
  cross.area = 838,
  category = c("PATTY", "m2"),
  fill = c(col_patty, col_m2),
  alpha = c(1, 1),
  cat.col = c(col_patty, col_m2),
  lty = "solid",
  cex = 1.2,
  cat.cex = 1.2,
  main = sprintf("%s  (|PATTY|=%d, |m2|=%d, overlap=%d)", id, A, B, OV),
  ind = FALSE
)
grid.draw(vp)
dev.off()


for (i in seq_len(nrow(df))) {
  id <- df$ID[i]
  Auniq <- df$PATTYuniq[i]
  Buniq <- df$m2uniq[i]
  OV    <- df$OV[i]
  
  A <- Auniq + OV
  B <- Buniq + OV
  
  png(file.path(out_dir, paste0(id, "_venn.png")), width = 1200, height = 1000, res = 200)
  grid.newpage()
  vp <- draw.pairwise.venn(
    area1 = A,
    area2 = B,
    cross.area = OV,
    category = c("PATTY", "m2"),
    fill = c(col_patty, col_m2),
    alpha = c(0.55, 0.55),
    cat.col = c(col_patty, col_m2),
    lty = "solid",
    cex = 1.2,
    cat.cex = 1.2,
    main = sprintf("%s  (|PATTY|=%d, |m2|=%d, overlap=%d)", id, A, B, OV),
    ind = FALSE
  )
  grid.draw(vp)
  dev.off()
}

