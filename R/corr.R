# plot correlations of real data

rm(list = ls())

sim_data_dir <- "~/Documents/MasterProject/R/illustrations/"

X <- t(as.matrix(read.table(file.path(sim_data_dir, "X.txt"), header = FALSE)))
colnames(X) <- paste0("snp_", 1:ncol(X))

cX <- abs(cor(X))

require(gplots)

hmcols<-colorRampPalette(c("grey98", "black"))(256)
par(cex.main = 0.8, cex.lab = 0.8)

heatmap.2(cX, dendrogram = "none", col = hmcols, Rowv = F,
          Colv = F, main = paste0("SNP correlation pattern"),
          density.info = "none", trace = "none",
          lhei = c(1, 4),
          lwid = c(0.9, 4),
          key = TRUE,
          keysize = 2.5,
          key.xlab = "",
          key.title = "",
          labRow = F, labCol = F)# margins = c(5, 5))
