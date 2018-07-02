options(stringsAsFactors = F)

nature_fpkm <- read.csv(file = "Supplementary Table 3 RefSeq.Gene.Expression.Single.Cells.fixDate.csv", row.names = 1, check.names = F)
nature_annot <- reshape2::colsplit(colnames(nature_fpkm), " #", c("CellType", "CellOrder"))
rownames(nature_annot) <- colnames(nature_fpkm)
nature_annot$Batch <- "Nature"

annot <- nature_annot
#annot <- subset(nature_annot, CellType %in% c("EC", "T1 pre-HSC", "T2 pre-HSC", "E12 HSC", "E14 HSC", "Adult HSC", "PreHSC"))
fpkm <- nature_fpkm[,rownames(annot)]

# What's Principal Component Analysis?
library(ggplot2)
library(RColorBrewer)

data <- log2(fpkm+1)[names(which(rowSums(fpkm > 1) > 3)), ]
data.sd <- apply(data, 1, sd)
data.gene <- names(sort(data.sd, decreasing = T)[1:1000])

# Dimension Reduction - PCA
res <- prcomp(data[data.gene,], scale. = T)

ggplot(cbind(annot[rownames(res$rotation),], res$rotation)) + 
  geom_point(aes(PC1, PC2, color = CellType, shape= Batch), size = 3) + 
  scale_color_manual(values = brewer.pal(12, "Paired")) + 
  theme_bw()

gridExtra::grid.arrange(grobs = lapply(c("Pecam1","Procr", "Itga2b","Spn","Ptprc","Kit","Slamf1", "Hoxa5", "Kdr","Hlf", 
                                         "Nrp2", "Nr2f2","Runx1", "Gfi1","Myb","Hmmr", "Mki67", "Plk1"), function(x){
  ggplot(cbind(annot[rownames(res$rotation),], res$rotation, t(data[,rownames(res$rotation)]))) + 
    geom_point(aes_string("PC1", "PC2", color = x, shape= "Batch"), size = 3) + 
    scale_colour_gradientn(colours = c("grey", "yellow", "red")) +
    ggtitle(x) + 
    theme_classic() +
    theme(plot.title = element_text(face = "bold.italic", size = 16, hjust = 0.5))
}))

# Hierarchical Clustering
hc <- hclust(dist(t(data[data.gene,])), method = "ward.D")
plot(hc)
plot(hc, hang = -1)

pheatmap::pheatmap(data[data.gene[1:100],], show_colnames = F, annotation_col = annot)

# DEG
degRes <- list()
for(i in sort(unique(annot$CellType))){
  typeIndicator <- (annot$CellType == i) + 0
  degRes[[i]] <- apply(data, 1, function(x) {
    cor(typeIndicator, x)
    })
}
tmp <- as.data.frame(degRes)
tmp1 <- apply(tmp, 2, function(x){
  rownames(tmp)[order(x, decreasing = T)[1:10]]
})
pheatmap::pheatmap(data[as.vector(tmp1), order(annot$CellType)], show_colnames = F, annotation_col = annot, 
                   cluster_cols = F, cluster_rows = F, scale = "row", )

#################################################################
plot(hc, col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", 
     col.axis = "#F38630", lwd = 3, lty = 3, sub = "", hang = -1, axes = FALSE)
axis(side = 2, at = seq(0, 1400, 100), col = "#F38630", labels = FALSE, 
     lwd = 2)
# add text in margin
mtext(seq(0, 1400, 100), side = 2, at = seq(0, 1400, 100), line = 1, 
      col = "#A38630", las = 2)

hcd = as.dendrogram(hc)
plot(hcd)
plot(hcd, type = "triangle")
cut(hcd, h = 400)
plot(cut(hcd, h = 400)$upper, main = "Upper tree of cut at h=400")
plot(cut(hcd, h = 400)$lower[[5]], main = "Second branch of lower tree with cut at h=400")

#
labelColors = dichromat::colorschemes$Categorical.12[1:9]
clusMember = cutree(hc, k = 9)

colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    #print(n)
    #print(a)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol, pch = c(NA, 15), col = labCol)
  }
  n
}
clusDendro = dendrapply(hcd, colLab)
plot(clusDendro, main = "Cool Dendrogram")

library(ape)
plot(as.phylo(hc), type = "unroot", show.node.label = T, show.tip.label = T, tip.color = labelColors[clusMember])
