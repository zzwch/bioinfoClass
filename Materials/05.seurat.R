setwd("/Users/hazard/Downloads/rstudio-export/")
options(stringsAsFactors = F)
rm(list = ls())
geneset <- list()
geneset$g1s <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2",
                 "Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2",
                 "Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7",
                 "Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1",
                 "Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1",
                 "Chaf1b","Brip1","E2f8")
geneset$g2m <- c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80",
                 "Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a",
                 "Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e",
                 "Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk",
                 "Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8",
                 "Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5",
                 "Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")
### data manipulation
#tmp <- read.delim(file = "refGene.UMI.count.xls", header = T, sep = "\t")
load(file = "05.seurat.Rdata")
tmp <- reshape2::colsplit(colnames(umi_all), pattern = "_sc", names = c("Library", "Barcode"))
tmp1 <- reshape2::colsplit(tmp$Library, "_", c("Location","Batch","Lane"))
umi_annot <- cbind(tmp1[,c("Location","Batch")], tmp[,"Barcode",drop=F])
rownames(umi_annot) <- colnames(umi_all)
#annot <- subset(umi_annot, Barcode %in% paste0("sc", setdiff(1:50, c(23,34))))
umi_annot$Location <- Hmisc::capitalize(gsub(pattern = '\\.[123]', replacement = "",umi_annot$Location))
annot <- umi_annot
umi <- umi_all[,rownames(annot)]
### librarry
library(Seurat)
library(RColorBrewer)
### colors
pal_rainbow <- c("red","orange","yellow","green","cyan","blue","purple")
bertie.color <- c("#D0B38A", "#A3D171", "#533E88", "#7957A3", "#000000", "#E63325", "#0A8041", "#C5208E", "#3CBDED", "#3B55A3", "#D691BE", "#D23E28", "#6474B6", "#4288C8", "#80A469", "#FFCF3F", "#FBCC9F", "#F06360", "#437BB2", "#A43E40", "#206767", "#779E43", "#258950", "#F7D059", "#ED803B") 
pal_dolphin <- c("#FF00AE", "#A020F1", "#000000", "#0403E5", "#FF8C01", "#8B0101", "#007502", "#FE0000", "#FFFF01", "#FF99CB", "#4A95FB", "#61FE69", "#9A7A01", "#017F8B", "#05FDFF", "grey")
pal_cluster <- c("#FF00AE", "#A020F1", "#000000", "#0403E5", "#FF8C01", 
                 "#8B0101", "#007502", "#FE0000", "#FFFF01", "#FF99CB", 
                 "#4A95FB", "#61FE69", "#9A7A01", "#017F8B", "#05FDFF",
                 "#D0B38A", "#533E88", "#D23E28", "#80A469", "#F06360")
pal_location <- brewer.pal(12, name = "Paired")
pal_stage <- c("indianred", "steelblue", "blue3")
### Create Seurat Object
pbmc <- CreateSeuratObject(raw.data = umi, meta.data = annot,
                           min.cells = 0, min.genes =-1, is.expr = 0, 
                           normalization.method = "LogNormalize", scale.factor = 1e5, do.scale = T, do.center = T, 
                           names.field = 1)
#pbmc <- NormalizeData(pbmc,scale.factor = 1e5)
# pbmc <- ScaleData(pbmc)
## Cell QC
VlnPlot(pbmc, c("nGene","nUMI"), group.by = "Location")
pbmc <- FilterCells(pbmc, subset.names = "nGene", low.thresholds = 2000)
pbmc <- FilterCells(pbmc, subset.names = "nUMI", low.thresholds = 20000, high.thresholds = 1e6)

VlnPlot(pbmc, c("Pecam1","Ptprc","Cdh5","Cd34"), nCol = 2, x.lab.rot = T, group.by = "Location", cols.use = pal_location)
## HVGs
pbmc <- FindVariableGenes(pbmc, x.low.cutoff = 0.8, x.high.cutoff = 8, y.cutoff = 1, y.high.cutoff = Inf)
length(pbmc@var.genes)
## PCA
pbmc <- RunPCA(object = pbmc, pcs.compute = 50)
PCAPlot(pbmc, pt.size= 3, cols.use = pal_location, group.by = "Location", dim.1 = 3, dim.2 = 4)
## confounding - cell cycle
pbmc <- CellCycleScoring(pbmc, g2m.genes = geneset$g2m, s.genes = geneset$g1s)
FeaturePlot(pbmc, features.plot = c("S.Score", "G2M.Score", "Mki67", "Hmmr"), reduction.use = "pca")
PCAPlot(pbmc, pt.size= 3, cols.use = pal_rainbow, group.by = "Phase", dim.1 = 1, dim.2 =2)
PCHeatmap(pbmc, pc.use = 9)
## tSNE
PCElbowPlot(pbmc, num.pc = 50)
pbmc <- RunTSNE(pbmc, reduction.use = "pca", dims.use = 1:10)
pbmc <- RunTSNE(pbmc, reduction.use = "pca", dims.use = setdiff(1:10,9))
pbmc@data <- as.matrix(pbmc@data)
pbmc <- ScaleData(pbmc, vars.to.regress = c("nUMI"), do.cpp = T)
TSNEPlot(pbmc, pt.size = 3, group.by = "Location",colors.use = pal_location)
## DimPlot
## Clustering
## DEGs
## Heatmap