library(clusterProfiler)
# bitr - ID translator
# GO analysis
# KEGG analysis
# Disease analysis
# Reactome pathway analysis
# DAVID functional analysis
# Universal enrichment analysis - using MsigDB gene set collections
# Visualization

load("human.3batch.clean.pbmc.v3.Rdata")
c5 <- subset(degs, cluster == "C5_Lymphatic")$gene
c5.fc <- subset(degs, cluster == "C5_Lymphatic")$avg_logFC

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
c5.idmap <- bitr(geneID = c5, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
c5.id <- c5.idmap$ENTREZID[match(c5, c5.idmap$SYMBOL)]

ggo <- groupGO(gene = c5, 
               OrgDb = org.Hs.eg.db, 
               keytype = "SYMBOL", 
               ont = "CC", 
               level = 5, 
               readable = F)
ego <- enrichGO(gene = c5, 
                OrgDb = org.Hs.eg.db, 
                keyType = "SYMBOL", 
                ont = "BP", # "BP", "CC", "MF", "ALL" 
                pvalueCutoff = 0.05, 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                minGSSize = 10, 
                maxGSSize = 500, 
                readable = F, 
                pool = F)
# ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)
# dropGO(x = ego, level = 5, term = NULL)
ego1 <- simplify(x = ego)
ego2 <- gofilter(x = ego, level = 4)
barplot(ego2, 
        colorBy = "p.adjust", showCategory = 10, font.size = 10, title = "clusterProfiler results") # see help
dotplot(ego, 
        colorBy = "p.adjust", showCategory = 10, font.size = 10, title = "clusterProfiler results") # see help
enrichMap(ego)
cnetplot(ego, categorySize = "pvalue", foldChange = setNames(c5.fc, c5), font.size = 1)
plotGOgraph(ego)

# using GSEA 
ego3 <- gseGO(geneList = sort(setNames(c5.fc, nm = c5), decreasing = T),
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              nPerm = 1000,
              minGSSize = 10,
              maxGSSize = 500,
              pvalueCutoff = 0.05, keyType = "SYMBOL", exponent = 2, pAdjustMethod = "BH")
dotplot(ego3)

# enrich KEGG
kk <- enrichKEGG(gene = c5.id, organism = "hsa", keyType = "kegg")
kk2 <- gseKEGG(geneList = sort(setNames(c5.fc, nm = c5.id), decreasing = T),
               organism = "hsa", 
               nPerm = 1000,
               minGSSize = 10,
               pvalueCutoff = 0.05, verbose = F)
gseaplot(kk2, geneSetID = "hsa04144")
browseKEGG(kk, 'hsa04144')
library("pathview")
hsa04110 <- pathview(gene.data  = sort(setNames(c5.fc, nm = c5.id), decreasing = T),
                     pathway.id = "hsa04144",
                     species    = "hsa",
                     limit      = list(gene=2, cpd=1))

# other databases
enrichDO() # Disease Ontology (DO) Semantic and Enrichment analysis
gseDO()
enrichNCG() # Network of Cancer Gene
gseNCG()
enrichDGN() # Disease Gene Network
gseDGN()
enrichPathway() # ReactomePA, require(ReactomePA)
gsePathway()
enrichDAVID() # using DAVIDWebService

# Using MSigDB gene set collections
enricher(gene, TERM2GENE = read.gmt('file to gmt'))
GSEA(geneList, TERM2GENE = read.gmt('file to gmt'))
seq2gene() #require(ChIPseeker)

# compare clusters
clusterGenes <- list()
for(i in unique(sort(degs$cluster))){
  ids <- bitr(geneID = subset(degs, cluster == i)$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = T)
  clusterGenes <- c(clusterGenes, list(ids$ENTREZID))
}
names(clusterGenes) <- unique(sort(degs$cluster))

ck <- compareCluster(geneCluster = clusterGenes, fun = "enrichKEGG")
dotplot(ck)
