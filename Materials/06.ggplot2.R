setwd("/data2/lzc/test")
options(stringsAsFactors = F)
rm(list = ls())
load("mouse.ec.Rdata")


library(ggplot2)
# using data.frame
ggData <- cbind(pbmc@meta.data, pbmc@dr$tsne@cell.embeddings, t(as.matrix(pbmc@data)))

# http://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html
# scatter plot
ggplot() + geom_point(mapping = aes(x = nUMI, y = nGene, color = Batch, size = nUMI), data = ggData)
ggplot(ggData) + geom_point(mapping = aes(x = nUMI, y = nGene))
ggplot(data = ggData, mapping = aes(x = nUMI, y = nGene)) + geom_point() + geom_density_2d()

# scatter plot using different color and shape
ggplot(data = ggData, mapping = aes(x = nUMI, y = nGene)) + geom_point(aes(color = Cluster, shape = nGene, size = nUMI))
ggplot(data = ggData, mapping = aes(x = nUMI, y = nGene)) + 
  geom_point(aes(color = Cluster, shape = Cluster, size = nUMI, alpha = nGene), color = "yellow")
# labs, themes, scale_
ggplot(data = ggData, mapping = aes(x = nUMI, y = nGene)) + 
  geom_point(aes(color = Cluster, size = nUMI)) +
  labs(x = "point_x", y = "point_y", title = "point____") + 
  scale_color_manual(values = c("red", "yellow", "blue", "green")) + 
  #scale_size_continuous(trans = "log10") 
  #scale_shape_manual(values = c(16, 17))
  theme(axis.text.x = element_text(color = "red", angle = 45, hjust = 1),
        panel.grid.major = element_line(colour = "grey"), panel.grid.minor = element_line(colour = "grey90"))
# geom_histogram
ggplot(ggData) + geom_histogram(aes(x = nUMI, fill = Cluster))
# geom_density_2d
ggplot(data = ggData, mapping = aes(x = nUMI, y = nGene)) + geom_point() + geom_density_2d(aes(color = ..level..))

# scatter plot with smooth line
ggplot(ggData, aes(nUMI, nGene)) + geom_point() + geom_smooth(se = F)
ggplot(ggData, aes(nUMI, nGene)) + geom_point(aes(color = Cluster)) + geom_smooth(aes(group = Cluster, color = Cluster), se = F) + 
  facet_wrap(facets = ~Cluster, nrow = 1, scales = "free") 

# scatter plot with encircle
# install.packages("ggalt")
ggplot(data = ggData, mapping = aes(x = tSNE_1, y = tSNE_2)) + 
  geom_point(aes(color = Cluster)) + 
  geom_encircle()   # encircle

# lollipop: segment, point, text
ggData2 <- as.data.frame(as.matrix(table(ggData$Cluster, ggData$Location)))
ggplot(ggData2) + 
  geom_point(aes(x = 1:nrow(ggData2), y = Freq), size = 5) + 
  geom_segment(aes(x = 1:nrow(ggData2), y =0,yend = Freq, xend = 1:nrow(ggData2)))+
  geom_text(aes(x = 1:nrow(ggData2), y = Freq, label = Freq), color = "red")
# boxplot, violin
boxplot(formula = nGene ~ Cluster, data = ggData, outline = F, col = c("red","blue","yellow","green"))
ggplot(ggData, mapping = aes(Cluster, nGene, fill = Cluster)) + 
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, fill = "white")
# bar
table(ggData$Cluster, ggData$Location)
ggplot(ggData) + 
  geom_bar(aes(Cluster, fill = Location, group = Location), stat = "count", position = "stack") +
  theme(axis.text.x.top = element_text(angle = 45))
# facet_warp
# ggplot list
library(gridExtra)
grid.arrange()
ggGrobs <- lapply(c("Bmx","Runx1","Pecam1"), function(gene){
  ggplot(ggData) + geom_point(aes_string("tSNE_1", "tSNE_2", color = gene)) + 
    scale_color_gradientn(colours = c("red","orange","yellow","green","cyan","blue", "purple"))
})

grid.arrange(grobs = ggGrobs, nrow = 1)
