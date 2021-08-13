#load libraries
library(flowCore)
library(uwot)
library(matrixStats)
library(reshape2)
library(ggplot2)
library(viridis)
library(FlowSOM)
library(ConsensusClusterPlus)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(shiny)
#make sure you also have run the flow functions .R file before any code or it won't work
#load, compensate, and format flow files
fcs.file.paths <- list.files("./data_modified/200909_D21_Invitro_ECR/subtypes/", full.names = T, pattern = ".fcs")
input.flowset <- read.flowSet(fcs.file.paths, transformation = F, truncate_max_range = F)
input.flowset <- fset.compensate(input.flowset)
elgcl <- elgcl.pars(input.flowset)
input.flowset <- fset.transform(input.flowset, elgcl)
input.flowset <- markernames.flow(input.flowset)
dat <- as.data.frame(fsApply(input.flowset, exprs))
head(dat)
#choose dimentions to run umap
dims <- c("CCL1", "IFNg", "Perforin", "TNF", "GzmA", "CD107a", "IL2")
dims.for.umap <- dims
#run umap. This took 1:01 minutes second on a computer with 16 cores (this code tells your computer to use all but one of its cores)
umap.dims <- umap(dat[, dims.for.umap], scale = TRUE, parallel::detectCores()-1, verbose = TRUE)
#format data
colnames(umap.dims) <- paste("umap", seq(ncol(umap.dims)), sep = ".")
mdat <- data.frame(sample = sampleNames(input.flowset),
                   phenotype.subset = sub(".fcs", "", sapply(sampleNames(input.flowset), function(i) strsplit(i, "_")[[1]][5])),
                   row.names = NULL)
mdat$phenotype.subset <- factor(mdat$phenotype.subset, levels = c("CD49a", "DP", "DN", "CD103"))
dat.umap.samples <- data.frame(dat,
                               umap.dims, 
                               sample = rep(sampleNames(input.flowset), fsApply(input.flowset, nrow)),
                               stim = "CD3",
                               tissue = "Lung")
dat.umap.samples <- merge(dat.umap.samples, mdat, by = "sample")
#create a directory to put results in
#dir.create("./data_results/umap/200909", recursive = T)
#save the file to be loaded later
#saveRDS(dat.umap.samples, file.path("./data_results/umap/200909", paste0(paste("CD8_SUBSETS_umap", Sys.Date(), sep = "_"), ".rds")))
#read a saved RDS file so the time consuming step of running umap doesn't have to be rerun.
dat.umap.samples <- readRDS("./data_results/umap/200909/CD8_SUBSETS_umap_2020-09-14.rds")
#prepare data for plotting
#select dimensions to be used for showing expression of each marker, though note that CD49a and CD103, 
#as shown above, were not used in the generation of the umap
dims <- c("CCL1", "IFNg", "Perforin", "TNF", "GzmA", "CD107a", "IL2",'CD49a','CD103')
dims.for.umap <- dims
dat.umap <- data.frame(dat.umap.samples[, dims.for.umap],
                       dat.umap.samples[, c("umap.1", "umap.2")])

#scale the data
dat.umap[, dims.for.umap] <- range.scaled(dat.umap[, dims.for.umap])
#melt the data for plotting
dat.melt <- melt(dat.umap, id = c("umap.1", "umap.2"), value.name = "expression", variable.name = "marker")
#plot the umap
umap.dims.plot3 <- ggplot(dat.melt, aes(umap.1, umap.2, z = expression)) +
  stat_summary_hex(bins = 100) +
  scale_fill_viridis(option = "plasma", name = "expression") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom") +
  facet_wrap(~marker,ncol = 2) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face='bold',color='black',size='11'),
        panel.background = element_blank())
umap.dims.plot3
#save the plot
ggsave(file.path("./data_results/umap/200909/", paste0("CD8_DIMS_umap_", Sys.Date(), ".pdf")), plot = umap.dims.plot3, height = 8, width = 8, units = "in")
#plot cell count faceted by subset
dat.umap <- dat.umap.samples[, c("umap.1", "umap.2", "phenotype.subset")]
events <- data.frame(setNames(data.frame(table(dat.umap$phenotype.subset)), nm = c("phenotype.subset", "total")),
                     umap.1 = max(dat.umap$umap.1)*.9,
                     umap.2 = max(dat.umap$umap.2))

umap.subset.plot <- ggplot(dat.umap, aes(umap.1, umap.2)) +
  geom_hex(aes(colour = ..count..), bins = 100) +
  scale_fill_viridis() +
  scale_color_viridis() +
  theme_dark() +
  geom_text(data = events, mapping = aes(label = prettyNum(total, big.mark = ","))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  facet_wrap(~phenotype.subset)

umap.subset.plot
#save plot
ggsave(file.path("./data_results/umap/200909/", paste0("CD8_SUBSETS_umap_", Sys.Date(), ".pdf")), plot = umap.subset.plot, height = 6, width = 8, units = "in")

#cluster with flowsom: get data ready
fsom <- ReadInput(input.flowset, scale = T)
fsom$elgcl <- elgcl
fsom$markers <- colnames(input.flowset)
fsom$seed <- 20040501
fsom$dims.used <- grep(paste0(dims, collapse = "|"), fsom$markers, value = T)
set.seed(fsom$seed)
#choose the number of nodes, which is the product of xdim and ydim below. These nodes are then osrted into clusters below.
fsom <- BuildSOM(fsom, colsToUse = sort(match(fsom$dims.used, fsom$markers)), xdim = 10, ydim = 10)
fsom <- BuildMST(fsom)
PlotStars(fsom)
#this number is the actual number of clusters
fsom$nbclust <- 18
#create directory to save stuff
dir.create("./data_results/fsoms/temp/Clustering Consensus", recursive = T)
#cluster the data
fsom$clustering <- ConsensusClusterPlus(t(fsom$map$codes), 
                                        maxK = fsom$nbclust, 
                                        reps = 100, 
                                        pItem = 1, 
                                        pFeature = 1, 
                                        title = "./data_results/fsoms/temp/Clustering Consensus", 
                                        plot = "png",
                                        clusterAlg = "hc", 
                                        innerLinkage = "average", 
                                        finalLinkage = "average",
                                        distance = "euclidean", 
                                        seed = 4242)
#prepare plots with new cluster data
fsom[["metaclustering"]] <- as.factor(fsom$clustering[[fsom$nbclust]]$consensusClass)
PlotStars(fsom, backgroundValues = fsom$metaclustering, legend = FALSE)
PlotNumbers(fsom, backgroundValues = fsom$metaclustering)
#use this to explore the data more fully. Note that only six markers can be viewed at once, so one of these genes must be commented out
fsom.explorer(fsom, 
              #"IL2", 
              "IFNg", 
              "TNF", 
              "CD107a", 
              "Perforin", 
              "GzmA", 
              "Ccl1",
              20000)
#pick dimensions to use in heatmap
dims <- c( "IL2", 
           "IFNg", 
           "TNF", 
           "CD107a", 
           "Perforin", 
           "GzmA", 
           "CCL1",
           "CD49a",
           "CD103")
#format data for heatmap plotting
fsom$dims.used <- grep(paste0(dims, collapse = "|"), fsom$markers, value = T)
cluster.medians <- as.data.frame(fsom.cluster.medians(fsom, dat.type = 'map', ranged = TRUE))
annotation_row <- cluster.medians["cluster"]
annotation_colors <- list(cluster = setNames(colorRampPalette(viridis(9))(nlevels(fsom$metaclustering)), levels(annotation_row$cluster)))
col.index <- fsom$markers[fsom$map$colsUsed]
cluster_rows <- hclust(dist(cluster.medians[, -grep("cluster", colnames(cluster.medians))][, col.index], method = "euclidean"), method = "average")
clustering_table <- as.numeric(table(fsom$metaclustering[fsom$map$mapping[, 1]]))
clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)

cluster.medians$clusternumber<-paste0('Cluster ',1:18)
labels_row <- paste0(cluster.medians$clusternumber#, " (", clustering_prop, "%)"#, " (", clustering_table, ")"
)
heatmapfeed<-cluster.medians[,dims]
library(tidyverse)
#plot heatmaps
pheatmap(cluster.medians[, dims], color = plasma(50),
         number_color = "black",
         fontsize = 8, fontsize_number = 6, 
         annotation_colors = annotation_colors)
clusterheatmap<-pheatmap(heatmapfeed, color = plasma(50),
                         cluster_rows = cluster_rows, labels_row =labels_row,
                         display_numbers = F, number_color = "black",
                         fontsize = 11, fontsize_number = 6)
clusterheatmap
#save plot
ggsave(file.path("./data_results/fsoms/", paste0("CD8_heatmap_", Sys.Date(), ".png")), plot = clusterheatmap,dpi=600)
#make directory to save stuff
dir.create("./data_results/fsoms/200909", recursive = T)
#save as a pdf
pdf(file.path("./data_results/fsoms/200909", paste0(paste("CD8_SUBSETS", Sys.Date(), "HeatMap", sep = "_"), ".pdf")))
fsom.heatmap(fsom, dat.type = "codes", dims = fsom$dims.used, clust.freq = T, clust.count = F, display_numbers = F, color.cut = 10)
dev.off()
#save this fsom as a RDS file so code doesn't need to be rerun
#saveRDS(fsom, file.path("./data_results/fsoms/200909", paste0(paste("CD8_SUBSETS", Sys.Date(), "fsom", sep = "_"), ".rds")))
#load RDS files
dat.umap.samples <- readRDS("./data_results/umap/200909/CD8_SUBSETS_umap_2020-09-14.rds")
fsom <- readRDS("./data_results/fsoms/200909/CD8_SUBSETS_2020-09-14_fsom.rds")
#format data and plot umap colored by cluster
fsom$node.counts <- matrix(table(fsom$map$mapping[,1 ]), ncol = 1, dimnames = list(NULL, "node.counts"))
fsom$umap.codes <- matrix(umap(X = fsom$map$codes, scale = TRUE, verbose = T), ncol = 2, dimnames = list(NULL, c("umap.1", "umap.2")))
umap.dat <- data.frame(range.scaled(fsom$map$codes), fsom$umap.codes, fsom$node.counts, cluster = fsom$metaclustering)
umap.dat.melt <- melt(umap.dat, id.vars = c("umap.1", "umap.2", "node.counts", "cluster"), value.name = "Expression", variable.name = "Marker")
fsom$umap.plots <- list(umap.dims = umap.dat.melt %>% ggplot(aes(umap.1, umap.2, color = Expression, fill = Expression)) + 
                          #geom_point(aes(size = node.counts), alpha = 0.7) +
                          geom_point(shape = 21, color = "black", stroke = 0.25, alpha = 0.7, aes(size = node.counts)) +
                          #scale_color_viridis(option = "plasma") +
                          scale_fill_viridis(option = "plasma") +
                          scale_size_area(max_size = 10) +
                          theme_dark() +
                          theme(axis.text = element_blank(),
                                axis.ticks = element_blank(),
                                panel.grid = element_blank(),
                                axis.title = element_blank(),
                                strip.text = element_text(size = 14)) +
                          facet_wrap(~Marker),
                        umap.clusters = umap.dat.melt %>% ggplot(aes(umap.1, umap.2, color = cluster, fill = cluster)) + 
                          #geom_point(aes(size = node.counts)) + 
                          geom_point(shape = 21, color = "black", stroke = 0.5, alpha = 0.7, aes(size = node.counts)) +
                          geom_polygon(data = plyr::ddply(umap.dat, "cluster", function(umap.dat) umap.dat[chull(umap.dat$umap.1, umap.dat$umap.2), ]), 
                                       alpha = 0.4, linetype = 0) +
                          scale_size_area(max_size = 10) +
                          #scale_color_manual(values = fsom$annotation$cluster.color) +
                          #scale_fill_manual(values = fsom$annotation$cluster.color) +
                          theme_void() +
                          theme(legend.position = c(.7,.7),
                                legend.text = element_text(size = 14)) +
                          guides(size = F)
)
dat <- data.frame(dat.umap.samples,
                  node = fsom$map$mapping[, 1],
                  cluster = fsom$metaclustering[fsom$map$mapping[,1]]
)
umapcoloredbycluster<-ggplot(dat, aes(umap.1, umap.2, color = cluster)) + geom_point(size = 1.5,alpha=.1)
umapcoloredbycluster
#count cells in each cluster. 
clusters.tabulated <- lapply(setNames(nm = unique(as.character(dat$sample))), function(i) table(dat$cluster[grep(i, dat$sample)]))
cluster.frame <- as.data.frame((bind_rows(clusters.tabulated)))
colnames(cluster.frame) <- paste("cluster", 1:ncol(cluster.frame), sep = ".")
rownames(cluster.frame) <- names(clusters.tabulated)
#make a directory to save this
#dir.create("./data_results/cluster_counts/200909", recursive = TRUE)
#save these cluster counts
write.csv(cluster.frame, file = file.path("./data_results/cluster_counts/200909/", paste0(paste("CD8_SUBSETS", Sys.Date(), "cluster_counts",sep = "_"), ".csv")))
#read cluster frame back in
cluster.frame <- read.csv("./data_results/cluster_counts/200909/CD8_SUBSETS_2020-09-14_cluster_counts.csv", row.names = 1)
#get frequencies on a per cluster basis
cluster.frame.freq <- as.data.frame(prop.table(as.matrix(cluster.frame), 1) * 100)
cluster.frame.freq$sample <- rownames(cluster.frame.freq)
cluster.dat <- merge(cluster.frame.freq, mdat, by = "sample")
#set color scheme
fourcolors<-c('#8744f6','#414141','#acfafa','#f29737')
#make dataframe to plot frequencies of each subset within each cluster
meltclusterdat<-melt(cluster.dat)
colnames(meltclusterdat)<- c('sample','phenotype.subset','cluster','value')

clusterlabels<-paste0('cluster.',1:18)
newclusterlabels<-paste0('Cluster ',1:18)
meltclusterdat$newcluster <- factor(meltclusterdat$cluster, levels = clusterlabels,labels = newclusterlabels)
meltclusterdat$subset <- factor(meltclusterdat$phenotype.subset, levels = c( "DP","CD49a","CD103","DN"))
#plot frequencies of each subset within each cluster
clusterboxplot<-ggplot(meltclusterdat, aes(x = subset, y = value, fill = subset)) + 
  geom_boxplot(outlier.shape = NA)+scale_fill_manual(values=fourcolors) +
  facet_wrap(~newcluster, scales = "free_y",ncol=3) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=9,color = 'black',face='bold'),
        axis.ticks.x = element_blank())+
  theme(strip.text.x = element_text(size = 11, color = "black", face = "bold"))+
  theme(axis.title.y=element_blank())+
  theme(legend.position = "none")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=3)+
  scale_color_manual(values=c('black','black','black','black'))+
  theme(strip.background = element_blank())
clusterboxplot
ggsave(file.path("./data_results/fsoms/", paste0("CD8_cluster_boxplots.dotplot.fromnate.mike_heightwidth_with_labels", Sys.Date(), ".png")), plot = clusterboxplot,dpi=700,height=6,width=4)
#save csv of data used 
write.csv(meltclusterdat,paste0('meltclusterdat.',Sys.Date(),'.csv'))
##FIN