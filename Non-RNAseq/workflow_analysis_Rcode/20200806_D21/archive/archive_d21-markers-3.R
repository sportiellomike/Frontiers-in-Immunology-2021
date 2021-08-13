#load libraries
library(flowCore)
library(FlowSOM)
library(ConsensusClusterPlus)
library(dplyr)
library(ggplot2)
library(viridis)
library(matrixStats)
library(RColorBrewer)
library(pheatmap)
library(shiny)
library(reshape2)

#load the files, compensate, transform, and format
input.flowset <- read.flowSet(list.files("./", pattern = ".fcs", full.names = T), transformation =F)

comp.matrix.edit <- as.matrix(read.csv("matrix.12.9.20.csv", row.names = 1, check.names = F))
all(sapply(colnames(comp.matrix.edit), function(i) unlist(strsplit(i, " :: "))[[1]]) == colnames(input.flowset[[1]]@description$`$SPILLOVER`))
colnames(comp.matrix.edit) <- colnames(input.flowset[[1]]@description$`$SPILLOVER`)
rownames(comp.matrix.edit) <- colnames(comp.matrix.edit)
input.flowset <- fsApply(input.flowset, function(i) compensate(i, comp.matrix.edit))
elgcl <- elgcl.pars(input.flowset)
input.flowset <- fset.transform(input.flowset, elgcl)
input.flowset <- markernames.flow(input.flowset)
#select which dims to cluster on
dims <- c(
  #'CD103',
  #'CD49a',
  'KLRG1',
  'CX3CR1',
  'CCR7',
  'CD69',
  'CD62L'
)
fsom <- ReadInput(input.flowset, scale = T)
fsom$elgcl <- elgcl
fsom$markers <- colnames(input.flowset)
#set seed and format
fsom$seed <- 585414
fsom$dims.used <- grep(paste0(dims, collapse = "|"), fsom$markers, value = T)
#saving the RDS so it can be quickloaded using readRDS later. 
#Of note, the filename loaded must be changed for your use
#saveRDS(fsom, file.path("./data_results/fsoms/", paste0("CD8_fsom_", Sys.Date(),'.rds')))
fsom<-readRDS('data_results/fsoms/CD8_fsom_2020-12-16.rds')
set.seed(fsom$seed)
fsom <- BuildSOM(fsom, colsToUse = sort(match(fsom$dims.used, fsom$markers)), xdim = 12, ydim = 12)
fsom <- BuildMST(fsom)

PlotStars(fsom)
PlotStars(fsom, view = 'grid')

#set the number of clusters for flowSOM to use
fsom$nbclust <- 20
#dir.create("./data_results/fsoms/temp/Clustering Consensus", recursive = T)
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
fsom[["metaclustering"]] <- as.factor(fsom$clustering[[fsom$nbclust]]$consensusClass)
PlotStars(fsom, backgroundValues = fsom$metaclustering, legend = FALSE)
PlotNumbers(fsom, backgroundValues = fsom$metaclustering)




#this will let you explore the data more thoroughly. 
#note you MUST comment out one of the markers as the function will only take 6 at a time. 
#Once launched, you can switch out whatever 6 you want.
fsom.explorer(fsom, #these determine the windows in the explorer.
             # "KLRG1", 
              'CD103',
              "CX3CR1",
              'CD49a',
              'CCR7',
              'CD69',
              'CD62L',
              20000)



fsom$dims.used <- c("KLRG1", 
                    "CX3CR1",
                    'CCR7',
                    'CD69',
                    'CD62L',
                    'CD103',
                    'CD49a')
fsom$dims.used <- grep(paste0(dims, collapse = "|"), fsom$markers, value = T)

#take a look at the heatmap
#fsom.heatmap(fsom, dat.type = "codes", fsom$dims.used)

dat.umap.samples <- readRDS("./data_results/umap/CD8_umap_2020-12-15.rds")
dat <- data.frame(dat.umap.samples,
                  node = fsom$map$mapping[, 1],
                  cluster = fsom$metaclustering[fsom$map$mapping[,1]]
)
#subset on samples that are CD44+
dat<-subset(dat,activation=='cd44')

#plot, colored by cluster. Note that this will take a few minutes.
ggplot(dat, aes(umap.1, umap.2, color = cluster)) + geom_point(alpha=.1) #+ scale_color_viridis(discrete=T)
#plot, faceted by cluster. Note that this will take a few minutes.
umapclustergg<-ggplot(dat, aes(umap.1, umap.2, color = cluster)) + 
  geom_point(alpha=.03) +
  facet_wrap(~cluster) 
umapclustergg
#save plot
ggsave(file.path("./data_results/umap/", paste0("CD8_umap_cluster20_hex_", Sys.Date(), ".png")), plot = umapclustergg,dpi=600)

#prepare plot of umap of fsom nodes
nodeumap<-setNames(data.frame(umap(fsom$map$codes)), nm=c('umap.1','umap.2'))
datnode<-data.frame(fsom$map$codes, nodeumap,cluster=fsom$metaclustering)
colnames(datnode)
#plot umap of fsom nodes
nodeumapclustergg<-ggplot(datnode, aes(umap.1, umap.2, color = cluster,size=4)) + 
  geom_point(alpha=.5) 
nodeumapclustergg

#prepare umap of nodes colored according to expression level of faceted markers
datnodemelt<-melt(datnode,id.vars = c('umap.1','umap.2','cluster'))
nodeumapclusterggmelted<-ggplot(datnodemelt, aes(umap.1, umap.2, color = value)) + 
  geom_point(alpha=.5)+
  facet_wrap(~variable)+scale_color_viridis(option='plasma')+
  theme_dark()+
  theme(panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())
nodeumapclusterggmelted

#prepare cluster frequency tables
table(dat$cluster[grep(unique(as.character(dat$sample))[1], dat$sample)])
clusters.tabulated <- lapply(setNames(nm = unique(as.character(dat$sample))), function(i) table(dat$cluster[grep(i, dat$sample)]))
cluster.frame <- as.data.frame((bind_rows(clusters.tabulated)))
colnames(cluster.frame) <- paste("cluster", 1:ncol(cluster.frame), sep = ".")
cluster.frame.freq <- as.data.frame(prop.table(as.matrix(cluster.frame), 1) * 100)
cluster.frame.freq$sample <- names(clusters.tabulated)
cluster.frame.1E6 <- as.data.frame(t(apply(cluster.frame, 1, function(x) x/sum(x)*1E6)))
meltclusterframefreq<-melt(cluster.frame.freq,id.vars = 'sample')
colnames(meltclusterframefreq)
colnames(meltclusterframefreq)<-c('sample','cluster','value')
colnames(meltclusterframefreq)
#rename samples
meltclusterframefreq$sample <- factor(meltclusterframefreq$sample, levels = c(
  'BAL1',
  'BAL2',
  'BAL3',
  'BAL4',
  'BAL5',
  'B1',
  'B2',
  'B3',
  'B4',
  'B5',
  'L1',
  'L2',
  'L3',
  'L4',
  'L5',
  'MLN1',
  'MLN2',
  'MLN3',
  'MLN4',
  'MLN5',
  'S1',
  'S2',
  'S3',
  'S4',
  'S5'
),
labels = c(
  'BAL 1',
  'BAL 2',
  'BAL 3',
  'BAL 4',
  'BAL 5',
  'Blood 1',
  'Blood 2',
  'Blood 3',
  'Blood 4',
  'Blood 5',
  'Lung 1',
  'Lung 2',
  'Lung 3',
  'Lung 4',
  'Lung 5',
  'MLN 1',
  'MLN 2',
  'MLN 3',
  'MLN 4',
  'MLN 5',
  'Spleen 1',
  'Spleen 2',
  'Spleen 3',
  'Spleen 4',
  'Spleen 5'
))

#rename clusters
clusterlabels<-paste0('cluster.',1:20)
newclusterlabels<-paste0('Cluster ',1:20)
meltclusterframefreq$cluster <- factor(meltclusterframefreq$cluster, levels = clusterlabels,labels = newclusterlabels)
#rename column names
colnames(meltclusterframefreq)<-c('sample','Cluster','value')
#make a cluster frequency plot that shows frequency of clusters faceted by sample colored by cluster
#this demonstrates that each biological replicate looks very similar and serves as a quality check. 
#You'll note that similar frequencies of the same clusters appear in the same organs. 
#Said another way: All the BAL samples look like the other BAL samples, 
#all the Blood samples look like the other Blood samples, etc.
meltclusterframefreqplot<-ggplot(meltclusterframefreq, 
                                 aes(x = Cluster, 
                                     y = value, 
                                     color = Cluster
                                 )) + 
  geom_point() + 
  scale_color_viridis(discrete=T,option = 'plasma')+
  facet_wrap(~sample, scales = "free_y") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(strip.text.x = element_text(size = 11, color = "black", face = "bold"))+
  theme(legend.text = element_text(face='bold',size=11),legend.title = element_text(face = 'bold',size=11))
meltclusterframefreqplot
#save plot
ggsave(file.path("./data_results/fsoms/", 
                 paste0("CD8_umap_meltclusterframefreq_", Sys.Date(), ".png")), 
       plot = meltclusterframefreqplot,dpi=700)

#set dims for heatmap
dims <- c("KLRG1", 
          'CD103',
          "CX3CR1",
          'CD49a',
          'CCR7',
          'CD69',
          'CD62L')
#prepare data used to read into heatmap plots
fsom$dims.used <- grep(paste0(dims, collapse = "|"), fsom$markers, value = T)
cluster.medians <- as.data.frame(fsom.cluster.medians(fsom, dat.type = 'map', ranged = TRUE))
annotation_row <- cluster.medians["cluster"]
annotation_colors <- list(cluster = setNames(colorRampPalette(viridis(9))(nlevels(fsom$metaclustering)), levels(annotation_row$cluster)))

col.index <- fsom$markers[fsom$map$colsUsed]
cluster_rows <- hclust(dist(cluster.medians[, -grep("cluster", colnames(cluster.medians))][, col.index], method = "euclidean"), method = "average")
clustering_table <- as.numeric(table(fsom$metaclustering[fsom$map$mapping[, 1]]))
clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)

cluster.medians$clusternumber<-paste0('Cluster ',1:20)
labels_row <- paste0(cluster.medians$clusternumber#, " (", clustering_prop, "%)"#, " (", clustering_table, ")"
)

#plot heatmaps
library(tidyverse)
#heatmap
pheatmap(cluster.medians[, dims], color = plasma(50),
         number_color = "black",
         fontsize = 8, fontsize_number = 6, 
         annotation_colors = annotation_colors)
#simplified heatmap
clusterheatmap<-pheatmap(cluster.medians[, dims], color = plasma(50),
                         cluster_rows = cluster_rows, labels_row =labels_row,
                         display_numbers = F,
                         fontsize = 11, fontsize_number = 6)
#save heatmap
ggsave(file.path("./data_results/fsoms/", 
                 paste0("CD8_heatmap_", Sys.Date(), ".png")), 
       plot = clusterheatmap,dpi=700)

#create column for mouse number
mouse<-paste0('mouse ',rep(1:5,5))
cluster.frame.freq$mouse<-mouse
#create column for clustersubset
clustersubset<-paste0(dat$subsetindex,'.cluster.',dat$cluster)
dat$clustersubset<-clustersubset

#prepare data for cluster boxplots that show what frequency of each subset is in each cluster
#start with dp subset, will do the other three subsets directly after
dat2<-subset(dat,subsetindex=='DP')
clusters.tabulated <- lapply(setNames(nm = unique(as.character(dat2$sample))), function(i) table(dat2$cluster[grep(i, dat2$sample)]))
cluster.frame <- as.data.frame((bind_rows(clusters.tabulated)))
colnames(cluster.frame) <- paste("cluster", 1:ncol(cluster.frame), sep = ".")
cluster.frame.freq <- as.data.frame(prop.table(as.matrix(cluster.frame), 1) * 100)
cluster.frame.freq$sample <- names(clusters.tabulated)
cluster.frame.1E6 <- as.data.frame(t(apply(cluster.frame, 1, function(x) x/sum(x)*1E6)))
meltclusterframefreq<-melt(cluster.frame.freq,id.vars = 'sample')
colnames(meltclusterframefreq)
colnames(meltclusterframefreq)<-c('sample','cluster','value')
colnames(meltclusterframefreq)
meltclusterframefreq$sample <- factor(meltclusterframefreq$sample, levels = c(
  'BAL1',
  'BAL2',
  'BAL3',
  'BAL4',
  'BAL5',
  'B1',
  'B2',
  'B3',
  'B4',
  'B5',
  'L1',
  'L2',
  'L3',
  'L4',
  'L5',
  'MLN1',
  'MLN2',
  'MLN3',
  'MLN4',
  'MLN5',
  'S1',
  'S2',
  'S3',
  'S4',
  'S5'
),
labels = c(
  'BAL 1',
  'BAL 2',
  'BAL 3',
  'BAL 4',
  'BAL 5',
  'Blood 1',
  'Blood 2',
  'Blood 3',
  'Blood 4',
  'Blood 5',
  'Lung 1',
  'Lung 2',
  'Lung 3',
  'Lung 4',
  'Lung 5',
  'MLN 1',
  'MLN 2',
  'MLN 3',
  'MLN 4',
  'MLN 5',
  'Spleen 1',
  'Spleen 2',
  'Spleen 3',
  'Spleen 4',
  'Spleen 5'
))
meltclusterframefreq$subset<-'DP'
meltdp<-meltclusterframefreq
#Okay now we're on the CD49a single positive subset
dat2<-subset(dat,subsetindex=='CD49a')
clusters.tabulated <- lapply(setNames(nm = unique(as.character(dat2$sample))), function(i) table(dat2$cluster[grep(i, dat2$sample)]))
cluster.frame <- as.data.frame((bind_rows(clusters.tabulated)))
colnames(cluster.frame) <- paste("cluster", 1:ncol(cluster.frame), sep = ".")
cluster.frame.freq <- as.data.frame(prop.table(as.matrix(cluster.frame), 1) * 100)
cluster.frame.freq$sample <- names(clusters.tabulated)
cluster.frame.1E6 <- as.data.frame(t(apply(cluster.frame, 1, function(x) x/sum(x)*1E6)))
meltclusterframefreq<-melt(cluster.frame.freq,id.vars = 'sample')
colnames(meltclusterframefreq)
colnames(meltclusterframefreq)<-c('sample','cluster','value')
colnames(meltclusterframefreq)
meltclusterframefreq$sample <- factor(meltclusterframefreq$sample, levels = c(
  'BAL1',
  'BAL2',
  'BAL3',
  'BAL4',
  'BAL5',
  'B1',
  'B2',
  'B3',
  'B4',
  'B5',
  'L1',
  'L2',
  'L3',
  'L4',
  'L5',
  'MLN1',
  'MLN2',
  'MLN3',
  'MLN4',
  'MLN5',
  'S1',
  'S2',
  'S3',
  'S4',
  'S5'
),
labels = c(
  'BAL 1',
  'BAL 2',
  'BAL 3',
  'BAL 4',
  'BAL 5',
  'Blood 1',
  'Blood 2',
  'Blood 3',
  'Blood 4',
  'Blood 5',
  'Lung 1',
  'Lung 2',
  'Lung 3',
  'Lung 4',
  'Lung 5',
  'MLN 1',
  'MLN 2',
  'MLN 3',
  'MLN 4',
  'MLN 5',
  'Spleen 1',
  'Spleen 2',
  'Spleen 3',
  'Spleen 4',
  'Spleen 5'
))
meltclusterframefreq$subset<-'CD49a'
meltcd49a<-meltclusterframefreq
#now we are onto DN subset
dat2<-subset(dat,subsetindex=='DN')
clusters.tabulated <- lapply(setNames(nm = unique(as.character(dat2$sample))), function(i) table(dat2$cluster[grep(i, dat2$sample)]))
cluster.frame <- as.data.frame((bind_rows(clusters.tabulated)))
colnames(cluster.frame) <- paste("cluster", 1:ncol(cluster.frame), sep = ".")
cluster.frame.freq <- as.data.frame(prop.table(as.matrix(cluster.frame), 1) * 100)
cluster.frame.freq$sample <- names(clusters.tabulated)
cluster.frame.1E6 <- as.data.frame(t(apply(cluster.frame, 1, function(x) x/sum(x)*1E6)))
meltclusterframefreq<-melt(cluster.frame.freq,id.vars = 'sample')
colnames(meltclusterframefreq)
colnames(meltclusterframefreq)<-c('sample','cluster','value')
colnames(meltclusterframefreq)
meltclusterframefreq$sample <- factor(meltclusterframefreq$sample, levels = c(
  'BAL1',
  'BAL2',
  'BAL3',
  'BAL4',
  'BAL5',
  'B1',
  'B2',
  'B3',
  'B4',
  'B5',
  'L1',
  'L2',
  'L3',
  'L4',
  'L5',
  'MLN1',
  'MLN2',
  'MLN3',
  'MLN4',
  'MLN5',
  'S1',
  'S2',
  'S3',
  'S4',
  'S5'
),
labels = c(
  'BAL 1',
  'BAL 2',
  'BAL 3',
  'BAL 4',
  'BAL 5',
  'Blood 1',
  'Blood 2',
  'Blood 3',
  'Blood 4',
  'Blood 5',
  'Lung 1',
  'Lung 2',
  'Lung 3',
  'Lung 4',
  'Lung 5',
  'MLN 1',
  'MLN 2',
  'MLN 3',
  'MLN 4',
  'MLN 5',
  'Spleen 1',
  'Spleen 2',
  'Spleen 3',
  'Spleen 4',
  'Spleen 5'
))
meltclusterframefreq$subset<-'DN'
meltdn<-meltclusterframefreq
#finally, we're doing the CD103 single positive subset
dat2<-subset(dat,subsetindex=='CD103')
clusters.tabulated <- lapply(setNames(nm = unique(as.character(dat2$sample))), function(i) table(dat2$cluster[grep(i, dat2$sample)]))
cluster.frame <- as.data.frame((bind_rows(clusters.tabulated)))
colnames(cluster.frame) <- paste("cluster", 1:ncol(cluster.frame), sep = ".")
cluster.frame.freq <- as.data.frame(prop.table(as.matrix(cluster.frame), 1) * 100)
cluster.frame.freq$sample <- names(clusters.tabulated)
cluster.frame.1E6 <- as.data.frame(t(apply(cluster.frame, 1, function(x) x/sum(x)*1E6)))
meltclusterframefreq<-melt(cluster.frame.freq,id.vars = 'sample')
colnames(meltclusterframefreq)
colnames(meltclusterframefreq)<-c('sample','cluster','value')
colnames(meltclusterframefreq)
meltclusterframefreq$sample <- factor(meltclusterframefreq$sample, levels = c(
  'BAL1',
  'BAL2',
  'BAL3',
  'BAL4',
  'BAL5',
  'B1',
  'B2',
  'B3',
  'B4',
  'B5',
  'L1',
  'L2',
  'L3',
  'L4',
  'L5',
  'MLN1',
  'MLN2',
  'MLN3',
  'MLN4',
  'MLN5',
  'S1',
  'S2',
  'S3',
  'S4',
  'S5'
),
labels = c(
  'BAL 1',
  'BAL 2',
  'BAL 3',
  'BAL 4',
  'BAL 5',
  'Blood 1',
  'Blood 2',
  'Blood 3',
  'Blood 4',
  'Blood 5',
  'Lung 1',
  'Lung 2',
  'Lung 3',
  'Lung 4',
  'Lung 5',
  'MLN 1',
  'MLN 2',
  'MLN 3',
  'MLN 4',
  'MLN 5',
  'Spleen 1',
  'Spleen 2',
  'Spleen 3',
  'Spleen 4',
  'Spleen 5'
))
meltclusterframefreq$subset<-'CD103'
meltcdcd103<-meltclusterframefreq
#now bind all four of the subsets together into "testmerge" dataframe
testmerge<-rbind(meltdp,meltdn,meltcd49a,meltcdcd103)
#now subset them based on sample
lung1<-subset(testmerge,sample=='Lung 1')
lung2<-subset(testmerge,sample=='Lung 2')
lung3<-subset(testmerge,sample=='Lung 3')
lung4<-subset(testmerge,sample=='Lung 4')
lung5<-subset(testmerge,sample=='Lung 5')
#bind those together by lung
lung<-rbind(lung1,lung2,lung3,lung4,lung5)
#create a new column for subset called "newsubset" as a quality check to ensure the levels are as expected
lung$newsubset <- factor(lung$subset, levels = c('DP','CD49a','CD103','DN'))
#rename cluster labels with a capital "C"
clusterlabels<-paste0('cluster.',1:20)
newclusterlabels<-paste0('Cluster ',1:20)
#create new column for cluster called "newcluster" as a quality check to ensure all levels are as expected
lung$newcluster<-factor(lung$cluster, levels = clusterlabels,labels = newclusterlabels)
#create a palette for use, one color per subset
fourcolors<-c('#8744f6','#414141','#acfafa','#f29737')
#make a plot
clusterboxplot<-ggplot(lung, aes(x = newsubset, y = value,   fill = newsubset)) + 
  geom_boxplot(outlier.shape = NA)+scale_fill_manual(values=fourcolors) +
  facet_wrap(~newcluster, scales = "free_y",ncol=4) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=11,color = 'black',face='bold'),
        axis.ticks.x = element_blank())+
  theme(strip.text.x = element_text(size = 11, color = "black", face = "bold"))+
  theme(axis.title.y=element_blank())+
  theme(legend.position = "none")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2)+
  scale_color_manual(values=c('black','black','black','black'))+
  theme(strip.background = element_blank())
clusterboxplot
#save plot
ggsave(file.path("./data_results/fsoms/", paste0("CD8_cluster_boxplots.dotplot_", Sys.Date(), ".png")), plot = clusterboxplot,dpi=700,height = 5,width=5)
