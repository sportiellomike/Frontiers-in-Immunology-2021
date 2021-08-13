##using 'datMERGE_UMAP', generate a FlowSOM clustering using select phenotypic markers
##various plots of clustering results
##includes components for Figure 9

##load libraries
library(FlowSOM)
library(ggplot2)

#source user-created functions
source("./workflow_analysis_Rcode/functions_depot/functions_NGL_clustering.R")
source("./workflow_analysis_Rcode/functions_depot/functions_NGL_plots.R")

if(!dir.exists("./data_results/20200806_D21/figure_9_components")){
  dir.create("./data_results/20200806_D21/figure_9_components", recursive = T)
}
##read in merged data; this .rds is a list
dat.mdat <- readRDS("./data_results/20200806_D21/20200806_D21_datMERGE_UMAP.rds")
class(dat.mdat)
sapply(dat.mdat, is.data.frame)#`$dat.mdat` = flow data (compensated, transformed), merged with meta-data
dat <- dat.mdat$dat.mdat
head(dat)
##for FlowSOM clustering, use the same dimensions as were used to generate the UMAP embedding
dims.of.interest <- dat.mdat$umap$dimensions;print(dims.of.interest)
##create FlowSOM object (list) and generate SOMs/cluster
args(fsom_scratch_build)
fsom <- fsom_scratch_build(dat)#create FlowSOM object; input data stored in '$data'
fsom$mdat <- dat.mdat$mdat#store the meta-data
fsom$elgcl <- dat.mdat$elgcl#store the .fcs bi-exponential transformation parameters; can be used to 'map' additional files if ever needed
fsom$seed.val <- dat.mdat$umap$seed.val#store the previously determined random-seed vaule
fsom$dims.used <- dat.mdat$umap$dimensions#store the dimensions used for clustering
set.seed(fsom$seed.val)
fsom <- BuildSOM(fsom, colsToUse = sort(match(fsom$dims.used, colnames(fsom$data))), xdim = 12, ydim = 12)#create the self-organizing maps
fsom <- BuildMST(fsom)#create a minimal-spanning tree 

PlotStars(fsom, range = "one")#visualization

##set the number of clusters; decided by the user; iterative process to determine optimum number
##SOMs are organized into clusters (metaclusters)
fsom$nbclust <- 20#number of clusters
cc.dir <- "./data_results/20200806_D21/fsoms/Clustering Consensus/"
if(!dir.exists(cc.dir)){
  dir.create(cc.dir, recursive = T)
}
fsom$clustering <- ConsensusClusterPlus::ConsensusClusterPlus(t(fsom$map$codes), 
                                                              maxK = fsom$nbclust, 
                                                              reps = 100, 
                                                              pItem = 1, 
                                                              pFeature = 1, 
                                                              title = cc.dir, 
                                                              plot = "png",
                                                              clusterAlg = "hc", 
                                                              innerLinkage = "average", 
                                                              finalLinkage = "average",
                                                              distance = "euclidean", 
                                                              seed = 4242)
fsom[["metaclustering"]] <- as.factor(fsom$clustering[[fsom$nbclust]]$consensusClass)#create a cluster column (factor)
PlotStars(fsom, range = "all", backgroundValues = fsom$metaclustering, legend = FALSE)

##heatmap (median expression per marker, per cluster)
args(heatmap.cluster)#fsom object as input to generate a heatmap; performs range scaling (lowest 1% set to 0, highest 1% set to 1)
heatmap.cluster(fsom, type = "codes", range.scale = FALSE)
heatmap.cluster(fsom, type = "codes", color = "viridis")
fsom$cluster.heatmap <- heatmap.cluster(fsom, type = "codes", color = "viridis")#store the heatmap
ggsave("./data_results/20200806_D21/fsoms/fsom_clusters_heat.png", plot = fsom$cluster.heatmap, width = 6, height = 8, dpi = 600)
ggsave("./data_results/20200806_D21/figure_9_components/figure_9_C_fsom_clusters_heat.png", 
       plot = fsom$cluster.heatmap, width = 6, height = 8, dpi = 600)
##
args(clusters.tabulated.list)#fsom object as input; generates tabulated cluster counts for each sample;merges with 'fsom$mdat'
fsom$clusters.tab <- clusters.tabulated.list(fsom)
print(lapply(fsom$clusters.tab, head))#list that contains counts, proportion (100%), and per million
##save .RDS
saveRDS(fsom, file.path("./data_results/20200806_D21/", paste0(paste("fsom", Sys.Date(), sep = "_"), ".rds")))
####
##restart here if clearing out R
fsom <- readRDS("./data_results/20200806_D21/fsom_2021-01-21.rds")
dat.mdat <- readRDS("./data_results/20200806_D21/20200806_D21_datMERGE_UMAP.rds")
dat <- dat.mdat$dat.mdat
####
args(umap.fsom.codes.plot)
umap.fsom.codes.plot(fsom, view = "marker")
umap.fsom.codes.plot(fsom, view = "cluster")

fsom$umap.codes.plots <- list(markers = umap.fsom.codes.plot(fsom, view = "marker"),
                              clusters = umap.fsom.codes.plot(fsom, view = "cluster")
)
fsom$umap.codes.plots#this list can be output as saved plots
##plot UMAP embedding (all cells), colored by FlowSOM cluster
dat.sub <- data.frame(dat.mdat$dat.mdat[, c("umap.1", "umap.2")], 
                  cluster = fsom$metaclustering[fsom$map$mapping[, 1]])[sample(1:length(fsom$map$mapping[, 1]), 100000), ]
head(dat.sub)
library(viridis)

cluster.labs <- c(paste0('Cluster ',1:20))
names(cluster.labs) <- c(1:20)
umap.cells.clusters <- ggplot(data = dat.sub, aes(umap.1, umap.2, fill = cluster)) +
  geom_hex(data = dat.sub[, c("umap.1", "umap.2")], fill = "gray", bins = 100, alpha = 0.5) +
  geom_hex(bins = 100) +
  theme_void() +
  facet_wrap(~cluster,
             ncol=4#,
            # labeller = labeller(cluster=cluster.labs)
             ) +
  guides(fill = FALSE)+theme(strip.text = element_text(size=12,face = 'bold')) +  scale_fill_viridis(discrete = T,option = 'plasma')
umap.cells.clusters#takes ~20 seconds to generate a viewable plot
ggsave("./data_results/20200806_D21/fsoms/UMAP_cells_FlowSOMclusters.height6.justnumbers.png", plot = umap.cells.clusters, width = 3.75, height = 6, dpi = 600)#once satisfied, write to a meaningful file name
##QC plot; biological replicates (5) per organ (CD44+) showing cluster frequency/proportion
attach(fsom$clusters.tab$count)
fsom$clusters.tab$count$organ.alias.replicate <- as.factor(paste0(organ.alias, replicate))
detach(fsom$clusters.tab$count)
args(prop.table.intersect)
dat.sub <- prop.table.intersect(input.data = fsom$clusters.tab$count,
                                input.data.factor = "organ.alias.replicate",
                                intersect.cols = c("organ.alias.replicate", "activation.state"),
                                intersect.val = "CD44+")
str(dat.sub)
dat.sub$organ <- as.factor(sub("[0-9]", "", dat.sub$sample))
dat.sub.melt <- reshape2::melt(dat.sub, id.vars = c("sample", "organ"), variable.name = "cluster", value.name = "proportion")
str(dat.sub.melt)


organ.labs <- c("Blood", "BAL",'Lung','MLN','Spleen')
names(organ.labs) <- c("B", "BAL",'L','MLN','S')

clusters.by.organ.replicate.plot <-ggplot(dat.sub.melt, aes(x = cluster, y = proportion, color = cluster)) + 
  geom_boxplot() +
  geom_point() + 
  viridis::scale_color_viridis(discrete = T, option = 'plasma') +
  facet_wrap(~organ, scales = "free_y",
             ncol = 1,
             labeller = labeller(organ=organ.labs)) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(face = 'bold',color='black'),
        axis.text.y = element_text(face = 'bold',color='black'),
        axis.ticks.x=element_blank()) +
  theme(strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.background = element_blank()) +
  theme(legend.text = element_text(face='bold', size = 11),
        legend.title = element_text(face = 'bold', size = 11))
clusters.by.organ.replicate.plot#save plot if needed
ggsave("./data_results/20200806_D21/fsoms/clusters.by.organ.replicate.plot.1col.png", plot = clusters.by.organ.replicate.plot, width = 3.75, height = 6, dpi = 600)

##Generate a UMAP embedding using cluster proportion within organ.replicate (CD44+)
set.seed(fsom$seed.val)
umap.organ.replicate <- uwot::umap(X = dat.sub[, sapply(dat.sub, is.numeric)])
dat.tmp <- cbind(dat.sub, setNames(data.frame(umap.organ.replicate), nm = paste("umap", 1:2, sep = ".")))
umap.organ.replicate.plot <- ggplot(dat.tmp, aes(x = umap.1, y = umap.2, color = organ)) + geom_point(size = 4, alpha = 0.5)
umap.organ.replicate.plot#save plot if needed
##plot - per organ - of each cluster, showing % of the 4 phenotypes
##refactor phenotype for plot order
attach(fsom$clusters.tab$count)
fsom$clusters.tab$count$organ.replicate.phenotype <- as.factor(paste(organ, replicate, phenotype, sep = "."))
detach(fsom$clusters.tab$count)
args(prop.table.intersect)
dat.sub <- prop.table.intersect(input.data = fsom$clusters.tab$count,
                                input.data.factor = "organ.replicate.phenotype",
                                intersect.cols = c("organ.replicate.phenotype", "activation.state"),
                                intersect.val = "CD44+")
str(dat.sub)
dat.sub <- cbind(dat.sub, 
                 lapply(tidyr::separate(data.frame(sample = dat.sub$sample), 
                                        sample, 
                                        into = c("organ", "replicate", "phenotype"), sep = "\\."), as.factor))
str(dat.sub)
dat.sub$phenotype <- factor(dat.sub$phenotype, levels = c('DP', 'CD49a', 'CD103', 'DN'))
unique(dat.sub$organ)
##split based organ
dat.sub.split <- lapply(split(dat.sub, dat.sub$organ), droplevels)
str(dat.sub.split$`IV_CD45- Lung`)
##output the splits as .csv; can calculate statistics using external programs if desired (Excel, Prism)
sapply(names(dat.sub.split), function(i){
  f.name <- file.path("./data_results/20200806_D21/fsoms/",
                      paste0(paste("fsom_clusters_proportion_phenotype", i, Sys.Date(), sep = "_"), ".csv")
  )
  write.csv(dat.sub.split[[i]], file = f.name, row.names = F)
})
##melt the splits
dat.sub.split.melt <- lapply(dat.sub.split, function(i){
  reshape2::melt(i, id.vars = c("sample", "organ", "replicate", "phenotype"), 
                 variable.name = "cluster", 
                 value.name = "proportion")
})
str(dat.sub.split.melt$`IV_CD45- Lung`)
#create a color-palette; one color per phenotype
four.colors.phenotype <- setNames(c('#8744f6','#414141','#acfafa','#f29737'), nm = levels(dat.sub$phenotype))
#make plots for all split types
p.list <- sapply(names(dat.sub.split.melt), function(i){
  ggplot(dat.sub.split.melt[[i]], aes(x = phenotype, y = proportion, fill = phenotype)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 2) +
    #scale_color_manual(values = four.colors.phenotype) +
    scale_fill_manual(values = four.colors.phenotype) +
    facet_wrap(~cluster, 
               scales = "free_y", 
               ncol = 4,
               labeller = labeller(cluster = setNames(paste("Cluster", 1:length(levels(dat.sub.split.melt[[i]]$cluster))), 
                                                      nm = 1:length(levels(dat.sub.split.melt[[i]]$cluster))))
    ) +
    theme(
      axis.title = element_blank()
      , axis.text.x = element_blank()
      , axis.text.y = element_text(size=11,color = 'black',face='bold')
      , axis.ticks.x = element_blank()
      , strip.text.x = element_text(size = 11, color = "black", face = "bold")
      , strip.background = element_blank()
      , legend.position = "none"
      , axis.title.y = element_blank()
    ) +
    labs(title = i)
}, simplify = FALSE, USE.NAMES = T)
p.list$`IV_CD45- Lung`
ggsave("./data_results/20200806_D21/figure_9_components/figure_9_D_cluster_boxplots_LungTissue.png", 
       plot = p.list$`IV_CD45- Lung` + theme(title = element_blank()),
       dpi = 600, width = 5, height =5)
##Generate statistical comparisons
library(rstatix)
comparisons<- list(c('DP','DN'),c('DP','CD49a'),c('DP','CD103'),c('CD49a','DN'),c('CD49a','CD103'),c('CD103','DN'))
dat.stats <- dat.sub.split$`IV_CD45- Lung`
dat.stats[,1:20] <- 1+dat.stats[1:20]
dat.stats[,1:20] <- log(dat.stats[1:20], 2)
colnames(dat.stats)[sapply(dat.stats, is.numeric)] <- paste("cluster", which(sapply(dat.stats, is.numeric)), sep = ".")
str(dat.stats)
head(dat.stats)

stats.results_ttest <- setNames(vector("list", length = length(grep("cluster", colnames(dat.stats)))), nm = grep("cluster", colnames(dat.stats), value = T))
str(stats.results_ttest)

for(i in grep("cluster", colnames(dat.stats), value = T)){
  stats.tmp <- dat.stats[, c(i, "phenotype", "replicate")]
  ttest <- pairwise_t_test(data = stats.tmp,
                           formula = as.formula(paste0(i," ~ phenotype")),
                           paired = TRUE,
                           comparisons = comparisons
  )
  stats.results_ttest [[i]] <- ttest
}
stats.table <- dplyr::bind_rows(stats.results_ttest)
write.csv(stats.table, "./data_results/20200806_D21/figure_9_components/figure_9_D_STATS_cluster_boxplots_LungTissue.allsixcomparisons.csv", row.names = F)
write.csv(stats.table, "./data_results/20200806_D21/figure_9_components/figure_9_D_STATS_cluster_boxplots_LungTissue.csv", row.names = F)
##
stats.results_anova <- setNames(vector("list", length = length(grep("cluster", colnames(dat.stats)))), nm = grep("cluster", colnames(dat.stats), value = T))
str(stats.results_anova)
for(i in grep("cluster", colnames(dat.stats), value = T)){
  stats.tmp <- dat.stats[, c(i, "phenotype", "replicate")]
  res.aov <- anova_test(data = stats.tmp, 
                        dv = str2lang(i), 
                        wid = replicate, 
                        within = phenotype)
  stats.results_anova[[i]] <- res.aov
}
res.aov <- anova_test(data = dat.s.c, 
                               dv = cluster.1, 
                               wid = replicate, 
                               within = phenotype)
test<-sapply(stats.results_anova, function(i) get_anova_table(i))
##more comparisons
args(fsom.stats)
stats.results <- sapply(c("aov", "dunnett"), function(i){
  fsom.stats(input.count.data = subset(fsom$clusters.tab$count, organ == "IV_CD45- Lung"),
             combine.for.new.factor = c("organ", "phenotype", "replicate"),
             intersect.col = "activation.state",
             intersect.val = "CD44+",
             new.factor.order = list(c("phenotype"), c('DP', 'CD49a', 'CD103', 'DN')),
             split.name = "organ",
             grouping.name = "phenotype",
             stats.test = i)
})

sapply(names(stats.results), function(i){
  f.name <- file.path("./data_results/20200806_D21/figure_9_components/", 
                      paste0(paste("figure_9_D_STATS", i, "cluster_boxplots", sep = "_"), 
                             ".csv"))
  write.csv(stats.results[[i]][length(stats.results[[i]])], f.name, row.names = T)
})
##Generate a UMAP embedding using cluster proportion withing organ.replicate.phenotype
set.seed(fsom$seed.val)
umap.organ.replicate.phenotype <- uwot::umap(X = dat.sub[, sapply(dat.sub, is.numeric)])
dat.tmp <- cbind(dat.sub, setNames(data.frame(umap.organ.replicate.phenotype), nm = paste("umap", 1:2, sep = ".")))
umap.organ.replicate.phenotype.plot <- ggplot(dat.tmp, aes(x = umap.1, y = umap.2, color = factor(paste(organ, phenotype, sep = ".")))) + geom_point(size = 4, alpha = 0.5)
umap.organ.replicate.phenotype.plot
##
##
