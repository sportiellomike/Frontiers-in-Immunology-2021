#load libraries
library(flowCore)
library(dplyr)
library(uwot)
library(ggplot2)
library(viridis)
library(matrixStats)
library(svglite)
library(extrafont)
#load fonts for figure creation if necessary
#font_import()
loadfonts(device = "win")
source("./workflow_analysis_Rcode/functions_depot/flow-functions-1-functions_flow.mike.12.17.2020.R")
#load files into flowset
exported.fcs <- list.files("./data_source/20200806_D21/exported_FCS/", pattern = ".fcs", full.names = T)
input.flowset <- read.flowSet(exported.fcs, transformation =F, truncate_max_range = F)
sampsub<-sub('export_Specimen_001_','',sampleNames(input.flowset))
#create indices of which files were exported as CD44+,tetramer+,vasculature associated cells, and tissue associated cells
cd44index<-grep('CD44+',sampsub)
cd44index
tetindex<-grep('tet',sampsub)
tetindex
vascindex<-grep('vasc',sampsub)
vascindex
tissueindex<-grep('tiss',sampsub)
tissueindex
bothindex<-grep('both',sampsub)
bothindex

length(cd44index)
length(tetindex)
length(vascindex)
length(tissueindex)
length(bothindex)

cd44.tissue.index<-intersect(cd44index,tissueindex)
cd44.vasc.index<-intersect(cd44index,vascindex)
tet.tissue.index<-intersect(tetindex,tissueindex)
tet.vasc.index<-intersect(tetindex,vascindex)
cd44.both.index<-intersect(cd44index,bothindex)
tet.both.index<-intersect(tetindex,bothindex)

length(cd44.tissue.index)
length(cd44.vasc.index)
length(tet.tissue.index)
length(tet.vasc.index)
length(cd44.both.index)
length(cd44.both.index)

spaceactivation<-vector(mode = 'character',length = 240)
spaceactivation[cd44.tissue.index]<-'cd44_tissue'
spaceactivation[cd44.vasc.index]<-'cd44_vasc'
spaceactivation[tet.tissue.index]<-'tet_tissue'
spaceactivation[tet.vasc.index]<-'tet_vasc'
spaceactivation[cd44.both.index]<-'cd44_both'
spaceactivation[tet.both.index]<-'tet_both'
spaceactivation
##
length(cd44index)
length(tetindex)
activation<-vector(mode = 'character',length = 240)
activation[cd44index]<-'cd44'
activation[tetindex]<-'tet'
activation

#create indices for what organ the samples come from
bloodindex<-grep('B1|B2|B3|B4|B5',sampsub)
bloodindex
balindex<-grep('BAL',sampsub)
balindex
mlnindex<-grep('MLN',sampsub)
mlnindex
spleenindex<-grep('S1|S2|S3|S4|S5',sampsub)
spleenindex
lungvascindex<-grep('vasc',sampsub)
lungvascindex
lungtissindex<-grep('tiss',sampsub)
lungtissindex

length(bloodindex)
length(balindex)
length(mlnindex)
length(spleenindex)
length(lungvascindex)
length(lungtissindex)



organindex<-vector(mode = 'character',length = 240)
organindex[bloodindex]<-'Blood'
organindex[balindex]<-'BAL'
organindex[mlnindex]<-'MLN'
organindex[spleenindex]<-'Spleen'
organindex[lungvascindex]<-'IV-CD45+ Lung'
organindex[lungtissindex]<-'IV-CD45- Lung'
organindex

#create index for integrin subsets (DP,DN, CD49a single positive, CD103 single positive)
cd49aindex<-grep('Q1',sampsub)
cd49aindex
cd103index<-grep('Q3',sampsub)
cd103index
dpindex<-grep('Q2',sampsub)
dpindex
dnindex<-grep('Q4',sampsub)
dnindex

length(cd49aindex)
length(cd103index)
length(dpindex)
length(dnindex)



subsetindex<-vector(mode = 'character',length = 240)
subsetindex[cd49aindex]<-'CD49a'
subsetindex[cd103index]<-'CD103'
subsetindex[dpindex]<-'DP'
subsetindex[dnindex]<-'DN'
subsetindex


#compensate and format
comp.matrix.edit <- as.matrix(read.csv("./data_source/20200806_D21/matrix.12.9.20.csv", row.names = 1, check.names = F))
all(sapply(colnames(comp.matrix.edit), function(i) unlist(strsplit(i, " :: "))[[1]]) == colnames(input.flowset[[1]]@description$`$SPILLOVER`))
colnames(comp.matrix.edit) <- colnames(input.flowset[[1]]@description$`$SPILLOVER`)
rownames(comp.matrix.edit) <- colnames(comp.matrix.edit)
input.flowset <- fsApply(input.flowset, function(i) compensate(i, comp.matrix.edit))

#remove the parameters for future clustering that were originally used to export fcs files
dimsusedduringexport <- c('CD8a','TCRB','L_D','CD44','CD45_IV','NP_PA tetramer','CD4')
input.flowset <- input.flowset[, -which(input.flowset[[1]]@parameters@data$desc %in% dimsusedduringexport)]

#transform
elgcl <- elgcl.pars(input.flowset)
names(elgcl@transforms) <- sub("/", "_", names(elgcl@transforms))
elgcl@transforms <- sapply(elgcl@transforms, function(i){ 
  i@output <- sub("/", "_", i@output)
  i@input <- sub("/", "_", i@input)
  i
})
elgcl@transforms <- elgcl@transforms[which(names(elgcl@transforms) %in% grep("FSC|SSC|Time", colnames(input.flowset), invert = T, value = T))]
input.flowset <- fset.transform(input.flowset, elgcl)
input.flowset <- markernames.flow(input.flowset)
dat <- input.flowset %>% fsApply(exprs) %>% as.data.frame()
head(dat)
#create columns in dataframe that correspond to factors and indices created above
dat$spaceactivation<-as.factor(rep(spaceactivation,fsApply(input.flowset,nrow)))
dat$organindex<-as.factor(rep(organindex,fsApply(input.flowset,nrow)))
dat$activation<-as.factor(rep(activation,fsApply(input.flowset,nrow)))
dat$subsetindex<-as.factor(rep(subsetindex,fsApply(input.flowset,nrow)))
str(dat)

#decide which dims to use to make umap. 
#These are the available options (except for the ones removed above that were used before fcs export)
#CD103 and CD49a were not used in clustering, as they are commented out, 
#but one may uncomment them if they would like them used
dims <- c(
  #'CD103',
  #'CD49a',
  'KLRG1',
  'CX3CR1',
  'CCR7',
  'CD69',
  'CD62L'
  )
dims.for.umap <- dims
#The area code of Rochester, NY, where this work was performed is 585. The area code for my hometown of Milwaukee where I was born and raised is 414. 
#This seemed as good a seed as any other. To perfectly reproduce our results, the seed must be identical; 
#however, to very closely approximate our results the seed should have no impact 
set.seed(585414)

#####THE STEPS BELOW IS ABSOLUTELY ESSENTIAL AND MUST BE 'UNCOMMENTED' OUT FOR THE BELOW STEPS TO WORK#####
#the step below is the formation of the UMAP. As it's a very time-consuming process, I commented it out. 
#However, this is absolutely essential for the steps we preformed below. 
#We saved the result as an RDS, then read that in which saves time (look belove for 'saveRDS' and later 'readrds').


#umap.dims <- umap(dat[, dims.for.umap], n_threads = parallel::detectCores()-1, verbose = TRUE)
#colnames(umap.dims) <- paste("umap", seq(ncol(umap.dims)), sep = ".")
#dat.umap.samples <- data.frame(dat,
#                               umap.dims, 
#                               sample = rep(fsApply(input.flowset, function(i) i@description$`TUBE NAME`), fsApply(input.flowset, nrow)),
#                               mouse = rep(fsApply(input.flowset, function(frame) frame@description$mouse), fsApply(input.flowset, nrow)),
#                               organ = rep(fsApply(input.flowset, function(frame) frame@description$organ), fsApply(input.flowset, nrow))
#                               )
#dat.umap.samples$condition <- as.factor(substr(dat.umap.samples$sample,1,1))
#dir.create("./data_results/umap", recursive = T)
#saveRDS(dat.umap.samples, file.path("./data_results/umap", paste0(paste("CD8_umap", Sys.Date(), sep = "_"), ".rds")))

dat.umap.samples <- readRDS("./data_results/umap/CD8_umap_2020-12-15.rds")
#for the next bit, we only want to include samples that are CD44+ and IV-CD45-
dat.umap.samples<-subset(dat.umap.samples, activation=='cd44')
dat.umap.samples<-subset(dat.umap.samples, organindex=='IV-CD45- Lung')

#addding back in CD49a and CD103 so that umaps can show expression, 
#though, to be clear, the umap was already made without taking their expression into account
dims <- c("CD49a", 
          'CD103',
          "CD69",
          'CD62L',
          'CCR7',
          'CX3CR1',
          'KLRG1')
dims.for.umap <- dims
dat.umap <- data.frame(dat.umap.samples[, dims.for.umap],
                       dat.umap.samples[, c("umap.1", "umap.2")])
dat.umap[, dims.for.umap] <- range.scaled(dat.umap[, dims.for.umap])
dat.melt <- reshape2::melt(dat.umap, id = c("umap.1", "umap.2"))

#umap plot with relative expression of each marker
umap.dims.plot2 <- ggplot(dat.melt, aes(umap.1, umap.2, z = value,color=value)) +
  stat_summary_hex(bins=100) +
  theme(text=element_text(family="Arial",face='bold'))+
  theme(panel.border = element_blank())+
  theme(panel.background = element_blank())+
  scale_fill_viridis(option = "plasma") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
  theme(legend.text = element_text(colour="black", size=11, face="bold"))+
  theme(legend.title = element_text(face = "bold",size=11))+
  facet_wrap(~variable,ncol=2)+
  theme(strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.background = element_blank())
umap.dims.plot2

#saving in different file types
ggsave(file.path("./data_results/umap/", paste0("CD8_umap_hex_", Sys.Date(), ".png")), plot = umap.dims.plot2,dpi=600)
ggsave(file.path("./data_results/umap/", paste0("CD8_umap_hex_", Sys.Date(), ".svg")), plot = umap.dims.plot2,dpi=600)
ggsave(file.path("./data_results/umap/", paste0("CD8_umap_hex_", Sys.Date(), ".tiff")), plot = umap.dims.plot2,dpi=600)

#make plot per organ
dat.umap.samples$organindex <- factor(dat.umap.samples$organindex, levels = c('IV-CD45- Lung','BAL','IV-CD45+ Lung','Blood','MLN','Spleen'),
                                      labels = c('Lung Tissue','Airway','Lung Vasculature','Blood','MLN','Spleen'))

umap.organ<-ggplot(dat.umap.samples, aes(umap.1, umap.2)) + geom_hex(bins=75) + facet_wrap(~organindex,ncol = 2)+
  theme_dark() + 
  theme(panel.border = element_blank())+
  theme(panel.background = element_blank())+
  scale_fill_viridis(option = "plasma") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
  theme(legend.text = element_text(colour="black", size=11, face="bold"))+
  theme(legend.title = element_text(face = "bold"))+
  theme(strip.text.x = element_text(size = 11, color = "black", face = "bold"))+
  theme(strip.background = element_blank())
umap.organ

#make plot per mouse
dat.umap.samples$mouse <- factor(dat.umap.samples$mouse, levels = c('1','2','3','4','5'),
                                 labels = c('Mouse 1','Mouse 2','Mouse 3', 'Mouse 4', 'Mouse 5'))
umap.mouse<-ggplot(dat.umap.samples, aes(umap.1, umap.2)) + geom_hex(bins=75) + facet_wrap(~mouse) +
  theme_dark() + 
  theme(panel.border = element_blank())+
  theme(panel.background = element_blank())+
  scale_fill_viridis(option = "plasma") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
       # legend.position = "bottom")
  theme(legend.text = element_text(colour="black", size=9, face="bold"))+
  theme(legend.title = element_text(face = "bold"))+
  theme(strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
umap.mouse

#make plot per sample
umap.sample<-ggplot(dat.umap.samples, aes(umap.1, umap.2)) + geom_hex(bins=75) + facet_wrap(~sample) +
  theme_dark() + 
  theme(panel.border = element_blank())+
  theme(panel.background = element_blank())+
  scale_fill_viridis(option = "plasma") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
     #   legend.position = "bottom") +
  theme(legend.text = element_text(colour="black", size=9, face="bold"))+
  theme(legend.title = element_text(face = "bold"))+
  theme(strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
umap.sample

#make plot per spaceactivation (index created above for both the type of marker for antigen experimence (tetramer+ or CD44+) and space (tissue or vasculature)
dat.umap.samples$spaceactivation <- factor(dat.umap.samples$spaceactivation, levels = c('cd44_tissue','cd44_vasc','tet_tissue','tet_vasc','cd44_both','tet_both'),
                                 labels = c('CD44+ IV-CD45-', 'CD44+ IV-CD45+','Tetramer+ IV-CD45-','Tetramer+ IV-CD45+','CD44+','Tetramer+'))
umap.spaceactivation<-ggplot(dat.umap.samples, aes(umap.1, umap.2)) + geom_hex(bins=75) + facet_wrap(~spaceactivation) +
  theme_dark() + 
  theme(panel.border = element_blank())+
  theme(panel.background = element_blank())+
  scale_fill_viridis(option = "plasma") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
  theme(legend.text = element_text(colour="black", size=9, face="bold"))+
  theme(legend.title = element_text(face = "bold"))+
  theme(strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
umap.spaceactivation

#save plots
ggsave(file.path("./data_results/umap/", paste0("CD8_umap_hex_organ_ncol2_", Sys.Date(), ".png")), plot = umap.organ,dpi=600)
ggsave(file.path("./data_results/umap/", paste0("CD8_umap_hex_mouse_", Sys.Date(), ".png")), plot = umap.mouse,dpi=600)
ggsave(file.path("./data_results/umap/", paste0("CD8_umap_hex_sample_", Sys.Date(), ".png")), plot = umap.sample,dpi=600)
ggsave(file.path("./data_results/umap/", paste0("CD8_umap_hex_spaceactivation_", Sys.Date(), ".png")), plot = umap.spaceactivation,dpi=600)


#load RDS so we can resubset it
dat.umap.samples <- readRDS("./data_results/umap/CD8_umap_2020-12-15.rds")
#cd44ind<-grep('CD44',dat.umap.samples)
#remember that "cd44_both" are samples that were exported that were CD44+ from both tissue 
#spaces (ie IV-CD45 positivity was not factored in) and does not include samples from the lung
dat.umap.samples<-subset(dat.umap.samples, spaceactivation == "cd44_both")


#addding back in cd49a and cd103 so that umaps can show expression, 
#though to be clear the umap was already made without taking their expression into account
dims <- c("CD49a", 
          'CD103',
          "CD69",
          'CD62L',
          'CCR7',
          'CX3CR1',
          'KLRG1')
dims.for.umap <- dims
dat.umap <- data.frame(dat.umap.samples[, dims.for.umap],
                       dat.umap.samples[, c("umap.1", "umap.2")])
dat.umap[, dims.for.umap] <- range.scaled(dat.umap[, dims.for.umap])
dat.melt <- reshape2::melt(dat.umap, id = c("umap.1", "umap.2"))
umap.dims.plot <- ggplot(dat.melt, aes(umap.1, umap.2, z = value)) +
  stat_summary_hex(bins = 100) +
  theme(strip.background =element_rect(fill="black"))+ 
  theme(panel.border = element_blank())+
  theme(panel.background = element_blank())+
  scale_fill_viridis(option = "plasma") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
  theme(legend.text = element_text(colour="black", size=11, face="bold"))+
  theme(legend.title = element_text(face = "bold"))+
  facet_wrap(~variable)+
  theme(strip.text.x = element_text(size = 11, color = "black", face = "bold"))+
  theme(strip.background = element_blank())

umap.dims.plot
ggsave(file.path("./data_results/umap/", paste0("CD8_umap_hex_cd44.tissue_cd44.tissue_", Sys.Date(), ".png")), plot = umap.dims.plot)

#make organ plot
dat.umap.samples$organ <- factor(dat.umap.samples$organ, levels = c('bal','blood','lung','mln','spleen'),
                                 labels = c('BAL','Blood','Lung','MLN','Spleen'))
umap.organ<-ggplot(dat.umap.samples, aes(umap.1, umap.2)) + geom_hex(bins=75) + facet_wrap(~organindex)+
  theme_dark() + 
  theme(panel.border = element_blank())+
  theme(panel.background = element_blank())+
  scale_fill_viridis(option = "plasma") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
  theme(legend.text = element_text(colour="black", size=9, face="bold"))+
  theme(legend.title = element_text(face = "bold"))+
  theme(strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
umap.organ

#make mouse plot
dat.umap.samples$mouse <- factor(dat.umap.samples$mouse, levels = c('1','2','3','4','5'),
                                 labels = c('Mouse 1','Mouse 2','Mouse 3', 'Mouse 4', 'Mouse 5'))
umap.mouse<-ggplot(dat.umap.samples, aes(umap.1, umap.2)) + geom_hex(bins=75) + facet_wrap(~mouse) +
  theme_dark() + 
  theme(panel.border = element_blank())+
  theme(panel.background = element_blank())+
  scale_fill_viridis(option = "plasma") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
  theme(legend.text = element_text(colour="black", size=9, face="bold"))+
  theme(legend.title = element_text(face = "bold"))+
  theme(strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
umap.mouse

#make sample plot
umap.sample<-ggplot(dat.umap.samples, aes(umap.1, umap.2)) + geom_hex(bins=75) + facet_wrap(~sample) +
  theme_dark() + 
  theme(panel.border = element_blank())+
  theme(panel.background = element_blank())+
  scale_fill_viridis(option = "plasma") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
  theme(legend.text = element_text(colour="black", size=9, face="bold"))+
  theme(legend.title = element_text(face = "bold"))+
  theme(strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
umap.sample

#make space activation plot
dat.umap.samples$spaceactivation <- factor(dat.umap.samples$spaceactivation, levels = c('cd44_tissue','cd44_vasc','tet_tissue','tet_vasc'),
                                           labels = c('CD44+ IV-CD45-', 'CD44+ IV-CD45+','Tetramer+ IV-CD45-','Tetramer+ IV-CD45+'))
umap.spaceactivation<-ggplot(dat.umap.samples, aes(umap.1, umap.2)) + geom_hex(bins=75) + facet_wrap(~spaceactivation) +
  theme_dark() + 
  theme(panel.border = element_blank())+
  theme(panel.background = element_blank())+
  scale_fill_viridis(option = "plasma") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
  theme(legend.text = element_text(colour="black", size=9, face="bold"))+
  theme(legend.title = element_text(face = "bold"))+
  theme(strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
umap.spaceactivation

#save files
ggsave(file.path("./data_results/umap/", paste0("CD8_umap_hex_cd44.tissue_organ_", Sys.Date(), ".png")), plot = umap.organ,dpi=600)
ggsave(file.path("./data_results/umap/", paste0("CD8_umap_hex_cd44.tissue_mouse_", Sys.Date(), ".png")), plot = umap.mouse,dpi=600)
ggsave(file.path("./data_results/umap/", paste0("CD8_umap_hex_cd44.tissue_sample_", Sys.Date(), ".png")), plot = umap.sample,dpi=600)
ggsave(file.path("./data_results/umap/", paste0("CD8_umap_hex_cd44.tissue_spaceactivation_", Sys.Date(), ".png")), plot = umap.spaceactivation,dpi=600)
