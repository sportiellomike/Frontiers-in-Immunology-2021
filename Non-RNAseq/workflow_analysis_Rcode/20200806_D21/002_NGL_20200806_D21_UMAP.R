##using 'datMERGE', generate a UMAP embedding using select phenotypic markers
##'faceted' plotting of subsets using factored columns
##includes components for Figure 9

##load libraries
library(ggplot2)
library(scales)

library(svglite)
library(extrafont)
#load fonts for figure creation if necessary; run 'font_import()' once (takes a few minutes)
#font_import()
loadfonts(device = "win")

#source user-created functions
source("./workflow_analysis_Rcode/functions_depot/functions_NGL_flow.R")
source("./workflow_analysis_Rcode/functions_depot/functions_NGL_clustering.R")

if(!dir.exists("./data_results/20200806_D21/figure_9_components")){
  dir.create("./data_results/20200806_D21/figure_9_components", recursive = T)
}
##read in merged data; this .rds is a list
dat.mdat <- readRDS("./data_results/20200806_D21/20200806_D21_datMERGE.rds")
class(dat.mdat)
sapply(dat.mdat, is.data.frame)#`$dat.mdat` = flow data (compensated, transformed), merged with meta-data
dat <- dat.mdat$dat.mdat
##quality-control plot checks
args(mdat.plots.merged.df)
mdat.plots.merged.df(dat, "location", c("CD8a", "CD45_IV")) + geom_hex(bins = 100)#double-check correct association of actual data with 'location'; checks out
mdat.plots.merged.df(dat, "phenotype", c("CD103", "CD49a")) + geom_hex(bins = 100)#double-check correct association of actual data with 'phenotype'; checks out
mdat.plots.merged.df(dat, "phenotype", c("CD49a", "CD69")) + geom_hex(bins = 100)#double-check correct association of actual data with 'phenotype'; checks out
mdat.plots.merged.df(dat, "phenotype", c("KLRG1", "CD69")) + geom_hex(bins = 100)#double-check correct association of actual data with 'phenotype'; checks out
mdat.plots.merged.df(dat, "organ", c("CD103", "CD49a")) + geom_hex(bins = 100)#double-check correct association of actual data with 'organ';checks out
mdat.plots.merged.df(dat, "organ", c("CD49a", "CD69")) + geom_hex(bins = 100)
##some additional exploration of flow data as associated with the 'IV_CD45- Lung' subset
mdat.plots.merged.df(subset(dat, organ == 'IV_CD45- Lung'), "organ", c("CD49a", "CD69")) + geom_hex(bins = 100)
ggplot(subset(dat, organ == 'IV_CD45- Lung'), aes(x = CD49a)) + geom_density() + geom_vline(xintercept = 2.5, color = "red")#CD49a+ ~2.5
ggplot(subset(dat, organ == 'IV_CD45- Lung'), aes(x = CD69)) + geom_density() + geom_vline(xintercept = 1.0, color = "red")#CD69+ ~1
mdat.plots.merged.df(subset(dat, organ == 'IV_CD45- Lung'), "organ", c("CD49a", "CD69")) + geom_hex(bins = 100) +
  geom_vline(xintercept = 2.5, color = "red") +
  geom_hline(yintercept = 1.0, color = "red")
mdat.plots.merged.df(subset(dat, organ == 'IV_CD45- Lung' & CD49a >= 2.5 & CD69 >= 1), "organ", c("CD103", "CD49a")) + geom_hex(bins = 100)
rm(dat)
##decide which dims to use for generating UMAP embedding; 
##available options:  c('KLRG1', 'CX3CR1', 'CCR7', 'CD69', 'CD62L'); other/remaining dimensions were used during export; 
##CD103 and CD49a were used during FlowJO export to define 4 'phenotypes'; quadrant gate: 4 'phenotypes'; 
##Since CD103 and CD49a define the 'phenotypes', they will be included in the UMAP embedding
dims.of.interest <- c(
  'CD49a', 
  'CD103',
  'KLRG1',
  'CX3CR1',
  'CCR7',
  'CD69',
  'CD62L'
)
##The area code of Rochester, NY - where this work was performed - is 585;
##The area code of Milwaukee, WI - where Mike was born born and raised - is 414; his mom got scared after his little fight and sent him off to Rochester; 
##585414 seemed as good a seed as any other...to perfectly reproduce the following results, the seed (along with R/package versions) must be identical; 
##however, to very closely approximate our results the seed should have no impact 
dat.mdat$umap$seed.val <- 585414#append to list
dat.mdat$umap$dimensions <- dims.of.interest#append to list
#####THE STEPS BELOW ARE ABSOLUTELY ESSENTIAL AND MUST BE 'UNCOMMENTED' FOR THE FOLLOWING FIGURE(S) GENERATION TO WORK#####
##the step below is the formation of the UMAP; be aware that it's a time-consuming process (embedding of 422,220 events)
##The final embedding (2-dimensions/columns) was saved as an .RDS ('saveRDS()'; compressed .R object), then read back in ('readRDS()') - as needed - to preclude re-running UMAP.
##In addition to saving the UMAP results, a merged data frame (input.flowset expression matrix (compensated and transformed) + meta-data factors + UMAP dimensions) was also saved (.RDS)
set.seed(dat.mdat$umap$seed.val)
umap.dims <- uwot::umap(X = dat.mdat$dat.mdat[, dat.mdat$umap$dimensions],
                        nn_method = "annoy",
                        n_threads = parallel::detectCores()-1,
                        verbose = TRUE)
##7 minutes to run using 11 cores/threads (Windows PC; i7-8700 CPU @3.2 GHz)
##5 minutes to run using 3 cores/threads (Mac OS; i5 CPU @3.2 GHz);...interesting time compared to more powerful PC...
ggplot(setNames(data.frame(umap.dims), nm = paste("umap", seq(ncol(umap.dims)), sep = ".")), aes(x = umap.1, y = umap.2)) + 
  geom_hex(bins = 200)#quick results/plot check
##
if(!dir.exists("./data_results/20200806_D21/umap")){
  dir.create("./data_results/20200806_D21/umap")
}
dat.mdat$umap$umap.dims <- umap.dims
dat.mdat$dat.mdat <- cbind(dat.mdat$dat.mdat,
                           setNames(data.frame(umap.dims), nm = paste("umap", seq(ncol(umap.dims)), sep = ".")))
saveRDS(dat.mdat, "./data_results/20200806_D21/20200806_D21_datMERGE_UMAP.rds")
head(dat.mdat$dat.mdat)
##R can be closed/cleared out at this point; 
##if restarting, re-establish workflow by loading appropriate libraries/functions (top of script) and reading in merged data
dat <- readRDS("./data_results/20200806_D21/20200806_D21_datMERGE_UMAP.rds")
class(dat)
sapply(dat, is.data.frame)#`$dat.mdat` = flow data (compensated, transformed), merged with meta-data and UMAP dimensions
dat <- dat$dat.mdat
##use range scaling to set values between 0-to-1 for marker expression;
##slightly different scaling effects depending on if you range scale the entire data set, or range scale only the subset to be plotted
scale.dims <- c(
  'CD49a', 
  'CD103',
  "CD69",
  'CD62L',
  'CCR7',
  'CX3CR1',
  'KLRG1'
)#same dimensions as was used to generate UMAP embedding
summary(dat[, scale.dims])
args(range.scaled)
dat <- range.scaled(dat, scale.dims)
summary(dat[, scale.dims])
##
source("./workflow_analysis_Rcode/functions_depot/functions_NGL_plots.R")
plot.dims <- c(scale.dims, paste("umap", 1:2, sep = "."))
##'UMAP Marker plots' plot list; will eventually write to file
plots.to.save_UMAP.markers <- list()
##plot UMAP embedding showing individual marker expression
plots.to.save_UMAP.markers$ALLdata <- umap.markers.plot(dat, plot.dims)
##show only 'CD44+' ('activation.state') and 'IV_CD45- Lung (tissue)' (organ) for the following plot
plots.to.save_UMAP.markers$CD44p.Lung.tissue <- umap.markers.plot(subset(dat, activation.state == 'CD44+' & organ =='IV_CD45- Lung'), plot.dims) + 
  stat_summary_hex(bins = 100)
##saving plots with different file types; 
##the quality/appearance of the generated plots can be influenced by operating system graphics devices; Windows vs MacOS sometimes shows differences
out.dir <- ("./data_results/20200806_D21/umap/")
#file.types <- c(".png", ".svg", ".tiff")#'.svg' and '.tiff' file-types are quite large;dropping those as a type; uncomment if needed
file.types <- c(".png"
                # ,".svg"
                # ,".tiff"
)
lapply(names(plots.to.save_UMAP.markers), function(j){
  sapply(file.types, function(i){
    f.name <- paste0(paste("UMAP_Markers", j, Sys.Date(), sep = "_"), i)
    ggsave(file.path(out.dir, f.name), plot = plots.to.save_UMAP.markers[[j]] +theme(plot.title =element_text(size=9, face='bold')), dpi = 600,        width = 3.75,        height = 4.4)
  })
})
ggsave("./data_results/20200806_D21/figure_9_components/figure_9_B_UMAP_markers_LungTissue_DATA.png",
       plot = plots.to.save_UMAP.markers$CD44p.Lung.tissue + theme(strip.text.x = element_text( size = 14,face='bold'),legend.position = "none"),
       dpi = 600,
       width = 3.75,
       height = 4.4)
library(ggpubr)
ggsave("./data_results/20200806_D21/figure_9_components/figure_9_B_UMAP_markers_LungTissue_LEGEND.png",
       plot = as_ggplot(get_legend(plots.to.save_UMAP.markers$CD44p.Lung.tissue + theme(legend.title = element_blank(), 
                                                                                        legend.direction = "horizontal", 
                                                                                        legend.text = element_text(angle = 300, vjust = 1.15), 
                                                                                        legend.text.align = 0,
                                                                                        plot.title =element_text(size=9, face='bold')))),
       dpi = 600,
       width = 3.75,
       height = 4.4)
##'UMAP Subset plots' plot list; will eventually write to file
plots.to.save_UMAP.subsets <- list()
##plot UMAP embedding for each organ
dat$organ <- factor(dat$organ, 
                    levels = c('IV_CD45- Lung', 'BAL', 'IV_CD45+ Lung', 'Blood', 'MLN', 'Spleen'),
                    labels = c('Lung Tissue', 'Airway', 'Lung Vasculature', 'Blood', 'MLN', 'Spleen'))#re-order factor for plotting

umap.count.plot(dat, "organ", squish.max.val = 250)
plots.to.save_UMAP.subsets$organ <- umap.count.plot(dat, "organ", squish.max.val = 250)
ggsave("./data_results/20200806_D21/figure_9_components/figure_9_A_UMAP_counts_organs_DATA.png",
       plot = umap.count.plot(dat, "organ", squish.max.val = 250, dummy.facet = TRUE) + theme(strip.text.x = element_text( size = 14,face='bold'),legend.position = "none"),
       dpi = 600,
       width = 3.75,
       height = 4.4)
library(ggpubr)
ggsave("./data_results/20200806_D21/figure_9_components/figure_9_A_UMAP_counts_organs_LEGEND.png",
       plot = as_ggplot(get_legend(umap.count.plot(dat, "organ", squish.max.val = 250, dummy.facet = TRUE) + theme(strip.text.x = element_text( size = 14,face='bold'),
                                                                                                                   legend.title = element_blank(), 
                                                                                                                   legend.direction = "horizontal", 
                                                                                                                   legend.text = element_text(angle = 300, vjust = 1.15), 
                                                                                                                   legend.text.align = 0))),
       dpi = 600,
       width = 3.75,
       height = 4.4)
#plot UMAP embedding for each mouse/replicate
dat$replicate <- factor(dat$replicate,
                        levels = c(1:5),
                        labels = paste("Mouse", 1:5, sep = "."))

umap.count.plot(dat, "replicate", squish.max.val = 200) + geom_hex(bins = 55)
plots.to.save_UMAP.subsets$replicate = umap.count.plot(dat, "replicate", squish.max.val = 200) + geom_hex(bins = 55)
##plot UMAP embedding per 'activation.state.location'
dat$activation.state.location <- factor(dat$activation.state.location,
                                        levels = c('CD44+_tissue', 'CD44+_vasc', 'tet+_tissue', 'tet+_vasc', 'CD44+_both', 'tet+_both'),
                                        labels = c('CD44+ IV_CD45-', 'CD44+ IV_CD45+', 'Tetramer+ IV_CD45-', 'Tetramer+ IV_CD45+', 'CD44+', 'Tetramer+'))

umap.count.plot(dat, "activation.state.location", squish.max.val = 400)
plots.to.save_UMAP.subsets$activation.state.location <- umap.count.plot(dat, "activation.state.location", squish.max.val = 400)
##saving above plot list with different file types; 
##the quality/appearance of the generated plots can be influenced by operating system graphics devices; Windows vs MacOS sometimes shows differences
out.dir <- ("./data_results/20200806_D21/umap/")
#file.types <- c(".png", ".svg", ".tiff")#'.svg' and '.tiff' file-types are quite large;dropping those as a type; uncomment if needed
file.types <- c(".png"
                # ,".svg"
                # ,".tiff"
)
lapply(names(plots.to.save_UMAP.subsets), function(j){
  sapply(file.types, function(i){
    f.name <- paste0(paste("UMAP_counts", j, Sys.Date(), sep = "_"), i)
    ggsave(file.path(out.dir, f.name), plot = plots.to.save_UMAP.subsets[[j]]+theme(plot.title =element_text(size=9, face='bold',color='red')), dpi = 600,        width = 3.75,        height = 4.4)
  })
})
##
##

