##generate a 'meaningful names' data.frame from existing .fcs file names; 
##to be merged with exported .fcs data following compensation and transformation;

##source user-created functions
source("./workflow_analysis_Rcode/functions_depot/functions_NGL_utility.R")
#list of exported (FlowJo export) .FCS files
exported.fcs <- list.files("./data_source/20200806_D21/exported_FCS/", pattern = ".fcs", full.names = T)
length(exported.fcs)#240 files; 
#4 tissues (Spleen, MLN, BAL, Blood) * 5 replicates each * 4 phenotypes with 2 'activation states': (4*5)*(4*2) = 160
#1 tissue (Lung) * 5 replicates each * 4 phenotypes with 2 'activation states' with 2 'locations': (1*5)*(4*2*2) = 80

#construct data.frame based on naming convention;
#element 1: organ+replicate
#element 2: activation state
#element 3: phenotype (quadrant gate (used as export gate))
#   element 3:  Q1 = CD103- CD49a+; alias = CD49a
#               Q2 = CD103+ CD49a+; alias = DP (double-positive)
#               Q3 = CD103+ CD49a-; alias = CD103
#               Q4 = CD103- CD49a-; alias = DN (double-negative)
#element 4: location (from tissue or vasculature; alias is usually "IV_CD45- Lung" (tissue) or "IV_CD45+ Lung" (vasculature))

args(meaningful.names.from.fcs.paths)#user-created (NGL) function
mdat <- meaningful.names.from.fcs.paths(fcs.file.paths = exported.fcs,
                                        split.by = '_',
                                        split.names = c("organ", "activation.state", "phenotype", "location")
)
mdat$location <- sub("tiss", "tissue", mdat$location)
mdat <- tidyr::separate(mdat, organ, c("organ", "replicate"), sep = "(?<=[A-Za-z])(?=[0-9])")#split organ into organ + replicate
head(mdat)

args(meaningful.names.check)#user-created (NGL) function
meaningful.names.check(meaningful.names.result = mdat, 
                       expected.numbers.vector = c(240,5,5,2,4,3)
)

mdat$phenotype.alias <- mdat$phenotype#preserve original phenotype column as an 'alias'
mdat$phenotype <- sapply(sapply(strsplit(mdat$phenotype, " "), function(j) getElement(j, 1)), function(i){
  if(i == "Q1"){
    i <- "CD49a"
  }else if(i == "Q2"){
    i <- "DP"
  }else if(i == "Q3"){
    i <- "CD103"
  }else if(i == "Q4"){
    i <- "DN"
  }
}, USE.NAMES = F)#create new phenotype names
head(mdat)

mdat$organ.alias <- mdat$organ#preserve original organ column as an 'alias'
mdat$organ <- sub("S$", "Spleen", sub("B$", "Blood", mdat$organ))
mdat$organ[which(mdat$location == "tissue")] <- "IV_CD45- Lung"
mdat$organ[which(mdat$location == "vasc")] <- "IV_CD45+ Lung"
unique(mdat$organ)

mdat$organ.replicate <- paste(mdat$organ, mdat$replicate, sep = ".")
head(mdat)

mdat$activation.state.location <- paste(mdat$activation.state, mdat$location, sep = "_")
head(mdat)
mdat$sample.og <- basename(exported.fcs)#preserve orginal .fcs filename

write.csv(mdat, "./data_results/20200806_D21/20200806_D21_meaningful_names.csv", row.names = F)
mdat.read <- read.csv("./data_results/20200806_D21/20200806_D21_meaningful_names.csv", stringsAsFactors = FALSE)#no factor levels
str(mdat.read)#all character vectors except for replicate (integer)
all(mdat == mdat.read)#checks out
mdat.read <- read.csv("./data_results/20200806_D21/20200806_D21_meaningful_names.csv", stringsAsFactors = TRUE, 
                      colClasses = c("sample" = "character",
                                     "sample.og" = "character",
                                     "replicate" = "factor"))#with factor levels for all but sample/sample.og
str(mdat.read)
all(mdat == mdat.read)#checks out
mdat <- mdat.read
rm(mdat.read)
##read, compensate, transform exported .fcs data;
##merge with mdat

##load libraries
library(flowCore)
library(ggplot2)
#source user-created functions
source("./workflow_analysis_Rcode/functions_depot/functions_NGL_flow.R")
##load exported (FlowJo export) .FCS files as a flowSet
exported.fcs <- list.files("./data_source/20200806_D21/exported_FCS/", pattern = ".fcs", full.names = T)
length(exported.fcs) == length(mdat$sample.og);print(paste(length(exported.fcs), length(mdat$sample.og)))
all(basename(exported.fcs) == mdat$sample.og)#
##
input.flowset <- read.flowSet(exported.fcs, transformation =FALSE, truncate_max_range = FALSE)
##compensate using a modified compensation matrix; the following comp matrix (.csv) was generated/exported in FlowJO (modified existing DIVA matrix)
comp.matrix.edit <- as.matrix(read.csv("./data_source/20200806_D21/matrix.12.9.20.csv", row.names = 1, check.names = F))
all(sapply(colnames(comp.matrix.edit), function(i) unlist(strsplit(i, " :: "))[[1]]) == colnames(input.flowset[[1]]@description$`$SPILLOVER`))
colnames(comp.matrix.edit) <- colnames(input.flowset[[1]]@description$`$SPILLOVER`)
rownames(comp.matrix.edit) <- colnames(comp.matrix.edit)
input.flowset <- fsApply(input.flowset, function(i) compensate(i, comp.matrix.edit))
##transform using estimate logicle/bi-exponential; need to estimate 'width (w)' parameter first, for each channel
elgcl <- elgcl.parms(input.flowset)
input.flowset <- fsApply(input.flowset, function(i) transform(i, elgcl))
##create dat for storing data, meta-data, modified compensation matrix, transformation parameters, and original column/marker  names
dat <- list(dat.mdat = NULL,
            mdat = mdat,
            comp.matrix.edit = comp.matrix.edit,
            elgcl = elgcl,
            columns.og = colnames(input.flowset),
            markers.og = markernames(input.flowset))
##modify colnames; for easier plotting/readability
colnames(input.flowset) <- markers.rename(input.flowset)
##quality-control plot checks
args(mdat.plots.input.flowset)
mdat.plots.input.flowset(input.flowset, mdat, "location", c("CD8a", "CD45_IV")) + geom_hex(bins = 100)#double-check correct association of actual data with exported location; checks out
mdat.plots.input.flowset(input.flowset, mdat, "phenotype", c("CD103", "CD49a")) + geom_hex(bins = 100)#double-check correct association of actual data with exported phenotype; checks out
mdat.plots.input.flowset(input.flowset, mdat, "phenotype", c("CD49a", "CD69")) + geom_hex(bins = 100)#double-check correct association of actual data with exported phenotype; checks out
mdat.plots.input.flowset(input.flowset, mdat, "organ", c("CD103", "CD49a")) + geom_hex(bins = 100)#double-check correct association of actual data with exported organ;checks out
##input.flowset expression matrix can be merged with mdat; 
##input.flowset expression matrix (dat) needs a 'merge by' column
merge.col <- colnames(mdat)[which(sapply(mdat, function(i) all(i == sampleNames(input.flowset))))]#name of column in mdat that matches input.flowset sample names
dat$dat.mdat = setNames(data.frame(fsApply(input.flowset, exprs), rep(sampleNames(input.flowset), fsApply(input.flowset, nrow))),
                        nm = c(colnames(input.flowset), merge.col))#two exported samples had zero events and are dropped during data.frame creation
cat("These samples have zero events (dropped):\n", paste0(sampleNames(input.flowset)[which(fsApply(input.flowset, nrow) == 0)], sep = '\n'))#dropped samples
dat$dat.mdat <- plyr::join(dat$dat.mdat, mdat, by = "sample.og")
mdat.plots.input.flowset(input.flowset, mdat, "organ", c("CD103", "CD49a")) + geom_hex(bins = 100)#compare facet: 'IV_CD45- Lung'
ggplot(subset(dat$dat.mdat, dat$dat.mdat$organ == 'IV_CD45- Lung'), aes(x = CD103, y = CD49a)) + geom_hex(bins = 200)#double-check merge results;same as above mdat.plot facet (pre-merge)
##save merged data as a .RDS object; two exported samples had zero events and are dropped during the merge; 
##merged data will be read back in ('readRDS()') for downstream analysis
saveRDS(dat, "./data_results/20200806_D21/20200806_D21_datMERGE.rds")
##
##