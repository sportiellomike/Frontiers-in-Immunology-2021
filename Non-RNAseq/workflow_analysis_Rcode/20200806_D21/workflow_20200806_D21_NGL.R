##001: generating a 'mdat' (meta-data) data.frame of meaningful names from FlowJO exported .fcs files
##FACS DIVA/FlowJo-exported .fcs files will have extraneous information in their names if using default values
##trim out these bits; split name into various factored columns; test for uniqueness/conformity
##resulting 'mdat' will be merged with input flow data for downstream analysis/plots
##read, compensate, transform, and merge exported .fcs data with mdat; store as .RDS for downstream use
analysis.step <- "001"
initials <- "NGL"
experiment.name <- basename(grep("20200806", list.dirs("./data_source", recursive = F), value = T))
this.step.script.name <- "dat_mdat_MERGE_DF"

script <- file.path("./workflow_analysis_Rcode/", 
                    experiment.name, 
                    paste0(paste(analysis.step, initials, experiment.name, this.step.script.name, sep = "_"), ".R")
);print(script)
file.edit(script)
##002: generating a UMAP using saved 'datMERGE' (exported .fcs files (compensated, transformed); merged with mdat);
##After UMAP is generated, 'datMERGE' is appended with UMAP dimensions and the list is updated/saved (.RDS)
##Various plots
analysis.step <- "002"
initials <- "NGL"
experiment.name <- basename(grep("20200806", list.dirs("./data_source", recursive = F), value = T))
this.step.script.name <- "UMAP"

script <- file.path("./workflow_analysis_Rcode/", 
                    experiment.name, 
                    paste0(paste(analysis.step, initials, experiment.name, this.step.script.name, sep = "_"), ".R")
);print(script)
file.edit(script)
##003: generating a FlowSOM clustering using saved '.datMERGE_UMAP' input data (compensated, transformed .fcs data);
##After FlowSOM object creation and clustering, 'datMERGE_UMAP' is appended with FlowSOM results and the list is updated/saved (.RDS)
##Various plots
analysis.step <- "003"
initials <- "NGL"
experiment.name <- basename(grep("20200806", list.dirs("./data_source", recursive = F), value = T))
this.step.script.name <- "FSOM"

script <- file.path("./workflow_analysis_Rcode/", 
                    experiment.name, 
                    paste0(paste(analysis.step, initials, experiment.name, this.step.script.name, sep = "_"), ".R")
);print(script)
file.edit(script)
##
##