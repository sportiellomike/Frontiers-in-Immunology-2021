##August 6th, 2020: D21 post-infection looking at various tissues/organs; standard T-cell memory markers; no treatments
##source files were pre-gated in FlowJo software and exported as subsets; live-singlets/CD8+/IV-CD45(+ or -)/CD44+ (or Tetramer +)
##    Figure 9
##September 9th, 2020: D21 post-infection; lung tissue; 4 phenotypes, 2 locations (IV); IV- (tissue) only; 6-hour restim (aCD3-aCD28)
##    Figure 3
source.experiments.names <- list.dirs("./data_source", recursive = F);print(source.experiments.names)

##WORKFLOW for August 6th, 2020 experiment; uses yyyymmdd naming convention
experiment.date <- "20200806"
experiment.name <- basename(grep(experiment.date, list.dirs("./data_source", recursive = F), value = T));print(experiment.name)
file.edit(file.path("./workflow_analysis_Rcode/", 
                    experiment.name, 
                    paste0(paste("workflow", experiment.name, "NGL", sep = "_"), ".R")))
##

