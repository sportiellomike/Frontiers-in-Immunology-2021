##create top-level folder
##top-level script;
##.Rproj (R Project at top-level for use of relative pathing)
## auto-fill can be used (...TAB)
if(any(grepl(".Rproj", list.files("./")))){
  message(".Rproj file exists at the top-level")
  message(list.files("./", pattern = ".Rproj"))
  message(list.files("./", pattern = ".Rproj", full.names = T))
}else{
  message("Create .Rproj file at top-level; associate with this existing directory")
}
##source data goes into './data_source/'; add source data via OS file management (Windows file explorer; MacOS Finder)
##experiment name should be unique: './data_source/yyyymmdd_UNIQUEINFORMATION_...'
if(!dir.exists("./data_source")){
  dir.create("./data_source")
}else{
  message("./data_source already exists")
}
##check source data/experiment names
source.experiments.names <- list.dirs("./data_source", recursive = F); print(source.experiments.names)
##generate derivative directories as they relate to established workflow;
##'./data_results'; './workflow_analysis_Rcode'; "./..."
derivative.directories <- c("data_results",
                            "workflow_analysis_Rcode")
sapply(derivative.directories, function(i){
  sapply(basename(source.experiments.names), function(j){
    if(!dir.exists(file.path(i,j))){
      dir.create(file.path(i,j), recursive = T)
    }else{
      message(paste(file.path(i,j), "already exists"))
    }
  })
})

if(!dir.exists("./workflow_analysis_Rcode/functions_depot")){
  dir.create("./workflow_analysis_Rcode/functions_depot")
}else{
  message("./workflow_analysis_Rcode/functions_depot already exists")
}
##
file.edit("./workflow_analysis_Rcode/workflow_NGL.R")
