##function to generate a sample name data.frame based on established naming convention of .fcs files
##function assumes files were exported from FlowJO (appended with 'export' prefix)
meaningful.names.from.fcs.paths <- function(fcs.file.paths, split.by, split.names){
  meaningful.names <- sub('export', '' , basename(unlist(fcs.file.paths)))#subs out FlowJO default export prefix
  meaningful.names <- sub('_Specimen_001_', '' , meaningful.names)#subs out default FACS DIVA name
  meaningful.names <- sub('_[0-9]+_', '_', meaningful.names)#subs out default FACS DIVA three-digit acquisition order
  meaningful.names <- sub('.fcs', '', meaningful.names)#subs out the .fcs suffix
  
  split.vals <- strsplit(meaningful.names, split.by)#if name is still meaningful/unique, splitting by underscore ('_') will usually result in meaningful elements
  
  if(unique(sapply(split.vals, length)) == length(split.names)){
    m.names <- cbind(sample = meaningful.names,
                     setNames(data.frame(sapply(seq(length(split.names)), 
                                                function(x) sapply(split.vals, getElement, x))), 
                              nm = split.names)
    )
  }else{
    stop("Check naming convention")
  }
  m.names
}
##function to check the result of 'meaningful.names.from.fcs.paths'
##user has to know the expected values for each column;if using a naming-convention, these values should be known
meaningful.names.check <- function(meaningful.names.result, expected.numbers.vector){
  c.names <- colnames(meaningful.names.result)
  c.numbers <- setNames(expected.numbers.vector, c.names)
  if(length(c.names) != length(expected.numbers.vector)){
    stop("mis-match between number of columns and expected numbers")
  }
  c.names.length <- sapply(c.names, function(i) length(unique(meaningful.names.result[,i])))
  if(all(c.names.length == c.numbers)){
    sapply(c.names, function(i){
      paste(paste0(paste("Expected number of", paste0(i, '(s)')), ':'), c.numbers[i])
    })
  }
}
##
##