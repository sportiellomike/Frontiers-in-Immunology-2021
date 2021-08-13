##
##function to calculate bi-exponential transformation parameters (parms) for each fluor channel
elgcl.parms <- function(compensated.set) {
  dat <- fsApply(compensated.set, exprs)#expression values for all flowset files
  channels.fluor <- names(markernames(compensated.set))#fluor channels;to be transformed
  dat <- dat[, channels.fluor]#fluors only;drop scatter, time
  elgcl <- estimateLogicle(new("flowFrame", dat), channels = channels.fluor)#estimates the transformation, per channel, based on input data (dat)
  elgcl
}
##
##function to rename colnames/markernames of a flowset;requires naming convention: "antigen SPACE fluor"
markers.rename <- function(input.flowset){
  m <- colnames(input.flowset)
  m[colnames(input.flowset) %in% names(markernames(input.flowset))] <- markernames(input.flowset)#replace PMT/channel with marker name
  m <- sub("-", ".", m)#sub out '-'
  m <- sapply(strsplit(m, " "), getElement, 1)#get only marker name;drop the fluor;follows established naming convention
  m
}
##
##function for generating ggplot objects based on factored meta-data columns and bi-variate marker pairs; add a plotting geometry; uses flowsets
mdat.plots.input.flowset <- function(input.flowset, mdat, mdat.factor.name, plot.pairs){
  if(!any(sapply(mdat, function(i) all(i == sampleNames(input.flowset))))){
    stop("No match between mdat and input.flowset sampleNames (by 'basename(fcs.file.path)')")
  }
  if(!is.factor(mdat[, mdat.factor.name])){
    stop("Need a factored column")
  }
  
  dat <- as.data.frame(fsApply(input.flowset, exprs))[, plot.pairs]
  dat[mdat.factor.name] <- rep(mdat[, mdat.factor.name], fsApply(input.flowset, nrow))
  
  require(ggplot2)
  ggplot(dat, aes_string(x = plot.pairs[[1]], y = plot.pairs[[2]])) + 
    facet_wrap(c(mdat.factor.name), nrow = ceiling(length(levels(mdat[, mdat.factor.name]))/2))
}
##
##function for generating ggplot objects based on factored meta-data columns and bi-variate marker pairs; add a plotting geometry; uses data.frames
mdat.plots.merged.df <- function(input.data.merged, mdat.factor.name, plot.pairs){

  if(!is.factor(input.data.merged[, mdat.factor.name])){
    stop("Need a factored column")
  }
  
  dat <- input.data.merged[, c(plot.pairs, mdat.factor.name)]
  
  require(ggplot2)
  ggplot(dat, aes_string(x = plot.pairs[[1]], y = plot.pairs[[2]])) + 
    facet_wrap(c(mdat.factor.name), nrow = ceiling(length(levels(mdat[, mdat.factor.name]))/2))
}
##
##