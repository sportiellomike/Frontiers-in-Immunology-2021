##
##
fsom.cluster.medians <- function(fsom, type = c("codes", "map")) {
 
  dat <- switch(type,
                codes = data.frame(fsom$map$codes, cluster = fsom$metaclustering),
                map = data.frame(fsom$map$medianValues, cluster = fsom$metaclustering))

  dat <- dat %>% dplyr::group_by(cluster) %>% dplyr::summarise_all(list(median))
  dat <- dat[, c("cluster", fsom$dims.used)]
  
}
##
##
range.scaled <- function(dat, cols, probs = c(0.01, 0.99)) {
  tmp <- as.matrix(dat[, cols])
  ranges <- matrixStats::colQuantiles(tmp, probs = probs)
  tmp <- t((t(tmp) - ranges[, 1]) / (ranges[, 2] - ranges[, 1]))
  tmp[tmp < 0] <- 0
  tmp[tmp > 1] <- 1
  dat[, cols] <- tmp
  dat
}
##
##
heatmap.cluster <- function(fsom, range.scale = TRUE, color = c("base", "viridis"), ...){#... to pass in updated quantiles to range.scaled function ('probs = c(#,#)')

  dat <- fsom.cluster.medians(fsom, ...)
  
  if(range.scale == TRUE){
    dat <- range.scaled(dat[, !colnames(dat) %in% "cluster"])
  }else{
    dat <- dat[, !colnames(dat) %in% "cluster"]
  }
  
  if(color == "base"){
    pheatmap::pheatmap(dat,
                       number_color = "black", fontsize = 11, fontsize_number = 6,
                       cluster_rows = FALSE,
                       angle_col = "45")
  }else if(color == "viridis"){
    color_heat <- viridis::viridis(option = "plasma", 50)
    pheatmap::pheatmap(dat, color = color_heat,
                       number_color = "black", fontsize = 11, fontsize_number = 6,
                       cluster_rows = FALSE,
                       angle_col = "45")
  }
}
##
##
#function to scratch build a FlowSOM object; 
#FlowSOM::ReadInput can be used but that function grows the input matrix using rbind for every frame/.fcs file in the flowset
#FlowSOM::ReadInput is unusable for large files/many files, so scratch build instead
fsom_scratch_build <- function(input.data, scale = FALSE, scale.cols = NULL){
  fsom <- list(data = NULL,
               pattern = ".fcs", 
               compensate = FALSE, 
               spillover = NULL, 
               transform = FALSE, 
               toTransform = NULL, 
               transformFunction = NULL, 
               scale = FALSE,
               silent = FALSE)
  class(fsom) <- "FlowSOM"
  
  if(class(input.data) == "flowSet"){
    message("Building from a flowSet")
    
    if(is.null(fsom$prettyColnames)){
      n <- input.flowset[[1]]@parameters@data[, "name"]
      d <- input.flowset[[1]]@parameters@data[, "desc"]
      d[is.na(d)] <- n[is.na(d)]
      fsom$prettyColnames <- paste0(d, " <", n, ">")
      names(fsom$prettyColnames) <- colnames(input.flowset)
    }
    
    fsom$metaData <- vector("list", length(sampleNames(input.flowset)))
    names(fsom$metaData) <- sampleNames(input.flowset)
    for(i in names(fsom$metaData)){
      fsom$metaData[[i]] <- c(1,nrow(input.flowset[[i]]))
    }
    rolling.total <- cumsum(unlist(sapply(fsom$metaData, function(i) i[2])))
    start.total <- unlist(sapply(fsom$metaData, function(i) i[2]))
    for(i in seq(length(fsom$metaData))){
      fsom$metaData[[i]] <- unname(c((rolling.total[i] - start.total[i])+1, rolling.total[i]))
    }
    
    if(scale == TRUE){
      message("Scaling scatter and fluor columns")
      fsom$scale <- TRUE
      #fsom$data <- fsApply(input.flowset, scale, use.exprs = TRUE)
      fsom$data <- fsApply(input.flowset, exprs)
      scale.index <- c(which(grepl("FSC|SSC", colnames(fsom$data))), which(colnames(fsom$data) %in% names(markernames(input.flowset))))
      fsom$data[, scale.index] <- scale(fsom$data[, scale.index])
    }else{
      fsom$data <- fsApply(input.flowset, exprs)
    }
    
  }else if(class(input.data) == "data.frame"){
    message("Building from a data.frame")
    
    if(is.null(fsom$prettyColnames)){
      c.names <- colnames(input.data)[sapply(input.data, is.numeric)]
      c.names <- c.names[!c.names %in% c("umap.1", "umap.2")]
      fsom$prettyColnames <- c.names
    }
    
    sample.fcs.col <- colnames(input.data)[sapply(input.data[1, ], function(i) grepl(".fcs", i))]#find column in input.data that contains a list of input files (.fcs)
    fcs.row.count <- table(input.data[, sample.fcs.col])
    fcs.row.end <- cumsum(fcs.row.count)
    fcs.row.start <- (fcs.row.end - fcs.row.count)+1

    fsom$metaData <- vector("list", length(fcs.row.count))
    names(fsom$metaData) <- names(fcs.row.count)
    for(i in names(fcs.row.count)){
      fsom$metaData[[i]] <- c(fcs.row.start[[i]], fcs.row.end[[i]])
    }
    
    message("Dropping non-numeric columns from input data:\n", 
            paste(colnames(dat)[!sapply(dat, is.numeric)], collapse = ", "))
    
    if(scale == TRUE){
      message("Scaling scatter and fluor columns")
      fsom$scale <- TRUE
      fsom$data <- as.matrix(input.data[, sapply(input.data, is.numeric)])
      scale.index <- scale.cols
      fsom$data[, scale.index] <- scale(fsom$data[, scale.index])
    }else{
      fsom$data <- as.matrix(input.data[, sapply(input.data, is.numeric)])
    }
    
  }
  
  fsom
}
##
##
##
##
sample_cluster.counts_fsom <- function(fsom){
  cluster.vec <- fsom$metaclustering[fsom$map$mapping[, 1]]
  if(length(cluster.vec) != max(unlist(fsom$metaData))){
    stop("cluster.vec does not match length of metaData")
  }
  cluster.counts <- sapply(names(fsom$metaData), function(i){
    cluster.vec[fsom$metaData[[i]][1]:fsom$metaData[[i]][2]]
  })
  
  cluster.counts <- matrix(t(sapply(cluster.counts, table)), 
                           nrow = length(cluster.counts), 
                           dimnames = list(names(cluster.counts),
                                           seq(length(levels(fsom$metaclustering)))
                           )
  )
}
##
##
clusters.tabulated.list <- function(fsom, sub.name.terms = NULL){
  c.counts <- sample_cluster.counts_fsom(fsom)
  
  if(!is.null(sub.name.terms)){
    samples = gsub(sub.name.terms, "", rownames(c.counts))
  }else{
    samples = rownames(c.counts)
  }
  
  clusters.tabulated <- list(count = cbind(sample = samples, as.data.frame(c.counts, row.names = F)),
                             per.million = cbind(sample = samples,
                                                 as.data.frame(t(apply(c.counts, 1, function(x) x/sum(x)*1E6)), row.names = F)),
                             proportion = cbind(sample = samples,
                                                as.data.frame(prop.table(c.counts, 1)*100, row.names = F))
  )
  
  if(!is.null(fsom$mdat)){
    if(!any(sapply(fsom$mdat, function(i) any(i %in% names(fsom$metaData))))){
      stop("No match between fsom$mdat and 'names(fsom$metaData)' (derived from: 'sampleNames(input.flowset)';'basename(fcs.file.path)')")
    }
    merge.col.name <- names(which(sapply(fsom$mdat, function(i) any(i %in% names(fsom$metaData)))))
    clusters.tabulated <- lapply(clusters.tabulated, function(i){
      colnames(i)[grep("sample", colnames(i))] <- merge.col.name
      i
    })
    clusters.tabulated <- sapply(clusters.tabulated, function(i) plyr::join(i, fsom$mdat, by = merge.col.name), simplify = F)
  }
  clusters.tabulated
}
##
##
colSums.intersect <- function(input.data, intersect.cols = c(col1, col2), intersect.vals = c(val1, val2)){
  colSums(input.data[intersect(which(input.data[[intersect.cols[1]]] == intersect.vals[1]), 
                               which(input.data[[intersect.cols[2]]] == intersect.vals[2])), which(sapply(input.data, is.numeric))])
}
##
##
prop.table.intersect <- function(input.data, input.data.factor, intersect.cols = c(col1, col2), intersect.val){
  dat.sub <- t(sapply(levels(input.data[[input.data.factor]]), function(i){
    colSums.intersect(input.data, intersect.cols = intersect.cols, intersect.vals = c(i, intersect.val))
  }))
  if(all(round(rowSums(prop.table(dat.sub, 1)*100), 10) == 100)){
    dat.sub <- prop.table(dat.sub, 1)*100
    if(any(class(dat.sub) != "data.frame")){
      dat.sub <- as.data.frame(dat.sub)
      dat.sub$sample <- as.factor(rownames(dat.sub))
      return(dat.sub)
    }
  }else{
    stop("Hmm...")
  }
}

##
##function to generate various statistical outputs; specific to 20200806_D21 dataset; room for improvement here (generalize/streamline)...but it works
fsom.stats <- function(input.count.data, combine.for.new.factor, intersect.col, intersect.val,
                       new.factor.order = list(factor.name, new.order), split.name, 
                       grouping.name, stats.test){
  # if(is.null(fsom$clusters.tab)){
  #   stop("Need cluster count data, stored at '$clusters.tab'")
  #   if(is.null(fsom$clusters.tab$count)){
  #     stop("Need cluster count data, stored at '$clusters.tab$count'")
  #   }
  # }
  
  new.factor.name <- paste0(combine.for.new.factor, collapse = ".")
  
  dat <- input.count.data
  dat[[new.factor.name]] <- as.factor(apply(dat[, combine.for.new.factor], 1, 
                                            function(i) paste0(i, collapse = ".")))
  dat <- prop.table.intersect(input.data = dat,
                              input.data.factor = new.factor.name,
                              intersect.cols = c(new.factor.name, intersect.col),
                              intersect.val = intersect.val)
  
  dat <- cbind(dat, 
               lapply(tidyr::separate(data.frame(sample = dat$sample), 
                                      sample, 
                                      into = c(combine.for.new.factor), sep = "\\."), as.factor)
  )
  
  dat[[new.factor.order[[1]]]] <- factor(dat[[new.factor.order[[1]]]], levels = new.factor.order[[2]])
  
  colnames(dat)[sapply(dat, is.numeric)] <- paste("cluster", which(sapply(dat, is.numeric)), sep = ".")
  
  dat.split <- lapply(split(dat, dat[[split.name]]), droplevels)
  
  if(stats.test == "aov"){
    aov.results <- lapply(dat.split, function(i){
      stats.results.list <- setNames(vector("list", length = (length(grep("cluster", colnames(i)))+1)), 
                                     nm = c(grep("cluster", colnames(i), value = T), "aov.summary.results"))
      for(j in grep("cluster", colnames(i), value = T)){
        stats.tmp <- i[, c(j, "phenotype", "replicate")]
        res.aov <- aov(as.formula(as.formula(paste(j, "~", grouping.name, sep = " "))),
                       data = stats.tmp
        )
        stats.results.list[[j]] <- res.aov
      }
      stats.results.list$aov.summary.results <- as.data.frame(t(sapply(stats.results.list[1:length(stats.results.list)-1], 
                                                                       function(i) unlist(summary(i)))))
      stats.results.list$aov.summary.results$sig.diff <- stats.results.list$aov.summary.results$`Pr(>F)1` <= 0.05
      stats.results.list
    })
    aov.results
  }else if(stats.test == "dunnett"){
    dunnett.results <- lapply(dat.split, function(i){
      stats.results.list <- setNames(vector("list", length = (length(grep("cluster", colnames(i)))+1)), 
                                     nm = c(grep("cluster", colnames(i), value = T), "dunnett.summary.results"))
      for(j in grep("cluster", colnames(i), value = T)){
        stats.tmp <- i[, c(j, "phenotype", "replicate")]
        res.DunnettTest <- DescTools::DunnettTest(x = stats.tmp[[j]],
                                                  g = stats.tmp[[grouping.name]],
                                                  control = 'DN'
        )
        stats.results.list[[j]] <- res.DunnettTest
      }
      dunnett.summary.results <- do.call("rbind", lapply(stats.results.list[1:(length(stats.results.list)-1)], function(i) i[["DN"]]))
      dunnett.summary.table <- matrix(c(dunnett.summary.results, 
                                        rep(1:(length(stats.results.list)-1), each = 3),
                                        as.integer(dunnett.summary.results[, "pval"] <= 0.05)), 
                                      ncol = 6, 
                                      dimnames = list(NULL, 
                                                      c(colnames(dunnett.summary.results), "cluster", "sig.diff"))
      )
      rownames(dunnett.summary.table) <- rownames(dunnett.summary.results)
      stats.results.list$dunnett.summary.results <- dunnett.summary.table
      stats.results.list
    })
    dunnett.results
  }else if(stats.test == "wilcox"){
   wilcox.results <- lapply(dat.split, function(i){
     stats.results.list <- setNames(vector("list", length = (length(grep("cluster", colnames(i)))+1)), 
                                    nm = c(grep("cluster", colnames(i), value = T), "wilcox.summary.results"))
     for(j in grep("cluster", colnames(i), value = T)){
       stats.tmp <- i[, c(j, "phenotype", "replicate")]
       res.wilcox <- pairwise.wilcox.test(x = stats.tmp[[j]],
                                          g = stats.tmp[[grouping.name]],
                                          p.adjust.method = "BH",
                                          paired = TRUE
       )
       stats.results.list[[j]] <- res.wilcox
     }
   })
  }
}
##
##