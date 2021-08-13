##
##
fset.compensate <- function(flowSet.samples){fsApply(flowSet.samples, function(frame) {
  ## compensate each 'frame'(sample) using stored compensation matrix (FCS header; keyword = $SPILL or SPILLOVER)
  print(paste0("compensating ", frame@description$`$FIL`))
  comp <- keyword(frame)[grep("SPILL", names(keyword(frame)), ignore.case = T)][[1]]
  frame_comped <- compensate(frame, comp)
  frame_comped
})
}
##
##
elgcl.pars <- function(compensated.set) {
  dat <- fsApply(compensated.set, exprs)
  channels.fluor <- as.vector(grep("FSC|SSC|Time", colnames(dat), ignore.case = T, value = T, invert = T))
  dat <- dat[, channels.fluor]
  elgcl <- estimateLogicle(new("flowFrame", dat), channels = channels.fluor)
  elgcl
}
##
##
fset.transform <- function(compensated.set, elgcl) {
  ## transform fcs data using ^ estimated logicile
  fsApply(compensated.set, function(frame) {
    print(paste0("transforming ", frame@description$`$FIL`))
    frame_trans <- transform(frame, elgcl)
    frame_trans
  })
}
##
##
markernames.flow <- function(fcs) {
  if(class(fcs) == "flowFrame") {
    colnames(fcs)[which(!is.na(fcs@parameters@data$desc))] <- markernames(fcs)
  } else if(class(fcs) == "flowSet") {
    colnames(fcs)[which(!is.na(fcs[[1]]@parameters@data$desc))] <- markernames(fcs)
  }
  colnames(fcs) <- sub("-", ".", colnames(fcs))
  colnames(fcs) <- sub(" ", ".", colnames(fcs))
  #colnames(fcs) <- sapply(colnames(fcs), function(i) strsplit(i, " ")[[1]][1])
  fcs
}
##
##
range.scaled <- function(dat, cols, probs = c(0.2, 0.99)) {
  tmp <- as.matrix(dat[, cols])
  ranges <- matrixStats::colQuantiles(tmp, probs = probs)
  if(ncol(tmp) == 1){
    tmp <- (tmp - ranges[1])/(ranges[2] - ranges[1])
  }else{
    tmp <- t((t(tmp) - ranges[, 1]) / (ranges[, 2] - ranges[, 1]))
  }
  tmp[tmp < 0] <- 0
  tmp[tmp > 1] <- 1
  dat[, cols] <- tmp
  dat
}
##
##
PlotStars.edit <- function (median.vals, fsom, markers = fsom$map$colsUsed, view = "MST", colorPalette = grDevices::colorRampPalette(c("#00007F", 
                                                                                                                                       "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", 
                                                                                                                                       "red", "#7F0000")), starBg = "white", backgroundValues = NULL, 
                            backgroundColor = function(n) {
                              grDevices::rainbow(n, alpha = 0.3)
                            }, backgroundLim = NULL, backgroundBreaks = NULL, backgroundSize = NULL, 
                            thresholds = NULL, legend = TRUE, query = NULL, main = "") 
{
  add.vertex.shape("star", clip = igraph.shape.noclip, plot = FlowSOM:::mystar, 
                   parameters = list(vertex.data = NULL, vertex.cP = colorPalette, 
                                     vertex.scale = TRUE, vertex.bg = starBg))
  if (is.null(thresholds)) {
    data <- median.vals[, markers, drop = FALSE]
    scale <- TRUE
  }
  else {
    if (fsom$transform) {
      warning("Thresholds should be given in the transformed space")
    }
    if (!is.null(fsom$scaled.center)) {
      thresholds <- scale(t(thresholds), center = fsom$scaled.center[markers], 
                          scale = fsom$scaled.scale[markers])
    }
    data <- t(sapply(seq_len(fsom$map$nNodes), function(i) {
      res = NULL
      for (m in seq_along(markers)) {
        res = c(res, sum(subset(fsom$data, fsom$map$mapping[, 
                                                            1] == i)[, markers[m]] > thresholds[m])/sum(fsom$map$mapping[, 
                                                                                                                         1] == i))
      }
      res
    }))
    scale <- FALSE
  }
  switch(view, MST = {
    layout <- fsom$MST$l
    lty <- 1
  }, grid = {
    layout <- as.matrix(fsom$map$grid)
    lty <- 0
  }, tSNE = {
    layout <- fsom$MST$l2
    lty <- 0
  }, stop("The view should be MST, grid or tSNE. tSNE will only work\n                   if you specified this when building the MST."))
  if (!is.null(backgroundValues)) {
    background <- computeBackgroundColor(backgroundValues, 
                                         backgroundColor, backgroundLim, backgroundBreaks)
    if (is.null(backgroundSize)) {
      backgroundSize <- fsom$MST$size
      backgroundSize[backgroundSize == 0] <- 3
    }
  }
  else {
    background <- NULL
  }
  oldpar <- graphics::par(no.readonly = TRUE)
  graphics::par(mar = c(1, 1, 1, 1))
  if (legend) {
    if (!is.null(backgroundValues)) {
      graphics::layout(matrix(c(1, 3, 2, 3), 2, 2, byrow = TRUE), 
                       widths = c(1, 2), heights = c(1))
    }
    else {
      graphics::layout(matrix(c(1, 2), 1, 2, byrow = TRUE), 
                       widths = c(1, 2), heights = c(1))
    }
    if (is.null(query)) {
      FlowSOM:::plotStarLegend(fsom$prettyColnames[markers], colorPalette(ncol(data)))
    }
    else {
      plotStarQuery(fsom$prettyColnames[markers], values = query == 
                      "high", colorPalette(ncol(data)))
    }
    if (!is.null(backgroundValues)) {
      PlotBackgroundLegend(backgroundValues, background)
    }
  }
  igraph::plot.igraph(fsom$MST$g, vertex.shape = "star", vertex.label = NA, 
                      vertex.size = fsom$MST$size, vertex.data = data, vertex.cP = colorPalette(ncol(data)), 
                      vertex.scale = scale, layout = layout, edge.lty = lty, 
                      mark.groups = background$groups, mark.col = background$col[background$values], 
                      mark.border = background$col[background$values], mark.expand = backgroundSize, 
                      main = main)
  graphics::par(oldpar)
  graphics::layout(1)
}
##
##
fsom.dat <- function(fsom) {
  if("node|cluster" %in% colnames(fsom$data)){
    fsom$data[, -grep("node|cluster", colnames(fsom$data))]
  }
  dat <- data.frame(fsom$data,
                    node = as.factor(fsom$map$mapping[, 1]), 
                    cluster = as.factor(fsom$metaclustering[fsom$map$mapping[, 1]])
  )
  names(dat)[1:length(fsom$markers)] <- fsom$markers
  dat
}
##
##
fsom.ridges <- function(fsom, cluster.select = NULL) {
  dat <- reshape2::melt(fsom.dat(fsom)[, c(fsom$dims.used, "cluster")],
                        id.vars = "cluster", value.name = "expression",
                        variable.name = "marker")
  dat$reference <- "no"
  dat.ref <- dat
  dat.ref$cluster <- "reference"
  dat.ref$reference <- "yes"
  
  if(!is.null(cluster.select)) {
    dat_plot <- rbind(dat[grep(paste0("^", cluster.select, "$", collapse = "|"), dat$cluster), ], dat.ref)
  } else {
    dat_plot <- rbind(dat, dat.ref)
  }
  
  ggplot() +
    ggridges::geom_density_ridges(data = dat_plot, aes(x = expression, y = cluster,
                                             color = reference, fill = reference), alpha = 0.3) +
    facet_wrap( ~ marker, scales = "free_x", nrow = 3) +
    ggridges::theme_ridges() +
    theme(axis.text = element_text(size = 7),
          strip.text = element_text(size = 7), legend.position = "none")
}
##
##
fsom.explorer <- function(input.fsom, name.X1, name.Y1, name.X2, name.Y2, name.X3, name.Y3, points, ...) {
  
  dat <- fsom.dat(input.fsom)
  colnames(dat)[grep("cluster", colnames(dat), ignore.case = T)] <- "cluster"
  colnames(dat)[grep("node|nodes", colnames(dat), ignore.case = T)] <- "node"
  
  heatmap <- fsom.heatmap(input.fsom, dat.type = "codes", fsom$dims.used, ...)
  
  shinyApp(
    ui = fluidPage(
      
      titlePanel("FlowSOM Explorer"),
      
      sidebarLayout(
        sidebarPanel(
          
          sliderInput('sampleSize', 'Sample Size', min=0, max=1000000,
                      value=points),
          
          selectInput('x', 'X.1', names(dat), name.X1),
          selectInput('y', 'Y.1', names(dat), name.Y1),
          
          selectInput('x2', 'X.2', names(dat), name.X2),
          selectInput('y2', 'Y.2', names(dat), name.Y2),
          
          selectInput('x3', 'X.3', names(dat), name.X3),
          selectInput('y3', 'Y.3', names(dat), name.Y3),
          
          selectInput('color', 'Color', c('None', names(dat))),
          
          if(length(grep("node", colnames(dat), ignore.case = T)) == 1) {
            selectInput('node', 'node', c('None', sort(unique(dat[, grep("node", colnames(dat), ignore.case = T)]))))
          } else {
            selectInput('node', 'node', c('None'))
          },
          
          if(length(grep("cluster", colnames(dat), ignore.case = T)) == 1) {
            selectInput('cluster', 'cluster', c('None', sort(unique(dat[, grep("cluster", colnames(dat), ignore.case = T)]))))
          } else {
            selectInput('cluster', 'cluster', c('None'))
          }
        ),
        
        mainPanel(
          fluidRow(
            splitLayout(cellWidths = c("50%", "50%"), plotOutput('fcs.plot1'), plotOutput('heatmap'))
          ),
          fluidRow(
            splitLayout(cellWidths = c("50%", "50%"), plotOutput('fcs.plot2'), plotOutput('fcs.plot3'))
          )
        )
      )),
    server = function(input, output) {
      
      dataset <- reactive({
        dat[sample(nrow(dat), input$sampleSize), ]
      })
      
      output$fcs.plot1 <- renderPlot({
        
        p <- ggplot(dataset(), aes_string(x=input$x, y=input$y)) + 
          geom_point(size = .5) +
          theme(axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18))
        
        if (input$color != 'None')
          p <- p + aes_string(color=input$color) + scale_color_viridis(option = "plasma")
        
        if (input$cluster != 'None')
          p <- p + geom_point(data = subset(dataset(), cluster == input$cluster), shape = 21, size = 3, stroke = 1, color = "red")
        
        if (input$node != 'None')
          p <- p + geom_point(data = subset(dataset(), node == input$node), shape = 21, size = 3, stroke = 1, color = "blue")
        
        print(p)
        
      })
      
      output$fcs.plot2 <- renderPlot({
        
        p <- ggplot(dataset(), aes_string(x=input$x2, y=input$y2)) + 
          geom_point(size = .5) +
          theme(axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18))
        
        if (input$color != 'None')
          p <- p + aes_string(color=input$color) + scale_color_viridis(option = "plasma")
        
        if (input$cluster != 'None')
          p <- p + geom_point(data = subset(dataset(), cluster == input$cluster), shape = 21, size = 3, stroke = 1, color = "red")
        
        if (input$node != 'None')
          p <- p + geom_point(data = subset(dataset(), node == input$node), shape = 21, size = 3, stroke = 1, color = "blue")
        
        print(p)
      })
      
      output$heatmap <- renderPlot({heatmap})
      
      output$fcs.plot3 <- renderPlot({
        
        p <- ggplot(dataset(), aes_string(x=input$x3, y=input$y3)) + 
          geom_point(size = .5) +
          theme(axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18))
        
        if (input$color != 'None')
          p <- p + aes_string(color=input$color) + scale_color_viridis(option = "plasma")
        
        if (input$cluster != 'None')
          p <- p + geom_point(data = subset(dataset(), cluster == input$cluster), shape = 21, size = 3, stroke = 1, color = "red")
        
        if (input$node != 'None')
          p <- p + geom_point(data = subset(dataset(), node == input$node), shape = 21, size = 3, stroke = 1, color = "blue")
        
        print(p)
      })
      
    },
    options = list()
  )
}
##
##
fsom.heatmap <- function(fsom, 
                         dat.type, 
                         dims, 
                         row_order = NULL,
                         clust.freq = T, 
                         clust.count = F,
                         display_numbers = FALSE, 
                         color.cut = 100, 
                         color.option = "brewer", ...) {
  
  if (!color.option %in% c("brewer", "viridis")) {
    stop("Use 'color.option' = 'brewer' or 'viridis'")
  }
  
  cluster.medians <- as.data.frame(fsom.cluster.medians(fsom, dat.type = dat.type, ranged = TRUE))
  rownames(cluster.medians) <- rownames(cluster.medians) # some weird bug here, so have to 'rename' rownames
  
  ## Calculate cluster frequencies
  clustering_table <- as.numeric(table(fsom$metaclustering[fsom$map$mapping[, 1]]))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  col.index <- fsom$markers[fsom$map$colsUsed]
  cluster_rows <- hclust(dist(cluster.medians[, -grep("cluster", colnames(cluster.medians))][, col.index], method = "euclidean"), method = "average")

  # custom row ordering
  if(!is.null(row_order)){
    cluster.medians <- cluster.medians[row_order, ]
    clustering_table <- clustering_table[row_order]
    clustering_prop <- clustering_prop[row_order]
    cluster_rows <- FALSE
  }
  
  # Colors for the heatmap
  if(color.option == "brewer"){
    color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(color.cut)
  } else {
    color_heat <- viridis(color.cut)
  }
  
  # appened cluster frequencies and or total counts
  if(clust.freq&clust.count){
    labels_row <- paste0(cluster.medians$cluster, " (", clustering_prop, " %)", " (", clustering_table, ")")
  } else if(clust.freq == TRUE&clust.count == FALSE) {
    labels_row <- paste0(cluster.medians$cluster, " (", clustering_prop , " %)")
  } else if(clust.freq == FALSE&clust.count == TRUE) {
    labels_row <- paste0(cluster.medians$cluster, " (", clustering_table , ")")
  } else {
    labels_row <- cluster.medians$cluster
  }
  
  # Annotation for the original clusters
  annotation_row <- cluster.medians["cluster"]
  annotation_colors <- list(cluster = setNames(colorRampPalette(viridis(9))(nlevels(fsom$metaclustering)), levels(annotation_row$cluster)))
  
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  
  pheatmap(cluster.medians[, dims], color = plasma(50),
           cluster_rows = cluster_rows, labels_row = labels_row,
           display_numbers = display_numbers, number_color = "black",
           fontsize = 8, fontsize_number = 6,  legend_breaks = legend_breaks,
           annotation_row = annotation_row, annotation_colors = annotation_colors,
           ...)
  
}
##
##
fsom.cluster.medians <- function(fsom, dat.type = c("map", "data", "codes"), ranged = TRUE) {
  
  if(dat.type == "map") {
    if(ranged){
      tmp <- range.scaled(fsom$map$medianValues)
    } else {
      tmp <- fsom$map$medianValues
    }
    colnames(tmp) <- fsom$markers
    tmp <- data.frame(tmp, cluster = fsom$metaclustering)
    tmp
  } else if(dat.type == "data") {
    if(ranged) {
      tmp <- range.scaled(fsom$data)
    } else {
      tmp <- fsom$data
    }
    colnames(tmp) <- fsom$markers
    tmp <- data.frame(tmp, cluster = fsom$metaclustering[fsom$map$mapping[, 1]])
    tmp
  } else if(dat.type == "codes") {
    if(ranged) {
      tmp <- range.scaled(fsom$map$codes)
    } else {
      tmp <- fsom$map$codes
    }
    #colnames(tmp) <- fsom$markers
    tmp <- data.frame(tmp, cluster = fsom$metaclustering)
    tmp
  }
  
  tmp <- tmp %>% group_by(cluster) %>% summarize_all(list(median))
  tmp
}
##
##
nodes.to.meta_fSOM <- function(fSOM.object, nodes) {
  levels(fSOM.object$metaclustering) <- c(levels(fSOM.object$metaclustering), length(levels(fSOM.object$metaclustering))+1)
  fSOM.object$metaclustering[nodes] <- length(levels(fSOM.object$metaclustering))
  fSOM.object$metaclustering <- factor(fSOM.object$metaclustering)
  levels(fSOM.object$metaclustering) <- c(1:length(levels(fSOM.object$metaclustering)))
  return(fSOM.object)
}
##
##


# use this function to make row or column names bold
# parameters:
#   mat: the matrix passed to pheatmap
#   rc_fun: either rownames or colnames
#   rc_names: vector of names that should appear in boldface
make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}

##
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
