##
##function to 'melt' input data (data.frame) by select markers ('plot.dims'); faceted plot of UMAP colored by level of marker expression
umap.markers.plot <- function(input.data, plot.dims){
  
  dat.melt <- reshape2::melt(input.data[, plot.dims],
                             id = paste("umap", 1:2, sep = "."),
                             variable.name = "Marker",
                             value.name = "Expression")
  require(ggplot2)
  ggplot(dat.melt, aes(x = umap.1, y = umap.2, z = Expression, color = Expression)) +
    stat_summary_hex(bins = 200) +
    theme(text=element_text(family="Arial",face='bold'))+
    theme(panel.border = element_blank())+
    theme(panel.background = element_blank())+
    viridis::scale_fill_viridis(option = "plasma") +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())+
    theme(legend.text = element_text(colour="black", size=9, face="bold"))+
    theme(legend.title = element_text(face = "bold",size=9))+
    facet_wrap(~Marker, ncol = 2)+
    theme(strip.text.x = element_text(size = 9, color = "black", face = "bold"),
          strip.background = element_blank())
}
##
##function to generate a faceted ggplot (using a factored column) showing count data; includes a value 'squish.max.val' to cap the top end of the color scale
umap.count.plot <- function(input.data, factor.for.facet, squish.max.val = NULL, dummy.facet = FALSE){
  
  require(ggplot2)
  require(scales)
  
  if(dummy.facet == TRUE){
    input.data[[factor.for.facet]] <- factor(input.data[[factor.for.facet]], levels = c(levels(input.data[[factor.for.facet]]), "blank"))
  }

  p <- ggplot(input.data, aes(umap.1, umap.2)) + 
    geom_hex(bins = 75) + 
    facet_wrap(factor.for.facet, ncol = 2, drop = F) +
    theme_dark() + 
    theme(panel.border = element_blank()) +
    theme(panel.background = element_blank()) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())+
    theme(legend.text = element_text(colour="black", size=9, face="bold"))+
    theme(legend.title = element_text(face = "bold"))+
    theme(strip.text.x = element_text(size = 9, color = "black", face = "bold"))+
    theme(strip.background = element_blank())
  
  if(is.null(squish.max.val)){
    p + viridis::scale_fill_viridis(option = "plasma")
  }else{
    
    if((squish.max.val/100)%%1 == 0 & (squish.max.val/100) >= 3){
      break.steps <- (squish.max.val/100)
      break.labels <- c(0, cumsum(rep(100, break.steps)))
      break.labels <- c(break.labels[1:(length(break.labels)-1)], paste0(break.labels[length(break.labels)], "+"))
    }else if((squish.max.val/50)%%1 == 0 & (squish.max.val/50) >= 4){
      break.steps <- (squish.max.val/50)
      break.labels <- c(0, cumsum(rep(50, break.steps)))
      break.labels <- c(break.labels[1:(length(break.labels)-1)], paste0(break.labels[length(break.labels)], "+"))
    }
    
    p + viridis::scale_fill_viridis(option = "plasma", limits = c(0, squish.max.val), oob = squish,
                                    labels = break.labels)
  }
}
##
##function to 'melt' FlowSOM input data; faceted plot of a UMAP (FlowSOM 'codes'), colored by level of marker expression or cluster with circle/point size related to node counts
umap.fsom.codes.plot <- function(fsom, range.scale = TRUE, view = c("marker", "cluster"), circle.size = 5){
  
  set.seed(fsom$seed.val)#initialize random seed for generating a UMAP of the codes/SOMs
  umap.codes <- matrix(uwot::umap(X = fsom$map$codes, 
                                  n_threads = parallel::detectCores()-1, 
                                  n_neighbors = 20,
                                  spread = 1,
                                  min_dist = 0.05,
                                  repulsion_strength = 0.1,
                                  verbose = F), 
                       ncol = 2, 
                       dimnames = list(NULL, c("umap.1", "umap.2"))
  )
  
  #node.counts <- tabulate(fsom$map$mapping[, 1])
  #cluster.counts <- tabulate(fsom$metaclustering[fsom$map$mapping[, 1]])
  #cluster.nodes_per <- tabulate(fsom$metaclustering)
  if(range.scale == TRUE){
    codes <- range.scaled(fsom$map$codes)
  }else{
    codes <- fsom$map$codes
  }
  
  dat <- data.frame(codes,
                    umap.codes, 
                    node.counts = tabulate(fsom$map$mapping[, 1]), 
                    cluster = fsom$metaclustering)
  dat.melt <- reshape2::melt(dat, 
                             id.vars = c("umap.1", "umap.2", "node.counts", "cluster"), 
                             value.name = "Expression", 
                             variable.name = "Marker")
  
  if(view == "marker"){
    ggplot(dat.melt, aes(umap.1, umap.2, color = Expression, fill = Expression)) + 
      geom_point(shape = 21, color = "black", stroke = 0.25, alpha = 0.7, aes(size = node.counts), show.legend = T) +
      viridis::scale_fill_viridis(option = "plasma") +
      scale_size_area(max_size = circle.size) +
      theme_dark() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_blank(),
            strip.text = element_text(size = 9)) +
      guides(size = F) +
      facet_wrap(~Marker)
  }else if(view == "cluster"){
    ggplot(dat.melt, aes(umap.1, umap.2, color = cluster, fill = cluster)) + 
      geom_polygon(data = plyr::ddply(dat, "cluster",
                                      function(dat) dat[chull(dat$umap.1, dat$umap.2), ]),
                   alpha = 0.4, linetype = 0,
                   show.legend = T) +
      scale_size_area(max_size = circle.size) +
      geom_point(shape = 21, color = "black", stroke = 0.5, alpha = 0.7, aes(size = node.counts), show.legend = T) +
      theme_void() +
      theme(legend.position = 'none',
            legend.text = element_text(size = 9))
  }
}
##
##function to plot existing UMAP embedding (all cells), colored and faceted by Flowcluster
# umap.cells.cluster <- function(input.umap.cells, fsom)
# dat <- data.frame(dat.mdat$dat.mdat[, c("umap.1", "umap.2")], cluster = fsom$metaclustering[fsom$map$mapping[, 1]])[sample(1:length(fsom$map$mapping[, 1]), 100000), ]
# 
# umap.cells.clusters <- ggplot(data = dat, aes(umap.1, umap.2, fill = cluster)) +
#   geom_hex(data = dat[, c("umap.1", "umap.2")], fill = "gray", bins = 100, alpha = 0.5) +
#   geom_hex(bins = 100) +
#   theme_void() +
#   facet_wrap(~cluster) +
#   guides(fill = FALSE)
