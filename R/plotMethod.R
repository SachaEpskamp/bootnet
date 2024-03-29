

plot.bootnet <- function(
  x, # bootnet object,
  statistics, # "edge" for normal bootstrap, c("strength","closeness","betweenness") for node and person
  plot, # Two types: difference and area. Show area for edges, difference for centralities
  graph,
  CIstyle = c("quantiles","SE"), #c("default","SE","quantiles"),
  rank = FALSE,
  # CIwidth = c("95%","99%","90%","75%"),
  sampleColor = "darkred",
  samplelwd = 1,
  meanColor = "black",
  meanlwd = 0.5,
  bootColor = "black",
  bootAlpha = 0.01,
  bootlwd = 0.9,
  areaAlpha = 0.1,
  order = c("id","sample","mean"),
  decreasing = TRUE,
  perNode = FALSE,
  # quantile = 2.5,
  legendNcol = 2, # Only for perNode plots.
  labels=TRUE,
  legend = TRUE,
  subsetRange = c(100,0),
  area = !perNode,
  alpha = 0.05, # Only used if plot = "difference"
  # bonferroni = FALSE,
  onlyNonZero = FALSE, # Set to TRUE to only show edges that are non-zero in the sample network
  differenceShowValue, # Show values in difference plots?
  differenceEdgeColor = TRUE, # Show blocks of edges as colors according to standard plot.
  verbose = TRUE,
  panels = TRUE, # defaults to TRUE if
  split0 = FALSE, # Set to TRUE to show # times as 0 and CI of non-zero seperate
  prop0 = ifelse(split0, TRUE, FALSE),
  prop0_cex = 1,
  prop0_alpha = 0.8,
  prop0_minAlpha = 0.25,
  subset,
  # subtitle,
  ...
){
  bonferroni <- FALSE
  # if (panels == "default"){
  #   panels <- length(statistics) > 1
  # }

  # Check for multiple graphs
  if (length(unique(x$sampleTable$graph)) > 1){
    if (missing(graph)){
      stop("Argument 'graph' can not be missing when multiple graphs have been estimated.")
    } else {
      x$sampleTable <- x$sampleTable[x$sampleTable$graph %in% graph,]
      x$bootTable <- x$bootTable[x$bootTable$graph %in% graph,]
    }
    # if (missing(subtitle)){
    #   subtitle <- graph
    # }
  }
  # if (missing(subtitle)){
  #   subtitle <- ""
  # }

  if (missing(statistics)){
    if (! x$type %in% c("person","node")){
      statistics <- "edge"
    } else {
      statistics <-  c("strength","outStrength","inStrength") #c("strength","outStrength","inStrength","closeness","betweenness")
    }
    statistics <- statistics[statistics %in% x$sampleTable$type]
  }
  if (any(statistics == "all")){
    if (! x$type %in% c("person","node")){
      statistics <- sort(unique(x$bootTable$type[x$bootTable$node2 != ""]))
    } else {
      statistics <-  sort(unique(x$bootTable$type[x$bootTable$node2 == ""]))
    }
    statistics <- statistics[statistics %in% x$sampleTable$type]
  }

  if (length(statistics) == 0){
    stop("No statistics to be plotted; use the 'statistics' argument.")
  }

  # Change first letter of statistics to lowercase:
  substr(statistics,0,1) <- tolower(substr(statistics,0,1))

  if (missing(plot)){
    if (x$type %in% c("person","node")){
      plot <- "area"
    } else {
      if (all(statistics %in% c("strength","closeness","betweenness","expectedInfluence","inStrength","outStrength","inExpectedInfluence","outExpectedInfluence","hybrid","rspbc"))){
        plot <- "difference"
      } else {
        plot <- "area"
      }
    }
  }

  if (plot == "difference" & x$type %in% c("person","node")){
    stop("plot = 'difference' is not supported for subset bootstrap.")
  }

  if (plot != "difference" & onlyNonZero){
    stop("'onlyNonZero' only supported for plot = 'difference'")
  }

  if (alpha != 0.05 & plot != "difference"){
    stop("'alpha' argument only used when plot = 'difference'")
  }
  if (isTRUE(bonferroni) & plot != "difference"){
    stop("'bonferroni' argument only used when plot = 'difference'")
  }



  # plot <- match.arg(plot)
  order <- match.arg(order)
  CIstyle <- match.arg(CIstyle)
  # CIwidth <- match.arg(CIwidth)

  if (CIstyle=="default"){
    if (rank){
      CIstyle <- "quantiles"
    } else {
      if(x$type %in% c("person","node")){
        CIstyle <- "quantiles"
      } else {
        # CIstyle <- ifelse(statistics %in% c("closeness","strength"),"SE","quantile")
        CIstyle <- "quantile"
      }
    }

  } else {
    if (x$type=="node" & any(CIstyle == "SE")){
      stop("'SE' style confidence intervals not supported for node dropping.")
      CIstyle <- "quantile"
    }
  }

  # if (! x$type %in% c("person","node")){
  #   CIstyle <- rep(CIstyle,length=length(statistics))
  # }

  if (any(statistics%in%c("strength", "closeness", "betweenness","expectedInfluence")) & any(statistics%in%c("edge","distance"))){
    stop("Plotting both centrality CIs and edge/distance CIs together is not supported.")
  }

  #   if (!quantile %in% c(2.5,1,5,25,50)){
  #     stop("Only quantiles 1, 2.5, 5, 25 and 50 are supported.")
  #   }

  ### Nodewise or personwise plots:
  if (x$type %in% c("person","node")){
    # Summarize:

    # Check if sample values have variance:
    removeStats <- numeric(0)
    for (i in seq_along(statistics)){
      if (sd(x$sampleTable$value[x$sampleTable$type == statistics[i]],na.rm=TRUE) == 0){
        warning(paste0("Statistic ",statistics[i]," does not contain any variance and is therefore not shown."))
        removeStats <- c(removeStats,i)
      }
    }
    if (length(removeStats)>0){
      statistics <- statistics[-removeStats]
    }

    if (perNode){
      x$bootTable <- rbind(x$bootTable,x$sampleTable)

    }
    Sum <- summary(x, statistic=statistics,perNode=perNode,rank=rank)

    if (x$type == "node"){
      Sum <- Sum[Sum$nNode <= (max(subsetRange)/100)*ncol(x$sample$graph),]
      Sum <- Sum[Sum$nNode >= (min(subsetRange)/100)*ncol(x$sample$graph),]
    } else {
      Sum <- Sum[Sum$nPerson <= (max(subsetRange)/100)*x$sample$nPerson,]
      Sum <- Sum[Sum$nPerson >= (min(subsetRange)/100)*x$sample$nPerson,]
    }


    if (CIstyle == "SE"){
      minArea <- "CIlower"
      maxArea <- "CIupper"
    } else {
      minArea <- "q2.5"
      maxArea <- "q97.5"
    }
    meanVar <- "mean"

    # Split0?
    if (split0){
      minArea <- paste0(minArea,"_non0")
      maxArea <- paste0(maxArea,"_non0")
      meanVar <- paste0(meanVar,"_non0")
    }

    if (plot == "area"){

      if (perNode){

        if (x$type == "node"){
          g <- ggplot(Sum, aes_string(x = 'nNode', y = meanVar, group = 'id', colour = 'id',ymin = minArea, ymax = maxArea, fill = "id"))

          if (panels){
            g <- g + facet_grid(type ~ ., scales = "free")
          }


          if (area){
            g <- g + geom_ribbon(colour = NA, alpha = areaAlpha)
          }

          g <- g +
            geom_line( lwd = samplelwd) + geom_point() +
            theme_bw() +
            xlab("Sampled nodes") + ylab("") +
            guides(fill=guide_legend(ncol=legendNcol),colour=guide_legend(ncol=legendNcol)) +
            scale_x_reverse(breaks = seq(0.9,0.1,by=-0.1) * ncol(x$sample$graph), labels=c(paste0(seq(90,10,by=-10),"%")))


        } else {

          g <- ggplot(Sum, aes_string(x = 'nPerson', y = meanVar, group = 'id', colour = 'id',ymin = minArea, ymax = maxArea, fill = "id"))

          if (isTRUE(panels)){
            g <- g + facet_grid(type ~ ., scales = "free")
          }


          if (area){
            g <- g + geom_ribbon(colour = NA, alpha = areaAlpha)
          }


          g <- g +
            geom_line( lwd = samplelwd) + geom_point() +theme_bw() +
            xlab("Sampled cases") + ylab("") +
            guides(fill=guide_legend(ncol=legendNcol),colour=guide_legend(ncol=legendNcol)) +
            scale_x_reverse(breaks = seq(0.9,0.1,by=-0.1) * x$sample$nPerson, labels=c(paste0(seq(90,10,by=-10),"%")))

        }


        if (identical(FALSE,legend)){
          g <- g + theme(legend.position = "none")
        }
        if (identical(FALSE,labels)){
          g <- g + theme(axis.text.y = element_blank())

        }
        return(g)

      } else {

        Sum <- Sum %>% filter(!is.na(mean))
        if (x$type == "node"){
          g <- ggplot(Sum, aes_string(x = 'nNode', y = meanVar, group = 'type', colour = 'type',ymin = minArea, ymax = maxArea, fill = "type"))
          if (area){
            g <- g + geom_ribbon(colour = NA, alpha = areaAlpha)
          }

          g <- g +
            geom_line( lwd = samplelwd) + geom_point() +
            theme_bw() +
            xlab("Sampled nodes") + ylab("Average correlation with original sample")+
            ylim(-1,1) + geom_hline(yintercept = 0) +
            scale_x_reverse(breaks = seq(0.9,0.1,by=-0.1) * ncol(x$sample$graph), labels=c(paste0(seq(90,10,by=-10),"%")))

        } else {

          g <- ggplot(Sum, aes_string(x = 'nPerson', y = meanVar, group = 'type', colour = 'type',ymin = minArea, ymax = maxArea, fill = "type"))
          if (area){
            g <- g + geom_ribbon(colour = NA, alpha = areaAlpha)
          }

          g <- g +
            geom_line( lwd = samplelwd) + geom_point() +
            theme_bw() +
            xlab("Sampled cases") + ylab("Average correlation with original sample")+
            ylim(-1,1) + geom_hline(yintercept = 0) +
            scale_x_reverse(breaks = seq(0.9,0.1,by=-0.1) *  x$sample$nPerson, labels=c(paste0(seq(90,10,by=-10),"%")),
                            limits = rev(range( Sum$nPerson))) + # c(0.9,0.1) *  x$sample$nPerson)
          scale_colour_discrete("") + scale_fill_discrete("") +
            theme(legend.position = "top")
        }


        if (identical(FALSE,legend)){
          g <- g + theme(legend.position = "none")
        }
        if (identical(FALSE,labels)){
          g <- g + theme(axis.text.y = element_blank())

        }
        return(g)
      }
    } else if (plot == "interval") {

      if (perNode){

        if (x$type == "node"){
          g <- ggplot(Sum, aes_string(x = 'nNode', y = meanVar, group = 'id', colour = 'id',ymin = minArea, ymax = maxArea, fill = "id")) +
            geom_errorbar(position =  position_dodge(width = 0.4)) +
            geom_point(position =  position_dodge(width = 0.4)) +
            geom_line(position =  position_dodge(width = 0.4)) +
            theme_bw() +
            xlab("Sampled nodes") + ylab("") +
            guides(fill=guide_legend(ncol=legendNcol),colour=guide_legend(ncol=legendNcol)) +
            scale_x_reverse(breaks = seq(0.9,0.1,by=-0.1) * ncol(x$sample$graph), labels=c(paste0(seq(90,10,by=-10),"%")),
                            limits = c( ncol(x$sample$graph)-1, 1))

          if (panels){
            g <- g + facet_grid(type ~ ., scales = "free")
          }
        } else {

          g <- ggplot(Sum, aes_string(x = 'nPeople', y = meanVar, group = 'id', colour = 'id',ymin = minArea, ymax = maxArea, fill = "id")) +
            geom_errorbar(position =  position_dodge(width = 0.4)) +
            geom_point(position =  position_dodge(width = 0.4)) +
            geom_line(position =  position_dodge(width = 0.4)) +
            theme_bw() +
            xlab("Sampled cases") + ylab("") +
            guides(fill=guide_legend(ncol=legendNcol),colour=guide_legend(ncol=legendNcol)) +
            scale_x_reverse(breaks = seq(0.9,0.1,by=-0.1) *  x$sample$nPerson, labels=c(paste0(seq(90,10,by=-10),"%")),
                            limits = c( ncol(x$sample$graph)-1, 1))

          if (panels){
            g <- g +   facet_grid(type ~ ., scales = "free")
          }

        }

        if (identical(FALSE,legend)){
          g <- g + theme(legend.position = "none")
        }
        if (identical(FALSE,labels)){
          g <- g + theme(axis.text.y = element_blank())

        }
        return(g)

      } else {

        if (x$type == "node"){
          g <- ggplot(Sum, aes_string(x = 'nNode', y = meanVar, group = 'type', colour = 'type',ymin = minArea, ymax = maxArea, fill = "type")) +
            geom_errorbar(position =  position_dodge(width = 0.4)) +
            geom_point(position =  position_dodge(width = 0.4)) +
            geom_line(position =  position_dodge(width = 0.4)) +
            theme_bw() +
            geom_hline(yintercept = 0) +
            xlab("Sampled nodes") + ylab("Average correlation with original sample") +
            scale_x_reverse(breaks = seq(0.9,0.1,by=-0.1) * ncol(x$sample$graph), labels=c(paste0(seq(90,10,by=-10),"%")),
                            limits = c( ncol(x$sample$graph)-1, 1)) +
            ylim(-1,1)
        } else {

          g <- ggplot(Sum, aes_string(x = 'nPeople', y = meanVar, group = 'type', colour = 'type',ymin = minArea, ymax = maxArea, fill = "type")) +
            geom_errorbar(position =  position_dodge(width = 0.4)) +
            geom_point(position =  position_dodge(width = 0.4)) +
            geom_line(position =  position_dodge(width = 0.4)) +
            theme_bw() +
            geom_hline(yintercept = 0) +
            xlab("Sampled cases") + ylab("Average correlation with original sample") +
            scale_x_reverse(breaks = seq(0.9,0.1,by=-0.1) * ncol(x$sample$graph), labels=c(paste0(seq(90,10,by=-10),"%")),
                            limits = c( ncol(x$sample$graph)-1, 1)) +
            ylim(-1,1)

        }

        if (identical(FALSE,legend)){
          g <- g + theme(legend.position = "none")
        }
        if (identical(FALSE,labels)){
          g <- g + theme(axis.text.y = element_blank())

        }
        return(g)
      }


    } else {
      stop("'line' plot not supported for node-wise bootstraps")
    }
  }



  ### DIFFERENCE PLOTS ####
  if (plot == "difference"){

    if (any(statistics %in% c("strength","betweenness","closeness","expectedInfluence",
                              "outStrength","outExpectedInfluence","inStrength","inExpectedInfluence","rspbc","hybrid", "eigenvector",
                              "bridgeStrength", "bridgeCloseness", "bridgeBetweenness",
                              "bridgeExpectedInfluence")) & any(statistics %in% c("edge","distance"))){
      stop("'difference' plot can not be made for centrality index and edge weights/distances at the same time.")
    }

    if (missing(differenceShowValue)){
      differenceShowValue <- any(statistics %in% c("strength","betweenness","closeness","expectedInfluence",
                                                   "outStrength","outExpectedInfluence","inStrength","inExpectedInfluence","rspbc","hybrid", "eigenvector",
                                                   "bridgeStrength", "bridgeCloseness", "bridgeBetweenness",
                                                   "bridgeExpectedInfluence"))
    }

    cent <- x$bootTable %>% filter(type %in% statistics) %>% dplyr::select(name,id,value,type)

    # Only non-zero in sample edges:
    if (onlyNonZero){
      if (!all(statistics == "edge")){
        stop("'onlyNonZero' only supported for statistics = 'edge'")
      }
      include <-     unique(x$sampleTable$id[x$sampleTable$type %in% statistics & x$sampleTable$value != 0])
      cent <- cent %>% filter(id %in% include)
    } else {
      include <-     unique(x$sampleTable$id[x$sampleTable$type %in% statistics])
      cent <- cent %>% filter(id %in% include)
    }

    # Set alpha level:
    if (bonferroni){
      nInclude <- length(include)
      alpha <- alpha / (nInclude*(nInclude-1)/2)
      if (verbose) message(paste0("Significance level (alpha) set to: ",format(signif(alpha,2),scientific = FALSE)))
    }

    if (verbose){
      exp <- expAlpha(alpha,length(x$boots))
      if (verbose) message(paste0("Expected significance level given number of bootstrap samples is approximately: ",format(signif(exp,2),scientific = FALSE)))
    }

    fullTable <- expand.grid(name = unique(cent$name),id1=unique(cent$id),id2=unique(cent$id),type = unique(cent$type),
                             stringsAsFactors = FALSE)

    Quantiles <- fullTable %>%
      left_join(dplyr::select(cent,name,id1=id,value1=value,type),by=c("name","id1","type")) %>%
      left_join(dplyr::select(cent,name,id2=id,value2=value,type),by=c("name","id2","type"))  %>%
      group_by(id1,id2,type) %>%
      summarize(lower = quantile(value2-value1,alpha/2),upper = quantile(value2-value1,1-alpha/2)) %>%
      mutate(contain0 = 0 >= lower & 0 <= upper)

    #bootmean:
    bootMeans <- x$bootTable %>% filter(type %in% statistics) %>% rename(id1=id) %>%
      group_by(id1,type) %>% summarize(mean = mean(value,na.rm=TRUE))

    sample <- x$sampleTable %>% filter(type %in% statistics) %>% dplyr::select(id1=id,value,type) %>%
      left_join(bootMeans,by=c("id1","type"))

    # Now for every node: minimal node equal to....
    DF <-  sample %>% group_by(type) %>%
      mutate(rank = order(order(value,mean))) %>% arrange(rank)

    DF2 <- DF %>% filter(type == statistics[[1]])

    if (onlyNonZero){
      if (!all(statistics == "edge")){
        stop("'onlyNonZero' only supported for statistics = 'edge'")
      }
      include <-     x$sampleTable$id[x$sampleTable$type %in% statistics & x$sampleTable$value != 0]
      DF2 <- DF2 %>% filter(id1 %in% include)
      DF <- DF %>% filter(id1 %in% include)
    }

    if (order == "sample"){
      levels <- DF2$id1[order(DF2$value,decreasing = !decreasing)]
    } else if (order == "mean"){
      levels <- DF2$id1[order(DF2$rank, decreasing = !decreasing)]
    } else  if (order == "id"){
      levels <- gtools::mixedsort(unique(sample$id1))
    }

    Quantiles$id1 <- factor(Quantiles$id1,levels=levels)
    Quantiles$id2 <- factor(Quantiles$id2,levels=levels)
    Quantiles$fill <- ifelse(Quantiles$id1 == Quantiles$id2, "same",
                             ifelse(Quantiles$contain0,"nonsig","sig"))
    DF$id2 <- DF$id1
    DF$id1 <- factor(DF$id1,levels=levels)
    DF$id2 <- factor(DF$id2,levels=levels)
    DF$label <- format(signif(DF$value,2),scientific = FALSE)
    DF$fill <- "same"

    lab <- statistics
    Quantiles$type <- factor(Quantiles$type, levels = statistics)
    DF$type <- factor(DF$type, levels = statistics)
    substr(lab,1,1) <- toupper(substr(lab,1,1))

    # If blocks are to be colored, change "same" to ID:
    if (differenceEdgeColor){
      Quantiles$fill <- ifelse(Quantiles$fill == "same",
                               as.character(Quantiles$id1),
                               Quantiles$fill)

      # Run qgraph to obtain colors:
      if (missing(graph)){
          samgraph <- plot(x$sample,DoNotPlot=TRUE)
      } else {
          samgraph <- plot(x$sample,DoNotPlot=TRUE,graph=graph)
      }
      Edgelist <- data.frame(
        from =  x$sample$labels[samgraph$Edgelist$from],
        to =  x$sample$labels[samgraph$Edgelist$to],
        col = samgraph$graphAttributes$Edges$color,
        stringsAsFactors=FALSE
      )

      if (any(samgraph$Edgelist$directed)){
          Edgelist$id <- paste0(Edgelist$from,"->",Edgelist$to)
      } else {
          Edgelist$id <- paste0(Edgelist$from,"--",Edgelist$to)
      }

      fullEdgelist <- data.frame(id=include,stringsAsFactors=FALSE) %>%
        dplyr::left_join(dplyr::select(Edgelist,.data[["id"]],.data[["col"]]),by = "id")
      fullEdgelist$col[is.na(fullEdgelist$col)] <- "white"

      colorValues <- fullEdgelist$col
      names(colorValues) <- fullEdgelist$id
      colorValues <- c(colorValues,"same" = "white","nonsig" = "lightgray","sig" = "black")

    } else {
      # Color values:
      colorValues <- c("same" = "white","nonsig" = "lightgray","sig" = "black")
    }

    g <- ggplot(Quantiles,aes(x=id1,y=id2,fill=fill)) +
      geom_tile(colour = 'white') + xlab("") + ylab("") +
      scale_fill_manual(values = colorValues) +
      theme(legend.position="none")

    if (panels){
      g <- g + facet_grid(~ type)
    }


    if (differenceShowValue){
      g <- g +  geom_text(data=DF,aes(label = label))
    }



    base_size <- 9
    g <- g + theme_grey(base_size = base_size) + labs(x = "",
                                                      y = "") + scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) + theme(legend.position = "none",
                                                 axis.ticks = element_blank(), axis.text.x = element_text(size = base_size *
                                                                                                            0.8, angle = 270, hjust = 0, colour = "grey50"))

    if (!labels){

      g <- g +   theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())
    }
    return(g)

    #
    #     cent <- x$bootTable %>% filter(type %in% statistics) %>% dplyr::select(name,node1,value,type)
    #     fullTable <- expand.grid(name = unique(cent$name),node1=unique(cent$node1),node2=unique(cent$node1),type = unique(cent$type),
    #                              stringsAsFactors = FALSE)
    #
    #     Quantiles <- fullTable %>%
    #       left_join(dplyr::select(cent,name,node1,value1=value,type),by=c("name","node1","type")) %>%
    #       left_join(dplyr::select(cent,name,node2=node1,value2=value,type),by=c("name","node2","type"))  %>%
    #       group_by(node1,node2,type) %>%
    #       summarize(lower = quantile(value2-value1,0.025),upper = quantile(value2-value1,0.975)) %>%
    #       mutate(contain0 = 0 >= lower & 0 <= upper)
    #
    #     #bootmean:
    #     bootMeans <- x$bootTable %>% filter(type %in% statistics) %>%
    #       group_by(node1,type) %>% summarize(mean = mean(value,na.rm=TRUE))
    #
    #     sample <- x$sampleTable %>% filter(type %in% statistics) %>% dplyr::select(node1,value,type) %>%
    #       left_join(bootMeans,by=c("node1","type"))
    #
    #     # Now for every node: minimal node equal to....
    #     DF <- Quantiles %>% left_join(dplyr::select(sample,node2=node1,value,type), by = c("node2","type")) %>%
    #       group_by(node1,type) %>%
    #       summarize(
    #         minNode = node2[contain0][which.min(value[contain0])],
    #         maxNode = node2[contain0][which.max(value[contain0])]
    #       ) %>% left_join(sample,by=c("node1","type")) %>% ungroup %>%
    #       mutate(valueMin = value[match(minNode,node1)],
    #              valueMax = value[match(maxNode,node1)],
    #              rank = order(order(value,mean))) %>% arrange(rank)
    #
    #     DF2 <- DF %>% filter(type == statistics[[1]])
    #
    #     if (order == "sample"){
    #       levels <- DF2$node1[order(DF2$value,decreasing = decreasing)]
    #     } else if (order == "mean"){
    #       levels <- DF2$node1[order(DF2$rank, decreasing=decreasing)]
    #     } else  if (order == "id"){
    #       levels <- x$sample$labels
    #     }
    #
    #     Quantiles$node1 <- factor(Quantiles$node1,levels=levels)
    #     Quantiles$node2 <- factor(Quantiles$node2,levels=levels)
    #     Quantiles$fill <- ifelse(Quantiles$node1 == Quantiles$node2, "same",
    #                              ifelse(Quantiles$contain0,"nonsig","sig"))
    #     DF$node2 <- DF$node1
    #     DF$node1 <- factor(DF$node1,levels=levels)
    #     DF$node2 <- factor(DF$node2,levels=levels)
    #     DF$label <- as.character(round(DF$value,2))
    #     DF$fill <- "same"
    #
    #     lab <- statistics
    #     Quantiles$type <- factor(Quantiles$type, levels = statistics)
    #     DF$type <- factor(DF$type, levels = statistics)
    #     substr(lab,1,1) <- toupper(substr(lab,1,1))
    #
    #     g <- ggplot(Quantiles,aes(x=node1,y=node2,fill=fill)) +
    #       geom_tile(colour = 'white') + xlab("") + ylab("") +
    #       scale_fill_manual(values = c("same" = "white","nonsig" = "lightgray","sig" = "black")) +
    #       geom_text(data=DF,aes(label = label))+ theme(legend.position="none") +
    #       facet_grid(~ type)
    #
    #
    #     base_size <- 9
    #     g <- g + theme_grey(base_size = base_size) + labs(x = "",
    #                                                       y = "") + scale_x_discrete(expand = c(0, 0)) +
    #       scale_y_discrete(expand = c(0, 0)) + theme(legend.position = "none",
    #                                                  axis.ticks = element_blank(), axis.text.x = element_text(size = base_size *
    #                                                                                                             0.8, angle = 270, hjust = 0, colour = "grey50"))
    #
    #     return(g)

  } else {

    # which to choose?
    if (CIstyle == "SE"){
      minArea <- "CIlower"
      maxArea <- "CIupper"
    } else {
      minArea <- "q2.5"
      maxArea <- "q97.5"
    }
    meanVar <- "mean"

    # Split0?
    if (split0){
      minArea <- paste0(minArea,"_non0")
      maxArea <- paste0(maxArea,"_non0")
      meanVar <- paste0(meanVar,"_non0")
    }


    if (any(statistics %in% c("strength","closeness","betweenness"))){
      warning("Bootstrapping CIs on centrality indices is NOT consistent. Interpret these plots with care.")
    }

    # Start plot:
    if (plot[[1]]=="line"){

      if (split0){
        stop("'split0' is not yet supported for plot = 'line'")
      }

      sampleTable <- x[['sampleTable']] %>% dplyr::filter(.data[['type']] %in% statistics) %>% dplyr::mutate(type = factor(.data[['type']], levels = statistics))
      bootTable <- x[['bootTable']] %>% dplyr::filter(.data[['type']] %in% statistics) %>% dplyr::mutate(type = factor(.data[['type']], levels = statistics))

      ### Ordering:
      if (order[[1]]=="id"){
        sampleTable$order <- match(as.character(sampleTable$id),gtools::mixedsort(as.character(sampleTable$id)))
      } else if (order[[1]]%in%c("sample","mean")){
        # Summarize first:
        summarySample <- sampleTable %>% dplyr::group_by(.data[['id']]) %>% dplyr::summarize( value = .data[['value']][.data[['type']]==statistics[[1]]])
        summaryBoot <- bootTable %>% dplyr::group_by(.data[['id']]) %>% dplyr::summarize( mean = mean(.data[['value']]))
        summary <- left_join(summarySample,summaryBoot,by="id")
        if (order[[1]] == "sample"){
          summary$order <- order(order(summary$value,summary$mean))
        } else {
          # summary$order <- dplyr::min_rank(summary$mean)
          summary$order <-  order(order(summary$mean))
        }
        sampleTable <- sampleTable %>% dplyr::left_join(dplyr::select(summary,.data[['id']],.data[['order']]), by = "id")
      } else stop(paste("'order'",order[[1]],"Not supported"))

      if (!decreasing){
        sampleTable$order <- -sampleTable$order
      }

      bootTable <- bootTable %>% dplyr::left_join(sampleTable %>% dplyr::select(.data[['order']],.data[['id']]), by = "id") %>%
        dplyr::arrange(dplyr::row_number(.data[['order']]))  %>%
        dplyr::mutate(id = gsub("^(E|N): ","",as.character(.data[['id']]))) %>%
        dplyr::mutate(
          id = factor(.data[['id']], levels = unique(.data[['id']]))
        )

      sampleTable <- sampleTable %>%
        dplyr::arrange(~dplyr::row_number(.data[['order']]))  %>%
        dplyr::mutate(id = gsub("^(E|N): ","",as.character(.data[['id']]))) %>%
        dplyr::mutate(
          id = factor(.data[['id']], levels = unique(.data[['id']]))
        )



      g <- ggplot(bootTable, aes_string(x = 'value', y = 'id', group = 'name')) +
        geom_path(alpha = bootAlpha, lwd = bootlwd) +
        geom_path(data = sampleTable, alpha=1, color = sampleColor, lwd = samplelwd) +
        # facet_grid(~ type, scales = "free") +
        theme_bw() +
        xlab("") +
        ylab("")

      if (panels){
        g <- g +  facet_grid(~ type, scales = "free")
      }


      if (identical(FALSE,legend)){
        g <- g + theme(legend.position = "none")
      }
      if (identical(FALSE,labels)){
        g <- g + theme(axis.text.y = element_blank())

      }

      return(g)
    } else if (plot[[1]] %in% c("interval","area")){
      # Compute summary stats:
      sumTable <- summary(x, statistics = statistics,rank=rank)  %>% ungroup %>% dplyr::mutate(type = factor(.data[['type']], levels = statistics))

      ### Ordering:
      if (order[[1]]=="id"){
        sumTable$order <- match(as.character(sumTable$id),gtools::mixedsort(as.character(sumTable$id)))
      } else if (order[[1]]%in%c("sample","mean")){
        # Summarize first:

        summary <- sumTable %>% dplyr::group_by(.data[['id']]) %>% dplyr::summarize(sample = .data[['sample']][.data[['type']]==statistics[[1]]], mean = mean(.data[[meanVar]],na.rm=TRUE))
        if (order[[1]]=="sample"){
          summary$order <- order(order(summary$sample,summary$mean))
        } else {
          summary$order <- dplyr::min_rank(summary$mean)
        }
        sumTable <- sumTable %>% dplyr::left_join(dplyr::select(summary,.data[['id']],.data[['order']]), by = "id")
      } else stop(paste("'order'",order[[1]],"Not supported"))

      if (!decreasing){
        sumTable$order <- -sumTable$order
      }

      # Reorder:
      sumTable <- sumTable %>%
        dplyr::arrange(dplyr::row_number(.data[['order']]))  %>%
        dplyr::mutate(id = gsub("^(E|N): ","",as.character(.data[['id']]))) %>%
        dplyr::mutate(
          id = factor(.data[['id']], levels = unique(.data[['id']]))
        )


      # Some fancy transformation:
      revTable <- function(x) x[nrow(x):1,]


      #     if (CIstyle == "SE"){
      #       minArea <- "CIlower"
      #       maxArea <- "CIupper"
      #     } else {
      #       minArea <- "q2.5"
      #       maxArea <- "q97.5"
      #     }


      sumTable$lbound <- sumTable[[minArea]]
      sumTable$ubound <- sumTable[[maxArea]]


      sumTable2 <- dplyr::bind_rows(
        sumTable %>% select(.data[['type']],.data[['id']],.data[['node1']],.data[['node2']],.data[['sample']],ci = .data[['lbound']], mean = .data[[meanVar]], .data[['prop0']]),
        revTable(sumTable %>% select(.data[['type']],.data[['id']],.data[['node1']],.data[['node2']],.data[['sample']],ci = .data[['ubound']], mean = .data[[meanVar]], .data[['prop0']]))
      )


      #     sumTable <- sumTable[gtools::mixedorder(sumTable$id),]
      #     sumTable$id <- factor(gsub("^(E|N): ","",as.character(sumTable$id)), levels = gsub("^(E|N): ","",unique(gtools::mixedsort(as.character(sumTable$id)))))

      if (plot == "area"){

        # Subset:
        if (!missing(subset)){
          sumTable2 <- sumTable2 %>% filter(.data[['id']] %in% subset | .data[['node1']] %in% subset | .data[['node2']] %in% subset)
          sumTable <- sumTable %>% filter(.data[['id']] %in% subset | .data[['node1']] %in% subset | .data[['node2']] %in% subset)

          sumTable2$id <- droplevels(sumTable2$id)
          sumTable$id <- droplevels(sumTable$id)

        }


        sumTable$numericID <- as.numeric(sumTable$id)
        sumTable2$numericID <- as.numeric(sumTable2$id)

        gathered_sumTable <- tidyr::gather(sumTable,"var","value",c("sample",meanVar))

        # add transparency:
        if (split0){
          gathered_sumTable$alpha <- ifelse(gathered_sumTable$var == meanVar,
                                            prop0_minAlpha + (1-prop0_minAlpha) * (1-gathered_sumTable$prop0),
                                            1)
        } else {
          gathered_sumTable$alpha <- 1
        }

        # Plot:
        g <- ggplot(gathered_sumTable, aes_string(x='value', y='numericID', colour = "var")) +
          geom_polygon(aes_string(x = "ci", y = "numericID"),fill = bootColor, colour = NA, alpha = areaAlpha, data = sumTable2) +
          geom_path(aes_string(x=meanVar,y="numericID"), colour = meanColor, lwd = meanlwd, data = sumTable) +
          geom_path(aes_string(x="sample",y="numericID"), colour = sampleColor, lwd = samplelwd, data = sumTable) +
          geom_point(aes_string(alpha = "alpha"),show.legend = c(alpha = FALSE, colour = TRUE)) +
          # geom_point(aes_string(x="mean",y="numericID"), colour = meanColor, data = sumTable) +
          # geom_point(aes_string(x="sample",y="numericID"), colour = sampleColor, data = sumTable) +
          theme_bw() +
          xlab("") +
          ylab("") +
          scale_y_continuous(breaks = seq(1:length(levels(sumTable$id))), labels = levels(sumTable$id), expand = par("pin")[1]/par("pin")[2] * 0.01 * c(1,1))+
          scale_color_manual("",values = c(meanColor,sampleColor), labels = c("Bootstrap mean","Sample")) +
          theme(legend.position = "top")
        #
        #
        # g <- ggplot(sumTable2, aes_string(x = "ci", y = "numericID")) +
        #   geom_polygon(fill = bootColor, colour = NA, alpha = areaAlpha) +
        #   geom_path(aes_string(x="mean",y="numericID"), colour = meanColor, lwd = meanlwd, data = sumTable) +
        #   geom_point(aes_string(x="mean",y="numericID"), colour = meanColor, data = sumTable) +
        #   geom_path(aes_string(x="sample",y="numericID"), colour = sampleColor, lwd = samplelwd, data = sumTable) +
        #   geom_point(aes_string(x="sample",y="numericID"), colour = sampleColor, data = sumTable) +
        #   theme_bw() +
        #   xlab("") +
        #   ylab("") +
        #   scale_y_continuous(breaks = seq(1:length(levels(sumTable$id))), labels = levels(sumTable$id))
        #
        if (isTRUE(panels)){
          g <- g + facet_grid(~ type, scales = "free")
        }

        if (identical(FALSE,legend)){
          g <- g + theme(legend.position = "none")
        }
        if (identical(FALSE,labels)){
          g <- g + theme(axis.text.y = element_blank())

        }


      } else {
        # Subset:
        if (!missing(subset)){
          sumTable <- sumTable %>% filter(.data[['id']] %in% subset | .data[['node1']] %in% subset | .data[['node2']] %in% subset)
          sumTable2 <- sumTable2 %>% filter(.data[['id']] %in% subset | .data[['node1']] %in% subset | .data[['node2']] %in% subset)

          # Relevel:
          sumTable2$id <- droplevels(sumTable2$id)
          sumTable$id <- droplevels(sumTable$id)
        }


        gathered_sumTable <- tidyr::gather(sumTable,"var","value",c("sample",meanVar))

        # add transparency:
        if (split0){
          gathered_sumTable$alpha <- ifelse(gathered_sumTable$var == meanVar,
             prop0_minAlpha + (1-prop0_minAlpha) * (1-gathered_sumTable$prop0),
             1)
          sumTable2$alpha <- prop0_minAlpha + (1-prop0_minAlpha) * (1-sumTable2$prop0)
        } else {
          gathered_sumTable$alpha <- 1
          sumTable2$alpha <- 1
        }



        g <- ggplot(gathered_sumTable, aes_string(x='value', y='id', group = 'id', colour = "var")) +
          geom_point(aes_string(alpha = 'alpha')) +
   geom_path(aes_string(y='id', group = 'id',x='ci',alpha = 'alpha'), colour = bootColor, data = sumTable2) +
          theme_bw() +
          xlab("") +
          ylab("") +
          scale_color_manual("",values = c(meanColor,sampleColor), labels = c("Bootstrap mean","Sample")) +
          theme(legend.position = "top") + scale_alpha(guide="none")

        # g <- ggplot(sumTable2, aes_string(x='sample', y='id', group = 'id')) +
        #   geom_path(aes_string(x='ci'), colour = bootColor) +
        #   geom_point(aes_string(x='mean'), colour = meanColor) +
        #   geom_point(colour = sampleColor) +
        #   theme_bw() +
        #   xlab("") +
        #   ylab("")

        if (isTRUE(panels)){
          g <- g + facet_grid(~ type, scales = "free")
        }

        if (identical(FALSE,legend)){
          g <- g + theme(legend.position = "none")
        }
        if (identical(FALSE,labels)){
          g <- g + theme(axis.text.y = element_blank())

        }


      }

      # Include boxes with prop0:

      if (prop0){

        g <- g + geom_label(aes(x=0,label=format(round(prop0, 2), nsmall = 2)), cex = prop0_cex * 2, data = sumTable,
                       label.padding = unit(0.1, "lines"),
                       label.size = 0.1, alpha = prop0_alpha,
                       colour = "black")
      }


      return(g)
    } else stop("Unsupported plot")
  }
}
