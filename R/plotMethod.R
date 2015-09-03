

plot.bootnet <- function(
  x, # bootnet object,
  statistics = c("strength", "closeness", "betweenness"),
  plot = c("area","line","interval"),
  sampleColor = "darkred",
  samplelwd = 1,
  bootColor = "black",
  bootAlpha = 0.01,
  bootlwd = 0.9,
  areaAlpha = 0.2,
  order = c("id","sample","mean"),
  decreasing = TRUE,
  perNode = FALSE,
  quantile = 2.5,
  legendNcol = 2, # Only for perNode plots.
  ...
){
  plot <- match.arg(plot)
  order <- match.arg(order)
  
  if (!quantile %in% c(2.5,1,5,25,50)){
    stop("Only quantiles 1, 2.5, 5, 25 and 50 are supported.")
  }
  
  ### Nodewise plots:
  if (x$type == "node"){
    # Summarize:
    Sum <- summary(x, statistic=statistics,perNode=perNode)
    minArea <- paste0("q",quantile)
    maxArea <- paste0("q",100-quantile)
    
    if (plot == "area"){
      
      if (perNode){
        g <- ggplot(Sum, aes_string(x = 'nNode', y = 'mean', group = 'id', colour = 'id',ymin = minArea, ymax = maxArea, fill = "id")) + 
          facet_grid(type ~ ., scales = "free") +
          geom_ribbon(colour = NA, alpha = areaAlpha) +
          geom_line( lwd = samplelwd) + geom_point() +
          theme_bw() + 
          xlab("# of nodes sampled") + ylab("") + 
          scale_x_reverse() + 
          guides(fill=guide_legend(ncol=legendNcol),colour=guide_legend(ncol=legendNcol)) 
        
        return(g)
        
      } else {
        g <- ggplot(Sum, aes_string(x = 'nNode', y = 'mean', group = 'type', colour = 'type',ymin = minArea, ymax = maxArea, fill = "type")) + 
          geom_ribbon(colour = NA, alpha = areaAlpha) +
          geom_line( lwd = samplelwd) + geom_point() +
          theme_bw() + 
          xlab("# of nodes sampled") + ylab("Average correlation with original sample") + 
          scale_x_reverse() + 
          ylim(0,1)
        
        return(g)
      } 
    } else if (plot == "interval") {
      
      if (perNode){
        
        g <- ggplot(Sum, aes_string(x = 'nNode', y = 'mean', group = 'id', colour = 'id',ymin = minArea, ymax = maxArea, fill = "id")) + 
          facet_grid(type ~ ., scales = "free") +
          geom_errorbar(position =  position_dodge(width = 0.4)) +
          geom_point(position =  position_dodge(width = 0.4)) +
          theme_bw() + 
          xlab("# of nodes sampled") + ylab("") + 
          scale_x_reverse() + 
          guides(fill=guide_legend(ncol=legendNcol),colour=guide_legend(ncol=legendNcol)) 
        
        
        return(g)
        
      } else {
        g <- ggplot(Sum, aes_string(x = 'nNode', y = 'mean', group = 'type', colour = 'type',ymin = minArea, ymax = maxArea, fill = "type")) + 
          geom_errorbar(position =  position_dodge(width = 0.4)) +
          geom_point(position =  position_dodge(width = 0.4)) +
          theme_bw() + 
          xlab("# of nodes sampled") + ylab("Average correlation with original sample") + 
          scale_x_reverse() + 
          ylim(0,1)
        
        return(g)
      }
      
      
    } else {
      stop("'line' plot not supported for node-wise bootstraps")
    }
  }
  
  # Start plot:
  if (plot[[1]]=="line"){
    sampleTable <- x[['sampleTable']] %>% dplyr::filter_(~type %in% statistics) %>% dplyr::mutate_(type = ~factor(type, levels = statistics))
    bootTable <- x[['bootTable']] %>% dplyr::filter_(~type %in% statistics) %>% dplyr::mutate_(type = ~factor(type, levels = statistics))
    
    ### Ordering:
    if (order[[1]]=="id"){
      sampleTable$order <- match(as.character(sampleTable$id),gtools::mixedsort(as.character(sampleTable$id)))
    } else if (order[[1]]%in%c("sample","mean")){
      # Summarize first:
      summary <- sampleTable %>% dplyr::group_by_(~id) %>% dplyr::summarize_(mean = ~mean(mean), value = ~value[type==statistics[[1]]])
      if (order[[1]] == "sample"){
        summary$order <- order(order(summary$value,summary$mean))
      } else {
        summary$order <- dplyr::min_rank(summary$mean)
      }
      sampleTable <- sampleTable %>% dplyr::left_join(dplyr::select_(summary,~id,~order), by = "id")
    } else stop(paste("'order'",order[[1]],"Not supported"))
    
    if (!decreasing){
      sampleTable$order <- -sampleTable$order
    }
    
    bootTable <- bootTable %>% dplyr::left_join(sampleTable %>% dplyr::select_(~order,~id), by = "id") %>%
      dplyr::arrange_(~dplyr::row_number(order))  %>%
      dplyr::mutate_(id = ~gsub("^(E|N): ","",as.character(id))) %>%
      dplyr::mutate_(
        id = ~factor(id, levels = unique(id))
      )
    
    sampleTable <- sampleTable %>%
      dplyr::arrange_(~dplyr::row_number(order))  %>%
      dplyr::mutate_(id = ~gsub("^(E|N): ","",as.character(id))) %>%
      dplyr::mutate_(
        id = ~factor(id, levels = unique(id))
      )
    
    g <- ggplot(bootTable, aes_string(x = 'value', y = 'id', group = 'name')) + 
      geom_path(alpha = bootAlpha, lwd = bootlwd) +
      geom_path(data = sampleTable, alpha=1, color = sampleColor, lwd = samplelwd) +
      facet_grid(~ type, scales = "free") +
      theme_bw() + 
      xlab("") +
      ylab("")
    
    return(g)
  } else if (plot[[1]] %in% c("interval","area")){
    # Compute summary stats:
    sumTable <- summary(x, statistics = statistics)  %>% ungroup %>% dplyr::mutate_(type = ~factor(type, levels = statistics))
    
    ### Ordering:
    if (order[[1]]=="id"){
      sumTable$order <- match(as.character(sumTable$id),gtools::mixedsort(as.character(sumTable$id)))
    } else if (order[[1]]%in%c("sample","mean")){
      # Summarize first:
      summary <- sumTable %>% dplyr::group_by_(~id) %>% dplyr::summarize_(sample = ~sample[type==statistics[[1]]], mean = ~mean(mean))
      if (order[[1]]=="sample"){
        summary$order <- order(order(summary$sample,summary$mean))
      } else {
        summary$order <- dplyr::min_rank(summary$mean)
      }
      sumTable <- sumTable %>% dplyr::left_join(dplyr::select_(summary,~id,~order), by = "id")
    } else stop(paste("'order'",order[[1]],"Not supported"))
    
    if (!decreasing){
      sumTable$order <- -sumTable$order
    }
    
    # Reorder:
    sumTable <- sumTable %>%
      dplyr::arrange_(~dplyr::row_number(order))  %>%
      dplyr::mutate_(id = ~gsub("^(E|N): ","",as.character(id))) %>%
      dplyr::mutate_(
        id = ~factor(id, levels = unique(id))
      )
    
    
    # Some fancy transformation:
    revTable <- function(x) x[nrow(x):1,]
    
    sumTable2 <- dplyr::rbind_list(
      sumTable %>% select_(~type,~id,~node1,~node2,~sample,ci = ~q2.5),
      revTable(sumTable %>% select_(~type,~id,~node1,~node2,~sample,ci = ~q97.5))
    )
    
    
    #     sumTable <- sumTable[gtools::mixedorder(sumTable$id),] 
    #     sumTable$id <- factor(gsub("^(E|N): ","",as.character(sumTable$id)), levels = gsub("^(E|N): ","",unique(gtools::mixedsort(as.character(sumTable$id)))))

    if (plot == "area"){
      
      
      sumTable$numericID <- as.numeric(sumTable$id)
      sumTable2$numericID <- as.numeric(sumTable2$id)
      
      g <- ggplot(sumTable2, aes_string(x = "ci", y = "numericID")) + 
        facet_grid(~ type, scales = "free") +
        geom_polygon(fill = bootColor, colour = NA, alpha = areaAlpha) +
        geom_path(aes_string(x="sample",y="numericID"), colour = sampleColor, lwd = samplelwd, data = sumTable) +
        geom_point(aes_string(x="sample",y="numericID"), colour = sampleColor, data = sumTable) +
        theme_bw() + 
        xlab("") +
        ylab("") + 
        scale_y_continuous(breaks = seq(1:length(levels(sumTable$id))), labels = levels(sumTable$id))
    
      
      return(g)
      
    } else {
      
      g <- ggplot(sumTable2, aes_string(x='sample', y='id', group = 'id')) + 
        geom_path(aes_string(x='ci'), colour = bootColor) +
        geom_point(colour = sampleColor) +
        facet_grid(~ type, scales = "free") +
        theme_bw() + 
        xlab("") +
        ylab("")
      
      return(g)
      
    }    
  } else stop("Unsupported plot")
  
}