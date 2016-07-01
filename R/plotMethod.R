

plot.bootnet <- function(
  x, # bootnet object,
  statistics, # "edge" for normal bootstrap, c("strength","closeness","betweenness") for node and person
  plot = c("area","interval","line"),
  CIstyle = c("default","SE","quantiles"),
  rank = FALSE,
  # CIwidth = c("95%","99%","90%","75%"),
  sampleColor = "darkred",
  samplelwd = 1,
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
  ...
){
  if (missing(statistics)){
    if (! x$type %in% c("person","node")){
      statistics <- "edge"
    } else {
      statistics <-  c("strength","closeness","betweenness") 
    }
    
    
  }
  plot <- match.arg(plot)
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
        CIstyle <- ifelse(statistics %in% c("closeness","strength"),"SE","quantile")
      }
    }

  } else {
    if (x$type=="node" & any(CIstyle == "SE")){
      stop("'SE' style confidence intervals not supported for node dropping.")
      CIstyle <- "quantile"
    }
  }
  
  if (! x$type %in% c("person","node")){
    CIstyle <- rep(CIstyle,length=length(statistics))
  }
  
  if (any(statistics%in%c("strength", "closeness", "betweenness")) & any(statistics%in%c("edge","distance"))){
    stop("Plotting both centrality CIs and edge/distance CIs together is not supported.")
  }
  
  #   if (!quantile %in% c(2.5,1,5,25,50)){
  #     stop("Only quantiles 1, 2.5, 5, 25 and 50 are supported.")
  #   }

  ### Nodewise or personwise plots:
  if (x$type %in% c("person","node")){
    # Summarize:
  
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
    
    if (plot == "area"){
      
      if (perNode){
      
        if (x$type == "node"){
          g <- ggplot(Sum, aes_string(x = 'nNode', y = 'mean', group = 'id', colour = 'id',ymin = minArea, ymax = maxArea, fill = "id")) + 
            facet_grid(type ~ ., scales = "free") 
          
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
          
          g <- ggplot(Sum, aes_string(x = 'nPerson', y = 'mean', group = 'id', colour = 'id',ymin = minArea, ymax = maxArea, fill = "id")) + 
            facet_grid(type ~ ., scales = "free")          
          if (area){
            g <- g + geom_ribbon(colour = NA, alpha = areaAlpha)
          }
      
          g <- g + 
            geom_line( lwd = samplelwd) + geom_point() +theme_bw() + 
            xlab("Sampled people") + ylab("") + 
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
          g <- ggplot(Sum, aes_string(x = 'nNode', y = 'mean', group = 'type', colour = 'type',ymin = minArea, ymax = maxArea, fill = "type"))         
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
 
          g <- ggplot(Sum, aes_string(x = 'nPerson', y = 'mean', group = 'type', colour = 'type',ymin = minArea, ymax = maxArea, fill = "type"))         
          if (area){
            g <- g + geom_ribbon(colour = NA, alpha = areaAlpha)
          }
          
          g <- g + 
            geom_line( lwd = samplelwd) + geom_point() +
            theme_bw() + 
            xlab("Sampled people") + ylab("Average correlation with original sample")+ 
            ylim(-1,1) + geom_hline(yintercept = 0) + 
            scale_x_reverse(breaks = seq(0.9,0.1,by=-0.1) *  x$sample$nPerson, labels=c(paste0(seq(90,10,by=-10),"%")))
          
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
          g <- ggplot(Sum, aes_string(x = 'nNode', y = 'mean', group = 'id', colour = 'id',ymin = minArea, ymax = maxArea, fill = "id")) + 
            facet_grid(type ~ ., scales = "free") +
            geom_errorbar(position =  position_dodge(width = 0.4)) +
            geom_point(position =  position_dodge(width = 0.4)) +
            geom_line(position =  position_dodge(width = 0.4)) +
            theme_bw() + 
            xlab("Sampled nodes") + ylab("") + 
            guides(fill=guide_legend(ncol=legendNcol),colour=guide_legend(ncol=legendNcol)) + 
            scale_x_reverse(breaks = seq(0.9,0.1,by=-0.1) * ncol(x$sample$graph), labels=c(paste0(seq(90,10,by=-10),"%")),
                            limits = c( ncol(x$sample$graph)-1, 1))
          
        } else {
          
          g <- ggplot(Sum, aes_string(x = 'nPeople', y = 'mean', group = 'id', colour = 'id',ymin = minArea, ymax = maxArea, fill = "id")) + 
            facet_grid(type ~ ., scales = "free") +
            geom_errorbar(position =  position_dodge(width = 0.4)) +
            geom_point(position =  position_dodge(width = 0.4)) +
            geom_line(position =  position_dodge(width = 0.4)) +
            theme_bw() + 
            xlab("Sampled people") + ylab("") + 
            guides(fill=guide_legend(ncol=legendNcol),colour=guide_legend(ncol=legendNcol)) + 
            scale_x_reverse(breaks = seq(0.9,0.1,by=-0.1) *  x$sample$nPerson, labels=c(paste0(seq(90,10,by=-10),"%")),
                            limits = c( ncol(x$sample$graph)-1, 1))
          
          
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
          g <- ggplot(Sum, aes_string(x = 'nNode', y = 'mean', group = 'type', colour = 'type',ymin = minArea, ymax = maxArea, fill = "type")) + 
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
          
          g <- ggplot(Sum, aes_string(x = 'nPeople', y = 'mean', group = 'type', colour = 'type',ymin = minArea, ymax = maxArea, fill = "type")) + 
            geom_errorbar(position =  position_dodge(width = 0.4)) +
            geom_point(position =  position_dodge(width = 0.4)) +
            geom_line(position =  position_dodge(width = 0.4)) +
            theme_bw() + 
            geom_hline(yintercept = 0) +
            xlab("Sampled people") + ylab("Average correlation with original sample") + 
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

  if (any(statistics %in% c("strength","closeness","betweenness"))){
    warning("Bootstrapping CIs on centrality indices is NOT consistent. Interpret these plots with care.")
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
      summarySample <- sampleTable %>% dplyr::group_by_(~id) %>% dplyr::summarize_( value = ~value[type==statistics[[1]]])
      summaryBoot <- bootTable %>% dplyr::group_by_(~id) %>% dplyr::summarize_( mean = ~mean(value))
      summary <- left_join(summarySample,summaryBoot,by="id")
      if (order[[1]] == "sample"){
        summary$order <- order(order(summary$value,summary$mean))
      } else {
        # summary$order <- dplyr::min_rank(summary$mean)
        summary$order <-  order(order(summary$mean))
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
    
    if (identical(FALSE,legend)){
      g <- g + theme(legend.position = "none")
    }
    if (identical(FALSE,labels)){
      g <- g + theme(axis.text.y = element_blank())
      
    }
    
    return(g)
  } else if (plot[[1]] %in% c("interval","area")){
    # Compute summary stats:
    sumTable <- summary(x, statistics = statistics,rank=rank)  %>% ungroup %>% dplyr::mutate_(type = ~factor(type, levels = statistics))
    
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
    

#     if (CIstyle == "SE"){
#       minArea <- "CIlower"
#       maxArea <- "CIupper"
#     } else {
#       minArea <- "q2.5"
#       maxArea <- "q97.5"  
#     }

    sumTable <- sumTable %>% mutate_(
      lbound = ~ifelse(CIstyle[match(type,statistics)] == "SE", CIlower,q2.5),
      ubound = ~ifelse(CIstyle[match(type,statistics)] == "SE", CIupper, q97.5)
    )
    
    
    sumTable2 <- dplyr::rbind_list(
      sumTable %>% select_(~type,~id,~node1,~node2,~sample,ci = ~lbound),
      revTable(sumTable %>% select_(~type,~id,~node1,~node2,~sample,ci = ~ubound))
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
      if (identical(FALSE,legend)){
        g <- g + theme(legend.position = "none")
      }
      if (identical(FALSE,labels)){
        g <- g + theme(axis.text.y = element_blank())
        
      }
      
      return(g)
      
    } else {
      
      g <- ggplot(sumTable2, aes_string(x='sample', y='id', group = 'id')) + 
        geom_path(aes_string(x='ci'), colour = bootColor) +
        geom_point(colour = sampleColor) +
        facet_grid(~ type, scales = "free") +
        theme_bw() + 
        xlab("") +
        ylab("")
      
      if (identical(FALSE,legend)){
        g <- g + theme(legend.position = "none")
      }
      if (identical(FALSE,labels)){
        g <- g + theme(axis.text.y = element_blank())
        
      }
      
      return(g)
      
    }    
  } else stop("Unsupported plot")
  
}