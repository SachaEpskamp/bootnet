expAlpha <- function(alpha, nBoots, reps = 1000) {
  c(sapply(alpha,function(a){
    sapply(nBoots,function(nb){
      mean(replicate(reps,quantile(runif(nb),a/2, type = 6))) + 
        (1 - mean(replicate(reps,quantile(runif(nb),1-a/2, type = 6))))
    })
  }))
}

differenceTest <- function(bootobject,x,y,measure = c("strength","closeness","betweenness"),alpha = 0.05,x2,y2, verbose=TRUE){
  
  if (!bootobject$type %in% c("nonparametric","parametric")){
    stop("Difference test requires type = 'nonparametric' or type = 'parametric'.")
  }
  
  stopifnot(class(bootobject) == "bootnet")
  
  if (verbose){
    exp <- expAlpha(alpha,length(bootobject$boots))
    message(paste0("Expected significance level given number of bootstrap samples is approximately: ",format(signif(exp,2),scientific = FALSE)))
  }
  
  if (any(measure %in% c("strength","betweenness","closeness")) & any(measure %in% c("edge","distance"))){
    stop("Difference test can not be made for centrality index and edge weights/distances at the same time.")
  }
  
  if (!missing(x2)){
    if (any(measure %in% c("edge","distance"))){
      opts <- paste0(c(x,x2),"--",c(x2,x))
      x <- opts[opts%in%bootobject$sampleTable$id][1]
    } else {
      warning("'x2' ignored")
    }
  }  
  if (!missing(y2)){
    if (any(measure %in% c("edge","distance"))){
      opts <- paste0(c(y,y2),"--",c(y2,y))
      y <- opts[opts%in%bootobject$sampleTable$id][1]
    } else {
      warning("'y2' ignored")
    }
  }  
  
  if (is.numeric(x)){
    if (any(measure %in% c("strength","betweenness","closeness"))){
      x <- bootobject$sample$labels[x]
    } else {
      stop("Numeric assignment not possible for edge or distance difference test")
    }
  }
  
  if (is.numeric(y)){
    if (any(measure %in% c("strength","betweenness","closeness"))){
      y <- bootobject$sample$labels[y]
    } else {
      stop("Numeric assignment not possible for edge or distance difference test")
    }
  }
 

  cent <- bootobject$bootTable %>% filter(type %in% measure) %>% dplyr::select(name,id1=id,value,type)
  
  if (!all(x %in% cent$id1)){
    stop("'x' is not a valid ID")
  }
  if (!all(y %in% cent$id1)){
    stop("'y' is not a valid ID")
  }
  
  fullTable <- expand.grid(name = unique(cent$name),id1=x,id2=y,type = measure,
                           stringsAsFactors = FALSE) 
  

  Quantiles <- fullTable %>% 
    dplyr::left_join(dplyr::select(cent,name,id1=id1,value1=value,type),by=c("name","id1","type")) %>% 
    dplyr::left_join(dplyr::select(cent,name,id2=id1,value2=value,type),by=c("name","id2","type"))  %>%
    dplyr::group_by(id1,id2,type) %>%
    dplyr::summarize(lower = quantile(value2-value1,alpha/2, type = 6),
              upper = quantile(value2-value1,1-alpha/2, type = 6)) %>%
    dplyr::mutate(contain0 = 0 >= lower & 0 <= upper) %>% 
    dplyr::mutate(significant = !contain0) %>%
    dplyr::select_("id1","id2","type","lower","upper","significant") %>%
    dplyr::rename(measure = type) %>% 
    as.data.frame
  
  #   Results <- list(
  #     node1 = Quantiles$node1,
  #     node2 = Quantiles$node2,
  #     measure = Quantiles$type,
  #     CIlower = Quantiles$lower,
  #     CIupper = Quantiles$upper
  #   )
  rownames(Quantiles) <- NULL
  
  
  return(Quantiles)
}

# 
# overlap <- function(x,measure = c("strength","closeness","betweenness"),
#                     order = c("value","order")){
#   
#   order <- match.arg(order)
#   measure <- match.arg(measure)
#   
#   cent <- x$bootTable %>% filter(type %in% measure) %>% dplyr::select(name,node1,value,type)
#   fullTable <- expand.grid(name = unique(cent$name),node1=unique(cent$node1),node2=unique(cent$node1),type = unique(cent$type),
#                            stringsAsFactors = FALSE) 
#   
#   Quantiles <- fullTable %>% 
#     left_join(dplyr::select(cent,name,node1,value1=value,type),by=c("name","node1","type")) %>% 
#     left_join(dplyr::select(cent,name,node2=node1,value2=value,type),by=c("name","node2","type"))  %>%
#     group_by(node1,node2,type) %>%
#     summarize(lower = quantile(value2-value1,0.025),upper = quantile(value2-value1,0.975)) %>%
#     mutate(contain0 = 0 >= lower & 0 <= upper)
#   
#   #bootmean:
#   bootMeans <- x$bootTable %>% filter(type %in% measure) %>%
#     group_by(node1,type) %>% summarize(mean = mean(value,na.rm=TRUE))
#   
#   sample <- x$sampleTable %>% filter(type %in% measure) %>% dplyr::select(node1,value,type) %>% 
#     left_join(bootMeans,by=c("node1","type"))
#   
#   # Now for every node: minimal node equal to....
#   DF <- Quantiles %>% left_join(dplyr::select(sample,node2=node1,value,type), by = c("node2","type")) %>%
#     group_by(node1,type) %>% 
#     summarize(
#       minNode = node2[contain0][which.min(value[contain0])],
#       maxNode = node2[contain0][which.max(value[contain0])]
#     ) %>% left_join(sample,by=c("node1","type")) %>% ungroup %>%
#     mutate(valueMin = value[match(minNode,node1)], 
#            valueMax = value[match(maxNode,node1)],
#            rank = order(order(value,mean))) %>% arrange(rank)
#   
#   if (order == "value"){
#     levels <- DF$node1[order(DF$value)]  
#   } else  if (order == "order"){
#     levels <- x$sample$labels
#   }
#   
#   Quantiles$node1 <- factor(Quantiles$node1,levels=levels)
#   Quantiles$node2 <- factor(Quantiles$node2,levels=levels)
#   Quantiles$fill <- ifelse(Quantiles$node1 == Quantiles$node2, "same",
#                            ifelse(Quantiles$contain0,"nonsig","sig"))
#   DF$node2 <- DF$node1
#   DF$node1 <- factor(DF$node1,levels=levels)
#   DF$node2 <- factor(DF$node2,levels=levels)
#   DF$label <- as.character(round(DF$value,2))
#   DF$fill <- "same"
#   
#   lab <- measure
#   substr(lab,1,1) <- toupper(substr(lab,1,1))
#   
#   g <- ggplot(Quantiles,aes(x=node1,y=node2,fill=fill)) + 
#     geom_tile(colour = 'white') + xlab("") + ylab("") + 
#     scale_fill_manual(values = c("same" = "white","nonsig" = "lightgray","sig" = "black")) + 
#     geom_text(data=DF,aes(label = label))+ theme(legend.position="none") + 
#     ggtitle(lab)
#   
#   base_size <- 9
#   g <- g + theme_grey(base_size = base_size) + labs(x = "",
#                                                     y = "") + scale_x_discrete(expand = c(0, 0)) +
#     scale_y_discrete(expand = c(0, 0)) + theme(legend.position = "none",
#                                                axis.ticks = element_blank(), axis.text.x = element_text(size = base_size *
#                                                                                                           0.8, angle = 330, hjust = 0, colour = "grey50"))
#   
#   return(g)
#   #   
#   #   plot(DF$value,DF$rank,type="o",xaxt="n",yaxt="n",bty="n",
#   #        xlab = measure, ylab = "node", pch = 17, cex = 2)
#   #   
#   #   for (i in seq_len(nrow(DF))){
#   #     
#   #     lines(c(DF$valueMin[i],DF$valueMax[i]),c(DF$rank[i],DF$rank[i]),lty=2)
#   #     lines(c(DF$valueMin[i],DF$valueMin[i]),c(DF$rank[i],DF$rank[DF$minNode[i]==DF$node1]),lty=2)
#   #     lines(c(DF$valueMax[i],DF$valueMax[i]),c(DF$rank[i],DF$rank[DF$maxNode[i]==DF$node1]),lty=2)
#   #     
#   #   }
#   #   
#   #   axis(1,at = round(seq(min(DF$value),max(DF$value),length=6),2))
#   #   axis(2,at = seq_len(nrow(DF)),labels = DF$node1,las=1)
#   #   
#   #   
#   #   Quantiles[Quantiles$node1 %in% c("V5","V6") & Quantiles$node2  %in% c("V5","V6"),]
#   #   
#   #   # ggplot(DF, aes(x=value,y=rank)) + 
#   #   #   geom_line() + geom_point() + 
#   #   #   theme_bw() + 
#   #   #   xlab(measure) + 
#   #   #   geom_line(aes(x=))
#   #   
#   #   invisible(NULL)
# }