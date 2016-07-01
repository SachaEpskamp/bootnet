comparePlot <- function(
  ..., # Bootnet objects
  statistics = c("strength", "closeness", "betweenness"),
  order = c("id","value"),
  decreasing = TRUE

){
  value <- NULL
  graph <- NULL
  
  order <- match.arg(order)
  
  dots <- list(...)
  
  if (is.null(names(dots))){
    names(dots) <- paste0("graph",seq_along(dots))
  }
  
  names(dots) <- ifelse(names(dots)=="",paste0("graph",seq_along(dots)), names(dots))
  
  if (!all(sapply(dots,is,"bootnet"))){
    stop("Only supply bootnet objects to comparePlot")
  }
  
#   tabs <- lapply(dots,function(x){
#     
#     sumTable <- summary(x, statistics = statistics)  %>% ungroup %>% dplyr::mutate_(type = ~factor(type, levels = statistics))
#     
#     ### Ordering:
#     if (order[[1]]=="id"){
#       sumTable$order <- match(as.character(sumTable$id),gtools::mixedsort(as.character(sumTable$id)))
#     } else if (order[[1]]=="value"){
#       # Summarize first:
#       summary <- sumTable %>% dplyr::group_by_(~id) %>% dplyr::summarize_(sample = ~sample[type==statistics[[1]]])
#       summary$order <- dplyr::min_rank(summary$sample)
#       sumTable <- sumTable %>% dplyr::left_join(dplyr::select_(summary,~id,~order), by = "id")
#     } else stop(paste("'order'",order[[1]],"Not supported"))
#     
#     if (!decreasing){
#       sumTable$order <- -sumTable$order
#     }
#     
#     # Reorder:
#     sumTable <- sumTable %>%
#       dplyr::arrange_(~dplyr::row_number(order))  %>%
#       dplyr::mutate_(id = ~gsub("^(E|N): ","",as.character(id))) %>%
#       dplyr::mutate_(
#         id = ~factor(id, levels = unique(id))
#       )
#     
#     
#     # Some fancy transformation:
#     revTable <- function(x) x[nrow(x):1,]
#     
#     sumTable2 <- dplyr::rbind_list(
#       sumTable %>% select_(~type,~id,~node1,~node2,~sample,ci = ~q2.5),
#       revTable(sumTable %>% select_(~type,~id,~node1,~node2,~sample,ci = ~q97.5))
#     )
#     return(sumTable2)
#   })
  
  


  # Add graph variables:
  for (i in seq_along(dots)){
    dots[[i]]$sampleTable$graph <- names(dots)[[i]]
    dots[[i]]$bootTable$graph <- names(dots)[[i]]
  }

  SampTab <- dplyr::rbind_all(lapply(dots,'[[','sampleTable')) %>%
        dplyr::filter_(~type %in% statistics)
  BootTab <- dplyr::rbind_all(lapply(dots,'[[','bootTable')) %>%
      dplyr::filter_(~type %in% statistics)

  if (order == "value"){

    sum <- BootTab %>% filter_(~type == statistics[[1]]) %>%
      group_by_(~id) %>%
      summarize_(mean = ~mean(value)) %>%
      arrange_(~mean)
    
    levs <- sum$id
    if (decreasing){
      levs <- rev(levs)
    }
    
    BootTab$id <- factor(BootTab$id, levels= levs)
    
    
  }
  
  g <- ggplot(BootTab, aes(x=id, y=value, fill = graph)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.5) + 
  facet_grid(type ~ ., scales = "free") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
  
  return(g)
}