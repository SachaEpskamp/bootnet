plot.bootnet <- function(
  x, # bootnet object,
  types = c("strength", "closeness", "betweenness"),
  plot = c("line","interval"),
  sampleColor = "darkred",
  samplelwd = 1.1,
  bootColor = "black",
  bootAlpha = 0.01,
  bootlwd = 1,
  order = c("id","value"),
  decreasing = TRUE,
  ...
  ){
  
  # Some comment
  # Start plot:
  if (plot[[1]]=="line"){
    sampleTable <- x[['sampleTable']] %>% dplyr::filter_(~type %in% types) %>% dplyr::mutate_(type = ~factor(type, levels = types))
    bootTable <- x[['bootTable']] %>% dplyr::filter_(~type %in% types) %>% dplyr::mutate_(type = ~factor(type, levels = types))
    
    ### Ordering:
    if (order[[1]]=="id"){
      sampleTable$order <- match(as.character(sampleTable$id),gtools::mixedsort(as.character(sampleTable$id)))
    } else if (order[[1]]=="value"){
      # Summarize first:
      summary <- sampleTable %>% dplyr::group_by_(~id) %>% dplyr::summarize_(value = ~value[type==types[[1]]])
      summary$order <- dplyr::min_rank(summary$value)
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
  } else if (plot[[1]]=="interval"){
    # Compute summary stats:

    sumTable <- summary(x, types = types)  %>% ungroup %>% dplyr::mutate_(type = ~factor(type, levels = types))
    
    ### Ordering:
    if (order[[1]]=="id"){
      sumTable$order <- match(as.character(sumTable$id),gtools::mixedsort(as.character(sumTable$id)))
    } else if (order[[1]]=="value"){
      # Summarize first:
      summary <- sumTable %>% dplyr::group_by_(~id) %>% dplyr::summarize_(sample = ~sample[type==types[[1]]])
      summary$order <- dplyr::min_rank(summary$sample)
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
    sumTable <- dplyr::rbind_list(
      sumTable %>% select_(~type,~id,~node1,~node2,~sample,ci = ~q2.5),
      sumTable %>% select_(~type,~id,~node1,~node2,~sample,ci = ~q97.5)
    )
    
 
#     sumTable <- sumTable[gtools::mixedorder(sumTable$id),] 
#     sumTable$id <- factor(gsub("^(E|N): ","",as.character(sumTable$id)), levels = gsub("^(E|N): ","",unique(gtools::mixedsort(as.character(sumTable$id)))))
    
    ggplot(sumTable, aes_string(x='sample', y='id', group = 'id')) + 
      geom_path(aes_string(x='ci'), colour = bootColor) +
      geom_point(colour = sampleColor) +
      facet_grid(~ type, scales = "free") +
      theme_bw() + 
      xlab("") +
      ylab("")
    
  } else stop("Unsupported plot")

}