plot.bootnet <- function(
  x, # bootnet object,
  types = c("strength", "closeness", "betweenness"),
  plot = c("line","interval"),
  sampleColor = "darkred",
  samplelwd = 1.5,
  bootColor = "black",
  bootAlpha = 0.01,
  bootlwd = 1
  ){

  # Start plot:
  if (plot[[1]]=="line"){
    sampleTable <- x[['sampleTable']] %>% dplyr::filter(type %in% types)
    bootTable <- x[['bootTable']] %>% dplyr::filter(type %in% types)
    
    # Reorder for nice paths:
    bootTable <- bootTable[gtools::mixedorder(bootTable$id),] 
    bootTable$id <- factor(gsub("^(E|N): ","",as.character(bootTable$id)), levels = gsub("^(E|N): ","",unique(gtools::mixedsort(as.character(bootTable$id)))))
    
    sampleTable <- sampleTable[gtools::mixedorder(sampleTable$id),] 
    sampleTable$id <- factor(gsub("^(E|N): ","",as.character(sampleTable$id)), levels = gsub("^(E|N): ","",unique(gtools::mixedsort(as.character(sampleTable$id)))))
    
    
    g <- ggplot(bootTable, aes(x = value, y = id, group = name)) + 
      geom_path(alpha = bootAlpha, lwd = bootlwd) +
      geom_path(data = sampleTable, alpha=1, color = sampleColor, lwd = samplelwd) +
      facet_grid(~ type, scales = "free") +
      theme_bw() + 
      xlab("") +
      ylab("")
    
    return(g)
  } else if (plot[[1]]=="interval"){
    # Compute summary stats:
    sumTable <- summary(x, types = types)
    
    # Some fancy transformation:
    sumTable <- dplyr::rbind_list(
      sumTable %>% select(type,id,node1,node2,sample,ci = q2.5),
      sumTable %>% select(type,id,node1,node2,sample,ci = q97.5)
    )
    
    # Reorder for consistancy:
    sumTable <- sumTable[gtools::mixedorder(sumTable$id),] 
    sumTable$id <- factor(gsub("^(E|N): ","",as.character(sumTable$id)), levels = gsub("^(E|N): ","",unique(gtools::mixedsort(as.character(sumTable$id)))))
    
    ggplot(sumTable, aes(x=sample, y=id, group = id)) + 
      geom_path(aes(x=ci), colour = bootColor) +
      geom_point(colour = sampleColor) +
      facet_grid(~ type, scales = "free") +
      theme_bw() + 
      xlab("") +
      ylab("")
    
  } else stop("Unsupported plot")

}