plot.bootnet <- function(
  x, # bootnet object,
  types = c("strength", "closeness", "betweenness"),
  sampleColor = "darkred",
  samplelwd = 1.5,
  bootColor = "black",
  bootAlpha = 0.01,
  bootlwd = 1
  ){
  sampleTable <- x[['sampleTable']] %>% dplyr::filter(type %in% types)
  bootTable <- x[['bootTable']] %>% dplyr::filter(type %in% types)
  
  # Reorder for nice paths:
  bootTable <- bootTable[gtools::mixedorder(bootTable$id),] 
  bootTable$id <- factor(gsub("^(E|N): ","",as.character(bootTable$id)), levels = gsub("^(E|N): ","",unique(gtools::mixedsort(as.character(bootTable$id)))))
  
  
  sampleTable <- sampleTable[gtools::mixedorder(sampleTable$id),] 
  sampleTable$id <- factor(gsub("^(E|N): ","",as.character(sampleTable$id)), levels = gsub("^(E|N): ","",unique(gtools::mixedsort(as.character(sampleTable$id)))))
  
  # Start plot:
  g <- ggplot(bootTable, aes(x = value, y = id, group = name)) + 
    geom_path(alpha = bootAlpha, lwd = bootlwd) +
    geom_path(data = sampleTable, alpha=1, color = sampleColor, lwd = samplelwd) +
    facet_grid(~ type, scales = "free") +
    theme_bw() + 
    xlab("") +
    ylab("")
    
  return(g)
}