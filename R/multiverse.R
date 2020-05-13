# Multiverse test:
multiverse <- function(x, labels = FALSE){
  
  # Obtain all networks:
  networks <- lapply(x$boots,'[[','graph')
  directed <- sapply(x$boots,'[[','directed')
  labs <- lapply(x$boots,'[[','labels')
  
  nFull <- x$sampleSize
  Ns <- sapply(x$boots,'[[','nPerson')
  
  # Resorder:
  networks <- networks[order(Ns,decreasing = TRUE)]
  directed <- directed[order(Ns,decreasing = TRUE)]
  labs <- labs[order(Ns,decreasing = TRUE)]
  Ns <- Ns[order(Ns,decreasing = TRUE)]
  
  # Make dfs:
  DFs <- lapply(seq_along(networks),function(i){
    if (directed[i]){
      graph <- networks[[i]]
      ind <- matrix(TRUE,nrow(graph),ncol(graph))
      weights <- graph[ind]
      edges <- paste0(labs[[i]][row(graph)[ind]]," -> ",labs[[i]][col(graph)[ind]])
    } else {
      graph <- networks[[i]]
      ind <- lower.tri(graph,diag=FALSE)
      weights <- graph[ind]
      edges <- paste0(labs[[i]][row(graph)[ind]]," -- ",labs[[i]][col(graph)[ind]])
    }
    
    df <- data.frame(
      id = i,
      weight=weights,
      edge = edges,
      n = Ns[i]
    )
    
    return(df)
  })
  
  # Make a big DF:
  bigDF <- do.call(rbind, DFs)
  
  # Proportion of data:
  bigDF$prop <- round(100 * bigDF$n / nFull,1)
  bigDF$prop <- factor(bigDF$prop, levels = sort(unique(bigDF$prop),decreasing = TRUE),
                       labels = paste0(sort(unique(bigDF$prop),decreasing = TRUE),"%"))
  
  # Max value:
  max <- max(abs(bigDF$weight),na.rm=TRUE)
  
  # Make the plot:
  p <- ggplot(bigDF, 
              aes_string(x="id",
                  y="edge",
                  fill="weight")) + 
    geom_tile() + 
    # geom_text(col = rgb(0.3,0.3,0.3), size = 1) +
    scale_fill_gradient2("",
                         low = "#BF0000", 
                         high = "#0000D5",
                         mid="white",limits = c(-max,max)) +
    theme_bw(base_size = 12) + 
    labs(x = "",y = "") + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(legend.position = "top", axis.ticks = element_blank(), 
          axis.text.y = element_blank(),
          legend.key.size = unit(0.01, "npc"),
          legend.key.width = unit(0.1,"npc") ) 
  
  if (x$type == "person"){
    p <- p + facet_grid( ~ prop, scales = "free")
  }
  
  if (labels){
    p <- p + theme(axis.text.x = element_text(size = 3))
  }
  
  return(p)
}