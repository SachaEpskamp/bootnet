# Compute tidy table from bootnetResult object:
# Result in data frame with entries:
# original (logical)
# name
# type
# node1
# node2
# value

statTable <- function(x,
                      name,
                      alpha = 1,
                      computeCentrality = TRUE,
                      statistics = c("edge","strength","closeness","betweenness"),
                      directed = FALSE,
                      bridgeArgs = list(),
                      includeDiagonal = FALSE,...){

  # If list, table for every graph!
  if (is.list(x$graph)){
    Tables <- list()
    for (i in seq_len(length(x$graph))){
      dummyobject <- x
      dummyobject$graph <- x$graph[[i]]
      dummyobject$directed <- x$directed[[i]]
      Tables[[i]] <- statTable(dummyobject,name=name,alpha=alpha,computeCentrality = computeCentrality,statistics=statistics,directed=dummyobject$directed,includeDiagonal=includeDiagonal)
      Tables[[i]]$graph <- names(x$graph)[[i]]
    }
    return(dplyr::bind_rows(Tables))
  }

  bridgeCentralityNames <- x[['labels']]

  # Statistics can be:
  # Change first letter of statistics to lowercase:
  substr(statistics,0,1) <- tolower(substr(statistics,0,1))
  validStatistics <-  c("intercept","edge","length","distance","closeness","betweenness","strength","expectedInfluence",
                        "outStrength","outExpectedInfluence","inStrength","inExpectedInfluence","rspbc","hybrid", "eigenvector",
                        "bridgeStrength", "bridgeCloseness", "bridgeBetweenness","bridgeInDegree","bridgeOutDegree",
                        "bridgeExpectedInfluence")
  if (!all(statistics %in% validStatistics)){
    stop(paste0("'statistics' must be one of: ",paste0("'",validStatistics,"'",collapse=", ")))
  }


  type <- NULL
  value <- NULL

  stopifnot(is(x, "bootnetResult"))
  tables <- list()
  if (is.null(x[['labels']])){
    x[['labels']] <- seq_len(ncol(x[['graph']]))
  }

  # edges:
  if (!directed){
    index <- upper.tri(x[['graph']], diag=FALSE)
    ind <- which(index, arr.ind=TRUE)
  } else {
    if (!includeDiagonal){
      index <- diag(ncol(x[['graph']]))!=1
    } else {
      index <- matrix(TRUE,ncol(x[['graph']]),ncol(x[['graph']]))
    }

    ind <- which(index, arr.ind=TRUE)
  }

  # Weights matrix:
  Wmat <- qgraph::getWmat(x)

  if ("edge" %in% statistics){
    tables$edges <- tibble::as_tibble(data.frame(
      name = name,
      type = "edge",
      node1 = x[['labels']][ind[,1]],
      node2 = x[['labels']][ind[,2]],
      value = Wmat[index],
      stringsAsFactors = FALSE
    ))
  }


  if ("length" %in% statistics){
    tables$length <- tibble::as_tibble(data.frame(
      name = name,
      type = "length",
      node1 = x[['labels']][ind[,1]],
      node2 = x[['labels']][ind[,2]],
      value = abs(1/abs(Wmat[index])),
      stringsAsFactors = FALSE
    ))
  }

  # Intercepts:
  if (!is.null(x[['intercepts']])){
    tables$intercepts <- tibble::as_tibble(data.frame(
      name = name,
      type = "intercept",
      node1 = x[['labels']],
      node2 = '',
      value = x[['intercepts']],
      stringsAsFactors = FALSE
    ))
  }

  if (computeCentrality){
    # Centrality analysis:
    if (all(x[['graph']]==0)){
      cent <- list(
        OutDegree = rep(0,ncol(x[['graph']])),
        InDegree = rep(0,ncol(x[['graph']])),
        Closeness = rep(0,ncol(x[['graph']])),
        Betweenness = rep(0,ncol(x[['graph']])),
        ShortestPathLengths = matrix(Inf,ncol(x[['graph']]),ncol(x[['graph']])),
        RSPBC = rep(0,ncol(x[['graph']])),
        Hybrid = rep(0,ncol(x[['graph']])),
        step1 = rep(0,ncol(x[['graph']])),
        expectedInfluence = rep(0,ncol(x[['graph']])),
        OutExpectedInfluence = rep(0,ncol(x[['graph']])),
        InExpectedInfluence = rep(0,ncol(x[['graph']])),
        bridgeStrength= rep(0,ncol(x[['graph']])),
        bridgeCloseness= rep(0,ncol(x[['graph']])),
        bridgeBetweenness= rep(0,ncol(x[['graph']])),
        bridgeExpectedInfluence= rep(0,ncol(x[['graph']]))
      )
      bridgeCentralityNames <- x[['labels']]
    } else {
      cent <- qgraph::centrality(Wmat, alpha = alpha, all.shortest.paths = FALSE)
      # EI <- expectedInf(Wmat, step="1")
      # names(EI) <- "expectedInfluence"
      bridgecen <- c("bridgeInDegree","bridgeOutDegree","bridgeStrength", "bridgeBetweenness", "bridgeCloseness", "bridgeExpectedInfluence")
      if(any(bridgecen %in% statistics)){
        bridgeArgs <- c(list(network=Wmat), bridgeArgs)
        if(is.null(bridgeArgs$communities)){
          warning("If bridge statistics are to be bootstrapped, the communities argument should be provided")
        }

        b <- do.call(networktools::bridge, args=bridgeArgs)

        # Rename:
        rename <- function(x,from,to){
          if (from %in% x){
            x[x==from] <- to
          }
          x
        }
        names(b) <- rename(names(b), "Bridge Indegree", "bridgeInDegree")
        names(b) <- rename(names(b), "Bridge Outdegree", "bridgeOutDegree")
        names(b) <- rename(names(b), "Bridge Strength", "bridgeStrength")
        names(b) <- rename(names(b), "Bridge Betweenness", "bridgeBetweenness")
        names(b) <- rename(names(b), "Bridge Closeness", "bridgeCloseness")
        names(b) <- rename(names(b), "Bridge Expected Influence (1-step)", "bridgeExpectedInfluence")
        names(b) <- rename(names(b), "Bridge Expected Influence (2-step)", "bridgeExpectedInfluence2step")

        b$communities <- NULL
      } else {
        b <- NULL
      }
      if(!is.null(bridgeArgs$useCommunities) && bridgeArgs$useCommunities[1] != "all"){
        b <- lapply(b, function(cen){cen[bridgeArgs$communities %in% bridgeArgs$useCommunities]})
        bridgeCentralityNames <- x[['labels']][bridgeArgs$communities %in% bridgeArgs$useCommunities]
      } else {
        bridgeCentralityNames <- x[['labels']]
      }
      cent <- c(cent, b)
    }

    # strength:
    if ("strength" %in% statistics & !directed){
      tables$strength <- tibble::as_tibble(data.frame(
        name = name,
        type = "strength",
        node1 = x[['labels']],
        node2 = '',
        value = cent[['OutDegree']],
        stringsAsFactors = FALSE
      ))
    }

    if ("outStrength" %in% statistics && directed){
      tables$outStrength <- tibble::as_tibble(data.frame(
        name = name,
        type = "outStrength",
        node1 = x[['labels']],
        node2 = '',
        value = cent[['OutDegree']],
        stringsAsFactors = FALSE
      ))
    }
    if ("inStrength" %in% statistics && directed){

      tables$inStrength <- tibble::as_tibble(data.frame(
        name = name,
        type = "inStrength",
        node1 = x[['labels']],
        node2 = '',
        value = cent[['InDegree']],
        stringsAsFactors = FALSE
      ))
    }

    # closeness:
    if ("closeness" %in% statistics){
      tables$closeness <- tibble::as_tibble(data.frame(
        name = name,
        type = "closeness",
        node1 = x[['labels']],
        node2 = '',
        value = cent[['Closeness']],
        stringsAsFactors = FALSE
      ))
    }


    # betweenness:
    if ("betweenness" %in% statistics){
      tables$betweenness <- tibble::as_tibble(data.frame(
        name = name,
        type = "betweenness",
        node1 = x[['labels']],
        node2 = '',
        value = cent[['Betweenness']],
        stringsAsFactors = FALSE
      ))
    }

    if ("distance" %in% statistics){
      tables$sp <- tibble::as_tibble(data.frame(
        name = name,
        type = "distance",
        node1 = x[['labels']][ind[,1]],
        node2 = x[['labels']][ind[,2]],
        value = cent[['ShortestPathLengths']][index],
        stringsAsFactors = FALSE
      ))
    }

    if ("expectedInfluence" %in% statistics && !directed){
      tables$expectedInfluence <- tibble::as_tibble(data.frame(
        name = name,
        type = "expectedInfluence",
        node1 = x[['labels']],
        node2 = '',
        value = cent[['OutExpectedInfluence']],
        stringsAsFactors = FALSE
      ))
    }

    # randomized shortest paths betweenness centrality:
    if ("rspbc" %in% statistics){
      tryrspbc <- try({
        tables$rspbc <- tibble::as_tibble(data.frame(
          name = name,
          type = "rspbc",
          node1 = x[['labels']],
          node2 = '',
          value = as.vector(NetworkToolbox::rspbc(abs(Wmat))),
          stringsAsFactors = FALSE
        ))
      })

      if (is(tryrspbc,"try-error")){
        tables$rspbc <- tibble::as_tibble(data.frame(
          name = name,
          type = "rspbc",
          node1 = x[['labels']],
          node2 = '',
          value = NA,
          stringsAsFactors = FALSE
        ))
      }
    }

    # hybrid:
    if ("hybrid" %in% statistics){

      tryhybrid <- try({
        tables$hybrid <- tibble::as_tibble(data.frame(
          name = name,
          type = "hybrid",
          node1 = x[['labels']],
          node2 = '',
          value = as.vector(NetworkToolbox::hybrid(abs(Wmat), BC = "random")),
          stringsAsFactors = FALSE
        ))
      })

      if (is(tryhybrid,"try-error")){
        tables$rspbc <- tibble::as_tibble(data.frame(
          name = name,
          type = "hybrid",
          node1 = x[['labels']],
          node2 = '',
          value = NA,
          stringsAsFactors = FALSE
        ))
      }

    }

    # eigenvector:
    if ("eigenvector" %in% statistics){

      tryeigenvector <- try({
        tables$eigenvector <- tibble::as_tibble(data.frame(
          name = name,
          type = "eigenvector",
          node1 = x[['labels']],
          node2 = '',
          value = as.vector(NetworkToolbox::eigenvector(Wmat)),
          stringsAsFactors = FALSE
        ))
      })

      if (is(tryeigenvector,"try-error")){
        tables$eigenvector <- tibble::as_tibble(data.frame(
          name = name,
          type = "eigenvector",
          node1 = x[['labels']],
          node2 = '',
          value = NA,
          stringsAsFactors = FALSE
        ))
      }

    }

    if ("outExpectedInfluence" %in% statistics && directed){

      tables$outExpectedInfluence <- tibble::as_tibble(data.frame(
        name = name,
        type = "outExpectedInfluence",
        node1 = x[['labels']],
        node2 = '',
        value = cent[['OutExpectedInfluence']],
        stringsAsFactors = FALSE
      ))
    }
    if ("inExpectedInfluence" %in% statistics && directed){

      tables$inExpectedInfluence <- tibble::as_tibble(data.frame(
        name = name,
        type = "inExpectedInfluence",
        node1 = x[['labels']],
        node2 = '',
        value = cent[['InExpectedInfluence']],
        stringsAsFactors = FALSE
      ))
    }


    # bridgeStrength:
    if ("bridgeStrength" %in% statistics){
      tables$bridgeStrength <- tibble::as_tibble(data.frame(
        name = name,
        type = "bridgeStrength",
        node1 = bridgeCentralityNames,
        node2 = '',
        value = cent[['bridgeStrength']],
        stringsAsFactors = FALSE
      ))
    }

    # bridgeCloseness:
    if ("bridgeCloseness" %in% statistics){
      tables$bridgeCloseness <- tibble::as_tibble(data.frame(
        name = name,
        type = "bridgeCloseness",
        node1 = bridgeCentralityNames,
        node2 = '',
        value = cent[['bridgeCloseness']],
        stringsAsFactors = FALSE
      ))
    }

    # bridgeBetweenness:
    if ("bridgeBetweenness" %in% statistics){
      tables$bridgeBetweenness <- tibble::as_tibble(data.frame(
        name = name,
        type = "bridgeBetweenness",
        node1 = bridgeCentralityNames,
        node2 = '',
        value = cent[['bridgeBetweenness']],
        stringsAsFactors = FALSE
      ))
    }

    # bridgeExpectedInfluence:
    if ("bridgeExpectedInfluence" %in% statistics){
      tables$bridgeExpectedInfluence <- tibble::as_tibble(data.frame(
        name = name,
        type = "bridgeExpectedInfluence",
        node1 = bridgeCentralityNames,
        node2 = '',
        value = cent[['bridgeExpectedInfluence']],
        stringsAsFactors = FALSE
      ))
    }

  }
  #   for (i in seq_along(tables)){
  #     tables[[i]]$id <- ifelse(tables[[i]]$node2=='',paste0("N: ",tables[[i]]$node1),paste0("E: ",tables[[i]]$node1, "--", tables[[i]]$node2))
  #   }

  for (i in seq_along(tables)){
    tables[[i]]$id <- ifelse(tables[[i]]$node2=='',tables[[i]]$node1,paste0(tables[[i]]$node1, ifelse(directed,"->","--"), tables[[i]]$node2))
  }

  tab <- dplyr::bind_rows(tables)
  tab$nNode <- x$nNode
  tab$nPerson <- x$nPerson

  # Compute rank:
  tab <- tab %>% group_by(.data[['type']]) %>%
    mutate(rank_avg = rank(value,ties.method = "average"),
           rank_min = rank(value,ties.method = "min"),
           rank_max = rank(value,ties.method = "max"))


  tab$graph <- "1"

  return(tab)
}
