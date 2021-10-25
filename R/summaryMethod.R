# Creates the summary table
summary.bootnet <- function(
  object, # bootnet object
  graph,
  statistics = c("edge", "intercept", "strength", "closeness", "betweenness","distance"), # stats to include in the table
  perNode = FALSE, # Set to true to investigate nodewise stabilty per node.
  rank = FALSE,
  tol = sqrt(.Machine$double.eps),
  ...
){
  if (length(unique(object$sampleTable$graph)) > 1 && missing(graph)){
    stop("Argument 'graph' can not be missing when multiple graphs have been estimated.")
  }
  if (!missing(graph)){
    object$sampleTable <- object$sampleTable[object$sampleTable$graph %in% graph,]
    object$bootTable <- object$bootTable[object$bootTable$graph %in% graph,]
  }
  
  
  naTo0 <- function(x){
    x[is.na(x)] <- 0
    x
  }
  
  if (rank){
    object$bootTable$value <- object$bootTable$rank_avg
    object$sampleTable$value <- object$sampleTable$rank_avg
    
    object$bootTable$value_min <- object$bootTable$rank_min
    object$sampleTable$rank_min <- object$sampleTable$rank_min
    
    object$bootTable$value_max <- object$bootTable$rank_max
    object$sampleTable$value_max <- object$sampleTable$rank_max
  } else {
    
    object$bootTable$value_min <- object$bootTable$value
    object$sampleTable$rank_min <- object$sampleTable$value
    
    object$bootTable$value_max <- object$bootTable$value
    object$sampleTable$value_max <- object$sampleTable$value
  }

  # Returns quantiles for type = "observation" and correlations with original for type = "node"
  if (!object$type %in% c("person","node")){

   
    if (object$type == "jackknife"){

      N <- object$sampleSize
      tab <- object$bootTable %>% 
        dplyr::filter(.data[['type']] %in% statistics) %>%
        dplyr::left_join(object$sampleTable %>% dplyr::select(.data[['type']],.data[['id']],.data[['node1']],.data[['node2']],sample = .data[['value']]), by=c("id","type","node1","node2")) %>%
        dplyr::mutate(PS = N*.data[['sample']] - (N-1)*.data[['value']]) %>%
        dplyr::group_by(.data[['type']], .data[['node1']], .data[['node2']], .data[['id']]) %>%
        dplyr::summarize(
          mean = mean(.data[['value']]),
          sample = mean(.data[['PS']]),
          var = (1/(N-1)) * sum((.data[['PS']] - .data[['value']])^2),
          CIlower = .data[['sample']] - 2 * sqrt(.data[['var']]/N),
          CIupper = .data[['sample']] + 2 * sqrt(.data[['var']]/N)
        )%>%
        dplyr::select(.data[['type']], .data[['id']], .data[['node1']], .data[['node2']], .data[['sample']], .data[['mean']], .data[['CIlower']], .data[['CIupper']])
      
    } else {

      tab <- object$bootTable %>% 
        dplyr::filter(.data[['type']] %in% statistics) %>%
        dplyr::group_by(.data[['type']], .data[['node1']], .data[['node2']], .data[['id']]) %>%
        dplyr::summarize(
          mean = mean(.data[['value']],na.rm=TRUE),
          var = var(.data[['value']],na.rm=TRUE),
          sd = sd(.data[['value']],na.rm=TRUE),
          prop0 = mean(abs(.data[['value']]) < tol) %>% naTo0,
          # q1 = ~quantile(value,1/100, na.rm = TRUE),
          q2.5 = quantile(.data[['value_min']], 2.5/100, na.rm = TRUE, type = 6) %>% naTo0,
          #                 q5 = ~quantile(value, 5/100, na.rm = TRUE),
          #                 q25 = ~quantile(value, 25/100, na.rm = TRUE),
          #                 q50 = ~quantile(value, 50/100, na.rm = TRUE),
          #                 q75 = ~quantile(value, 75/100, na.rm = TRUE),
          #                 q95 = ~quantile(value, 95/100, na.rm = TRUE),
          q97.5 = quantile(.data[['value_max']], 97.5/100, na.rm = TRUE, type = 6) %>% naTo0,
          q2.5_non0 = quantile(.data[['value_min']][!abs(.data[['value_min']]) < tol], 2.5/100, na.rm = TRUE, type = 6) %>% naTo0,
          mean_non0 = mean(.data[['value']][!abs(.data[['value']]) < tol],  na.rm = TRUE) %>% naTo0,
          q97.5_non0 = quantile(.data[['value_max']][!abs(.data[['value_max']]) < tol], 97.5/100, na.rm = TRUE, type = 6) %>% naTo0,
          var_non0 = var(.data[['value']][!abs(.data[['value']]) < tol],na.rm=TRUE) %>% naTo0,
          sd_non0 = sd(.data[['value']][!abs(.data[['value']]) < tol],na.rm=TRUE) %>% naTo0
          # q99 = ~quantile(value, 99/100, na.rm = TRUE)
        ) %>%
        dplyr::left_join(object$sampleTable %>% dplyr::select(.data[['type']],.data[['id']],.data[['node1']],.data[['node2']],sample = .data[['value']]), by=c("id","type","node1","node2"))   %>%
        dplyr::mutate(
          CIlower = .data[['sample']]-2*.data[['sd']], CIupper = .data[['sample']] + 2*.data[['sd']],
          CIlower_non0 = .data[['mean_non0']] - 2*.data[['sd_non0']], CIupper_non0 = .data[['mean_non0']] + 2*.data[['sd_non0']]) %>%
        dplyr::select(.data[['type']], .data[['id']], .data[['node1']], .data[['node2']], .data[['sample']], .data[['mean']], .data[['sd']], .data[['CIlower']], .data[['CIupper']],
                       .data[['q2.5']], .data[['q97.5']], .data[['q2.5_non0']],  .data[['mean_non0']],  .data[['q97.5_non0']],  .data[['var_non0']],  .data[['sd_non0']], .data[['prop0']])
    }
    
    
    
  } else {

    # Nodewise
    tab <- object$bootTable %>% 
      dplyr::filter(.data[['type']] %in% statistics) %>% 
      dplyr::left_join(object$sampleTable %>% dplyr::select(.data[['type']],.data[['id']],.data[['node1']],.data[['node2']],sample = .data[['value']]), by=c("id","type","node1","node2"))
    
    if (perNode){
      tab <- tab %>% group_by(.data[['id']], .data[['type']], .data[['nNode']], .data[['nPerson']])  %>%
        dplyr::summarize(
          mean = mean(.data[['value']],na.rm=TRUE),
          var = var(.data[['value']],na.rm=TRUE),
          sd = sd(.data[['value']],na.rm=TRUE),
          q1 = quantile(.data[['value']],1/100, na.rm = TRUE, type = 6) %>% naTo0,
          q2.5 = quantile(.data[['value']], 2.5/100, na.rm = TRUE, type = 6) %>% naTo0,
          q5 = quantile(.data[['value']], 5/100, na.rm = TRUE, type = 6) %>% naTo0,
          q25 = quantile(.data[['value']], 25/100, na.rm = TRUE, type = 6) %>% naTo0,
          q50 = quantile(.data[['value']], 50/100, na.rm = TRUE, type = 6) %>% naTo0,
          q75 = quantile(.data[['value']], 75/100, na.rm = TRUE, type = 6) %>% naTo0,
          q95 = quantile(.data[['value']], 95/100, na.rm = TRUE, type = 6) %>% naTo0,
          q97.5 = quantile(.data[['value']], 97.5/100, na.rm = TRUE, type = 6) %>% naTo0,
          q99 = quantile(.data[['value']], 99/100, na.rm = TRUE, type = 6) %>% naTo0,
          prop0 = mean(abs(.data[['value']]) < tol),
          q2.5_non0 = quantile(.data[['value']][!abs(.data[['value']]) < tol], 2.5/100, na.rm = TRUE, type = 6) %>% naTo0,
          mean_non0 = mean(.data[['value']][!abs(.data[['value']]) < tol],  na.rm = TRUE) %>% naTo0,
          q97.5_non0 = quantile(.data[['value']][!abs(.data[['value']]) < tol], 97.5/100, na.rm = TRUE, type = 6) %>% naTo0,
          var_non0 = var(.data[['value']][!abs(.data[['value']]) < tol],na.rm=TRUE) %>% naTo0,
          sd_non0 = sd(.data[['value']][!abs(.data[['value']]) < tol],na.rm=TRUE) %>% naTo0
        ) %>% mutate(
          CIlower = .data[['mean']] - 2*.data[['sd']], CIupper = .data[['mean']] + 2*.data[['sd']],
          CIlower_non0 = .data[['mean_non0']] - 2*.data[['sd_non0']], CIupper_non0 = .data[['mean_non0']] + 2*.data[['sd_non0']]
          ) %>% arrange(.data[['nNode']],.data[['nPerson']])
      
    } else {

      tab <- tab %>% group_by(.data[['name']], .data[['type']], .data[['nNode']], .data[['nPerson']])  %>%
        summarize(cor = suppressWarnings(cor(.data[['value']],sample, use = "pairwise.complete.obs"))) %>%
        dplyr::group_by(.data[['nNode']], .data[['nPerson']], .data[['type']]) %>%
        dplyr::summarize(
          mean = mean(.data[['cor']],na.rm=TRUE),
          var = var(.data[['cor']],na.rm=TRUE),
          sd = sd(.data[['cor']],na.rm=TRUE),
          q1 = quantile(.data[['cor']],1/100, na.rm = TRUE, type = 6) %>% naTo0,
          q2.5 = quantile(.data[['cor']], 2.5/100, na.rm = TRUE, type = 6) %>% naTo0,
          q5 = quantile(.data[['cor']], 5/100, na.rm = TRUE, type = 6) %>% naTo0,
          q25 = quantile(.data[['cor']], 25/100, na.rm = TRUE, type = 6) %>% naTo0,
          q50 = quantile(.data[['cor']], 50/100, na.rm = TRUE, type = 6) %>% naTo0,
          q75 = quantile(.data[['cor']], 75/100, na.rm = TRUE, type = 6) %>% naTo0,
          q95 = quantile(.data[['cor']], 95/100, na.rm = TRUE, type = 6) %>% naTo0,
          q97.5 = quantile(.data[['cor']], 97.5/100, na.rm = TRUE, type = 6) %>% naTo0,
          q99 = quantile(.data[['cor']], 99/100, na.rm = TRUE, type = 6) %>% naTo0
        ) %>% arrange(.data[['nNode']], .data[['nPerson']])
      
    }
  }
  
  
  return(tab)
}


