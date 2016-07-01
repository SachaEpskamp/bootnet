# Creates the summary table
summary.bootnet <- function(
  object, # bootnet object
  statistics = c("edge", "intercept", "strength", "closeness", "betweenness","distance"), # stats to include in the table
  perNode = FALSE, # Set to true to investigate nodewise stabilty per node.
  rank = FALSE,
  ...
){
  
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
#     tab <- object$bootTable %>% 
#       dplyr::filter_(~type %in% statistics) %>%
#       dplyr::group_by_(~type, ~node1, ~node2, ~id) %>%
#       dplyr::summarize_(
#         mean = ~mean(value,na.rm=TRUE),
#         var = ~var(value,na.rm=TRUE),
#         sd = ~sd(value,na.rm=TRUE),
#         prop0 = ~mean(value == 0),
#         q1 = ~quantile(value,1/100, na.rm = TRUE),
#         q2.5 = ~quantile(value, 2.5/100, na.rm = TRUE),
#         q5 = ~quantile(value, 5/100, na.rm = TRUE),
#         q25 = ~quantile(value, 25/100, na.rm = TRUE),
#         q50 = ~quantile(value, 50/100, na.rm = TRUE),
#         q75 = ~quantile(value, 75/100, na.rm = TRUE),
#         q95 = ~quantile(value, 95/100, na.rm = TRUE),
#         q97.5 = ~quantile(value, 97.5/100, na.rm = TRUE),
#         q99 = ~quantile(value, 99/100, na.rm = TRUE)
#       ) %>%
#       dplyr::left_join(object$sampleTable %>% dplyr::select_(~type,~id,~node1,~node2,sample = ~value), by=c("id","type","node1","node2")) %>%
#       dplyr::select_(~type, ~id, ~node1, ~node2, ~sample, ~mean, ~var, ~q1, ~q2.5, ~q5, ~q25, ~q50, ~q75, ~q95, ~q97.5, ~q99)
#     

   
    if (object$type == "jackknife"){

      N <- object$sampleSize
      tab <- object$bootTable %>% 
        dplyr::filter_(~type %in% statistics) %>%
        dplyr::left_join(object$sampleTable %>% dplyr::select_(~type,~id,~node1,~node2,sample = ~value), by=c("id","type","node1","node2")) %>%
        dplyr::mutate_(PS = ~N*sample - (N-1)*value) %>%
        dplyr::group_by_(~type, ~node1, ~node2, ~id) %>%
        dplyr::summarize_(
          mean = ~mean(value),
          sample = ~mean(PS),
          var = ~(1/(N-1)) * sum((PS - value)^2),
          CIlower = ~sample - 2 * sqrt(var/N),
          CIupper = ~sample + 2 * sqrt(var/N)
        )%>%
        dplyr::select_(~type, ~id, ~node1, ~node2, ~sample, ~mean, ~CIlower, ~CIupper)
      
    } else {
      tab <- object$bootTable %>% 
        dplyr::filter_(~type %in% statistics) %>%
        dplyr::group_by_(~type, ~node1, ~node2, ~id) %>%
        dplyr::summarize_(
          mean = ~mean(value,na.rm=TRUE),
          var = ~var(value,na.rm=TRUE),
          sd = ~sd(value,na.rm=TRUE),
          prop0 = ~mean(value == 0),
          # q1 = ~quantile(value,1/100, na.rm = TRUE),
          q2.5 = ~quantile(value_min, 2.5/100, na.rm = TRUE),
          #                 q5 = ~quantile(value, 5/100, na.rm = TRUE),
          #                 q25 = ~quantile(value, 25/100, na.rm = TRUE),
          #                 q50 = ~quantile(value, 50/100, na.rm = TRUE),
          #                 q75 = ~quantile(value, 75/100, na.rm = TRUE),
          #                 q95 = ~quantile(value, 95/100, na.rm = TRUE),
          q97.5 = ~quantile(value_max, 97.5/100, na.rm = TRUE)
          # q99 = ~quantile(value, 99/100, na.rm = TRUE)
        ) %>%
        dplyr::left_join(object$sampleTable %>% dplyr::select_(~type,~id,~node1,~node2,sample = ~value), by=c("id","type","node1","node2"))   %>%
        dplyr::mutate_(CIlower = ~sample-2*sd, CIupper = ~sample + 2*sd) %>%
        dplyr::select_(~type, ~id, ~node1, ~node2, ~sample, ~mean, ~sd, ~CIlower, ~CIupper,
                       ~q2.5, ~q97.5)
    }
    
    
    
  } else {

    # Nodewise
    tab <- object$bootTable %>% 
      dplyr::filter_(~type %in% statistics) %>% 
      dplyr::left_join(object$sampleTable %>% dplyr::select_(~type,~id,~node1,~node2,sample = ~value), by=c("id","type","node1","node2"))
    
    if (perNode){
      tab <- tab %>% group_by_(~id, ~type, ~nNode, ~nPerson)  %>%
        dplyr::summarize_(
          mean = ~mean(value,na.rm=TRUE),
          var = ~var(value,na.rm=TRUE),
          sd = ~sd(value,na.rm=TRUE),
          q1 = ~quantile(value,1/100, na.rm = TRUE),
          q2.5 = ~quantile(value, 2.5/100, na.rm = TRUE),
          q5 = ~quantile(value, 5/100, na.rm = TRUE),
          q25 = ~quantile(value, 25/100, na.rm = TRUE),
          q50 = ~quantile(value, 50/100, na.rm = TRUE),
          q75 = ~quantile(value, 75/100, na.rm = TRUE),
          q95 = ~quantile(value, 95/100, na.rm = TRUE),
          q97.5 = ~quantile(value, 97.5/100, na.rm = TRUE),
          q99 = ~quantile(value, 99/100, na.rm = TRUE)
        ) %>% mutate_(CIlower = ~mean - 2*sd, CIupper = ~mean + 2*sd) %>% arrange_(~nNode,~nPerson)
      
    } else {
      tab <- tab %>% group_by_(~name, ~type, ~nNode, ~nPerson)  %>%
        summarize_(cor = ~ suppressWarnings(cor(value,sample, use = "pairwise.complete.obs"))) %>%
        dplyr::group_by_(~nNode, ~nPerson, ~type) %>%
        dplyr::summarize_(
          mean = ~mean(cor,na.rm=TRUE),
          var = ~var(cor,na.rm=TRUE),
          sd = ~sd(cor,na.rm=TRUE),
          q1 = ~quantile(cor,1/100, na.rm = TRUE),
          q2.5 = ~quantile(cor, 2.5/100, na.rm = TRUE),
          q5 = ~quantile(cor, 5/100, na.rm = TRUE),
          q25 = ~quantile(cor, 25/100, na.rm = TRUE),
          q50 = ~quantile(cor, 50/100, na.rm = TRUE),
          q75 = ~quantile(cor, 75/100, na.rm = TRUE),
          q95 = ~quantile(cor, 95/100, na.rm = TRUE),
          q97.5 = ~quantile(cor, 97.5/100, na.rm = TRUE),
          q99 = ~quantile(cor, 99/100, na.rm = TRUE)
        ) %>% arrange_(~nNode, ~nPerson)
      
    }
  }
  
  
  return(tab)
}


