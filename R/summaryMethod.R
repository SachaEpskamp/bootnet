# Creates the summary table
summary.bootnet <- function(
  object, # bootnet object
  types = c("edge", "intercept", "strength", "closeness", "betweenness","distance"), # stats to include in the table
  ...
  ){

  tab <- object$bootTable %>% 
    dplyr::filter_(~type %in% types) %>%
    dplyr::group_by_(~type, ~node1, ~node2, ~id) %>%
    dplyr::summarize_(
        mean = ~mean(value),
        var = ~var(value),
        sd = ~sd(value),
        prop0 = ~mean(value == 0),
        q1 = ~quantile(value,1/100),
        q2.5 = ~quantile(value, 2.5/100),
        q5 = ~quantile(value, 5/100),
        q25 = ~quantile(value, 25/100),
        q50 = ~quantile(value, 50/100),
        q75 = ~quantile(value, 75/100),
        q95 = ~quantile(value, 95/100),
        q97.5 = ~quantile(value, 97.5/100),
        q99 = ~quantile(value, 99/100)
      ) %>%
    dplyr::left_join(object$sampleTable %>% dplyr::select_(~type,~id,~node1,~node2,sample = ~value), by=c("id","type","node1","node2")) %>%
    dplyr::select_(~type, ~id, ~node1, ~node2, ~sample, ~mean, ~var, ~q1, ~q2.5, ~q5, ~q25, ~q50, ~q75, ~q95, ~q97.5, ~q99)
    
  return(tab)
}


