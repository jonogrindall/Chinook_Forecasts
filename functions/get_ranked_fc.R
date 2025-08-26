get_ranked_fc <- function(fc) {
  age_cols <- c("Age 2", "Age 3", "Age 4", "Age 5")
  
  age_rank_lists <- list()
  summary_lists <- list()
  
  for (stock in names(fc)) {
    stock_fc <- fc[[stock]]
    
    # Get rank matrix from retro
    ranks_bal <- stock_fc$retro.pm$retro.pm.all.varyrs
    ranked_bal <- rankModels(ranks_bal)
    
    # Build cumulative rank summary
    cum_ranked_bal <- ranked_bal$cumulativerank %>%
      tibble::rownames_to_column(var = "Model")
    
    # Build table with forecasts and cumulative rank sum
    forecast_sorted <- round(stock_fc$table.ptf) %>%
      tibble::rownames_to_column(var = "Model") %>%
      dplyr::left_join(cum_ranked_bal, by = "Model") %>%
      dplyr::arrange(rank.sum)
    
    existing_age_cols <- intersect(age_cols, names(forecast_sorted))
    forecast_ranksum <- forecast_sorted %>%
      dplyr::select(Model, all_of(existing_age_cols), rank.sum)
    
    # Extract rank.avg per age class
    ranked_ages <- ranked_bal[setdiff(names(ranked_bal),
                                      c("cumulativerank", "bestmodel", "cumulativerankSorted"))]
    
    age_ranked_bal <- purrr::map(
      ranked_ages,
      ~ .x %>%
        tibble::rownames_to_column(var = "Model") %>%
        dplyr::select(Model, rank.avg) %>%
        dplyr::arrange(rank.avg)
    )
    
    # Store all tied best models per age as collapsed strings
    best_models_table <- tibble(
      AgeClass = names(age_ranked_bal),
      Models = age_ranked_bal
    ) %>%
      rowwise() %>%
      mutate(
        Models = {
          top_models <- Models$Model[Models$rank.avg == min(Models$rank.avg)]
          
          # Determine youngest age dynamically
          age_classes <- names(age_ranked_bal)
          youngest_age <- age_classes[grepl("^Age", age_classes)] %>% sort() %>% first()
          
          # Apply Naive5 vs Sib logic only for youngest age
          if (AgeClass == youngest_age && "Naive5" %in% top_models) {
            top_models <- top_models[!grepl("^Sib", top_models)]
          }
          
          paste(top_models, collapse = ", ")
        }
      ) %>%
      ungroup() %>%
      as.data.frame()
    
    # Store separately
    age_rank_lists[[stock]] <- age_ranked_bal
    summary_lists[[stock]] <- list(
      ranksum = forecast_ranksum,
      bestmodel = best_models_table
    )
  }
  
  # Return as a list of two elements; summary last for convenience
  return(list(
    age_rank_avg = age_rank_lists,
    ranked_forecasts = summary_lists
  ))
}
