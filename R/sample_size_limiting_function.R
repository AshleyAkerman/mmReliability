
#' Title: Limit the sample sizes required to only the MDC values
#'
#' NOTES: Used for sample size tab in the supplementary material document
#'
#' @param stages_info vector of strings denoting which levels ("BL", "LH", "MH", and "HH")
#' @param MDC_main_table dataframe of the MDC estimates as standardised effect sizes
#' @param sample_size_table dataframe of the large sample size dataframe
#' @return dataframe of sample sizes required for each MDC estimate
#' @keywords reliability, intra-class correlation coefficient, coefficient of variation
#' @export
#'

sample_size_limiting_function <- function(stages_info,
                                          MDC_main_table,
                                          sample_size_table){

  # Make a dataframe with the MDCs presented as a standardised change score
  MDC_as_standardised <- cbind(MDC_main_table$Variable_name, MDC_main_table$Levels,
                               MDC_main_table %>% summarise(eff1 = (MDC70_comb_gp / sdDiff),
                                                                 eff2 = (MDC80_comb_gp / sdDiff),
                                                                 eff3 = (MDC90_comb_gp / sdDiff),
                                                                 eff4 = (MDC95_comb_gp / sdDiff),
                                                                 eff5 = (MDC99_comb_gp / sdDiff)))


  names(MDC_as_standardised) <- c("outcome", "stages", "MDC70", "MDC80", "MDC90", "MDC95", "MDC99")

  # New empty vector
  req_samples_reduced <- c()

  # For each variable and each stage, select only the values in the main sheet that are equal to the MDC values
  for (name in 1: length(variable_names)){

    for (stage in 1: length(stages_info)){

      # First holder df is just the data from the large sample size table, but subset to run through each outcome measure and level
       limited1 <- sample_size_table %>% dplyr::select(stages, outcome, effect_size1, within_sample_size, between_sample_size) %>%

         dplyr::filter(outcome == variable_names[name] & stages == stages_info[stage])

       # Second holder is the same as above, but from the MDC as standardised changes df
       limited2 <- MDC_as_standardised %>% dplyr::filter(outcome == variable_names[name] & stages == stages_info[stage])

       # Save the row outputs that are equal to: the corresponding row in the main req samples which match the MDC as standardised effect size values
       row_output <- which(limited1$effect_size1 == round(limited2$MDC70, digits = 2) |
                            limited1$effect_size1 == round(limited2$MDC80, digits = 2)   |
                            limited1$effect_size1 == round(limited2$MDC90, digits = 2)  |
                            limited1$effect_size1 == round(limited2$MDC95, digits = 2)  |
                            limited1$effect_size1 == round(limited2$MDC99, digits = 2) )

      req_samples_reduced <- rbind(req_samples_reduced, limited1[row_output, ])

    }

  }


  return(req_samples_reduced)

}



