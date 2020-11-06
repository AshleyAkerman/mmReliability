
#' Title: Extracting reliability Estimates from Raw Data
#'
#' This function will take a value of the number of variable required, and a long and wide format of the data frame
#' Requires the reliability_function from the Hopkins spreadsheet and also requires the psych package for ICC
#'
#' @param n a numeric values indicating the number of variables in the dataframe
#' @param long a dataframe in long-form
#' @param wide a dataframe in wide-form
#' @return list of calculated reliability statistics (hopkins) and icc (psych) and raw data for effect size calculations
#' @keywords reliability, intra-class correlation coefficient, coefficient of variation
#' @export
#'
reliability_extraction_function <- function(n, long, wide){

  # There are n variables, and a subject column, so need to generate the numbers for the variables
  # It will start at 2 for that reason, and range to the length + 1.
  variable_placement <- c(2:(n + 1))

  # List of dataframes to store outputs.
  # This will end up being length equal to the number of variables
  list_of_dataframes <- list()

  # List of dataframes to store outputs for effect sizes
  # This will end up being length equal to the number of variables
  effect_size_list <- list()

  # list for ICCs from the psych package
  # This will end up being length equal to the number of variables
  newICC_list <- list()

  # For each variable number...
  for (variable in 1: length(variable_placement)) {

    # If the length of the variables is only 1, then do not need to subset the dataframe
    if(n == 1) {

      df <- wide

    }

    # But if it is larger than 1 variable, need to adjust the dataframe
    if(n > 1) {

      # Adjusting data for independent analysis of time points
      # The variable + 1 is because the Subject column needs to be accounted for
      # The second part of the cbind takes the columns of the respective data frame, starting from the column with the first variable, and then each column corresponding to that variable
      df <- cbind(wide$Subject, wide[, c(seq(variable + 1, (length(long) - (ifelse(variable + 1 == max(variable_placement), 1, 2))) * 3,
                                             by = length(long) - 2))])

    }

    # Add just the data to a df for ICC calculation
    icc_data <- df[,-1]
    # Save the output
    icc_output <- psych::ICC(log(icc_data))
    # Extract just the ICC[3, 1] form.
    newICC_list[[variable]] <- icc_output$results[3,]

    # Manipulate the df from above using
    df <- df_manip(df)
    # Re-order so that it matches the requirement for the reliability function
    df <- df[, col_order]

    # Saving the raw data for the effect size calculations later
    effect_size_list[[variable]] <- df
    # Add to the list of dataframes
    list_of_dataframes[[variable]] <- reliability_function(df)

  }

  large_list <- append(list_of_dataframes, newICC_list)
  large_list <- append(large_list, effect_size_list)

  return(large_list)

}
