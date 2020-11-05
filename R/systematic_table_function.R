
#' Title: extract systematic differences for reliability data
#'
#' @param data_list list of data tables
#' @param NHS string input taking form of "Y" or "N"
#' @param LHS string input taking form of "Y" or "N"
#' @param MHS string input taking form of "Y" or "N"
#' @param HHS string input taking form of "Y" or "N"
#' @param levels numeric value indicating number of heating levels (2: NHS and HHS; 4: NHS, LHS, MHS, and HHS)
#' @return list of figures (1st minor index) followed by sample size requirements (2nd and 3rd index)
#' @keywords reliability, intra-class correlation coefficient, coefficient of variation
#' @export
#'

systematic_table_function <- function(data_list, NHS, LHS, MHS, HHS, levels) {

  # New variable for use later and output
  p_value_list <- c()

  # For every table in the list of tables....
  for (table in 1: length(data_list)) {

    # Generate a new table that is just the raw data from each table
    new_table <- data_list[[table]][1:4]

    # Simple linear model on the data for systematic differences
    new_table$Subject <- factor(new_table$Subject)  # change subject to factor
    new_table <- tidyr::gather(new_table, Trial, measurement, Trial1:Trial3, factor_key = TRUE)  # change from wide to long format

    levels(new_table$Trial) <- c(1, 2, 3)   # Rename and change Trial to a factor
    new_table$Trial <- factor(new_table$Trial)

    # Anova results from lme
    lme_results <- stats::lm(log(measurement) ~ Trial, data = stats::na.omit(new_table))
    p_value <- stats::anova(lme_results)$'Pr(>F)'[1] # Extract just the p-value

    # Add the result to the new list
    p_value_list <- rbind(p_value_list, p_value)


  }

  # Remove row values
  row.names(p_value_list) <- NULL

  # Each variable (n=4) will be repeated 4 times because there are 4 levels for each  (unless specified)
  Variable_name <- rep(seq(1, length(variable_names)), each = levels)

  # Need to determine the names based on what is provided in the function
  ## There are very few options for the manuscripts, so they are provided below

  # If there are all 4 levels - sudomotor and cardiovascular
  if(NHS == "Y" & LHS == "Y" & MHS == "Y" & HHS == "Y") {

    Comparison <- rep(c("Baseline", "Low Heat Strain", "Moderate Heat Strain", "High Heat Strain"), times = length(variable_names))

  }

  # If there are just 3 levels (i.e., no baseline) - sudomotor
  if(NHS == "N" & LHS == "Y" & MHS == "Y" & HHS == "Y") {

    Comparison <- rep(c("Low Heat Strain", "Moderate Heat Strain", "High Heat Strain"), times = length(variable_names))

  }

  # If there are just 3 levels (i.e., no high heat strain) - microvascular
  if(NHS == "Y" & LHS == "Y" & MHS == "Y" & HHS == "N") {

    Comparison <- rep(c("Baseline", "Low Heat Strain", "Moderate Heat Strain"), times = length(variable_names))

  }

  # If there are just 2 levels (i.e., no low or moderate) - cardiovascular squats
  if(NHS == "Y" & LHS == "N" & MHS == "N" & HHS == "Y") {

    Comparison <- rep(c("Baseline", "High Heat Strain"), times = length(variable_names))

  }

  # If there is just 1 level - sudomotor
  if(NHS == "Y" & LHS == "N" & MHS == "N" & HHS == "N") {

    Comparison <- rep(c("Baseline"), times = length(variable_names))

  }


  # Combine the data
  PvalueTable <- data.frame(Variable_name = Variable_name,
                            Comparison = Comparison,
                            P_value = round(p_value_list, digits = 3))

  # Change the values of each to the proper names
  PvalueTable$Variable_name <- rep(variable_names, each = levels)

  return(PvalueTable)

}



