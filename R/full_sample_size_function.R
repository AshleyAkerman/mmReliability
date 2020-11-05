
#' Title: Convert MDC to standardised effects size and extract sample size data and power plot
#'
#' NOTES: In the case of sudomotor, the one level will be designated via NHS
#' Use the plot_and_sizes function to generate the sample sizes required and the plots for those data (i.e., power curves)
#'
#' @param df dataframe
#' @param NHS string input taking form of "Y" or "N"
#' @param LHS string input taking form of "Y" or "N"
#' @param MHS string input taking form of "Y" or "N"
#' @param HHS string input taking form of "Y" or "N"
#' @param SESOI vector of numeric values defining the smallest effect size of interest (as %). Default is: 1, 3, 5, 10, 20, 25, 30%
#' @return list of figures (1st minor index) followed by sample size requirements (2nd and 3rd index)
#' @keywords reliability, intra-class correlation coefficient, coefficient of variation
#' @export
#'

full_sample_size_function <- function(df,
                                      NHS,
                                      LHS,
                                      MHS,
                                      HHS,
                                      SESOI = c(1, 3, 5, 10, 20, 25, 30)){

  # Need a blank list to start.

  # At the end there will be lists within lists.
  # The major index is the outcome measures
  # The minor indices is length equal to 3 * levels of heat strain
  # 1: Grob of within plot and between plot at a given level of heat strain
  # 2: sample sizes required for within subject design at a given level of heat strain
  # 3: sample sizes required for within subject design at a given level of heat strain

  sample_sizes <- list()

  for (name in 1: length(variable_names)){

    ## If NHS is a yes
    if(NHS == "Y") {

      # Extract the MDC for BL from each variable
      effects_of_MDC_BL <- df %>% dplyr::filter(Variable_name == variable_names[name] & Levels == "BL") %>%

        dplyr::summarise(eff1 = (MDC70_comb_gp / sdDiff),
                         eff2 = (MDC80_comb_gp / sdDiff),
                         eff3 = (MDC90_comb_gp / sdDiff),
                         eff4 = (MDC95_comb_gp / sdDiff),
                         eff5 = (MDC99_comb_gp / sdDiff),
                         eff6 = ((grandmean_sampleSize * (SESOI / 100)) / sdDiff))


      # Generate the effect sizes from values extracted above
      effect_sizes_BL <- round(sort(unique(unlist(effects_of_MDC_BL))), digits = 2)

      ## The plot_and_sizes_function:
      # Most studies would never have more than 50 participants, so stop there but generate values from 6 to 50, in increments of 1
      # Combine the sample sizes and the effect sizes, then run a power calculation using a t-test on each of those rows of the df
      # i.e., taking the sample size, and the effect size and then extracting the provided power
      # Plots the data with a vertical line at 0.8 power
      # Then iterates through the number of effect sizes desired
      # For each one filter just for that effect size and the smallest sample required for 80% power
      BL_results <- plot_and_sizes_function(effect_sizes_BL)

    }


    ## If LHS is a yes
    if(LHS == "Y") {

      # Extract the MDC for LH from each variable
      effects_of_MDC_LH <- df %>% dplyr::filter(Variable_name == variable_names[name] & Levels == "LH") %>%

        dplyr::summarise(eff1 = (MDC70_comb_gp / sdDiff),
                         eff2 = (MDC80_comb_gp / sdDiff),
                         eff3 = (MDC90_comb_gp / sdDiff),
                         eff4 = (MDC95_comb_gp / sdDiff),
                         eff5 = (MDC99_comb_gp / sdDiff),
                         eff6 = ((grandmean_sampleSize * (SESOI / 100)) / sdDiff))


      # Generate the effect sizes from values extracted above
      effect_sizes_LH <- round(sort(unique(unlist(effects_of_MDC_LH))), digits = 2)

      ## The plot_and_sizes_function:
      # Most studies would never have more than 50 participants, so stop there but generate values from 6 to 50, in increments of 1
      # Combine the sample sizes and the effect sizes, then run a power calculation using a t-test on each of those rows of the df
      # i.e., taking the sample size, and the effect size and then extracting the provided power
      # Plots the data with a vertical line at 0.8 power
      # Then iterates through the number of effect sizes desired
      # For each one filter just for that effect size and the smallest sample required for 80% power
      LH_results <- plot_and_sizes_function(effect_sizes_LH)

    }

    ## If MHS is a yes
    if(MHS == "Y") {

      # Extract the MDC for MH from each variable
      effects_of_MDC_MH <- df %>% dplyr::filter(Variable_name == variable_names[name] & Levels == "MH") %>%

        dplyr::summarise(eff1 = (MDC70_comb_gp / sdDiff),
                         eff2 = (MDC80_comb_gp / sdDiff),
                         eff3 = (MDC90_comb_gp / sdDiff),
                         eff4 = (MDC95_comb_gp / sdDiff),
                         eff5 = (MDC99_comb_gp / sdDiff),
                         eff6 = ((grandmean_sampleSize * (SESOI / 100)) / sdDiff))


      # Generate the effect sizes from values extracted above
      effect_sizes_MH <- round(sort(unique(unlist(effects_of_MDC_MH))), digits = 2)

      ## The plot_and_sizes_function:
      # Most studies would never have more than 50 participants, so stop there but generate values from 6 to 50, in increments of 1
      # Combine the sample sizes and the effect sizes, then run a power calculation using a t-test on each of those rows of the df
      # i.e., taking the sample size, and the effect size and then extracting the provided power
      # Plots the data with a vertical line at 0.8 power
      # Then iterates through the number of effect sizes desired
      # For each one filter just for that effect size and the smallest sample required for 80% power
      MH_results <- plot_and_sizes_function(effect_sizes_MH)

    }


    ## If HHS is a yes
    if(HHS == "Y") {

      # Extract the MDC for MH from each variable
      effects_of_MDC_HH <- df %>% dplyr::filter(Variable_name == variable_names[name] & Levels == "HH") %>%

        dplyr::summarise(eff1 = (MDC70_comb_gp / sdDiff),
                         eff2 = (MDC80_comb_gp / sdDiff),
                         eff3 = (MDC90_comb_gp / sdDiff),
                         eff4 = (MDC95_comb_gp / sdDiff),
                         eff5 = (MDC99_comb_gp / sdDiff),
                         eff6 = ((grandmean_sampleSize * (SESOI / 100)) / sdDiff))



      # Generate the effect sizes from values extracted above
      effect_sizes_HH <- round(sort(unique(unlist(effects_of_MDC_HH))), digits = 2)

      ## The plot_and_sizes_function:
      # Most studies would never have more than 50 participants, so stop there but generate values from 6 to 50, in increments of 1
      # Combine the sample sizes and the effect sizes, then run a power calculation using a t-test on each of those rows of the df
      # i.e., taking the sample size, and the effect size and then extracting the provided power
      # Plots the data with a vertical line at 0.8 power
      # Then iterates through the number of effect sizes desired
      # For each one filter just for that effect size and the smallest sample required for 80% power
      HH_results <- plot_and_sizes_function(effect_sizes_HH)

    }


    ## Final sample size df dependent on what was given in the function
    ## There are very few options for the manuscripts, so they are provided below

    # If there are all 4 levels - sudomotor and cardiovascular
    if(NHS == "Y" & LHS == "Y" & MHS == "Y" & HHS == "Y") {

      sample_sizes[[name]] <- c(BL_results, LH_results, MH_results, HH_results)

    }

    # If there are just 3 levels (i.e., no baseline) - sudomotor
    if(NHS == "N" & LHS == "Y" & MHS == "Y" & HHS == "Y") {

      sample_sizes[[name]] <- c(LH_results, MH_results, HH_results)

    }

    # If there are just 3 levels (i.e., no high heat strain) - microvascular
    if(NHS == "Y" & LHS == "Y" & MHS == "Y" & HHS == "N") {

      sample_sizes[[name]] <- c(BL_results, LH_results, MH_results)

    }

    # If there are just 2 levels (i.e., no low or moderate) - cardiovascular squats
    if(NHS == "Y" & LHS == "N" & MHS == "N" & HHS == "Y") {

      sample_sizes[[name]] <- c(BL_results, HH_results)

    }

    # If there is just 1 level - sudomotor
    if(NHS == "Y" & LHS == "N" & MHS == "N" & HHS == "N") {

      sample_sizes[[name]] <- c(BL_results)

    }

  }

  # Return the appropriate df
  return(sample_sizes)

}


