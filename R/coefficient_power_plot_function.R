
#' Title: Coefficient of variation and Power plots for different tools and effect sizes
#'
#' NOTES: Function to generate plots for sample size estimations based on the coefficient of variation of two tools given a similar mean estimate
#' NOTES: Or alternatively a power plot for sample size estimations based on effect size and power.
#'
#' @param MDC_input vector of values referring to the minimum detectable change estimates. Only required for power curve
#' @param grand_mean numeric value for the grand mean. If coefficient figure then this is the value consistent for both measurement tools. If a power curve, it is the grand mean of the one tool.
#' @param fig_type string input taking form of "coefficient" or "power". Coefficient is to compare sample sizes across CV values. Power is to compare sample sizes across power for given effect sizes.
#' @param SESOI vector of values for the smallest effect size of interest (default: 1, 3, 5, 10, 20, 25, 30%)
#' @param CoV vector of values for the coefficient of variation range important for the "coefficient" plot. Default is a sequence from 4 to 50 with increments of 2.
#' @param CV1 numeric value for the coefficient of variation. If coefficient figure then the 1st CV would be the poorer tool, if power figure then CV for the only tool (and CV2 == NA)
#' @param CV2 numeric value for the Coefficient of variation. If coefficient figure this is required as the better of the two tools.
#' @param all_percent boolean input. Default is FALSE. TRUE when the smallest effect size of interest and the MDC are the same unit (i.e., %), so don't need to transfer to absolute values.
#' @return figure either for coefficient of variation or power curve and estimating sample size.
#' @keywords reliability, intra-class correlation coefficient, coefficient of variation
#' @export
#'

coeff_power_plot_fun <- function(MDC_input = NA,
                                 grand_mean,
                                 fig_type,
                                 SESOI = c(1, 3, 5, 10, 20, 25, 30),
                                 CoV = seq(4, 50, 2),
                                 CV1 = NA,
                                 CV2 = NA,
                                 all_percent = FALSE) {

  if(fig_type == "coefficient"){

    if(all_percent == FALSE){

      # What are the SESOI as raw values
      effects_1 <- grand_mean * (SESOI / 100)

      # Convert to long format whereby each iteration of n samples with have each of the effect sizes next to it
      conditions_1 <- expand.grid(cov = CoV, effect_sizes = effects_1)

      conditions_1$sd_diff <- sqrt(2) * (((conditions_1$cov / 100)) * grand_mean)

      # First measurement tool
      output_1 <- c()

      # Loop over the rows and run a power calculation
      for (i in 1: nrow(conditions_1)){
        # Store the results of the sample size in the vector saved above
        output_1 <- c(output_1,

                      pwr::pwr.t.test(d = conditions_1[[i, 2]] / conditions_1[[i, 3]],
                                      sig.level = 0.05,
                                      power = 0.80,
                                      type = "paired")$n)

      }


      # Make a new df with the results and the corresponding conditions of the calculations
      power_curve_df_1 <- dplyr::bind_cols(
        conditions_1,
        samples_within = output_1)

      # Filter to samples < 200 people for ease of display
      power_curve_df_1 <- power_curve_df_1 %>% filter(samples_within < 200)

      # Change the values of the 'change' column
      # For every row of the dataframe...
      for (i in 1: nrow(power_curve_df_1)){
        # And every value in the Change scores (grand_mean * (SESOI / 100))
        for (j in 1: length(effects_1)){
          # Change the value in absolute units back to what it was in SESOI units
          if(power_curve_df_1[i, 2] == effects_1[j]){

            power_curve_df_1[i, 2] <- SESOI[j]

          }

        }

      }

      # Plot for first measurement tool
      plot_1 <- ggplot(power_curve_df_1,
                       aes(x = cov, y = samples_within)) +
        geom_line(aes(color = factor(effect_sizes)), size = 0.5) +
        geom_vline(xintercept = CV1, linetype = 2, colour = 'red') +
        geom_vline(xintercept = CV2, linetype = 2, colour = 'blue') +
        scale_x_continuous("Coefficient of Variation (%)", breaks = CoV, expand = c(0, 0)) +
        scale_y_continuous("Sample Size", breaks = seq(0, 200, 10), expand = c(0, 0)) +
        scale_color_viridis_d("Effect Size", guide = guide_legend(override.aes = list(size = 3,
                                                                                      alpha = 1) ) ) +
        theme_bw(base_size = 15)

      return(plot_1)

    }

    if(all_percent == TRUE) {

      # What are the SESOI as raw values
      effects_1 <- SESOI

      # Convert to long format whereby each iteration of n samples with have each of the effect sizes next to it
      conditions_1 <- expand.grid(cov = CoV, effect_sizes = effects_1)

      conditions_1$sd_diff <- sqrt(2) * (((conditions_1$cov / 100)) * grand_mean)

      # First measurement tool
      output_1 <- c()

      # Loop over the rows and run a power calculation
      for (i in 1: nrow(conditions_1)){
        # Store the results of the sample size in the vector saved above
        output_1 <- c(output_1,

                      pwr::pwr.t.test(d = conditions_1[[i, 2]] / conditions_1[[i, 3]],
                                      sig.level = 0.05,
                                      power = 0.80,
                                      type = "paired")$n)

      }


      # Make a new df with the results and the corresponding conditions of the calculations
      power_curve_df_1 <- dplyr::bind_cols(
        conditions_1,
        samples_within = output_1)

      # Filter to samples < 200 people for ease of display
      power_curve_df_1 <- power_curve_df_1 %>% filter(samples_within < 200)

      # Change the values of the 'change' column
      # For every row of the dataframe...
      for (i in 1: nrow(power_curve_df_1)){
        # And every value in the Change scores (grand_mean * (SESOI / 100))
        for (j in 1: length(effects_1)){
          # Change the value in absolute units back to what it was in SESOI units
          if(power_curve_df_1[i, 2] == effects_1[j]){

            power_curve_df_1[i, 2] <- SESOI[j]

          }

        }

      }

      # Plot for first measurement tool
      plot_2 <- ggplot(power_curve_df_1,
                       aes(x = cov, y = samples_within)) +
        geom_line(aes(color = factor(effect_sizes)), size = 0.5) +
        geom_vline(xintercept = CV1, linetype = 2, colour = 'red') +
        geom_vline(xintercept = CV2, linetype = 2, colour = 'blue') +
        scale_x_continuous("Coefficient of Variation (%)", breaks = CoV, expand = c(0, 0)) +
        scale_y_continuous("Sample Size", breaks = seq(0, 200, 10), expand = c(0, 0)) +
        scale_color_viridis_d("Effect Size", guide = guide_legend(override.aes = list(size = 3,
                                                                                      alpha = 1) ) ) +
        theme_bw(base_size = 15)

      return(plot_2)

    }



  }


  if(fig_type == "power"){

    if(all_percent == FALSE) {

      sd_diff <- sqrt(2) * ((CV1 / 100) * grand_mean)
      MDC_as_Change <- round(((MDC_input/grand_mean) * 100), 0)
      SESOI <- sort(c(SESOI, MDC_as_Change))
      effects_1 <- grand_mean * (SESOI / 100)
      effects_2 <- round(effects_1 / sd_diff, 2)


      # Generate new empty list
      output_list <- c()

      # Most studies would never have more than 50 participants, so stop there but generate values from 6 to 50, in increments of 1
      # This will be consistent for each of the iterations below
      n_samples <- seq(6, 100, 1)

      # Convert to long format whereby each iteration of n samples with have each of the effect sizes next to it
      conditions <- expand.grid(n = n_samples, effect_sizes = effects_2)

      # Then run a power calculation using a t-test on each of those rows
      # i.e., taking the sample size, and the effect size and then extracting the provided power
      power_curve_within <- sapply(seq_len(nrow(conditions)), function(i)
        pwr::pwr.t.test(n = conditions[[i , 1]],
                        d = conditions[[i, 2]],
                        sig.level = 0.05,
                        type = "paired")$power)

      # Bind the columns together so that the data frame above (i.e., 'conditions' has the achieved power next to it))
      power_curve_df <- dplyr::bind_cols(
        conditions,
        power_within = power_curve_within)

      # Plot the data with a vertical line at 0.8 power
      plot_3 <- ggplot2::ggplot(power_curve_df, ggplot2::aes(x = power_within, y = n)) +
        ggplot2::geom_line(ggplot2::aes(color = factor(effect_sizes)), size = 0.5) +
        ggplot2::geom_vline(xintercept = 0.8, linetype = 2, colour = 'gray30') +
        ggplot2::scale_x_continuous("Power", breaks = seq(0, 1, .1), expand = c(0, 0)) +
        ggplot2::scale_y_continuous("Sample Size", breaks = seq(6, 206, 10), expand = c(0, 0)) +
        scale_color_viridis_d("Effect Size", guide = guide_legend(override.aes = list(size = 3,
                                                                                      alpha = 1) ) ) +
        ggplot2::theme_bw(base_size = 12)

      return(plot_3)

    }

    if(all_percent == TRUE) {

      sd_diff <- sqrt(2) * ((CV1 / 100) * grand_mean)
      MDC_as_Change <- round(MDC_input, 0)
      SESOI <- sort(c(SESOI, MDC_as_Change))
      effects_1 <- round(SESOI / sd_diff, 2)


      # Generate new empty list
      output_list <- c()

      # Most studies would never have more than 50 participants, so stop there but generate values from 6 to 50, in increments of 1
      # This will be consistent for each of the iterations below
      n_samples <- seq(6, 100, 1)

      # Convert to long format whereby each iteration of n samples with have each of the effect sizes next to it
      conditions <- expand.grid(n = n_samples, effect_sizes = effects_1)

      # Then run a power calculation using a t-test on each of those rows
      # i.e., taking the sample size, and the effect size and then extracting the provided power
      power_curve_within <- sapply(seq_len(nrow(conditions)), function(i)
        pwr::pwr.t.test(n = conditions[[i , 1]],
                        d = conditions[[i, 2]],
                        sig.level = 0.05,
                        type = "paired")$power)

      # Bind the columns together so that the data frame above (i.e., 'conditions' has the achieved power next to it))
      power_curve_df <- dplyr::bind_cols(
        conditions,
        power_within = power_curve_within)

      # Plot the data with a vertical line at 0.8 power
      plot_4 <- ggplot2::ggplot(power_curve_df, ggplot2::aes(x = power_within, y = n)) +
        ggplot2::geom_line(ggplot2::aes(color = factor(effect_sizes)), size = 0.5) +
        ggplot2::geom_vline(xintercept = 0.8, linetype = 2, colour = 'gray30') +
        ggplot2::scale_x_continuous("Power", breaks = seq(0, 1, .1), expand = c(0, 0)) +
        ggplot2::scale_y_continuous("Sample Size", breaks = seq(6, 206, 10), expand = c(0, 0)) +
        scale_color_viridis_d("Effect Size", guide = guide_legend(override.aes = list(size = 3,
                                                                                      alpha = 1) ) ) +
        ggplot2::theme_bw(base_size = 12)

      return(plot_4)

    }

  }


}
