
#' Title: sample size calculation and power curve plot for Reliability Manuscripts
#'
#' @param input A vector of effect sizes is passed to the function (typically increments of 0.5)
#' @return list of required sample sizes and figures
#' @keywords reliability, intra-class correlation coefficient, coefficient of variation
#' @export
#'

plot_and_sizes_function <- function(input) {

  # Check it is a vector
  input <- as.vector(input)

  # Generate new empty list
  output_list <- c()

  # Most studies would never have more than 50 participants, so stop there but generate values from 6 to 50, in increments of 1
  # This will be consistent for each of the iterations below
  n_samples <- seq(5, 100, 1)

  # Convert to long format whereby each iteration of n samples with have each of the effect sizes next to it
  conditions <- expand.grid(n = n_samples, effect_sizes = input)

  # Then run a power calculation using a t-test on each of those rows
  # i.e., taking the sample size, and the effect size and then extracting the provided power
  power_curve_within <- sapply(seq_len(nrow(conditions)), function(i)
    pwr::pwr.t.test(n = conditions[[i , 1]],
               d = conditions[[i, 2]],
               sig.level = 0.05,
               type = "paired")$power)

  power_curve_between <- sapply(seq_len(nrow(conditions)), function(i)
    pwr::pwr.t.test(n = conditions[[i , 1]],
               d = conditions[[i, 2]],
               sig.level = 0.05,
               type = "two.sample")$power)

  # Bind tthe columns together so that the data frame above (i.e., 'conditions' has the achieved power next to it))
  power_curve_df <- dplyr::bind_cols(
    conditions,
    power_within = power_curve_within,
    power_between = power_curve_between)

  # Plot the data with a vertical line at 0.8 power
  plot_within <- ggplot2::ggplot(power_curve_df, ggplot2::aes(x = power_within, y = n)) +
    ggplot2::geom_line(ggplot2::aes(color = factor(effect_sizes)), size = 0.5) +
    ggplot2::geom_vline(xintercept = 0.8, linetype = 2, colour = 'gray30') +
    ggplot2::scale_x_continuous("Power", breaks = seq(0, 1, .1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous("Sample Size", breaks = seq(5, 100, 5), expand = c(0, 0)) +
    ggplot2::scale_color_viridis_d("Effect Size") +
    ggplot2::ggtitle("A. Power curve for paired samples") +
    #ggplot2::theme(plot.title = element_text(size = 12, face = "bold")) +
    ggplot2::theme_bw(base_size = 12)

  # Plot the data with a vertical line at 0.8 power
  plot_between <- ggplot2::ggplot(power_curve_df, ggplot2::aes(x = power_between, y = n)) +
    ggplot2::geom_line(ggplot2::aes(color = factor(effect_sizes)), size = 0.5) +
    ggplot2::geom_vline(xintercept = 0.8, linetype = 2, colour = 'gray30') +
    ggplot2::scale_x_continuous("Power", breaks = seq(0, 1, .1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous("Sample Size", breaks = seq(5, 100, 5), expand = c(0, 0)) +
    ggplot2::scale_color_viridis_d("Effect Size") +
    ggplot2::ggtitle("B. Power curve for independent samples") +
    #ggplot2::theme(plot.title = element_text(size = 12, face = "bold")) +
    ggplot2::theme_bw(base_size = 12)

  # Need to use Grob functions because this will mean that excess processing power is not taken up plotting each jointplot
  # If not fussed, then use grid.arrange(plot_within, plot_between, ncol = 2)
  # To save, it would then be ggsave(..)
  # To plot the grob, use grid::grid.draw(joint_plot)
  joint_plot <- gridExtra::arrangeGrob(plot_within, plot_between, nrow = 2) # changed from ncol = 2

  #grid::grid.draw(joint_plot)

  # Make a new vector to hold the required sample sizes
  sample_size_req_within <- c()

  # Iterate through the number of effect sizes desired
  # For each one filter just for that effect size and the smallest sample required for 80% power
  for (i in 1: length(input)) {
    all_powers <- power_curve_df %>% dplyr::filter(effect_sizes == input[i] & power_within >= 0.80)
    smallest_power <- all_powers[1,]
    sample_size_req_within <- c(sample_size_req_within, smallest_power[[1]])
  }

  # Bind it together so you now have the smallest sample size required for a given effect size, to achieve 80% power
  sample_size_req_within <- cbind(sample_size_req_within, input)

  # Make a new vector to hold the required sample sizes
  sample_size_req_between <- c()

  # Iterate through the number of effect sizes desired
  for (i in 1: length(input)) {
    all_powers <- power_curve_df %>% dplyr::filter(effect_sizes == input[i] & power_between >= 0.80)
    smallest_power <- all_powers[1,]
    sample_size_req_between <- c(sample_size_req_between, smallest_power[[1]])
  }

  # Bind it together so you now have the smallest sample size required for a given effect size, to achieve 80% power
  sample_size_req_between <- cbind(sample_size_req_between, input)


  output_list <- list(joint_plot, sample_size_req_within, sample_size_req_between)

  return(output_list)

}


