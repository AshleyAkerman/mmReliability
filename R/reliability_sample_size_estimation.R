
#' Title: Sample Size Estimations For Reliability Studies
#'
#' Function that takes a given width of a confidence interval (w), number of replicates (n), z score for confidence interval required (z), and planned ICC value (p_h)
#' Formulas given by Shoukri et al., 2004 (equation 7)
#' Shoukri MM, Asyali MH & Donner A (2004). Sample size requirements for the design of reliability study: review and new results. Statistical Methods in Medical Research 13, 251-271.
#'
#' For: w = 0.2, n = 3, z = 1.96, p_h = 0.8, type = "single, should equal k = 35, k1 = 36 to fit with Table 3 in Shoukri et al., 2004
#'
#' @param n numeric value for number of replicates
#' @param w numeric value or vector of values for width of the confidence interval
#' @param p_h numeric value or vector of values for expected/planned ICC
#' @param z numeric value indicating the z score value for a given confidence interval (e.g., 1.96 for 95% CI)
#' @param type string input taking form of "single" or "multiple".
#' @return dataframe consisting of two values (type = "single") for k and K1, or large dataframe with k and k1 values for a range of w and p_h values.
#' @keywords reliability, intra-class correlation coefficient
#' @export
#'


reliability_sample_size_estimation <- function(w,
                                               n,
                                               z,
                                               p_h,
                                               type){

  if(type == "single"){

    k <- (8 * (z ^ 2)) * ((1 - p_h) ^ 2) *
      (1 + ((n - 1) * p_h)) ^ 2 /
      ((w ^ 2) * n * (n - 1))

    return(data.frame(k = round(k),
                      k1 = round(k + 1)))

  }

  if(type == "multiple"){

    ## Varying ICC and width of confidence interval, each by unit increments of 0.01

    # This will expand the grid of values so that there is now a dataframe consisting of a row for each unique pair of ICC and width
    w_and_p_grid <- expand.grid(w = w, p_h = p_h)

    # Use sapply to loop over the dataframe, and for each row (i.e., unique ICC and width estimates), run the reliability function, taking those values as estimates
    # This would also be possible if you wanted to vary replicate number too. But for these purposes, we are sticking with 3.

    sample_sizes <- data.frame(k = NA, k1 = NA)

    for(i in 1: nrow(w_and_p_grid)){

      sample_sizes[i, 1] <-

        (8 * (z ^ 2)) * ((1 - p_h) ^ 2) *
        (1 + ((n - 1) * p_h)) ^ 2 /
        ((w ^ 2) * n * (n - 1))

      sample_sizes[i, 2] <- sample_sizes[i, 1] + 1


    }


    # Unlist the sample sizes and take every 2nd element corresponding to either the k or k1 value.
    # Save to a new dataframe
    result_df <- data.frame(w = w_and_p_grid$w,
                            p_h = w_and_p_grid$p_h,
                            k = unlist(sample_sizes)[seq(1, length(sample_sizes), 2)],
                            k1 = unlist(sample_sizes)[seq(2, length(sample_sizes), 2)])

    return(result_df)

  }




}

