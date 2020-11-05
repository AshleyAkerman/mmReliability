
#' Title: reliability statistics calculator (from Hopkins 2015)
#'
#' This function is used in the reliability extraction function to generate the estimates of the various reliability statistics
#' The output is a large table with various parameters in columns and each outcome measure in rows
#' Prior to its use, one must generate a long-form and wide-form df of the data, and a vector called variable_names that stores the variables of interest
#'
#' Analysis and calculations derived from pairwise excel sheet: https://www.sportsci.org/resource/stats/rerlycalc.html
#' Hopkins WG (2015). Spreadsheets for analysis of validity and reliability. Sportscience 19, 36-42 (sportsci.org/2015/ValidRely.htm)
#'
#' The excel sheet provides a readily accessible, widely-used resource for individuals running reliability analysis.
#' The intent of converting to R function was for ease (rather than running multiple outcomes over the same excel sheet), but also to allow for greater transparency with respect to the underlying math within  the excel formulas and conditional formatting.
#' This has also provided a useful teaching tool for walking through the steps of reliability analysis with students.
#'
#' Code will work for cases of missing values. e.g., if n participants have only two repeat measurements rather than three.
#' Overall purpose: To get the average reliability of the measurements, instead of averaging the reliability of consecutive trials, one must weight their (individual estimates) squares by the degrees of freedom and take the square root.
#'
#' There are instances where the CIs were not required for the manuscripts, so will be 'commented out'. If an individual wished to use these parameters, just remove the '#'
#' If working with a Mac, this would be achieved quickly by highlighting the required code and using cmd+shift+c
#' In all cases, useable data does not contain a space between the '#' and the variable/data
#' On the other hand, all titles contain a space between the '#' and the text
#'
#' NOTES:
#' The reliability function was written to only allow for 3 independent repeat measures (i.e., consistent with our design and the maximum often completed in the field).
#' The full R implementation, i.e., allowing for more than 3 replicates is almost complete and will be uploaded to OSF as a separate script when complete.
#'
#' The script also mixes base and tidyverse functionality. As with most of the R language, there are multiple ways to achieve the same outcome.
#' The choice of base or tidyverse was chosen based on ease of achieving a given outcome when writing the script.
#' Similarly, and often considered blasphemy in R, for/while/if loops were used in some instances where apply/lapply/sapply function are also possible.
#' Often this was for ease of teaching/demonstrating the calculations for an individual not accustomed to any program language.
#'
#' @param df dataframe input in wide form
#' @return dataframe of the calculated reliability statistics


#' @keywords reliability, intra-class correlation coefficient, coefficient of variation
#' @export

reliability_function <- function(df) {

  ##------##
  ## BASIC STATISTICS: SAMPLE SIZE AND DEGREES OF FREEDOM (individual and combined)
  ##------##

  # Derived values for means, sds, and sizes
  Df1means <- apply(df, 2, mean, na.rm = TRUE)
  Df1sds <- apply(df, 2, sd, na.rm = TRUE)

  # Sample sizes
  # How many values are not NA. Could also achieve with:   length(df$Trial1[!is.na(df$Trial1)])
  Trial1n <- sum((!is.na(df$Trial1)) == "TRUE")
  Trial2n <- sum((!is.na(df$Trial2)) == "TRUE")
  Trial3n <- sum((!is.na(df$Trial3)) == "TRUE")

  Trial1nLOG <- sum((!is.na(df$Trial1log)) == "TRUE")
  Trial2nLOG <- sum((!is.na(df$Trial2log)) == "TRUE")
  Trial3nLOG <- sum((!is.na(df$Trial3log)) == "TRUE")

  # Generate the sumproduct function from excel that determines how many different participants are contributing to the overall mean
  # This would work with an "apply" loop

  df2 <- df # saved as a new dataframe for ease of debugging

  df2[is.na(df2)] <- 0    # Adjust so that the sumproduct function takes a missing value as a zero rather than NA

  # 'For' loop to derive the samples contributing to overall mean
  output <- vector(length = nrow(df2))
  for (i in 1:nrow(df2)) {
    output[i] <- (df2$Trial1[i] * df2$Trial1[i] * Trial1n +
                    df2$Trial2[i] * df2$Trial2[i] * Trial2n +
                    df2$Trial3[i] * df2$Trial3[i] * Trial3n)
  }

  # 'For' loop to derive the samples contributing to overall mean (LOG data; this is unlikely to be different to calculation above, but provided for completeness)
  outputLOG <- vector(length = nrow(df2))
  for (i in 1:nrow(df2)) {
    outputLOG[i] <- (df2$Trial1log[i] * df2$Trial1log[i] * Trial1nLOG +
                       df2$Trial2log[i] * df2$Trial2log[i] * Trial2nLOG +
                       df2$Trial3log[i] * df2$Trial3log[i] * Trial3nLOG)
  }

  # Not the mean number of subjects, but number of indpendent samples making up the mean
  meanTrialn <- length(output[output > 0])
  meanTrialnLOG <- length(outputLOG[outputLOG > 0])

  # Since the calculations are derived from consecutive pairwise comparisons (i.e., 1vs2, and 2vs3), these comparison will be denoted throughout as '12' and '23' respectively.

  # Dervied values for sizes and df
  Trial12n <- sum((!is.na(df$RawDiff12)) == "TRUE")
  Trial23n <- sum((!is.na(df$RawDiff23)) == "TRUE")

  Trial12nLOG <- sum((!is.na(df$LogDiff12)) == "TRUE")
  Trial23nLOG <- sum((!is.na(df$LogDiff23)) == "TRUE")

  Trial1df <- Trial1n - 1
  Trial2df <- Trial2n - 1
  Trial3df <- Trial3n - 1

  Trial1dfLOG <- Trial1nLOG - 1
  Trial2dfLOG <- Trial2nLOG - 1
  Trial3dfLOG <- Trial3nLOG - 1

  Trial12df <- Trial12n - 1
  Trial23df <- Trial23n - 1

  Trial12dfLOG <- Trial12nLOG - 1
  Trial23dfLOG <- Trial23nLOG - 1

  ##------##
  # BASIC STATISTICS: Grand means, SDs, and statistics used for CI calculations (individual and combined)
  ##------##

  # Used for calculations of CI's and averaged reliability statistics
  GrandMean <- ((Df1means[[2]] * Trial1n) + (Df1means[[3]] * Trial2n) +
                  (Df1means[[4]] * Trial3n)) / (Trial1n + Trial2n + Trial3n)

  GrandMeanLOG <- ((Df1means[[7]] * Trial1n) + (Df1means[[8]] * Trial2n) +
                     (Df1means[[9]] * Trial3n)) / (Trial1n + Trial2n + Trial3n)

  GrandSD <- sqrt((Df1sds[[2]] * Df1sds[[2]] * Trial1df +
                     Df1sds[[3]] * Df1sds[[3]] * Trial2df +
                     Df1sds[[4]] * Df1sds[[4]] * Trial3df) /
                    (Trial1df + Trial2df + Trial3df))

  GrandSDLOG <- sqrt((Df1sds[[7]] * Df1sds[[7]] * Trial1df +
                        Df1sds[[8]] * Df1sds[[8]] * Trial2df +
                        Df1sds[[9]] * Df1sds[[9]] * Trial3df) /
                       (Trial1df + Trial2df + Trial3df))

  # Effective number of trials
  ENOT12 <- 1 + Trial12df/(Trial1n + Trial2n - Trial12n - 1)  # Effective number of trials for 1 vs 2
  ENOT23 <- 1 + Trial23df/(Trial2n + Trial3n - Trial23n - 1)  # Effective number of trials for 2 vs 3

  # 'Adjusted for independence of consecutive statistics and for any missing values'
  # The '3' is because there are only 3 trials. This would need to be adjusted for any other number of trials
  CombDF <- (1 - (0.22 * (Trial1n + Trial2n + Trial3n)) /
               (3 * meanTrialn)) * (Trial12df + Trial23df)

  CombDFLOG <- (1 - (0.22 * (Trial1nLOG + Trial2nLOG + Trial3nLOG)) /
                  (3 * meanTrialnLOG)) * (Trial12dfLOG + Trial23dfLOG)

  ENOTcomb <- 1 + CombDF/(meanTrialn - 1)

  # Sum the product of a matrix of:
  #      (SD of the change from T2 - T1) * (SD of the change from T2 - T1) * (df for T2-T1)
  #      (SD of the change from T3 - T1) * (SD of the change from T3 - T1) * (df for T3-T2)
  # Then divide by the sum of the df's

  SDchange <- sqrt((Df1sds[[5]] * Df1sds[[5]] * Trial12df + Df1sds[[6]] * Df1sds[[6]] * Trial23df) /
                     (Trial12df + Trial23df))

  SDchangeLOG <- sqrt((Df1sds[[10]] * Df1sds[[10]] * Trial12df + Df1sds[[11]] * Df1sds[[11]] * Trial23df) /
                        (Trial12df + Trial23df))

  SDchangeLOGasFactor <- exp((SDchangeLOG / 100))

  SDchangeLOGasCV <- (100 * (exp((SDchangeLOG / 100)))) - 100

  ##------------------------------------------------------------------##
  ##------------------------------------------------------------------##
  ## Measures of Reliability via Raw Values ##
  ##------------------------------------------------------------------##
  ##------------------------------------------------------------------##

  ##----##
  ##----##

  ## Change in Mean and 95% CI##

  # NOTE: The change in mean is divided by the "pure" between-subject SD (or: sqrt(SD ^ 2 - sd ^ 2))
  # SD is the mean SD for the trials included in the analysis (derived from averaging the variances, weighted by the degrees of freedom)
  # sd is the typical error

  ##----##
  ##----##

  # All 95% confidence intervals are denoted as 'lCL' and 'uCL' for lower and upper confidence limits respectively.

  # Need to take minus 1 and divide by 2 for the CL to get the same result at TINV in excel
  ChMean12lCL <- Df1means[[5]] - stats::qt(1-(1- CL/100)/2, Trial12df) * Df1sds[[5]] / sqrt(Trial12n)
  ChMean12uCL <- Df1means[[5]] + stats::qt(1-(1- CL/100)/2, Trial12df) * Df1sds[[5]] / sqrt(Trial12n)

  # Need to take minus 1 and divide by 2 for the CL to get the same result at TINV in excel
  ChMean23lCL <- Df1means[[6]] - stats::qt(1-(1- CL/100)/2, Trial23df) * Df1sds[[6]] / sqrt(Trial23n)
  ChMean23uCL <- Df1means[[6]] + stats::qt(1-(1- CL/100)/2, Trial23df) * Df1sds[[6]] / sqrt(Trial23n)

  ##----##
  ##----##

  ## Raw Typical Error and 95% CI

  # NOTE: The TE and CIs are divided by the "pure" between-subject SD (or: sqrt(SD ^ 2 - sd ^ 2))
  # SD is the mean SD for the trials included in the analysis (derived from averaging the variances, weighted by the degrees of freedom.  sd is the typical error)

  # Confidence intervals around a variation (or otherwise)
  # Need to use the chi-squared distribution to find the critical values associated with a probability of 0.025 (for one tail of a 95% CI) and 0.975 (for the other tail of a 95% CI)

  # The equation for finding the CI around a population SD is:
  # (n-1)s^2 / X2 (alpha/2) = Lower CI tail
  # (n-2)s^2 / X2 (1-(alpha/2)) = Upper CI tail
  # Where n-1 is the df, and X2 is the chi-squared distribution at either 0.025 (LCI) or 0.975 (UCI) to equal a 95% CL.

  ##----##
  ##----##

  TE12 <- Df1sds[[5]] / sqrt(2)
  # sqrt(df * TE^2/chiinv((1 - CL/100)/2, df)) - excel formula
  TE12lCL <- sqrt(((Trial12df * (TE12 ^ 2)) / stats::qchisq(p = 1-(1 - CL/100)/2, df = Trial12df)))
  # sqrt(df * TE^2/chiinv((1-(1 - CL/100))/2, df)) - excel formula
  TE12uCL <- sqrt((((Trial12df * (TE12 ^ 2)) / stats::qchisq(p = (1 - CL/100)/2, df = Trial12df))))


  TE23 <- Df1sds[[6]] / sqrt(2)
  # sqrt(df * TE^2/chiinv((1 - CL/100)/2, df)) - excel formula
  TE23lCL <- sqrt(((Trial23df * (TE23 ^ 2)) / stats::qchisq(p = 1-(1 - CL/100)/2, df = Trial23df)))
  # sqrt(df * TE^2/chiinv((1-(1 - CL/100))/2, df)) - excel formula
  TE23uCL <- sqrt((((Trial23df * (TE23 ^ 2)) / stats::qchisq(p = (1 - CL/100)/2, df = Trial23df))))

  # Combined calculation
  TEcomb <- sqrt((TE12 * TE12 * Trial12df + TE23 * TE23 * Trial23df) /
                   (Trial12df + Trial23df))

  # These estimations for the CI are different to the ones provided by the excel sheet (see link above)
  # This is due - I think - to the different chisq value provided by the CHIINV function (excel) and the qchisq function (R) when the degrees of freedom are non-integers
  # Currently the chisq inverse on excel (p = 0.025, df = 10.92) equals 20.4831774; R equals 21.80588
  # Currently the chisq invest on excel (p = 0.975, df = 10.92) equals 3.2469728; R equals 3.769503

  # The R answer is corroborated on: homepage.divms.uiowa.edu/~mbognar/applets/chisq.html
  # But, excel is corroborated on: fourmilab.ch/rrpkp/experriments/analysis/chiCalc.html
  # Therefore the CI estimations here will be slightly more liberal than the estimation from the excel sheet (i.e., tighter limits)

  TEcomblCL <- sqrt((CombDF * TEcomb ^ 2) /
                      stats::qchisq(p = ((1 - CL/100)/2),
                             df = CombDF, lower.tail = FALSE))

  TEcombuCL <- sqrt((CombDF * TEcomb ^ 2) /
                      stats::qchisq(p = (1 - (1 - CL/100)/2),
                             df = CombDF, lower.tail = FALSE))

  TE12MULTIPLYDIVIDED <- sqrt(TE12uCL/TE12lCL)
  TE23MULTIPLYDIVIDED <- sqrt(TE23uCL/TE23lCL)

  TE12BiasCF <- 1 + (1 / (4 * Trial12df))
  TE23BiasCF <- 1 + (1 / (4 * Trial23df))

  TEcombMULTIPLYDIVIDED <- sqrt(TEcombuCL/TEcomblCL)
  TEcombBiasCF <- 1 + (1 / (4 * CombDF))

  ##------------------------------------------------------------------##
  ##------------------------------------------------------------------##
  ## Measures of Reliability via Standardised Values (ICC > 0) ##

  # NOTE: If ICC is zero or negative then the standardized values should not be calculated because the pure SD is effectively zero
  ##------------------------------------------------------------------##
  ##------------------------------------------------------------------##

  # The equivalent of the "sumproduct" matrix in the excel sheet.
  # Shortened and saved here for convenience
  matrix12calc <- (Df1sds[[2]] * Df1sds[[2]] * Trial1df +
                     Df1sds[[3]] * Df1sds[[3]] * Trial2df)

  matrix23calc <- (Df1sds[[3]] * Df1sds[[3]] * Trial2df +
                     Df1sds[[4]] * Df1sds[[4]] * Trial3df)

  matrix12calcLOG <- (Df1sds[[7]] * Df1sds[[7]] * Trial1dfLOG +
                        Df1sds[[8]] * Df1sds[[8]] * Trial2dfLOG)

  matrix23calcLOG <- (Df1sds[[8]] * Df1sds[[8]] * Trial2dfLOG +
                        Df1sds[[9]] * Df1sds[[9]] * Trial3dfLOG)

  ##----##
  ##----##
  ## Pearson correlations Calculations # 1

  # NOTE: Corrected for small-sample bias using the factor at the bottom of the 2nd section below
  ##----##
  ##----##

  Pearson12 <- stats::cor(df$Trial1, df$Trial2, use = "complete.obs") *
    (1+ (1- (stats::cor(df$Trial1, df$Trial2, use = "complete.obs") ^ 2)) / 2 / (Trial12n - 3))

  Pearson23 <- stats::cor(df$Trial2, df$Trial3, use = "complete.obs") *
    (1+ (1- (stats::cor(df$Trial2, df$Trial3, use = "complete.obs") ^ 2)) / 2 / (Trial23n - 3))

  # Using fisher transformation And inverse of fisher transformation
  # See: https://www.personality-project.org/r/psych/help/fisherz.html

  FisherR12 <- psych::fisherz(Pearson12)
  FisherR23 <- psych::fisherz(Pearson23)
  FisherSE12 <- (Trial12n - 3)
  FisherSE23 <- (Trial23n - 3)

  FisherRcomb <- (FisherR12 * FisherSE12 + FisherR23 * FisherSE23) / (FisherSE12 + FisherSE23)

  PearsonComb <- psych::fisherz2r(FisherRcomb)

  ##----##
  ##----##
  ## ICC calculations

  # NOTE: The 1-sd ^ 2 / SD ^ 2, where sd is the TE and SD is the mean b/w subject SD in the two trials (derived by weighting the variances by their degrees of freedom)
  # NOTE: Also corrected for small-sample bias using the same factor as the Pearson.
  ##----##
  ##----##

  ICC12 <- (1 - TE12 ^ 2 /(matrix12calc/ (Trial1df + Trial2df))) *
    (1 + (1- (1 - TE12 ^ 2 /(matrix12calc/ (Trial1df + Trial2df))) ^ 2) /
       (Trial1n + Trial2n - Trial12n - 3))

  ICC23 <- (1 - TE23 ^ 2 /(matrix23calc/ (Trial2df + Trial3df))) *
    (1 + (1- (1 - TE23 ^ 2 /(matrix23calc/ (Trial2df + Trial3df))) ^ 2) /
       (Trial2n + Trial3n - Trial23n - 3))

  ICCcomb <- (1 - TEcomb ^ 2 / GrandSD ^ 2) *
    (1 + (1 - (1 - TEcomb ^ 2/ GrandSD ^ 2) ^ 2) /
       (meanTrialn - 3))

  # CI's using F values for ICC
  F12 <- 1 + (ICC12 * ENOT12) / (1 - ICC12)
  F23 <- 1 + (ICC23 * ENOT23) / (1 - ICC23)
  CombF <- 1 + (ICCcomb * ENOTcomb) / (1 - ICCcomb)

  # Numerator degrees of freedom for F-test
  numDf12 <- Trial1n + Trial2n - Trial12n - 1
  numDf23 <- Trial2n + Trial3n - Trial23n - 1
  numerDfcomb <- meanTrialn - 1

  # Denominator degreees of freedom for F-test
  denomDf12 <- Trial12df
  denomDf23 <- Trial23df
  denomDfcomb <- CombDF

  # F values
  fLower12 <- F12 / stats::qf(p = (1- CL/100)/2,
                       df1 = numDf12, df2 = denomDf12,
                       lower.tail = FALSE)

  fLower23 <- F23 / stats::qf(p = (1- CL/100)/2,
                       df1 = numDf23, df2 = denomDf23,
                       lower.tail = FALSE)

  fLowercomb  <- CombF / stats::qf(p = (1- CL/100)/2,
                            df1 = numerDfcomb, df2 = denomDfcomb,
                            lower.tail = FALSE)

  fUpper12 <- F12 * stats::qf(p = (1- CL/100)/2,
                       df1 = denomDf12, df2 = numDf12,
                       lower.tail = FALSE)

  fUpper23 <- F23 * stats::qf(p = (1- CL/100)/2,
                       df1 = denomDf23, df2 = numDf23,
                       lower.tail = FALSE)

  fUppercomb <- CombF * stats::qf(p = (1- CL/100)/2,
                           df1 = denomDfcomb, df2 = numerDfcomb,
                           lower.tail = FALSE)

  # Values for the combined upper and lower are different from excel sheet
  # This is likely due to a similar problem as the upper and lower CI from above (related to chiinv; see TE)

  ICCcomblCL <- (fLowercomb - 1) / ((fLowercomb + ENOTcomb) - 1)
  ICCcombuCL <- (fUppercomb - 1) / ((fUppercomb + ENOTcomb) - 1)

  ICC12lCL <- (fLower12 - 1) / ((fLower12 + ENOT12) - 1)
  ICC12uCL <- (fUpper12 - 1) / ((fUpper12 + ENOT12) - 1)

  ICC23lCL <- (fLower23 - 1) / ((fLower23 + ENOT23) - 1)
  ICC23uCL <- (fUpper23 - 1) / ((fUpper23 + ENOT23) - 1)

  ##----##
  ##----##
  ## Variance Calculations
  ##----##
  ##----##

  # Pure variance calculations and SE of PV
  PureVariance12 <- (matrix12calc / (Trial1df + Trial2df)) - (TE12 ^ 2)
  PureVariance23 <- (matrix23calc / (Trial2df + Trial3df)) - (TE23 ^ 2)
  PureVarianceComb <- (GrandSD ^ 2) - (TEcomb ^ 2)

  SEPureVariance12 <- sqrt(2 * ((((matrix12calc / (Trial1df + Trial2df)) ^ 2) / numDf12) -
                                  (TE12 ^ 4 / denomDf12)))

  SEPureVariance23 <- sqrt(2 * ((((matrix23calc / (Trial2df + Trial3df)) ^ 2) / numDf23) -
                                  (TE23 ^ 4 / denomDf23)))

  SEPureVarianceComb <- sqrt(((2 * ((GrandSD ^ 4) / numerDfcomb)) -
                                ((TEcomb ^ 4) / denomDfcomb)))

  SEPV12lCL <- PureVariance12 + (stats::qnorm((1- CL/100)/2) * SEPureVariance12)
  SEPV23lCL <- PureVariance23 + (stats::qnorm((1- CL/100)/2) * SEPureVariance23)
  SEPVcomblCL <- PureVarianceComb + (stats::qnorm((1- CL/100)/2) * SEPureVarianceComb)

  SEPV12uCL <- PureVariance12 - (stats::qnorm((1- CL/100)/2) * SEPureVariance12)
  SEPV23uCL <- PureVariance23 - (stats::qnorm((1- CL/100)/2) * SEPureVariance23)
  SEPVcombuCL <- PureVarianceComb - (stats::qnorm((1- CL/100)/2) * SEPureVarianceComb)

  ##----##
  ##----##
  ## Smallest effect from pure SD (RAW VALUES)

  # NOTE: Based on a cohen's d of 0.2, 0.5, and 0.8 - entered above
  # NOTE: Given that any number multiplied by itself equals a positive number, a negative number cannot have a sqrt.
  # In the instance below this is solved by providing the minus before the number (i.e., to make positive)

  # The code for providing the confidence intervals (for one of the effect sizes) is provided for clarity, but not run in this instance
  ##----##
  ##----##

  #SmallestEffpureSD12 <- eff_interest_sml * (sqrt(PureVariance12))
  #SmallestEffpureSD12lCL <- eff_interest_sml * ifelse(SEPV12lCL > 0, (sqrt(SEPV12lCL)), -(sqrt(-SEPV12lCL)))
  #SmallestEffpureSD12uCL <- eff_interest_sml * (sqrt(SEPV12uCL))

  #SmallestEffpureSD23 <- eff_interest_sml * (sqrt(PureVariance23))
  #SmallestEffpureSD23lCL <- eff_interest_sml * ifelse(SEPV23lCL > 0, (sqrt(SEPV23lCL)), -(sqrt(-SEPV23lCL)))
  #SmallestEffpureSD23uCL <- eff_interest_sml * (sqrt(SEPV23uCL))

  #SmallestEffpureSDcomb <- eff_interest_sml * (sqrt(PureVarianceComb))
  #SmallestEffpureSDcomblCL <- eff_interest_sml * ifelse(SEPVcomblCL > 0, (sqrt(SEPVcomblCL)), -(sqrt(-SEPVcomblCL)))
  #SmallestEffpureSDcombuCL <- eff_interest_sml * (sqrt(SEPVcombuCL))

  #MedEffectpureSD12 <- eff_interest_med * (sqrt(PureVariance12))
  #MedEffectpureSD23 <- eff_interest_med * (sqrt(PureVariance23))
  #MedEffectpureSDcomb <- eff_interest_med * (sqrt(PureVarianceComb))

  #LrgEffectpureSD12 <- eff_interest_lrg * (sqrt(PureVariance12))
  #LrgEffectpureSD23 <- eff_interest_lrg * (sqrt(PureVariance23))
  #LrgEffectpureSDcomb <- eff_interest_lrg * (sqrt(PureVarianceComb))

  ##----##
  ##----##
  ## Smallest effect from observed SD (RAW VALUES)

  # NOTE: Based on a cohen's d of 0.2, 0.5, and 0.8 - entered below

  # The code for providing the confidence intervals (for one of the effect sizes) is provided for clarity, but not run in this instance
  ##----##
  ##----##

  #SmallestEffobsSD12 <- eff_interest_sml * (sqrt(matrix12calc/ (Trial1df + Trial2df)))
  #SmallestEffobsSD12lCL <- sqrt((numDf12 * (SmallestEffobsSD12 ^ 2)) / qchisq(p = ((1 - CL/100)/2), df = numDf12, lower.tail = FALSE))
  #SmallestEffobsSD12uCL <- sqrt((numDf12 * (SmallestEffobsSD12 ^ 2)) / qchisq(p = (1 - (1 - CL/100)/2), df = numDf12, lower.tail = FALSE))

  #SmallestEffobsSD23 <- eff_interest_sml * (sqrt(matrix23calc/ (Trial2df + Trial3df)))
  #SmallestEffobsSD23lCL <- sqrt((numDf23 * (SmallestEffobsSD23 ^ 2)) / qchisq(p = ((1 - CL/100)/2), df = numDf23, lower.tail = FALSE))
  #SmallestEffobsSD23uCL <- sqrt((numDf23 * (SmallestEffobsSD23 ^ 2)) / qchisq(p = (1 - (1 - CL/100)/2), df = numDf23, lower.tail = FALSE))

  #SmallestEffobsSDcomb <- eff_interest_sml * GrandSD
  #SmallestEffobsSDcomblCL <- sqrt((numerDfcomb * (SmallestEffobsSDcomb ^ 2)) / qchisq(p = ((1 - CL/100)/2), df = numerDfcomb, lower.tail = FALSE))
  #SmallestEffobsSDcombuCL <- sqrt((numerDfcomb * (SmallestEffobsSDcomb ^ 2)) / qchisq(p = (1 - (1 - CL/100)/2), df = numerDfcomb, lower.tail = FALSE))

  #MedEffectobsSD12 <- eff_interest_med * (sqrt(matrix12calc/ (Trial1df + Trial2df)))
  #MedEffectobsSD23 <- eff_interest_med * (sqrt(matrix12calc/ (Trial2df + Trial3df)))
  #MedEffectobsSDcomb <- eff_interest_med * GrandSD

  #LrgEffectobsSD12 <- eff_interest_lrg * (sqrt(matrix12calc/ (Trial1df + Trial2df)))
  #LrgEffectobsSD23 <- eff_interest_lrg * (sqrt(matrix12calc/ (Trial2df + Trial3df)))
  #LrgEffectobsSDcomb <- eff_interest_lrg * GrandSD

  ##----##
  ##----##
  ## Pearson correlations Calculations # 2

  # NOTE: Corrected for small-sample bias using the factor at the bottom of this section
  ##----##

  Pearson12lCL <- (exp((2* ((0.5 * log((1 + Pearson12) / (1 - Pearson12))) - stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial12n - 3))))) - 1) /
    (exp((2* ((0.5 * log((1 + Pearson12) / (1 - Pearson12))) - stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial12n - 3))))) + 1)

  Pearson23lCL <- (exp((2* ((0.5 * log((1 + Pearson23) / (1 - Pearson23))) - stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial23n - 3))))) - 1) /
    (exp((2* ((0.5 * log((1 + Pearson23) / (1 - Pearson23))) - stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial23n - 3))))) + 1)

  Pearson12uCL <- (exp((2* ((0.5 * log((1 + Pearson12) / (1 - Pearson12))) + stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial12n - 3))))) - 1) /
    (exp((2* ((0.5 * log((1 + Pearson12) / (1 - Pearson12))) + stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial12n - 3))))) + 1)

  Pearson23uCL <- (exp((2* ((0.5 * log((1 + Pearson23) / (1 - Pearson23))) + stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial23n - 3))))) - 1) /
    (exp((2* ((0.5 * log((1 + Pearson23) / (1 - Pearson23))) + stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial23n - 3))))) + 1)

  Pearson12BiasCF <- 1 + ((1 - (stats::cor(df$Trial1, df$Trial2, use = "complete.obs") ^ 2)) / 2) / (Trial12n - 3)
  Pearson23BiasCF <- 1 + ((1 - (stats::cor(df$Trial2, df$Trial3, use = "complete.obs") ^ 2)) / 2) / (Trial23n - 3)

  ##----##
  ## Standardised values for Change in mean ##
  ##----##

  ChMean12standardized <- Df1means[[5]] / sqrt(((matrix12calc/ (Trial1df + Trial2df)) - (TE12 ^ 2)))
  ChMean12standardizedlCL <- ChMean12lCL / sqrt(((matrix12calc/ (Trial1df + Trial2df)) - (TE12 ^ 2)))
  ChMean12standardizeduCL <- ChMean12uCL / sqrt(((matrix12calc/ (Trial1df + Trial2df)) - (TE12 ^ 2)))
  ChMean12standardizedPLUSMINUS <- (ChMean12standardizeduCL - ChMean12standardizedlCL) / 2

  ChMean23standardized <- Df1means[[6]] / sqrt(((matrix23calc/ (Trial2df + Trial3df)) - (TE23 ^ 2)))
  ChMean23standardizedlCL <- ChMean23lCL  / sqrt(((matrix23calc/ (Trial2df + Trial3df)) - (TE23 ^ 2)))
  ChMean23standardizeduCL <- ChMean23uCL / sqrt(((matrix23calc/ (Trial2df + Trial3df)) - (TE23 ^ 2)))
  ChMean23standardizedPLUSMINUS <- (ChMean23standardizeduCL - ChMean23standardizedlCL) / 2

  ##----##
  ## Standardised values for Typical Error ##

  # NOTE: Standardized TE should be doubled to interpret its magnitude using the thresholds for change in the mean
  ##----##

  TE12standardized <- TE12 / sqrt(((matrix12calc/ (Trial1df + Trial2df)) - (TE12 ^ 2)))
  TE12standardizedlCL <- TE12lCL / sqrt(((matrix12calc/ (Trial1df + Trial2df)) - (TE12 ^ 2)))
  TE12standardizeduCL <- TE12uCL / sqrt(((matrix12calc/ (Trial1df + Trial2df)) - (TE12 ^ 2)))
  TE12standardizedMULTIPLYDIVIDED <- sqrt(TE12standardizeduCL/TE12standardizedlCL)
  TE12standardizedBiasCF <- TE12BiasCF

  TE23standardized <- TE23 / sqrt(((matrix23calc/ (Trial2df + Trial3df)) - (TE23 ^ 2)))
  TE23standardizedlCL <- TE23lCL / sqrt(((matrix23calc/ (Trial2df + Trial3df)) - (TE23 ^ 2)))
  TE23standardizeduCL <- TE23uCL / sqrt(((matrix23calc/ (Trial2df + Trial3df)) - (TE23 ^ 2)))
  TE23standardizedMULTIPLYDIVIDED <- sqrt(TE23standardizeduCL/TE23standardizedlCL)
  TE23standardizedBiasCF <- TE23BiasCF

  TEcombstandardized <- TEcomb/ (sqrt((GrandSD ^ 2)))
  TEcombstandardizedlCL <- TEcomblCL / (sqrt(GrandSD ^ 2))
  TEcombstandardizeduCL <- TEcombuCL / (sqrt(GrandSD ^ 2))
  TE23standardizedMULTIPLYDIVIDED <- sqrt(TEcombstandardizeduCL/TEcombstandardizedlCL)
  TE23standardizedBiasCF <- TEcombBiasCF

  ##------------------------------------------------------------------##
  ##------------------------------------------------------------------##

  ## Measures of Reliability via LOG-TRANSFORMED DATA ###

  # Use log transformation if:
  #   - data can only be positive and non-zero
  #   - between-subject SD or change-score SD is more than 20% of the mean
  #   - magnitude of the percent of the factor error is expected to be similar for all subjects

  ##------------------------------------------------------------------##
  ##------------------------------------------------------------------##

  # function to take backtransform each value and divide by 100
  BackTrans <- function(x) {
    exp((x/100))
  }

  # Only interested in the last 5 elements
  BackTransfMeans <- unlist(lapply(Df1means[7:11], BackTrans))
  SDasFactor <- unlist(lapply(Df1sds[7:11], BackTrans))

  # Function to convert the SD to a coefficient
  SDtoCV <- function(x) {
    (100 * x) - 100
  }

  SDasCV <- unlist(lapply(SDasFactor, SDtoCV))

  BackTransfGrandMeanLOG <- exp((GrandMeanLOG/100))
  BackTransfGrandSDLOG <- exp((GrandSDLOG/100))
  SDasCVcomb <- (100 * BackTransfGrandSDLOG) - 100

  # Percent change in mean
  ChMean12PercentLOG <- (100 * (exp((Df1means[[10]]/100)))) - 100
  ChMean23PercentLOG <- (100 * (exp((Df1means[[11]]/100)))) - 100

  # Change in mean and CI
  ChMean12LOGlCL <- Df1means[[10]] - stats::qt(1-(1- CL/100)/2, (Trial12nLOG - 1)) * Df1sds[[10]] / sqrt(Trial12nLOG)
  ChMean12LOGuCL <- Df1means[[10]] + stats::qt(1-(1- CL/100)/2, (Trial12nLOG - 1)) * Df1sds[[10]] / sqrt(Trial12nLOG)

  ChMean23LOGlCL <- Df1means[[11]] - stats::qt(1-(1- CL/100)/2, (Trial23nLOG - 1)) * Df1sds[[11]] / sqrt(Trial23nLOG)
  ChMean23LOGuCL <- Df1means[[11]] + stats::qt(1-(1- CL/100)/2, (Trial23nLOG - 1)) * Df1sds[[11]] / sqrt(Trial23nLOG)

  ##----##
  ##----##
  ## Typical Error
  ##----##
  ##----##

  TE12LOG <- Df1sds[[10]] / (sqrt(2))
  TE23LOG <- Df1sds[[11]] / (sqrt(2))

  TE12lCLLOG <- sqrt(((Trial12dfLOG * (TE12LOG ^ 2)) / stats::qchisq(p = 1-(1 - CL/100)/2, df = Trial12dfLOG)))
  TE12uCLLOG <- sqrt((((Trial12dfLOG * (TE12LOG ^ 2)) / stats::qchisq(p = (1 - CL/100)/2, df = Trial12dfLOG))))

  TE23lCLLOG <- sqrt(((Trial23dfLOG * (TE23LOG ^ 2)) / stats::qchisq(p = 1-(1 - CL/100)/2, df = Trial23dfLOG)))
  TE23uCLLOG <- sqrt((((Trial23dfLOG * (TE23LOG ^ 2)) / stats::qchisq(p = (1 - CL/100)/2, df = Trial23dfLOG))))

  TEcombLOG <- sqrt((TE12LOG * TE12LOG * Trial12dfLOG + TE23LOG * TE23LOG * Trial23dfLOG) /
                      (Trial12dfLOG + Trial23dfLOG))

  TEcomblCLLOG <- sqrt((CombDFLOG * TEcombLOG ^ 2) / stats::qchisq(p = ((1 - CL/100)/2), df = CombDFLOG, lower.tail = FALSE))

  TEcombuCLLOG <- sqrt((CombDFLOG * TEcombLOG ^ 2) / stats::qchisq(p = (1 - (1 - CL/100)/2), df = CombDFLOG, lower.tail = FALSE))

  # Bias correction factors
  TE12BiasCFLOG <- 1 + (1 / (4 * Trial2dfLOG))
  TE23BiasCFLOG <- 1 + (1 / (4 * Trial23dfLOG))

  TEcombBiasCFLOG <- 1 + (1 / (4 * CombDFLOG))

  ##----##
  ##----##
  ## Intraclass r and associated calculations
  ##----##
  ##----##

  Intraclassr12LOG <- (1 - (TE12LOG ^ 2) / (matrix12calcLOG / (Trial1dfLOG + Trial2dfLOG))) *
    (1 + ((1 - ((1 - (TE12LOG ^ 2) / (matrix12calcLOG / (Trial1dfLOG + Trial2dfLOG))) ^ 2))) / (Trial1nLOG + Trial2nLOG - Trial12nLOG - 3))

  Intraclassr23LOG <- (1 - (TE23LOG ^ 2) / (matrix23calcLOG / (Trial2dfLOG + Trial3dfLOG))) *
    (1 + ((1 - ((1 - (TE23LOG ^ 2) / (matrix23calcLOG / (Trial2dfLOG + Trial3dfLOG))) ^ 2))) / (Trial2nLOG + Trial3nLOG - Trial23nLOG - 3))

  IntraclassrCombLOG <- (1 - ((TEcombLOG ^ 2)) / (GrandSDLOG ^ 2)) * (1 + ((1 - ((1 - ((TEcombLOG ^ 2)) / (GrandSDLOG ^ 2)) ^ 2)) / (meanTrialnLOG - 3)))

  ##----##
  ##----##
  ## Effective number of trials
  ##----##
  ##----##

  ENOT12LOG <- 1 + Trial12dfLOG/(Trial1nLOG + Trial2nLOG - Trial12nLOG - 1)  # Effective number of trials for 1 vs 2
  ENOT23LOG <- 1 + Trial23dfLOG/(Trial2nLOG + Trial3nLOG - Trial23nLOG - 1)  # Effective number of trials for 2 vs 3
  ENOTcombLOG <- 1 + CombDFLOG/(meanTrialnLOG - 1)

  ##----##
  ##----##
  # F values
  ##----##
  ##----##

  F12LOG <- 1 + (Intraclassr12LOG * ENOT12LOG) / (1 - Intraclassr12LOG)
  F23LOG <- 1 + (Intraclassr23LOG * ENOT23LOG) / (1 - Intraclassr23LOG)
  CombFLOG <- 1 + (IntraclassrCombLOG * ENOTcombLOG) / (1 - IntraclassrCombLOG)

  ##----##
  ##----##
  ## F-test statistics
  ##----##
  ##----##

  # Numerator degrees of freedom for F-test
  numDf12LOG <- Trial12dfLOG
  numDf23LOG <- Trial2nLOG + Trial3nLOG - Trial23nLOG - 1
  numerDfcombLOG <- meanTrialnLOG - 1

  # Denominator degreees of freedom for F-test
  denomDf12LOG <- Trial12dfLOG
  denomDf23LOG <- Trial23dfLOG
  denomDfcombLOG <- CombDFLOG

  # F values (the combined CL are also not the same as the excel sheet due to different values in the FINV vs. qf when not using positive integers)
  fLower12LOG <- F12LOG / stats::qf(p = (1- CL/100)/2, df1 = numDf12LOG, df2 = denomDf12LOG, lower.tail = FALSE)
  fLower23LOG <- F23LOG / stats::qf(p = (1- CL/100)/2, df1 = numDf23LOG, df2 = denomDf23LOG, lower.tail = FALSE)
  fLowercombLOG  <- CombFLOG / stats::qf(p = (1- CL/100)/2, df1 = numerDfcombLOG, df2 = denomDfcombLOG, lower.tail = FALSE)

  fUpper12LOG <- F12LOG * stats::qf(p = (1- CL/100)/2, df1 = denomDf12LOG, df2 = numDf12LOG, lower.tail = FALSE)
  fUpper23LOG <- F23LOG * stats::qf(p = (1- CL/100)/2, df1 = denomDf23LOG, df2 = numDf23LOG, lower.tail = FALSE)
  fUppercombLOG <- CombFLOG * stats::qf(p = (1- CL/100)/2, df1 = denomDfcombLOG, df2 = numerDfcombLOG, lower.tail = FALSE)

  ##----##
  ##----##
  ## Intraclass r #2
  ##----##
  ##----##

  Intraclassr12lCLLOG <- (fLower12LOG - 1) / (fLower12LOG + ENOT12LOG - 1)
  Intraclassr12uCLLOG <- (fUpper12LOG - 1) / (fUpper12LOG + ENOT12LOG - 1)

  Intraclassr23lCLLOG <- (fLower23LOG - 1) / (fLower23LOG + ENOT23LOG - 1)
  Intraclassr23uCLLOG <- (fUpper23LOG - 1) / (fUpper23LOG + ENOT23LOG - 1)

  IntraclassrComblCLLOG <- (fLowercombLOG - 1) / (fLowercombLOG + ENOTcombLOG - 1)
  IntraclassrCombuCLLOG <- (fUppercombLOG - 1) / (fUppercombLOG + ENOTcombLOG - 1)

  ##----##
  ##----##
  ## Pearson r
  ##----##
  ##----##

  Pearson12LOG <- stats::cor(df$Trial1log, df$Trial2log, use = "complete.obs") * (1+ (1- (stats::cor(df$Trial1log, df$Trial2log, use = "complete.obs") ^ 2)) / 2 / (Trial12nLOG - 3))

  Pearson23LOG <- stats::cor(df$Trial2log, df$Trial3log, use = "complete.obs") * (1+ (1- (stats::cor(df$Trial2log, df$Trial3log, use = "complete.obs") ^ 2)) / 2 / (Trial23nLOG - 3))

  Pearson12lCLLOG <- (exp((2* ((0.5 * log((1 + Pearson12LOG) / (1 - Pearson12LOG))) - stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial12nLOG - 3))))) - 1) /
    (exp((2* ((0.5 * log((1 + Pearson12LOG) / (1 - Pearson12LOG))) - stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial12nLOG - 3))))) + 1)
  Pearson23lCLLOG <- (exp((2* ((0.5 * log((1 + Pearson23LOG) / (1 - Pearson23LOG))) - stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial23nLOG - 3))))) - 1) /
    (exp((2* ((0.5 * log((1 + Pearson23LOG) / (1 - Pearson23LOG))) - stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial23nLOG - 3))))) + 1)

  Pearson12uCLLOG <- (exp((2* ((0.5 * log((1 + Pearson12LOG) / (1 - Pearson12LOG))) + stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial12nLOG - 3))))) - 1) /
    (exp((2* ((0.5 * log((1 + Pearson12LOG) / (1 - Pearson12LOG))) + stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial12nLOG - 3))))) + 1)
  Pearson23uCLLOG <- (exp((2* ((0.5 * log((1 + Pearson23LOG) / (1 - Pearson23LOG))) + stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial23nLOG - 3))))) - 1) /
    (exp((2* ((0.5 * log((1 + Pearson23LOG) / (1 - Pearson23LOG))) + stats::qnorm((1 - ((1- CL/100)/2))) / sqrt((Trial23nLOG - 3))))) + 1)

  Pearson12BiasCFLOG <- 1 + ((1 - (stats::cor(df$Trial1log, df$Trial2log, use = "complete.obs") ^ 2)) / 2) / (Trial12nLOG - 3)
  Pearson23BiasCFLOG <- 1 + ((1 - (stats::cor(df$Trial2log, df$Trial3log, use = "complete.obs") ^ 2)) / 2) / (Trial23nLOG - 3)

  ##----##
  ##----##
  ## Fisher r
  ##----##
  ##----##

  FisherR12LOG <- psych::fisherz(Pearson12LOG)
  FisherR23LOG <- psych::fisherz(Pearson23LOG)
  FisherSE12LOG <- (Trial12nLOG - 3)
  FisherSE23LOG <- (Trial23nLOG - 3)

  FisherRcombLOG <- (FisherR12LOG * FisherSE12LOG + FisherR23LOG * FisherSE23LOG) / (FisherSE12LOG + FisherSE23LOG)

  PearsonCombLOG <- psych::fisherz2r(FisherRcombLOG)

  ##----##
  ##----##
  ## Pure Variance
  ##----##
  ##----##

  PureVariance12LOG <- (matrix12calcLOG / (Trial1dfLOG + Trial2dfLOG)) - (TE12LOG ^ 2)
  PureVariance23LOG <- (matrix23calcLOG / (Trial2dfLOG + Trial3dfLOG)) - (TE23LOG ^ 2)
  PureVarianceCombLOG <- (GrandSDLOG ^ 2) - (TEcombLOG ^ 2)

  SEPureVariance12LOG <- sqrt(2 * ((((matrix12calcLOG / (Trial1dfLOG + Trial2dfLOG)) ^ 2) / numDf12LOG) - (TE12LOG ^ 4 / denomDf12LOG)))
  SEPureVariance23LOG <- sqrt(2 * ((((matrix23calcLOG / (Trial2dfLOG + Trial3dfLOG)) ^ 2) / numDf23LOG) - (TE23LOG ^ 4 / denomDf23LOG)))
  SEPureVarianceCombLOG <- sqrt(((2 * ((GrandSDLOG ^ 4) / numerDfcombLOG)) - ((TEcombLOG ^ 4) / denomDfcombLOG)))

  SEPV12lCLLOG <- PureVariance12LOG + (stats::qnorm((1- CL/100)/2) * SEPureVariance12LOG)
  SEPV23lCLLOG <- PureVariance23LOG + (stats::qnorm((1- CL/100)/2) * SEPureVariance23LOG)
  SEPVcomblCLLOG <- PureVarianceCombLOG + (stats::qnorm((1- CL/100)/2) * SEPureVarianceCombLOG)

  SEPV12uCLLOG <- PureVariance12LOG - (stats::qnorm((1- CL/100)/2) * SEPureVariance12LOG)
  SEPV23uCLLOG <- PureVariance23LOG - (stats::qnorm((1- CL/100)/2) * SEPureVariance23LOG)
  SEPVcombuCLLOG <- PureVarianceCombLOG - (stats::qnorm((1- CL/100)/2) * SEPureVarianceCombLOG)

  ##----##
  ##----##

  ## Smallest effect from pure SD (LOG VALUES)

  # The minus value within the sqrt bracket is added because it is a negative number

  # NOTE: Based on a cohen's d of 0.2, 0.5, and 0.8 - entered below

  # The code for providing the confidence intervals (for one of the effect sizes) is provided for clarity, but not run in this instance

  ##----##
  ##----##

  #SmallestEffpureSD12LOG <- eff_interest_sml * (sqrt(PureVariance12LOG))
  #SmallestEffpureSD12lCLLOG <- eff_interest * ifelse(SEPV12lCLLOG > 0, (sqrt(SEPV12lCLLOG)), -(sqrt(-SEPV12lCLLOG)))
  #SmallestEffpureSD12uCLLOG <- eff_interest * (sqrt(SEPV12uCLLOG))

  #SmallestEffpureSD23LOG <- eff_interest_sml * (sqrt(PureVariance23LOG))
  #SmallestEffpureSD23lCLLOG <- eff_interest * ifelse(SEPV23lCLLOG > 0, (sqrt(SEPV23lCLLOG)), -(sqrt(-SEPV23lCLLOG)))
  #SmallestEffpureSD23uCLLOG <- eff_interest * (sqrt(SEPV23uCLLOG))

  #SmallestEffpureSDcombLOG <- eff_interest_sml * (sqrt(PureVarianceCombLOG))
  #SmallestEffpureSDcomblCLLOG <- eff_interest * ifelse(SEPVcomblCLLOG > 0, (sqrt(SEPVcomblCLLOG)), -(sqrt(-SEPVcomblCLLOG)))
  #SmallestEffpureSDcombuCLLOG <- eff_interest * (sqrt(SEPVcombuCLLOG))

  #MedEffectpureSD12LOG <- eff_interest_med * (sqrt(PureVariance12LOG))
  #MedEffectpureSD23LOG <- eff_interest_med * (sqrt(PureVariance23LOG))
  #MedEffectpureSDcombLOG <- eff_interest_med * (sqrt(PureVarianceCombLOG))

  #LrgEffectpureSD12LOG <- eff_interest_lrg * (sqrt(PureVariance12LOG))
  #LrgEffectpureSD23LOG <- eff_interest_lrg * (sqrt(PureVariance23LOG))
  #LrgEffectpureSDcombLOG <- eff_interest_lrg * (sqrt(PureVarianceCombLOG))

  ##----##
  ##----##

  ## Smallest effect from observed SD (LOG VALUES)

  # NOTE: Based on a cohen's d of 0.2, 0.5, and 0.8 - entered below

  # The code for providing the confidence intervals (for one of the effect sizes) is provided for clarity, but not run in this instance

  ##----##
  ##----##

  #SmallestEffobsSD12LOG <- eff_interest_sml * (sqrt(matrix12calcLOG/ (Trial1dfLOG + Trial2dfLOG)))
  #SmallestEffobsSD12lCLLOG <- sqrt((numDf12LOG * (SmallestEffobsSD12LOG ^ 2)) / qchisq(p = ((1 - CL/100)/2), df = numDf12LOG, lower.tail = FALSE))
  #SmallestEffobsSD12uCLLOG <- sqrt((numDf12LOG * (SmallestEffobsSD12LOG ^ 2)) / qchisq(p = (1 - (1 - CL/100)/2), df = numDf12LOG, lower.tail = FALSE))

  #SmallestEffobsSD23LOG <- eff_interest_sml * (sqrt(matrix23calcLOG/ (Trial2dfLOG + Trial3dfLOG)))
  #SmallestEffobsSD23lCLLOG <- sqrt((numDf23LOG * (SmallestEffobsSD23LOG ^ 2)) / qchisq(p = ((1 - CL/100)/2), df = numDf23LOG, lower.tail = FALSE))
  #SmallestEffobsSD23uCLLOG <- sqrt((numDf23LOG * (SmallestEffobsSD23LOG ^ 2)) / qchisq(p = (1 - (1 - CL/100)/2), df = numDf23LOG, lower.tail = FALSE))

  #SmallestEffobsSDcombLOG <- eff_interest_sml * GrandSDLOG
  #SmallestEffobsSDcomblCLLOG <- sqrt((numerDfcombLOG * (SmallestEffobsSDcombLOG ^ 2)) / qchisq(p = ((1 - CL/100)/2), df = numerDfcombLOG, lower.tail = FALSE))
  #SmallestEffobsSDcombuCLLOG <- sqrt((numerDfcombLOG * (SmallestEffobsSDcombLOG ^ 2)) / qchisq(p = (1 - (1 - CL/100)/2), df = numerDfcombLOG, lower.tail = FALSE))

  #MedEffectobsSD12LOG <- eff_interest_med * (sqrt(matrix12calcLOG/ (Trial1dfLOG + Trial2dfLOG)))
  #MedEffectobsSD23LOG <- eff_interest_med * (sqrt(matrix23calcLOG/ (Trial2dfLOG + Trial3dfLOG)))
  #MedEffectobsSDcombLOG <- eff_interest_med * GrandSDLOG

  #LrgEffectobsSD12LOG <- eff_interest_lrg * (sqrt(matrix12calcLOG/ (Trial1dfLOG + Trial2dfLOG)))
  #LrgEffectobsSD23LOG <- eff_interest_lrg * (sqrt(matrix23calcLOG/ (Trial2dfLOG + Trial3dfLOG)))
  #LrgEffectobsSDcombLOG <- eff_interest_lrg * GrandSDLOG

  ##----##
  ##----##

  # Standardized change in means

  ##----##
  ##----##

  StdizedChangeMean12LOG <- Df1means[[10]] / sqrt((matrix12calcLOG / (Trial1dfLOG + Trial2dfLOG)) - (TE12LOG ^ 2))
  StdizedChangeMean23LOG <- Df1means[[11]] / sqrt((matrix23calcLOG / (Trial2dfLOG + Trial3dfLOG)) - (TE23LOG ^ 2))

  StdizedChangeMean12lCLLOG <- ChMean12LOGlCL / sqrt((matrix12calcLOG / (Trial1dfLOG + Trial2dfLOG)) - (TE12LOG ^ 2))
  StdizedChangeMean23lCLLOG <- ChMean23LOGlCL / sqrt((matrix23calcLOG / (Trial2dfLOG + Trial3dfLOG)) - (TE23LOG ^ 2))

  StdizedChangeMean12uCLLOG <- ChMean12LOGuCL / sqrt((matrix12calcLOG / (Trial1dfLOG + Trial2dfLOG)) - (TE12LOG ^ 2))
  StdizedChangeMean23uCLLOG <- ChMean23LOGuCL / sqrt((matrix23calcLOG / (Trial2dfLOG + Trial3dfLOG)) - (TE23LOG ^ 2))

  ##----##
  ##----##

  # Standardized typical errors

  ##----##
  ##----##

  StdizedTE12LOG <- TE12LOG / sqrt((matrix12calcLOG / (Trial1dfLOG + Trial2dfLOG)) - (TE12LOG ^ 2))
  StdizedTE23LOG <- TE23LOG / sqrt((matrix23calcLOG / (Trial2dfLOG + Trial3dfLOG)) - (TE23LOG ^ 2))

  StdizedTE12lCLLOG <- TE12lCLLOG / sqrt((matrix12calcLOG / (Trial1dfLOG + Trial2dfLOG)) - (TE12LOG ^ 2))
  StdizedTE23lCLLOG <- TE23lCLLOG / sqrt((matrix23calcLOG / (Trial2dfLOG + Trial3dfLOG)) - (TE23LOG ^ 2))

  StdizedTE12uCLLOG <- TE12uCLLOG / sqrt((matrix12calcLOG / (Trial1dfLOG + Trial2dfLOG)) - (TE12LOG ^ 2))
  StdizedTE23uCLLOG <- TE23uCLLOG / sqrt((matrix23calcLOG / (Trial2dfLOG + Trial3dfLOG)) - (TE23LOG ^ 2))

  StdizedTEcombLOG <- TEcombLOG / sqrt(GrandSDLOG ^ 2)
  StdizedTEcomblCLLOG <- TEcomblCLLOG / sqrt(GrandSDLOG ^ 2)
  StdizedTEcombuCLLOG <- TEcombuCLLOG / sqrt(GrandSDLOG ^ 2)

  ##----##
  ## FACTORS:  ##

  # NOTE: Use when changes and mean and errors exceed 50%
  ##----##

  # Factor change in mean
  FactChMean12 <- exp((Df1means[[10]]/100))
  FactChMean12lCL <- exp((ChMean12LOGlCL/100))
  FactChMean12uCL <- exp((ChMean12LOGuCL/100))
  FactChMean12CLasMULTIPLYDIVIDED <- sqrt(FactChMean12uCL/FactChMean12lCL)

  FactChMean23 <- exp((Df1means[[11]]/100))
  FactChMean23lCL <- exp((ChMean23LOGlCL/100))
  FactChMean23uCL <- exp((ChMean23LOGuCL/100))
  FactChMean23CLasMULTIPLYDIVIDED <- sqrt(FactChMean23uCL/FactChMean23lCL)

  # Factor typical error
  FactTE12 <- exp(TE12LOG/100)
  FactTE12lCL <- exp(TE12lCLLOG/100)
  FactTE12uCL <- exp(TE12uCLLOG/100)
  FactTE12CLasMULTIPLYDIVIDED <- sqrt(FactTE12uCL/FactTE12lCL)
  FactTE12BiasCF <- (FactTE12 ^ TE12BiasCFLOG) / FactTE12

  FactTE23 <- exp(TE23LOG/100)
  FactTE23lCL <- exp(TE23lCLLOG/100)
  FactTE23uCL <- exp(TE23uCLLOG/100)
  FactTE23CLasMULTIPLYDIVIDED <- sqrt(FactTE23uCL/FactTE23lCL)
  FactTE23BiasCF <- (FactTE23 ^ TE23BiasCFLOG) / FactTE23

  FactTEcomb <- exp(TEcombLOG/100)
  FactTEcomblCL <- exp(TEcomblCLLOG/100)
  FactTEcombuCL <- exp(TEcombuCLLOG/100)
  FactTEcombCLasMULTIPLYDIVIDED <- sqrt(FactTEcombuCL/FactTEcomblCL)
  FactTEcombBiasCF <- (FactTEcomb ^ TEcombBiasCFLOG) / FactTEcomb

  ## Smallest effect Pure

  #FactSmallestEffpure12 <- exp(SmallestEffpureSD12LOG /100)
  #FactSmallestEffectPure12lCL <- exp(SmallestEffpureSD12lCLLOG /100)
  #FactSmallestEffectPure12uCL <- exp(SmallestEffpureSD12uCLLOG /100)

  #FactSmallestEffpure23 <- exp(SmallestEffpureSD23LOG /100)
  #FactSmallestEffectPure23lCL <- exp(SmallestEffpureSD23lCLLOG /100)
  #FactSmallestEffectPure23uCL <- exp(SmallestEffpureSD23uCLLOG /100)

  #FactSmallestEffpureComb <- exp(SmallestEffpureSDcombLOG /100)
  #FactSmallestEffectPureComblCL <- exp(SmallestEffpureSDcomblCLLOG /100)
  #FactSmallestEffectPureCombuCL <- exp(SmallestEffpureSDcombuCLLOG /100)

  #FactMedEffectpure12 <- exp(MedEffectpureSD12LOG /100)
  #FactMedEffectpure23 <- exp(MedEffectpureSD23LOG /100)
  #FactMedEffectpureComb <- exp(MedEffectpureSDcombLOG /100)

  #FactLrgEffectpure12 <- exp(LrgEffectpureSD12LOG /100)
  #FactLrgEffectpure23 <- exp(LrgEffectpureSD23LOG /100)
  #FactLrgEffectpureComb <- exp(LrgEffectpureSDcombLOG /100)

  ## Smallest effect Observed

  #FactSmallestEffobs12 <- exp(SmallestEffobsSD12LOG /100)
  #FactSmallestEffectObs12lCL <- exp(SmallestEffobsSD12lCLLOG /100)
  #FactSmallestEffectObs12uCL <- exp(SmallestEffobsSD12uCLLOG /100)

  #FactSmallestEffobs23 <- exp(SmallestEffobsSD23LOG /100)
  #FactSmallestEffectObs23lCL <- exp(SmallestEffobsSD23lCLLOG /100)
  #FactSmallestEffectObs23uCL <- exp(SmallestEffobsSD23uCLLOG /100)

  #FactSmallestEffobsComb <- exp(SmallestEffobsSDcombLOG /100)
  #FactSmallestEffectObsComblCL <- exp(SmallestEffobsSDcomblCLLOG /100)
  #FactSmallestEffectObsCombuCL <- exp(SmallestEffobsSDcombuCLLOG /100)

  #FactMedEffectobs12 <- exp(MedEffectobsSD12LOG /100)
  #FactMedEffectobs23 <- exp(MedEffectobsSD23LOG /100)
  #FactMedEffectobsComb <- exp(MedEffectobsSDcombLOG /100)

  #FactLrgEffectobs12 <- exp(LrgEffectobsSD12LOG /100)
  #FactLrgEffectobs23 <- exp(LrgEffectobsSD23LOG /100)
  #FactLrgEffectobsComb <- exp(LrgEffectobsSDcombLOG /100)

  ##----##
  ## PERCENTAGES:  ##

  # NOTE: Use when changes and mean and errors are less than 50%

  # 100 * log(value) is used for ease of interpretation as the data are aleady effectively a % when expressed as changes in the typical error.
  ##----##

  ## Percentage Change in mean
  PercentChMean12 <- 100 * exp((Df1means[[10]]/100)) - 100
  PercentChMean12lCL <- 100 * exp((ChMean12LOGlCL/100)) - 100
  PercentChMean12uCL <- 100 * exp((ChMean12LOGuCL/100)) - 100
  PercentChMean12CLplusminus <- (PercentChMean12uCL - PercentChMean12lCL) / 2

  PercentChMean23 <- 100 * exp((Df1means[[11]]/100)) - 100
  PercentChMean23lCL <- 100 * exp((ChMean23LOGlCL/100)) - 100
  PercentChMean23uCL <- 100 * exp((ChMean23LOGuCL/100)) - 100
  PercentChMean23CLplusminus <- (PercentChMean23uCL - PercentChMean23lCL) / 2

  PercentTE12 <- 100 * exp((TE12LOG / 100)) - 100
  PercentTE12lCL <- 100 * exp((TE12lCLLOG / 100)) - 100
  PercentTE12uCL <- 100 * exp((TE12uCLLOG / 100)) - 100
  PercentTE12CLasMULTIPLYDIVIDED <- sqrt(PercentTE12uCL / PercentTE12lCL)
  PercentTE12BiasCF <- (100 * (FactTE12 ^ TE12BiasCF) - 100) / ((100 * FactTE12) - 100)

  PercentTE23 <- 100 * exp((TE23LOG / 100)) - 100
  PercentTE23lCL <- 100 * exp((TE23lCLLOG / 100)) - 100
  PercentTE23uCL <- 100 * exp((TE23uCLLOG / 100)) - 100
  PercentTE23CLasMULTIPLYDIVIDED <- sqrt(PercentTE23uCL / PercentTE23lCL)
  PercentTE23BiasCF <- (100 * (FactTE23 ^ TE23BiasCF) - 100) / ((100 * FactTE23) - 100)

  PercentTEcomb <- 100 * exp((TEcombLOG / 100)) - 100
  PercentTEcomblCL <- 100 * exp((TEcomblCLLOG / 100)) - 100
  PercentTEcombuCL <- 100 * exp((TEcombuCLLOG / 100)) - 100
  PercentTEcombCLasMULTIPLYDIVIDED <- sqrt(PercentTEcombuCL / PercentTEcomblCL)
  PercentTEcombBiasCF <- (100 * (FactTEcomb ^ TEcombBiasCF) - 100) / ((100 * FactTEcomb) - 100)

  ## Smallest effect Pure

  #PercentSmallestEffpure12 <- (100 * FactSmallestEffpure12) - 100
  #PercentSmallestEffectPure12lCL <- (100 * FactSmallestEffectPure12lCL) - 100
  #PercentSmallestEffectPure12uCL <- (100 * FactSmallestEffectPure12uCL) - 100

  #PercentSmallestEffpure23 <- (100 * FactSmallestEffpure23) - 100
  #PercentSmallestEffectPure23lCL <- (100 * FactSmallestEffectPure23lCL) - 100
  #PercentSmallestEffectPure23uCL <- (100 * FactSmallestEffectPure23uCL) - 100

  #PercentSmallestEffpureComb <- (100 * FactSmallestEffpureComb) - 100
  #PercentSmallestEffectPureComblCL <- (100 * FactSmallestEffectPureComblCL) - 100
  #PercentSmallestEffectPureCombuCL <- (100 * FactSmallestEffectPureCombuCL) - 100

  #PercentMedEffectpure12 <- (100 * FactMedEffectpure12) - 100
  #PercentMedEffectpure23 <- (100 * FactMedEffectpure23) - 100
  #PercentMedEffectpureComb <- (100 * FactMedEffectpureComb) - 100

  #PercentLrgEffectpure12 <- (100 * FactLrgEffectpure12) - 100
  #PercentLrgEffectpure23 <- (100 * FactLrgEffectpure23) - 100
  #PercentLrgEffectpureComb <- (100 * FactLrgEffectpureComb) - 100

  ## Smallest effect Observed

  #PercentSmallestEffobs12 <- (100 * FactSmallestEffobs12) - 100
  #PercentSmallestEffectObs12lCL <- (100 * FactSmallestEffectObs12lCL) - 100
  #PercentSmallestEffectObs12uCL <- (100 * FactSmallestEffectObs12uCL) - 100

  #PercentSmallestEffobs23 <- (100 * FactSmallestEffobs23) - 100
  #PercentSmallestEffectObs23lCL <- (100 * FactSmallestEffectObs23lCL) - 100
  #PercentSmallestEffectObs23uCL <- (100 * FactSmallestEffectObs23uCL) - 100

  #PercentSmallestEffobsComb <- (100 * FactSmallestEffobsComb) - 100
  #PercentSmallestEffectObsComblCL <- (100 * FactSmallestEffectObsComblCL) - 100
  #PercentSmallestEffectObsCombuCL <- (100 * FactSmallestEffectObsCombuCL) - 100

  #PercentMedEffectobs12 <- (100 * FactMedEffectobs12) - 100
  #PercentMedEffectobs23 <- (100 * FactMedEffectobs23) - 100
  #PercentMedEffectobsComb <- (100 * FactMedEffectobsComb) - 100

  #PercentLrgEffectobs12 <- (100 * FactLrgEffectobs12) - 100
  #PercentLrgEffectobs23 <- (100 * FactLrgEffectobs23) - 100
  #PercentLrgEffectobsComb <- (100 * FactLrgEffectobsComb) - 100

  ##----##
  ## STANDARDIZED:  ##
  ##----##

  ## Standardized change in means (see code above for calculations)
  StdizedChangeMean12CLplusminus <- (StdizedChangeMean12uCLLOG - StdizedChangeMean12lCLLOG) / 2
  StdizedChangeMean23CLplusminus <- (StdizedChangeMean23uCLLOG - StdizedChangeMean23lCLLOG) / 2

  ## Standardized typical errors
  StdizedTE12CLplusminus <- sqrt((StdizedTE12uCLLOG / StdizedTE12lCLLOG))
  StdizedTE23CLplusminus <- sqrt((StdizedTE23uCLLOG / StdizedTE23lCLLOG))
  StdizedTEcombCLplusminus <- sqrt((StdizedTEcombuCLLOG / StdizedTEcomblCLLOG))

  ##----##
  ## CORRELATIONS:  ##
  ##----##

  ICC12LOG <- Intraclassr12LOG
  ICC12lCLLOG <- Intraclassr12lCLLOG
  ICC12uCLLOG <- Intraclassr12uCLLOG

  ICC23LOG <- Intraclassr23LOG
  ICC23lCLLOG <- Intraclassr23lCLLOG
  ICC23uCLLOG <- Intraclassr23uCLLOG

  ICCcombLOG <- IntraclassrCombLOG
  ICCcomblCLLOG <- IntraclassrComblCLLOG
  ICCcombuCLLOG <- IntraclassrCombuCLLOG

  ##---------------------------------##
  ## MINIMUM DETECTABLE CHANGE  ##

  # Note: From Beaton, 2000, Spine, and see Bland and Altman 1996, BMJ.
  # 1.96 is the 2-sided z value for the 95% CI; sqrt(2 or 3) accounts for the variance of 2 or 3 measurements
  # As a %, this equates to the other derivations
  ##---------------------------------##

  ## 70% confidence level

  MDC70_12 <- TE12 * (stats::qnorm(0.70 + (1 - 0.70) / 2)) * (sqrt(2))
  MDC70_23 <- TE23 * (stats::qnorm(0.70 + (1 - 0.70) / 2)) * (sqrt(2))
  MDC70_comb <- TEcomb * (stats::qnorm(0.70 + (1 - 0.70) / 2)) * (sqrt(3))

  MDC70_rel_12 <- (MDC70_12/GrandMean) * 100
  MDC70_rel_23 <- (MDC70_23/GrandMean) * 100
  MDC70_rel_comb <- (MDC70_comb/GrandMean) * 100

  MDC70_12LOG <- PercentTE12 * (stats::qnorm(0.70 + (1 - 0.70) / 2)) * (sqrt(2))
  MDC70_23LOG <- PercentTE23 * (stats::qnorm(0.70 + (1 - 0.70) / 2)) * (sqrt(2))
  MDC70_combLOG <- PercentTEcomb * (stats::qnorm(0.70 + (1 - 0.70) / 2)) * (sqrt(3))

  MDC70_rel_12LOG <- (MDC70_12LOG/GrandMeanLOG) * 100
  MDC70_rel_23LOG <- (MDC70_23LOG/GrandMeanLOG) * 100
  MDC70_rel_combLOG <- (MDC70_combLOG/GrandMeanLOG) * 100

  # Group level calculations - divided by sqrt of sample size. n as 4 for conservative estimate (i.e., any group more than 4)

  MDC70_12_gp <- (TE12 / sqrt(4)) * (stats::qnorm(0.70 + (1 - 0.70) / 2)) * (sqrt(2))
  MDC70_23_gp <- (TE23 / sqrt(4)) * (stats::qnorm(0.70 + (1 - 0.70) / 2)) * (sqrt(2))
  MDC70_comb_gp <- (TEcomb / sqrt(4)) * (stats::qnorm(0.70 + (1 - 0.70) / 2)) * (sqrt(3))

  MDC70_12LOG_gp <- (PercentTE12 / sqrt(4)) * (stats::qnorm(0.70 + (1 - 0.70) / 2)) * (sqrt(2))
  MDC70_23LOG_gp <- (PercentTE23 / sqrt(4)) * (stats::qnorm(0.70 + (1 - 0.70) / 2)) * (sqrt(2))
  MDC70_combLOG_gp <- (PercentTEcomb / sqrt(4)) * (stats::qnorm(0.70 + (1 - 0.70) / 2)) * (sqrt(3))

  ## 80% confidence level

  MDC80_12 <- TE12 * (stats::qnorm(0.80 + (1 - 0.80) / 2)) * (sqrt(2))
  MDC80_23 <- TE23 * (stats::qnorm(0.80 + (1 - 0.80) / 2)) * (sqrt(2))
  MDC80_comb <- TEcomb * (stats::qnorm(0.80 + (1 - 0.80) / 2)) * (sqrt(3))

  MDC80_rel_12 <- (MDC80_12/GrandMean) * 100
  MDC80_rel_23 <- (MDC80_23/GrandMean) * 100
  MDC80_rel_comb <- (MDC80_comb/GrandMean) * 100

  MDC80_12LOG <- PercentTE12 * (stats::qnorm(0.80 + (1 - 0.80) / 2)) * (sqrt(2))
  MDC80_23LOG <- PercentTE23 * (stats::qnorm(0.80 + (1 - 0.80) / 2)) * (sqrt(2))
  MDC80_combLOG <- PercentTEcomb * (stats::qnorm(0.80 + (1 - 0.80) / 2)) * (sqrt(3))

  MDC80_rel_12LOG <- (MDC80_12LOG/GrandMeanLOG) * 100
  MDC80_rel_23LOG <- (MDC80_23LOG/GrandMeanLOG) * 100
  MDC80_rel_combLOG <- (MDC80_combLOG/GrandMeanLOG) * 100

  # Group level calculations - divided by sqrt of sample size. n as 4 for conservative estimate (i.e., any group more than 4)

  MDC80_12_gp <- (TE12 / sqrt(4)) * (stats::qnorm(0.80 + (1 - 0.80) / 2)) * (sqrt(2))
  MDC80_23_gp <- (TE23 / sqrt(4)) * (stats::qnorm(0.80 + (1 - 0.80) / 2)) * (sqrt(2))
  MDC80_comb_gp <- (TEcomb / sqrt(4)) * (stats::qnorm(0.80 + (1 - 0.80) / 2)) * (sqrt(3))

  MDC80_12LOG_gp <- (PercentTE12 / sqrt(4)) * (stats::qnorm(0.80 + (1 - 0.80) / 2)) * (sqrt(2))
  MDC80_23LOG_gp <- (PercentTE23 / sqrt(4)) * (stats::qnorm(0.80 + (1 - 0.80) / 2)) * (sqrt(2))
  MDC80_combLOG_gp <- (PercentTEcomb / sqrt(4)) * (stats::qnorm(0.80 + (1 - 0.80) / 2)) * (sqrt(3))

  ## 90% confidence level

  MDC90_12 <- TE12 * (stats::qnorm(0.90 + (1 - 0.90) / 2)) * (sqrt(2))
  MDC90_23 <- TE23 * (stats::qnorm(0.90 + (1 - 0.90) / 2)) * (sqrt(2))
  MDC90_comb <- TEcomb * (stats::qnorm(0.90 + (1 - 0.90) / 2)) * (sqrt(3))

  MDC90_rel_12 <- (MDC90_12/GrandMean) * 100
  MDC90_rel_23 <- (MDC90_23/GrandMean) * 100
  MDC90_rel_comb <- (MDC90_comb/GrandMean) * 100

  MDC90_12LOG <- PercentTE12 * (stats::qnorm(0.90 + (1 - 0.90) / 2)) * (sqrt(2))
  MDC90_23LOG <- PercentTE23 * (stats::qnorm(0.90 + (1 - 0.90) / 2)) * (sqrt(2))
  MDC90_combLOG <- PercentTEcomb * (stats::qnorm(0.90 + (1 - 0.90) / 2)) * (sqrt(3))

  MDC90_rel_12LOG <- (MDC90_12LOG/GrandMeanLOG) * 100
  MDC90_rel_23LOG <- (MDC90_23LOG/GrandMeanLOG) * 100
  MDC90_rel_combLOG <- (MDC90_combLOG/GrandMeanLOG) * 100

  # Group level calculations - divided by sqrt of sample size. n as 4 for conservative estimate (i.e., any group more than 4)

  MDC90_12_gp <- (TE12 / sqrt(4)) * (stats::qnorm(0.90 + (1 - 0.90) / 2)) * (sqrt(2))
  MDC90_23_gp <- (TE23 / sqrt(4)) * (stats::qnorm(0.90 + (1 - 0.90) / 2)) * (sqrt(2))
  MDC90_comb_gp <- (TEcomb / sqrt(4)) * (stats::qnorm(0.90 + (1 - 0.90) / 2)) * (sqrt(3))

  MDC90_12LOG_gp <- (PercentTE12 / sqrt(4)) * (stats::qnorm(0.90 + (1 - 0.90) / 2)) * (sqrt(2))
  MDC90_23LOG_gp <- (PercentTE23 / sqrt(4)) * (stats::qnorm(0.90 + (1 - 0.90) / 2)) * (sqrt(2))
  MDC90_combLOG_gp <- (PercentTEcomb / sqrt(4)) * (stats::qnorm(0.90 + (1 - 0.90) / 2)) * (sqrt(3))

  ## 95% confidence level

  MDC95_12 <- TE12 * (stats::qnorm(0.95 + (1 - 0.95) / 2)) * (sqrt(2))
  MDC95_23 <- TE23 * (stats::qnorm(0.95 + (1 - 0.95) / 2)) * (sqrt(2))
  MDC95_comb <- TEcomb * (stats::qnorm(0.95 + (1 - 0.95) / 2)) * (sqrt(3))

  MDC95_rel_12 <- (MDC95_12/GrandMean) * 100
  MDC95_rel_23 <- (MDC95_23/GrandMean) * 100
  MDC95_rel_comb <- (MDC95_comb/GrandMean) * 100

  MDC95_12LOG <- PercentTE12 * (stats::qnorm(0.95 + (1 - 0.95) / 2)) * (sqrt(2))
  MDC95_23LOG <- PercentTE23 * (stats::qnorm(0.95 + (1 - 0.95) / 2)) * (sqrt(2))
  MDC95_combLOG <- PercentTEcomb * (stats::qnorm(0.95 + (1 - 0.95) / 2)) * (sqrt(3))

  MDC95_rel_12LOG <- (MDC95_12LOG/GrandMeanLOG) * 100
  MDC95_rel_23LOG <- (MDC95_23LOG/GrandMeanLOG) * 100
  MDC95_rel_combLOG <- (MDC95_combLOG/GrandMeanLOG) * 100

  # Group level calculations - divided by sqrt of sample size. n as 4 for conservative estimate (i.e., any group more than 4)

  MDC95_12_gp <- (TE12 / sqrt(4)) * (stats::qnorm(0.95 + (1 - 0.95) / 2)) * (sqrt(2))
  MDC95_23_gp <- (TE23 / sqrt(4)) * (stats::qnorm(0.95 + (1 - 0.95) / 2)) * (sqrt(2))
  MDC95_comb_gp <- (TEcomb / sqrt(4)) * (stats::qnorm(0.95 + (1 - 0.95) / 2)) * (sqrt(3))

  MDC95_12LOG_gp <- (PercentTE12 / sqrt(4)) * (stats::qnorm(0.95 + (1 - 0.95) / 2)) * (sqrt(2))
  MDC95_23LOG_gp <- (PercentTE23 / sqrt(4)) * (stats::qnorm(0.95 + (1 - 0.95) / 2)) * (sqrt(2))
  MDC95_combLOG_gp <- (PercentTEcomb / sqrt(4)) * (stats::qnorm(0.95 + (1 - 0.95) / 2)) * (sqrt(3))

  ## 99% confidence level

  MDC99_12 <- TE12 * (stats::qnorm(0.99 + (1 - 0.99) / 2)) * (sqrt(2))
  MDC99_23 <- TE23 * (stats::qnorm(0.99 + (1 - 0.99) / 2)) * (sqrt(2))
  MDC99_comb <- TEcomb * (stats::qnorm(0.99 + (1 - 0.99) / 2)) * (sqrt(3))

  MDC99_rel_12 <- (MDC99_12/GrandMean) * 100
  MDC99_rel_23 <- (MDC99_23/GrandMean) * 100
  MDC99_rel_comb <- (MDC99_comb/GrandMean) * 100

  MDC99_12LOG <- PercentTE12 * (stats::qnorm(0.99 + (1 - 0.99) / 2)) * (sqrt(2))
  MDC99_23LOG <- PercentTE23 * (stats::qnorm(0.99 + (1 - 0.99) / 2)) * (sqrt(2))
  MDC99_combLOG <- PercentTEcomb * (stats::qnorm(0.99 + (1 - 0.99) / 2)) * (sqrt(3))

  MDC99_rel_12LOG <- (MDC99_12LOG/GrandMeanLOG) * 100
  MDC99_rel_23LOG <- (MDC99_23LOG/GrandMeanLOG) * 100
  MDC99_rel_combLOG <- (MDC99_combLOG/GrandMeanLOG) * 100

  # Group level calculations - divided by sqrt of sample size. n as 4 for conservative estimate (i.e., any group more than 4)

  MDC99_12_gp <- (TE12 / sqrt(4)) * (stats::qnorm(0.99 + (1 - 0.99) / 2)) * (sqrt(2))
  MDC99_23_gp <- (TE23 / sqrt(4)) * (stats::qnorm(0.99 + (1 - 0.99) / 2)) * (sqrt(2))
  MDC99_comb_gp <- (TEcomb / sqrt(4)) * (stats::qnorm(0.99 + (1 - 0.99) / 2)) * (sqrt(3))

  MDC99_12LOG_gp <- (PercentTE12 / sqrt(4)) * (stats::qnorm(0.99 + (1 - 0.99) / 2)) * (sqrt(2))
  MDC99_23LOG_gp <- (PercentTE23 / sqrt(4)) * (stats::qnorm(0.99 + (1 - 0.99) / 2)) * (sqrt(2))
  MDC99_combLOG_gp <- (PercentTEcomb / sqrt(4)) * (stats::qnorm(0.99 + (1 - 0.99) / 2)) * (sqrt(3))

  ##---------------------------------##
  ##---------------------------------##
  ## TABLES DATA ##
  ##---------------------------------##
  ##---------------------------------##

  ## Raw Values Table

  ## NOTE: Have removed CI around MDC and effect sizes
  ## NOTE: Have also removed MDC and effect sizes for 12 and 23 and just retained the combined

  # Descriptive Stats
  Trial1 <- c(paste(round(Df1means[[2]], digits = 2), " (", round(Df1sds[[2]], digits = 2), ")"))
  Trial2 <- c(paste(round(Df1means[[3]], digits = 2), " (", round(Df1sds[[3]], digits = 2), ")"))
  Trial3 <- c(paste(round(Df1means[[4]], digits = 2), " (", round(Df1sds[[4]], digits = 2), ")"))
  Average <- c(paste(round(GrandMean, digits = 2),  " (", round(GrandSD, digits = 2), ")"))

  # All ICCs
  tableICC12 <- c(paste(round(ICC12, digits = 2), "[", round(ICC12lCL, digits = 2), ",", round(ICC12uCL, digits = 2), "]"))
  tableICC23 <- c(paste(round(ICC23, digits = 2), "[", round(ICC23lCL, digits = 2), ",", round(ICC23uCL, digits = 2), "]"))
  tableICCmean <- c(paste(round(ICCcomb, digits = 2), "[", round(ICCcomblCL, digits = 2), ",", round(ICCcombuCL, digits = 2), "]"))

  # All TE's
  tableTE12 <- c(paste(round(TE12, digits = 2), "[", round(TE12lCL, digits = 2), ",", round(TE12uCL, digits = 2), "]"))
  tableTE23 <- c(paste(round(TE23, digits = 2), "[", round(TE23lCL, digits = 2), ",", round(TE23uCL, digits = 2), "]"))
  tableTEmean <- c(paste(round(TEcomb, digits = 2), "[", round(TEcomblCL, digits = 2), ",", round(TEcombuCL, digits = 2), "]"))

  # All Pearson's
  tablePearson12 <- c(paste(round(Pearson12, digits = 2), "[", round(Pearson12lCL, digits = 2), ",", round(Pearson12uCL, digits = 2), "]"))
  tablePearson23 <- c(paste(round(Pearson23, digits = 2), "[", round(Pearson23lCL, digits = 2), ",", round(Pearson23uCL, digits = 2), "]"))
  tablePearsonmean <- c(paste(round(PearsonComb, digits = 2), "[", "NA", "]"))

  ## All effect sizes of interest from pure and observed SD

  #tableSmallestEffpure12 <- c(round(SmallestEffpureSD12, digits = 2))
  #tableSmallestEffpure23 <- c(round(SmallestEffpureSD23, digits = 2))
  #tableSmallestEffpureComb <- c(round(SmallestEffpureSDcomb, digits = 2))
  #tableSmallestEffobs12 <- c(round(SmallestEffobsSD12, digits = 2))
  #tableSmallestEffobs23 <- c(round(SmallestEffobsSD23, digits = 2))
  #tableSmallestEffobsComb <- c(round(SmallestEffobsSDcomb, digits = 2))

  #tableMedEffectpure12 <- c(round(MedEffectpureSD12, digits = 2))
  #tableMedEffectpure23 <- c(round(MedEffectpureSD23, digits = 2))
  #tableMedEffectpureComb <- c(round(MedEffectpureSDcomb, digits = 2))
  #tableMedEffectobs12 <- c(round(MedEffectobsSD12, digits = 2))
  #tableMedEffectobs23 <- c(round(MedEffectobsSD23, digits = 2))
  #tableMedEffectobsComb <- c(round(MedEffectobsSDcomb, digits = 2))

  #tableLrgEffectpure12 <- c(round(LrgEffectpureSD12, digits = 2))
  #tableLrgEffectpure23 <- c(round(LrgEffectpureSD23, digits = 2))
  #tableLrgEffectpureComb <- c(round(LrgEffectpureSDcomb, digits = 2))
  #tableLrgEffectobs12 <- c(round(LrgEffectobsSD12, digits = 2))
  #tableLrgEffectobs23 <- c(round(LrgEffectobsSD23, digits = 2))
  #tableLrgEffectobsComb <- c(round(LrgEffectobsSDcomb, digits = 2))

  ## Log Values Table

  # Descriptive Stats
  Trial1log <- c(paste(round(Df1means[[7]], digits = 1), " (", round(Df1sds[[7]], digits = 1), ")"))
  Trial2log <- c(paste(round(Df1means[[8]], digits = 1), " (", round(Df1sds[[8]], digits = 1), ")"))
  Trial3log <- c(paste(round(Df1means[[9]], digits = 1), " (", round(Df1sds[[9]], digits = 1), ")"))
  Averagelog <- c(paste(round(GrandMeanLOG, digits = 2), " (", round(GrandSDLOG, digits = 2), ")"))

  # All ICCs
  tableICC12log <- c(paste(round(ICC12LOG, digits = 2), "[", round(ICC12lCLLOG, digits = 2), ",", round(ICC12uCLLOG, digits = 2), "]"))
  tableICC23log <- c(paste(round(ICC23LOG, digits = 2), "[", round(ICC23lCLLOG, digits = 2), ",", round(ICC23uCLLOG, digits = 2), "]"))
  tableICCmeanlog <- c(paste(round(ICCcombLOG, digits = 2), "[", round(ICCcomblCLLOG, digits = 2), ",", round(ICCcombuCLLOG, digits = 2), "]"))

  # All TE's
  tableTE12log <- c(paste(round(PercentTE12, digits = 2), "[", round(PercentTE12lCL, digits = 2), ",", round(PercentTE12uCL, digits = 2), "]"))
  tableTE23log <- c(paste(round(PercentTE23, digits = 2), "[", round(PercentTE23lCL, digits = 2), ",", round(PercentTE23uCL, digits = 2), "]"))
  tableTEmeanlog <- c(paste(round(PercentTEcomb, digits = 2), "[", round(PercentTEcomblCL, digits = 2), ",", round(PercentTEcombuCL, digits = 2), "]"))

  # All Pearson's
  tablePearson12log <- c(paste(round(Pearson12LOG, digits = 2), "[", round(Pearson12lCLLOG, digits = 2), ",", round(Pearson12uCLLOG, digits = 2), "]"))
  tablePearson23log <- c(paste(round(Pearson23LOG, digits = 2), "[", round(Pearson23lCLLOG, digits = 2), ",", round(Pearson23uCLLOG, digits = 2), "]"))
  tablePearsonmeanlog <- c(paste(round(PearsonCombLOG, digits = 2), "[", "NA", "]"))

  ## All effect sizes of interest from pure and observed SD

  #tableSmallestEffpure12log <- c(round(PercentSmallestEffpure12, digits = 2))
  #tableSmallestEffpure23log <- c(round(PercentSmallestEffpure23, digits = 2))
  #tableSmallestEffpureComblog <- c(round(PercentSmallestEffpureComb, digits = 2))
  #tableSmallestEffobs12log <- c(round(PercentSmallestEffobs12, digits = 2))
  #tableSmallestEffobs23log <- c(round(PercentSmallestEffobs23, digits = 2))
  #tableSmallestEffobsComblog <- c(round(PercentSmallestEffobsComb, digits = 2))

  #tableMedEffectpure12log <- c(round(PercentMedEffectpure12, digits = 2))
  #tableMedEffectpure23log <- c(round(PercentMedEffectpure23, digits = 2))
  #tableMedEffectpureComblog <- c(round(PercentMedEffectpureComb, digits = 2))
  #tableMedEffectobs12log <- c(round(PercentMedEffectobs12, digits = 2))
  #tableMedEffectobs23log <- c(round(PercentMedEffectobs23, digits = 2))
  #tableMedEffectobsComblog <- c(round(PercentMedEffectobsComb, digits = 2))

  #tableLrgEffectpure12log <- c(round(PercentLrgEffectpure12, digits = 2))
  #tableLrgEffectpure23log <- c(round(PercentLrgEffectpure23, digits = 2))
  #tableLrgEffectpureComblog <- c(round(PercentLrgEffectpureComb, digits = 2))
  #tableLrgEffectobs12log <- c(round(PercentLrgEffectobs12, digits = 2))
  #tableLrgEffectobs23log <- c(round(PercentLrgEffectobs23, digits = 2))
  #tableLrgEffectobsComblog <- c(round(PercentLrgEffectobsComb, digits = 2))

  ##---------------------------------##
  ##---------------------------------##
  ## Final Values Table
  ##---------------------------------##
  ##---------------------------------##

  ValuesTable <- data.frame(Trial1 = c(Trial1, Trial1log), Trial2 = c(Trial2, Trial2log), Trial3 = c(Trial3, Trial3log), Average = c(Average, Averagelog),
                            Trial1n = c(Trial1n, Trial1nLOG), Trial2n = c(Trial2n, Trial2nLOG), Trial3n = c(Trial3n, Trial3nLOG),

                            ICC12 = c(tableICC12, tableICC12log), ICC23 = c(tableICC23, tableICC23log), ICCmean = c(tableICCmean, tableICCmeanlog),
                            TE12 = c(tableTE12, tableTE12log), TE23 = c(tableTE23, tableTE23log), TEmean = c(tableTEmean, tableTEmeanlog),

                            R12 = c(tablePearson12, tablePearson12log), R23 = c(tablePearson23, tablePearson23log), Rmean = c(tablePearsonmean, tablePearsonmeanlog),

                            #MDC70_12 = c(round(MDC70_12, digits = 2), round(MDC70_12LOG, digits = 2)),
                            #MDC70_23 = c(round(MDC70_23, digits = 2), round(MDC70_23LOG, digits = 2)),
                            MDC70_comb = c(round(MDC70_comb, digits = 2), round(MDC70_combLOG, digits = 2)),
                            #MDC70_rel_12 = c(round(MDC70_rel_12, digits = 2), round(MDC70_rel_12LOG, digits = 2)),
                            #MDC70_rel_23 = c(round(MDC70_rel_23, digits = 2), round(MDC70_rel_23LOG, digits = 2)),
                            #MDC70_rel_comb = c(round(MDC70_rel_comb, digits = 2), round(MDC70_rel_combLOG, digits = 2)),

                            #MDC80_12 = c(round(MDC80_12, digits = 2), round(MDC80_12LOG, digits = 2)),
                            #MDC80_23 = c(round(MDC80_23, digits = 2), round(MDC80_23LOG, digits = 2)),
                            MDC80_comb = c(round(MDC80_comb, digits = 2), round(MDC80_combLOG, digits = 2)),
                            #MDC80_rel_12 = c(round(MDC80_rel_12, digits = 2), round(MDC80_rel_12LOG, digits = 2)),
                            #MDC80_rel_23 = c(round(MDC80_rel_23, digits = 2), round(MDC80_rel_23LOG, digits = 2)),
                            #MDC80_rel_comb = c(round(MDC80_rel_comb, digits = 2), round(MDC80_rel_combLOG, digits = 2)),

                            #MDC90_12 = c(round(MDC90_12, digits = 2), round(MDC90_12LOG, digits = 2)),
                            #MDC90_23 = c(round(MDC90_23, digits = 2), round(MDC90_23LOG, digits = 2)),
                            MDC90_comb = c(round(MDC90_comb, digits = 2), round(MDC90_combLOG, digits = 2)),
                            #MDC90_rel_12 = c(round(MDC90_rel_12, digits = 2), round(MDC90_rel_12LOG, digits = 2)),
                            #MDC90_rel_23 = c(round(MDC90_rel_23, digits = 2), round(MDC90_rel_23LOG, digits = 2)),
                            #MDC90_rel_comb = c(round(MDC90_rel_comb, digits = 2), round(MDC90_rel_combLOG, digits = 2)),

                            #MDC95_12 = c(round(MDC95_12, digits = 2), round(MDC95_12LOG, digits = 2)),
                            #MDC95_23 = c(round(MDC95_23, digits = 2), round(MDC95_23LOG, digits = 2)),
                            MDC95_comb = c(round(MDC95_comb, digits = 2), round(MDC95_combLOG, digits = 2)),
                            #MDC95_rel_12 = c(round(MDC95_rel_12, digits = 2), round(MDC95_rel_12LOG, digits = 2)),
                            #MDC95_rel_23 = c(round(MDC95_rel_23, digits = 2), round(MDC95_rel_23LOG, digits = 2)),
                            #MDC95_rel_comb = c(round(MDC95_rel_comb, digits = 2), round(MDC95_rel_combLOG, digits = 2)),

                            #MDC99_12 = c(round(MDC99_12, digits = 2), round(MDC99_12LOG, digits = 2)),
                            #MDC99_23 = c(round(MDC99_23, digits = 2), round(MDC99_23LOG, digits = 2)),
                            MDC99_comb = c(round(MDC99_comb, digits = 2), round(MDC99_combLOG, digits = 2)),
                            #MDC99_rel_12 = c(round(MDC99_rel_12, digits = 2), round(MDC99_rel_12LOG, digits = 2)),
                            #MDC99_rel_23 = c(round(MDC99_rel_23, digits = 2), round(MDC99_rel_23LOG, digits = 2)),
                            #MDC99_rel_comb = c(round(MDC99_rel_comb, digits = 2), round(MDC99_rel_combLOG, digits = 2)),

                            MDC70_comb_gp = c(round(MDC70_comb_gp, digits = 4), round(MDC70_combLOG_gp, digits = 4)),
                            MDC80_comb_gp = c(round(MDC80_comb_gp, digits = 4), round(MDC80_combLOG_gp, digits = 4)),
                            MDC90_comb_gp = c(round(MDC90_comb_gp, digits = 4), round(MDC90_combLOG_gp, digits = 4)),
                            MDC95_comb_gp = c(round(MDC95_comb_gp, digits = 4), round(MDC95_combLOG_gp, digits = 4)),
                            MDC99_comb_gp = c(round(MDC99_comb_gp, digits = 4), round(MDC99_combLOG_gp, digits = 4)),

                            #SmallestEffpureSD12 = c(tableSmallestEffpure12, tableSmallestEffpure12log),
                            #SmallestEffpureSD23 = c(tableSmallestEffpure23, tableSmallestEffpure23log),
                            #SmallestEffpureSDcomb = c(tableSmallestEffpureComb, tableSmallestEffpureComblog),
                            #SmallestEffobsSD12 = c(tableSmallestEffobs12, tableSmallestEffobs12log),
                            #SmallestEffobsSD23 = c(tableSmallestEffobs23, tableSmallestEffobs23log),
                            #SmallestEffobsSDcomb = c(tableSmallestEffobsComb, tableSmallestEffobsComblog),

                            #MedEffectpureSD12 = c(tableMedEffectpure12, tableMedEffectpure12log),
                            #MedEffectpureSD23 = c(tableMedEffectpure23, tableMedEffectpure23log),
                            #MedEffectpureSDcomb = c(tableMedEffectpureComb, tableMedEffectpureComblog),
                            #MedEffectobsSD12 = c(tableMedEffectobs12, tableMedEffectobs12log),
                            #MedEffectobsSD23 = c(tableMedEffectobs23, tableMedEffectobs23log),
                            #MedEffectobsSDcomb = c(tableMedEffectobsComb, tableMedEffectobsComblog),

                            #LrgEffectpureSD12 = c(tableLrgEffectpure12, tableLrgEffectpure12log),
                            #LrgEffectpureSD23 = c(tableLrgEffectpure23, tableLrgEffectpure23log),
                            #LrgEffectpureSDcomb = c(tableLrgEffectpureComb, tableLrgEffectpureComblog),
                            #LrgEffectobsSD12 = c(tableLrgEffectobs12, tableLrgEffectobs12log),
                            #LrgEffectobsSD23 = c(tableLrgEffectobs23, tableLrgEffectobs23log),
                            #LrgEffectobsSDcomb = c(tableLrgEffectobsComb, tableLrgEffectobsComblog))

                            grandmean_sampleSize = round(GrandMean, digits = 4), # This is for sample size calculations

                            sdDiff = round(sqrt(2) * ((PercentTEcomb / 100) * GrandMean), digits = 4)) # This is for sample size calculations

  return(ValuesTable)

}
