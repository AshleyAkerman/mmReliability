#' Heart rate variability dataset for use with Myths and Methodologies Reliability Papers
#'
#' Average (mean) heart rate is designated as ‘AverageRate’, whereas the mean and median of the RR interval
#' are designated as ‘AverageRR’ and ‘MedianRR’ respectively.
#' SD of the heart rate is designated as ‘SDRate’ whereas SD of the RR interval is designated as “SDRR’
#' The coefficient of variation of the RR interval (i.e., ratio between SDRR and AverageRR) is designated
#' as ‘CVRR’ and reflects overall autonomic activity.
#' The SD of the successive RR interval differences is designated as ‘SDSD’ and reflects short-term variability.
#' The root mean square of the successive difference between normal heart beats is designated as ‘RMSSD’
#' and reflects overall parasympathetic activity.
#' The total power (in the frequency domain) of RR interval variability is the total variance and is designated
#'  as ‘TotalPower’
#' The power in the low frequency and high frequency bands (in the frequency domain) are designated as ‘LFms, LFrel, LFnu’
#' and ‘HFms, HFrel, HFnu’ and the ratio between these measures as ‘LF_HF’. When power in each frequency band is estimated
#' relative to TotalPower it is expressed as a percentage (‘rel’) and in normal units (‘nu’).
#' Nonlinear assessments of the SD of instantaneous beat-to-beat interval variability (derived from the Poincare plot) is
#' designated as ‘SD1’, and the continuous long-term RR interval variability (from the same plot) is designated as ‘SD2’.

#'
#' @docType data
#'
#' @usage data(heart_rate_variability_data)
#'
#' @keywords datasets

"heart_rate_variability_data"


