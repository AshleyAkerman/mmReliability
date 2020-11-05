
#' Title: Extracting standardised effect sizes from data
#'
#'
#' @param data_list list of data tables
#' @param comparison string taking form of either "consecutive" or "incremental
#' @param levels numeric value for the number of comparisons (e.g., 2 levels is PRE and POST; 4 levels is NHS, LHS, MHS, HHS)
#' @keywords reliability, intra-class correlation coefficient, coefficient of variation
#' @return output will be a table (dataframe) of the effect size calculations
#' @export


effectsize_table_function <- function(data_list, comparison, levels) {

  # If the levels = 2 then there is 1 comparison, if there are 4 then there are 3 comparisons
  if (levels == 2) {
    no_comparisons <- 1
  }

  if (levels == 4) {
    no_comparisons <- 3
  }

  # Sequences of the variables.
  # The length of data_list will be equal to the number of variables * number of levels
  # e.g., index 1 = var 1, level 1; index 2 = var 1, level 2; index 3 = var 1, level 3; index 4 = var 1, level 4; index 5 = var 2, level 1.
  # The sequence of variables will therefore give the index equal to the start of a new variable
  variable_sequences <- seq(1, length(data_list), levels)

  # Generate new variables to save for later for consecutive
  Stzd_Effect_L1_cons <- c()
  Stzd_Effect_L2_cons <- c()
  Stzd_Effect_L3_cons <- c()

  Raw_Effect_L1_cons <- c()
  Raw_Effect_L2_cons <- c()
  Raw_Effect_L3_cons <- c()

  Size_L1_cons <- c()
  Size_L2_cons <- c()
  Size_L3_cons <- c()

  SD_of_delta_L1_cons <- c()
  SD_of_delta_L2_cons <- c()
  SD_of_delta_L3_cons <- c()

  Cor_change_L1_cons <- c()
  Cor_change_L2_cons <- c()
  Cor_change_L3_cons <- c()


  # Generate new variables to save for later for incremental
  Stzd_Effect_L1_inc <- c()
  Stzd_Effect_L2_inc <- c()
  Stzd_Effect_L3_inc <- c()

  Raw_Effect_L1_inc <- c()
  Raw_Effect_L2_inc <- c()
  Raw_Effect_L3_inc <- c()

  Size_L1_inc <- c()
  Size_L2_inc <- c()
  Size_L3_inc <- c()

  SD_of_delta_L1_inc <- c()
  SD_of_delta_L2_inc <- c()
  SD_of_delta_L3_inc <- c()

  Cor_change_L1_inc <- c()
  Cor_change_L2_inc <- c()
  Cor_change_L3_inc <- c()

  ##-------------------------------------------------##
  ##-------------------------------------------------##

  ## If there is only 1 comparisons due to 2 levels...

  if(no_comparisons == 1){

    # For every index that is the start of a variable in the list of data frames...
    for (variable in 1: length(variable_sequences)){

      # New variables that will be cleared after looping over each index
      # Relate to different levels (i.e., lasers [microv.] or trials [other manuscripts])
      L1 <- c()
      L2 <- c()
      L3 <- c()

      # This will generate the index for each table that is related to one of the variables, i.e., 4 numbers starting with the first index that is stored in the variable above
      seqs <- seq(variable_sequences[variable], variable_sequences[variable] + 3)

      # Bind together the BL and HHS from each of the tables related to the first, second, and third trial respectively.
      L1 <- cbind(L1, data_list[[seqs[1]]][[2]],
                  data_list[[seqs[2]]][[2]])

      L2 <- cbind(L2, data_list[[seqs[1]]][[3]],
                  data_list[[seqs[2]]][[3]])

      L3 <- cbind(L3, data_list[[seqs[1]]][[4]],
                  data_list[[seqs[2]]][[4]])

      ##-------------------------------------------------##

      # NHS to HHS effect size Level (Laser or Trial) 1
      BLtoHHL1 <- effsize::cohen.d(L1[,2], L1[,1], paired = TRUE, na.rm = TRUE)
      BLtoHHL1_mean <- mean(L1[,2] - L1[,1], na.rm = TRUE)
      BLtoHHL1_sd <- sqrt(sum(((stats::na.omit(L1[, 2] - L1[, 1])) - BLtoHHL1_mean) ^ 2) /
                        (length(stats::na.omit(L1[, 2] - L1[, 1])) - 1))
      BLtoHHL1_cor <- stats::cor(L1[,1], L1[,2], use = "complete.obs")
      # NHS to HHS effect size Level (Laser or Trial) 2
      BLtoHHL2 <- effsize::cohen.d(L2[,2], L2[,1], paired = TRUE, na.rm = TRUE)
      BLtoHHL2_mean <- mean(L2[,2] - L2[,1], na.rm = TRUE)
      BLtoHHL2_sd <- sqrt(sum(((stats::na.omit(L2[, 2] - L2[, 1])) - BLtoHHL2_mean) ^ 2) /
                            (length(stats::na.omit(L2[, 2] - L2[, 1])) - 1))
      BLtoHHL2_cor <- stats::cor(L2[,1], L2[,2], use = "complete.obs")
      # NHS to HHS effect size Level (Laser or Trial) 3
      BLtoHHL3 <- effsize::cohen.d(L3[,2], L3[,1], paired = TRUE, na.rm = TRUE)
      BLtoHHL3_mean <- mean(L3[,2] - L3[,1], na.rm = TRUE)
      BLtoHHL3_sd <- sqrt(sum(((stats::na.omit(L3[, 2] - L3[, 1])) - BLtoHHL3_mean) ^ 2) /
                            (length(stats::na.omit(L3[, 2] - L3[, 1])) - 1))
      BLtoHHL3_cor <- stats::cor(L3[,1], L3[,2], use = "complete.obs")

      ##-------------------------------------------------##
      ##-------------------------------------------------##

      # What is the effect size and associated 95% CI - incremental
      Stzd_Effect_L1_inc <- rbind(Stzd_Effect_L1_inc,
                                  paste(round(BLtoHHL1$estimate, digits = 2), "[",
                                        round(BLtoHHL1$conf.int[1], digits = 2), ",",
                                        round(BLtoHHL1$conf.int[2], digits = 2) ,"]"))


      # What is the effect size and associated 95% CI - incremental
      Stzd_Effect_L2_inc <- rbind(Stzd_Effect_L2_inc,
                                  paste(round(BLtoHHL2$estimate, digits = 2), "[",
                                        round(BLtoHHL2$conf.int[1], digits = 2), ",",
                                        round(BLtoHHL2$conf.int[2], digits = 2) ,"]"))


      # What is the effect size and associated 95% CI - incremental
      Stzd_Effect_L3_inc <- rbind(Stzd_Effect_L3_inc,
                                  paste(round(BLtoHHL3$estimate, digits = 2), "[",
                                        round(BLtoHHL3$conf.int[1], digits = 2), ",",
                                        round(BLtoHHL3$conf.int[2], digits = 2) ,"]"))

      ##-------------------------------------------------##

      # What is the absolute effect for L1 - incremental?
      Raw_Effect_L1_inc <- rbind(Raw_Effect_L1_inc,
                                 round(BLtoHHL1_mean, digits = 2))


      # What is the absolute effect for L2 - incremental?
      Raw_Effect_L2_inc <- rbind(Raw_Effect_L2_inc,
                                 round(BLtoHHL2_mean, digits = 2))


      # What is the absolute effect for L3 - incremental?
      Raw_Effect_L3_inc <- rbind(Raw_Effect_L3_inc,
                                 round(BLtoHHL3_mean, digits = 2))

      ##-------------------------------------------------##

      # What is the SD of the change for L1 - incremental?
      SD_of_delta_L1_inc <- rbind(SD_of_delta_L1_inc,
                                  round(BLtoHHL1_sd, digits = 2))

      # What is the SD of the change for L2 - incremental?
      SD_of_delta_L2_inc <- rbind(SD_of_delta_L2_inc,
                                  round(BLtoHHL2_sd, digits = 2))

      # What is the SD of the change for L3 - incremental?
      SD_of_delta_L3_inc <- rbind(SD_of_delta_L3_inc,
                                  round(BLtoHHL3_sd, digits = 2))

      ##-------------------------------------------------##

      # How many samples after removing missing values for L1 - incremental?
      Size_L1_inc <- rbind(Size_L1_inc,
                           nrow(stats::na.omit(cbind(L1[,1], L1[,2]))))

      # How many samples after removing missing values for L2 - incremental?
      Size_L2_inc <- rbind(Size_L2_inc,
                           nrow(stats::na.omit(cbind(L2[,1], L2[,2]))))

      # How many samples after removing missing values for L3 - incremental?
      Size_L3_inc <- rbind(Size_L3_inc,
                           nrow(stats::na.omit(cbind(L3[,1], L3[,2]))))

      ##-------------------------------------------------##

      # Correlation between measures for L1 - incremental ?
      Cor_change_L1_inc <- rbind(Cor_change_L1_inc,
                                 round(BLtoHHL1_cor, digits = 2))

      # Correlation between measures for L2 - incremental ?
      Cor_change_L2_inc <- rbind(Cor_change_L2_inc,
                                 round(BLtoHHL2_cor, digits = 2))

      # Correlation between measures for L3 - incremental ?
      Cor_change_L3_inc <- rbind(Cor_change_L3_inc,
                                 round(BLtoHHL3_cor, digits = 2))


    }

    ##-------------------------------------------------##

    ## Generate the table dependent on the function input

    # Each variable (n=4) will be repeated 3 times because there are 3 comparisons
    Variable_name <- rep(seq(1, length(variable_names)), each = no_comparisons)

    # There are 3 comparisons for 4 variables - incremental
    Comparison_inc <- rep(c("NHS-HHS"), times = length(variable_names))

    EffectSizeTable <- data.frame(Variable_name = Variable_name,
                                  Comparison = Comparison_inc,

                                  Stzd_Effect_L1_inc = Stzd_Effect_L1_inc,
                                  Stzd_Effect_L2_inc = Stzd_Effect_L2_inc,
                                  Stzd_Effect_L3_inc = Stzd_Effect_L3_inc,

                                  Size_L1_inc = Size_L1_inc,
                                  Size_L2_inc = Size_L2_inc,
                                  Size_L3_inc = Size_L3_inc,

                                  Raw_Effect_L1_inc = Raw_Effect_L1_inc,
                                  Raw_Effect_L2_inc = Raw_Effect_L2_inc,
                                  Raw_Effect_L3_inc = Raw_Effect_L3_inc,

                                  SD_of_delta_L1_inc = SD_of_delta_L1_inc,
                                  SD_of_delta_L2_inc = SD_of_delta_L2_inc,
                                  SD_of_delta_L3_inc = SD_of_delta_L3_inc,

                                  Cor_change_L1_inc = Cor_change_L1_inc,
                                  Cor_change_L2_inc = Cor_change_L2_inc,
                                  Cor_change_L3_inc = Cor_change_L3_inc)


    # Change the values of each to the proper names
    EffectSizeTable$Variable_name <- rep(variable_names, each = no_comparisons)

  }

  ##-------------------------------------------------##
  ##-------------------------------------------------##


  ### 3 comparisons and therefore choice of incremental or consecutive


  ## If there are 3 comparisons due to 4 levels...
  if(no_comparisons == 3){

    # For every index that is the start of a variable in the list of data frames...
    for (variable in 1: length(variable_sequences)){

      # New variables that will be cleared after looping over each index
      # Relate to different levels (i.e., lasers [microv.] or trials [other manuscripts])
      L1 <- c()
      L2 <- c()
      L3 <- c()

      # This will generate the index for each table that is related to one of the variables, i.e., 4 numbers starting with the first index that is stored in the variable above
      seqs <- seq(variable_sequences[variable], variable_sequences[variable] + 3)

      # Bind together the BL, LHS, MHS, and HHS from each of the tables related to the first, second, and third trial respectively.
      L1 <- cbind(L1, data_list[[seqs[1]]][[2]],
                  data_list[[seqs[2]]][[2]],
                  data_list[[seqs[3]]][[2]],
                  data_list[[seqs[4]]][[2]])

      L2 <- cbind(L2, data_list[[seqs[1]]][[3]],
                  data_list[[seqs[2]]][[3]],
                  data_list[[seqs[3]]][[3]],
                  data_list[[seqs[4]]][[3]])

      L3 <- cbind(L3, data_list[[seqs[1]]][[4]],
                  data_list[[seqs[2]]][[4]],
                  data_list[[seqs[3]]][[4]],
                  data_list[[seqs[4]]][[4]])

      ##-------------------------------------------------##

      # NHS to LHS effect size Level (Laser or Trial) 1
      BLtoLHL1 <- effsize::cohen.d(L1[,2], L1[,1], paired = TRUE, na.rm = TRUE)
      BLtoLHL1_mean <- mean(L1[,2] - L1[,1], na.rm = TRUE)
      BLtoLHL1_sd <- sqrt(sum(((stats::na.omit(L1[, 2] - L1[, 1])) - BLtoLHL1_mean) ^ 2) /
                            (length(stats::na.omit(L1[, 2] - L1[, 1])) - 1))
      BLtoLHL1_cor <- stats::cor(L1[,1], L1[,2], use = "complete.obs")
      # NHS to LHS effect size Level (Laser or Trial) 2
      BLtoLHL2 <- effsize::cohen.d(L2[,2], L2[,1], paired = TRUE, na.rm = TRUE)
      BLtoLHL2_mean <- mean(L2[,2] - L2[,1], na.rm = TRUE)
      BLtoLHL2_sd <- sqrt(sum(((stats::na.omit(L2[, 2] - L2[, 1])) - BLtoLHL2_mean) ^ 2) /
                            (length(stats::na.omit(L2[, 2] - L2[, 1])) - 1))
      BLtoLHL2_cor <- stats::cor(L2[,1], L2[,2], use = "complete.obs")
      # NHS to LHS effect size Level (Laser or Trial) 3
      BLtoLHL3 <- effsize::cohen.d(L3[,2], L3[,1], paired = TRUE, na.rm = TRUE)
      BLtoLHL3_mean <- mean(L3[,2] - L3[,1], na.rm = TRUE)
      BLtoLHL3_sd <- sqrt(sum(((stats::na.omit(L3[, 2] - L3[, 1])) - BLtoLHL3_mean) ^ 2) /
                            (length(stats::na.omit(L3[, 2] - L3[, 1])) - 1))
      BLtoLHL3_cor <- stats::cor(L3[,1], L3[,2], use = "complete.obs")

      ##-------------------------------------------------##

      # NHS to MHS effect size Level (Laser or Trial) 1
      BLtoMHL1 <- effsize::cohen.d(L1[,3], L1[,1], paired = TRUE, na.rm = TRUE)
      BLtoMHL1_mean <- mean(L1[,3] - L1[,1], na.rm = TRUE)
      BLtoMHL1_sd <- sqrt(sum(((stats::na.omit(L1[, 3] - L1[, 1])) - BLtoMHL1_mean) ^ 2) /
                            (length(stats::na.omit(L1[, 3] - L1[, 1])) - 1))
      BLtoMHL1_cor <- stats::cor(L1[,1], L1[,3], use = "complete.obs")
      # NHS to MHS effect size Level (Laser or Trial) 2
      BLtoMHL2 <- effsize::cohen.d(L2[,3], L2[,1], paired = TRUE, na.rm = TRUE)
      BLtoMHL2_mean <- mean(L2[,3] - L2[,1], na.rm = TRUE)
      BLtoMHL2_sd <- sqrt(sum(((stats::na.omit(L2[, 3] - L2[, 1])) - BLtoMHL2_mean) ^ 2) /
                            (length(stats::na.omit(L2[, 3] - L2[, 1])) - 1))
      BLtoMHL2_cor <- stats::cor(L2[,1], L2[,3], use = "complete.obs")
      # NHS to MHS effect size Level (Laser or Trial) 3
      BLtoMHL3 <- effsize::cohen.d(L3[,3], L3[,1], paired = TRUE, na.rm = TRUE)
      BLtoMHL3_mean <- mean(L3[,3] - L3[,1], na.rm = TRUE)
      BLtoMHL3_sd <- sqrt(sum(((stats::na.omit(L3[, 3] - L3[, 1])) - BLtoMHL3_mean) ^ 2) /
                            (length(stats::na.omit(L3[, 3] - L3[, 1])) - 1))
      BLtoMHL3_cor <- stats::cor(L3[,1], L3[,3], use = "complete.obs")

      ##-------------------------------------------------##

      # LHS to MHS effect size Level (Laser or Trial) 1
      LHtoMHL1 <- effsize::cohen.d(L1[,3], L1[,2], paired = TRUE, na.rm = TRUE)
      LHtoMHL1_mean <- mean(L1[,3] - L1[,2], na.rm = TRUE)
      LHtoMHL1_sd <- sqrt(sum(((stats::na.omit(L1[, 3] - L1[, 2])) - LHtoMHL1_mean) ^ 2) /
                            (length(stats::na.omit(L1[, 3] - L1[, 2])) - 1))
      LHtoMHL1_cor <- stats::cor(L1[,2], L1[,3], use = "complete.obs")
      # LHS to MHS effect size Level (Laser or Trial) 2
      LHtoMHL2 <- effsize::cohen.d(L2[,3], L2[,2], paired = TRUE, na.rm = TRUE)
      LHtoMHL2_mean <- mean(L2[,3] - L2[,2], na.rm = TRUE)
      LHtoMHL2_sd <- sqrt(sum(((stats::na.omit(L2[, 3] - L2[, 2])) - LHtoMHL2_mean) ^ 2) /
                            (length(stats::na.omit(L2[, 3] - L2[, 2])) - 1))
      LHtoMHL2_cor <- stats::cor(L2[,2], L2[,3], use = "complete.obs")
      # LHS to MHS effect size Level (Laser or Trial) 3
      LHtoMHL3 <- effsize::cohen.d(L3[,3], L3[,2], paired = TRUE, na.rm = TRUE)
      LHtoMHL3_mean <- mean(L3[,3] - L3[,2], na.rm = TRUE)
      LHtoMHL3_sd <- sqrt(sum(((stats::na.omit(L3[, 3] - L3[, 2])) - LHtoMHL3_mean) ^ 2) /
                            (length(stats::na.omit(L3[, 3] - L3[, 2])) - 1))
      LHtoMHL3_cor <- stats::cor(L3[,2], L3[,3], use = "complete.obs")

      ##-------------------------------------------------##

      # NHS to HHS effect size Level (Laser or Trial) 1
      BLtoHHL1 <- effsize::cohen.d(L1[,4], L1[,1], paired = TRUE, na.rm = TRUE)
      BLtoHHL1_mean <- mean(L1[,4] - L1[,1], na.rm = TRUE)
      BLtoHHL1_sd <- sqrt(sum(((stats::na.omit(L1[, 4] - L1[, 1])) - BLtoHHL1_mean) ^ 2) /
                            (length(stats::na.omit(L1[, 4] - L1[, 1])) - 1))
      BLtoHHL1_cor <- stats::cor(L1[,1], L1[,4], use = "complete.obs")
      # NHS to HHS effect size Level (Laser or Trial) 2
      BLtoHHL2 <- effsize::cohen.d(L2[,4], L2[,1], paired = TRUE, na.rm = TRUE)
      BLtoHHL2_mean <- mean(L2[,4] - L2[,1], na.rm = TRUE)
      BLtoHHL2_sd <- sqrt(sum(((stats::na.omit(L2[, 4] - L2[, 1])) - BLtoHHL2_mean) ^ 2) /
                            (length(stats::na.omit(L2[, 4] - L2[, 1])) - 1))
      BLtoHHL2_cor <- stats::cor(L2[,1], L2[,4], use = "complete.obs")
      # NHS to HHS effect size Level (Laser or Trial) 3
      BLtoHHL3 <- effsize::cohen.d(L3[,4], L3[,1], paired = TRUE, na.rm = TRUE)
      BLtoHHL3_mean <- mean(L3[,4] - L3[,1], na.rm = TRUE)
      BLtoHHL3_sd <- sqrt(sum(((stats::na.omit(L3[, 4] - L3[, 1])) - BLtoHHL3_mean) ^ 2) /
                            (length(stats::na.omit(L3[, 4] - L3[, 1])) - 1))
      BLtoHHL3_cor <- stats::cor(L3[,1], L3[,4], use = "complete.obs")

      ##-------------------------------------------------##

      # MHS to HHS effect size Level (Laser or Trial) 1
      MHtoHHL1 <- effsize::cohen.d(L1[,4], L1[,3], paired = TRUE, na.rm = TRUE)
      MHtoHHL1_mean <- mean(L1[,4] - L1[,3], na.rm = TRUE)
      MHtoHHL1_sd <- sqrt(sum(((stats::na.omit(L1[, 4] - L1[, 3])) - MHtoHHL1_mean) ^ 2) /
                            (length(stats::na.omit(L1[, 4] - L1[, 3])) - 1))
      MHtoHHL1_cor <- stats::cor(L1[,3], L1[,4], use = "complete.obs")
      # MHS to HHS effect size Level (Laser or Trial) 2
      MHtoHHL2 <- effsize::cohen.d(L2[,4], L2[,3], paired = TRUE, na.rm = TRUE)
      MHtoHHL2_mean <- mean(L2[,4] - L2[,3], na.rm = TRUE)
      MHtoHHL2_sd <- sqrt(sum(((stats::na.omit(L2[, 4] - L2[, 3])) - MHtoHHL2_mean) ^ 2) /
                            (length(stats::na.omit(L2[, 4] - L2[, 3])) - 1))
      MHtoHHL2_cor <- stats::cor(L2[,3], L2[,4], use = "complete.obs")
      # MHS to HHS effect size Level (Laser or Trial) 3
      MHtoHHL3 <- effsize::cohen.d(L3[,4], L3[,3], paired = TRUE, na.rm = TRUE)
      MHtoHHL3_mean <- mean(L3[,4] - L3[,3], na.rm = TRUE)
      MHtoHHL3_sd <- sqrt(sum(((stats::na.omit(L3[, 4] - L3[, 3])) - MHtoHHL3_mean) ^ 2) /
                            (length(stats::na.omit(L3[, 4] - L3[, 3])) - 1))
      MHtoHHL3_cor <- stats::cor(L3[,3], L3[,4], use = "complete.obs")

      ##-------------------------------------------------##
      ##-------------------------------------------------##

      # What is the effect size for L1 and associated 95% CI - consecutive
      Stzd_Effect_L1_cons <- rbind(Stzd_Effect_L1_cons,
                                   paste(round(BLtoLHL1$estimate, digits = 2), "[",
                                         round(BLtoLHL1$conf.int[1], digits = 2), ",",
                                         round(BLtoLHL1$conf.int[2], digits = 2) ,"]"),
                                   paste(round(LHtoMHL1$estimate, digits = 2), "[",
                                         round(LHtoMHL1$conf.int[1], digits = 2), ",",
                                         round(LHtoMHL1$conf.int[2], digits = 2) ,"]"),
                                   paste(round(MHtoHHL1$estimate, digits = 2), "[",
                                         round(MHtoHHL1$conf.int[1], digits = 2), ",",
                                         round(MHtoHHL1$conf.int[2], digits = 2) ,"]"))

      # What is the effect size and associated 95% CI - incremental
      Stzd_Effect_L1_inc <- rbind(Stzd_Effect_L1_inc,
                                  paste(round(BLtoLHL1$estimate, digits = 2), "[",
                                        round(BLtoLHL1$conf.int[1], digits = 2), ",",
                                        round(BLtoLHL1$conf.int[2], digits = 2) ,"]"),
                                  paste(round(BLtoMHL1$estimate, digits = 2), "[",
                                        round(BLtoMHL1$conf.int[1], digits = 2), ",",
                                        round(BLtoMHL1$conf.int[2], digits = 2) ,"]"),
                                  paste(round(BLtoHHL1$estimate, digits = 2), "[",
                                        round(BLtoHHL1$conf.int[1], digits = 2), ",",
                                        round(BLtoHHL1$conf.int[2], digits = 2) ,"]"))

      ##-------------------------------------------------##

      # What is the effect size for L1 and associated 95% CI - consecutive
      Stzd_Effect_L2_cons <- rbind(Stzd_Effect_L2_cons,
                                   paste(round(BLtoLHL2$estimate, digits = 2), "[",
                                         round(BLtoLHL2$conf.int[1], digits = 2), ",",
                                         round(BLtoLHL2$conf.int[2], digits = 2) ,"]"),
                                   paste(round(LHtoMHL2$estimate, digits = 2), "[",
                                         round(LHtoMHL2$conf.int[1], digits = 2), ",",
                                         round(LHtoMHL2$conf.int[2], digits = 2) ,"]"),
                                   paste(round(MHtoHHL2$estimate, digits = 2), "[",
                                         round(MHtoHHL2$conf.int[1], digits = 2), ",",
                                         round(MHtoHHL2$conf.int[2], digits = 2) ,"]"))

      # What is the effect size and associated 95% CI - incremental
      Stzd_Effect_L2_inc <- rbind(Stzd_Effect_L2_inc,
                                  paste(round(BLtoLHL2$estimate, digits = 2), "[",
                                        round(BLtoLHL2$conf.int[1], digits = 2), ",",
                                        round(BLtoLHL2$conf.int[2], digits = 2) ,"]"),
                                  paste(round(BLtoMHL2$estimate, digits = 2), "[",
                                        round(BLtoMHL2$conf.int[1], digits = 2), ",",
                                        round(BLtoMHL2$conf.int[2], digits = 2) ,"]"),
                                  paste(round(BLtoHHL2$estimate, digits = 2), "[",
                                        round(BLtoHHL2$conf.int[1], digits = 2), ",",
                                        round(BLtoHHL2$conf.int[2], digits = 2) ,"]"))

      ##-------------------------------------------------##

      # What is the effect size for L1 and associated 95% CI - consecutive
      Stzd_Effect_L3_cons <- rbind(Stzd_Effect_L3_cons,
                                   paste(round(BLtoLHL3$estimate, digits = 2), "[",
                                         round(BLtoLHL3$conf.int[1], digits = 2), ",",
                                         round(BLtoLHL3$conf.int[2], digits = 2) ,"]"),
                                   paste(round(LHtoMHL3$estimate, digits = 2), "[",
                                         round(LHtoMHL3$conf.int[1], digits = 2), ",",
                                         round(LHtoMHL3$conf.int[2], digits = 2) ,"]"),
                                   paste(round(MHtoHHL3$estimate, digits = 2), "[",
                                         round(MHtoHHL3$conf.int[1], digits = 2), ",",
                                         round(MHtoHHL3$conf.int[2], digits = 2) ,"]"))

      # What is the effect size and associated 95% CI - incremental
      Stzd_Effect_L3_inc <- rbind(Stzd_Effect_L3_inc,
                                  paste(round(BLtoLHL3$estimate, digits = 2), "[",
                                        round(BLtoLHL3$conf.int[1], digits = 2), ",",
                                        round(BLtoLHL3$conf.int[2], digits = 2) ,"]"),
                                  paste(round(BLtoMHL3$estimate, digits = 2), "[",
                                        round(BLtoMHL3$conf.int[1], digits = 2), ",",
                                        round(BLtoMHL3$conf.int[2], digits = 2) ,"]"),
                                  paste(round(BLtoHHL3$estimate, digits = 2), "[",
                                        round(BLtoHHL3$conf.int[1], digits = 2), ",",
                                        round(BLtoHHL3$conf.int[2], digits = 2) ,"]"))

      ##-------------------------------------------------##

      # What is the absolute effect for L1 - consecutive?
      Raw_Effect_L1_cons <- rbind(Raw_Effect_L1_cons,
                                  round(BLtoLHL1_mean, digits = 2),
                                  round(LHtoMHL1_mean, digits = 2),
                                  round(MHtoHHL1_mean, digits = 2))

      # What is the absolute effect for L1 - incremental?
      Raw_Effect_L1_inc <- rbind(Raw_Effect_L1_inc,
                                 round(BLtoLHL1_mean, digits = 2),
                                 round(BLtoMHL1_mean, digits = 2),
                                 round(BLtoHHL1_mean, digits = 2))

      ##-------------------------------------------------##

      # What is the absolute effect for L2 - consecutive?
      Raw_Effect_L2_cons <- rbind(Raw_Effect_L2_cons,
                                  round(BLtoLHL2_mean, digits = 2),
                                  round(LHtoMHL2_mean, digits = 2),
                                  round(MHtoHHL2_mean, digits = 2))

      # What is the absolute effect for L2 - incremental?
      Raw_Effect_L2_inc <- rbind(Raw_Effect_L2_inc,
                                 round(BLtoLHL2_mean, digits = 2),
                                 round(BLtoMHL2_mean, digits = 2),
                                 round(BLtoHHL2_mean, digits = 2))

      ##-------------------------------------------------##

      # What is the absolute effect for L1 - consecutive?
      Raw_Effect_L3_cons <- rbind(Raw_Effect_L3_cons,
                                  round(BLtoLHL3_mean, digits = 2),
                                  round(LHtoMHL3_mean, digits = 2),
                                  round(MHtoHHL3_mean, digits = 2))

      # What is the absolute effect for L3 - incremental?
      Raw_Effect_L3_inc <- rbind(Raw_Effect_L3_inc,
                                 round(BLtoLHL3_mean, digits = 2),
                                 round(BLtoMHL3_mean, digits = 2),
                                 round(BLtoHHL3_mean, digits = 2))

      ##-------------------------------------------------##

      # What is the SD of the chaange for L1 - consecutive?
      SD_of_delta_L1_cons <- rbind(SD_of_delta_L1_cons,
                                   round(BLtoLHL1_sd, digits = 2),
                                   round(LHtoMHL1_sd, digits = 2),
                                   round(MHtoHHL1_sd, digits = 2))

      # What is the SD of the chaange for L1 - incremental?
      SD_of_delta_L1_inc <- rbind(SD_of_delta_L1_inc,
                                  round(BLtoLHL1_sd, digits = 2),
                                  round(BLtoMHL1_sd, digits = 2),
                                  round(BLtoHHL1_sd, digits = 2))

      ##-------------------------------------------------##

      # What is the SD of the chaange for L2 - consecutive?
      SD_of_delta_L2_cons <- rbind(SD_of_delta_L2_cons,
                                   round(BLtoLHL2_sd, digits = 2),
                                   round(LHtoMHL2_sd, digits = 2),
                                   round(MHtoHHL2_sd, digits = 2))

      # What is the SD of the chaange for L2 - incremental?
      SD_of_delta_L2_inc <- rbind(SD_of_delta_L2_inc,
                                  round(BLtoLHL2_sd, digits = 2),
                                  round(BLtoMHL2_sd, digits = 2),
                                  round(BLtoHHL2_sd, digits = 2))

      ##-------------------------------------------------##

      # What is the SD of the chaange for L3 - consecutive?
      SD_of_delta_L3_cons <- rbind(SD_of_delta_L3_cons,
                                   round(BLtoLHL3_sd, digits = 2),
                                   round(LHtoMHL3_sd, digits = 2),
                                   round(MHtoHHL3_sd, digits = 2))

      # What is the SD of the chaange for L3 - incremental?
      SD_of_delta_L3_inc <- rbind(SD_of_delta_L3_inc,
                                  round(BLtoLHL3_sd, digits = 2),
                                  round(BLtoMHL3_sd, digits = 2),
                                  round(BLtoHHL3_sd, digits = 2))

      ##-------------------------------------------------##

      # How many samples after removing missing values for L1 - consecutive?
      Size_L1_cons <- rbind(Size_L1_cons,
                            nrow(stats::na.omit(cbind(L1[,1], L1[,2]))),
                            nrow(stats::na.omit(cbind(L1[,2], L1[,3]))),
                            nrow(stats::na.omit(cbind(L1[,3], L1[,4]))))

      # How many samples after removing missing values for L1 - incremental?
      Size_L1_inc <- rbind(Size_L1_inc,
                           nrow(stats::na.omit(cbind(L1[,1], L1[,2]))),
                           nrow(stats::na.omit(cbind(L1[,1], L1[,3]))),
                           nrow(stats::na.omit(cbind(L1[,1], L1[,4]))))

      ##-------------------------------------------------##

      # How many samples after removing missing values for L2 - consecutive?
      Size_L2_cons <- rbind(Size_L2_cons,
                            nrow(stats::na.omit(cbind(L2[,1], L2[,2]))),
                            nrow(stats::na.omit(cbind(L2[,2], L2[,3]))),
                            nrow(stats::na.omit(cbind(L2[,3], L2[,4]))))

      # How many samples after removing missing values for L2 - incremental?
      Size_L2_inc <- rbind(Size_L2_inc,
                           nrow(stats::na.omit(cbind(L2[,1], L2[,2]))),
                           nrow(stats::na.omit(cbind(L2[,1], L2[,3]))),
                           nrow(stats::na.omit(cbind(L2[,1], L2[,4]))))

      ##-------------------------------------------------##

      # How many samples after removing missing values for L3 - consecutive?
      Size_L3_cons <- rbind(Size_L3_cons,
                            nrow(stats::na.omit(cbind(L3[,1], L3[,2]))),
                            nrow(stats::na.omit(cbind(L3[,2], L3[,3]))),
                            nrow(stats::na.omit(cbind(L3[,3], L3[,4]))))

      # How many samples after removing missing values for L3 - incremental?
      Size_L3_inc <- rbind(Size_L3_inc,
                           nrow(stats::na.omit(cbind(L3[,1], L3[,2]))),
                           nrow(stats::na.omit(cbind(L3[,1], L3[,3]))),
                           nrow(stats::na.omit(cbind(L3[,1], L3[,4]))))

      ##-------------------------------------------------##

      # Correlation between measures for L1 - consecutive ?
      Cor_change_L1_cons <- rbind(Cor_change_L1_cons,
                                  round(BLtoLHL1_cor, digits = 2),
                                  round(LHtoMHL1_cor, digits = 2),
                                  round(MHtoHHL1_cor, digits = 2))

      # Correlation between measures for L1 - incremental ?
      Cor_change_L1_inc <- rbind(Cor_change_L1_inc,
                                 round(BLtoLHL1_cor, digits = 2),
                                 round(BLtoMHL1_cor, digits = 2),
                                 round(BLtoHHL1_cor, digits = 2))

      ##-------------------------------------------------##

      # Correlation between measures for L2 - consecutive ?
      Cor_change_L2_cons <- rbind(Cor_change_L2_cons,
                                  round(BLtoLHL2_cor, digits = 2),
                                  round(LHtoMHL2_cor, digits = 2),
                                  round(MHtoHHL2_cor, digits = 2))

      # Correlation between measures for L2 - incremental ?
      Cor_change_L2_inc <- rbind(Cor_change_L2_inc,
                                 round(BLtoLHL2_cor, digits = 2),
                                 round(BLtoMHL2_cor, digits = 2),
                                 round(BLtoHHL2_cor, digits = 2))

      ##-------------------------------------------------##

      # Correlation between measures for L3 - consecutive ?
      Cor_change_L3_cons <- rbind(Cor_change_L3_cons,
                                  round(BLtoLHL3_cor, digits = 2),
                                  round(LHtoMHL3_cor, digits = 2),
                                  round(MHtoHHL3_cor, digits = 2))

      # Correlation between measures for L3 - incremental ?
      Cor_change_L3_inc <- rbind(Cor_change_L3_inc,
                                 round(BLtoLHL3_cor, digits = 2),
                                 round(BLtoMHL3_cor, digits = 2),
                                 round(BLtoHHL3_cor, digits = 2))

    }

    ##-------------------------------------------------##

    ## Generate the table dependent on the function input

    # Each variable (n=4) will be repeated 3 times because there are 3 comparisons
    Variable_name <- rep(seq(1, length(variable_names)), each = no_comparisons)

    # There are 3 comparisons for 4 variables - consecutive
    Comparison_cons <- rep(c("NHS-LHS", "LHS-MHS", "MHS-HHS"), times = length(variable_names))
    # There are 3 comparisons for 4 variables - incremental
    Comparison_inc <- rep(c("NHS-LHS", "NHS-MHS", "NHS-HHS"), times = length(variable_names))

    # Consecutive table
    if (comparison == "consecutive"){

      # Combine the data
      EffectSizeTable <- data.frame(Variable_name = Variable_name,
                                    Comparison = Comparison_cons,

                                    Stzd_Effect_L1_cons = Stzd_Effect_L1_cons,
                                    Stzd_Effect_L2_cons = Stzd_Effect_L2_cons,
                                    Stzd_Effect_L3_cons = Stzd_Effect_L3_cons,

                                    Size_L1_cons = Size_L1_cons,
                                    Size_L2_cons = Size_L2_cons,
                                    Size_L3_cons = Size_L3_cons,

                                    Raw_Effect_L1_cons = Raw_Effect_L1_cons,
                                    Raw_Effect_L2_cons = Raw_Effect_L2_cons,
                                    Raw_Effect_L3_cons = Raw_Effect_L3_cons,

                                    SD_of_delta_L1_cons = SD_of_delta_L1_cons,
                                    SD_of_delta_L2_cons = SD_of_delta_L2_cons,
                                    SD_of_delta_L3_cons = SD_of_delta_L3_cons,

                                    Cor_change_L1_cons = Cor_change_L1_cons,
                                    Cor_change_L2_cons = Cor_change_L2_cons,
                                    Cor_change_L3_cons = Cor_change_L3_cons)


      # Change the values of each to the proper names
      EffectSizeTable$Variable_name <- rep(variable_names, each = no_comparisons)

    }

    # Incremental table
    if (comparison == "incremental"){

      # Combine the data
      EffectSizeTable <- data.frame(Variable_name = Variable_name,
                                    Comparison = Comparison_inc,

                                    Stzd_Effect_L1_inc = Stzd_Effect_L1_inc,
                                    Stzd_Effect_L2_inc = Stzd_Effect_L2_inc,
                                    Stzd_Effect_L3_inc = Stzd_Effect_L3_inc,

                                    Size_L1_inc = Size_L1_inc,
                                    Size_L2_inc = Size_L2_inc,
                                    Size_L3_inc = Size_L3_inc,

                                    Raw_Effect_L1_inc = Raw_Effect_L1_inc,
                                    Raw_Effect_L2_inc = Raw_Effect_L2_inc,
                                    Raw_Effect_L3_inc = Raw_Effect_L3_inc,

                                    SD_of_delta_L1_inc = SD_of_delta_L1_inc,
                                    SD_of_delta_L2_inc = SD_of_delta_L2_inc,
                                    SD_of_delta_L3_inc = SD_of_delta_L3_inc,

                                    Cor_change_L1_inc = Cor_change_L1_inc,
                                    Cor_change_L2_inc = Cor_change_L2_inc,
                                    Cor_change_L3_inc = Cor_change_L3_inc)


      # Change the values of each to the proper names
      EffectSizeTable$Variable_name <- rep(variable_names, each = no_comparisons)

    }


  }


  return(EffectSizeTable)

}

