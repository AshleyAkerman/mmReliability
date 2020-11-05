
#' Title: Formatting sample size data and power plots
#'
#' NOTES:
#' In the case of sudomotor, the one level will be designated via NHS
#' In the case of cardiovascular it will be NHS and HHS
#'
#' @param df dataframe
#' @param NHS string input taking form of "Y" or "N"
#' @param LHS string input taking form of "Y" or "N"
#' @param MHS string input taking form of "Y" or "N"
#' @param HHS string input taking form of "Y" or "N"
#' @return a dataframe with the number of participants required for each sample size and each effect size
#' @keywords reliability, intra-class correlation coefficient, coefficient of variation
#' @export
#'

formatting_req_samples_function <- function(df, NHS, LHS, MHS, HHS){

  ## Now need to add the rows of the data together dependent on the function input
  # Few combinations in manuscripts, so again will only provide those options

  ##------------------------------##
  ##------------------------------##

  # If there are all 4 levels - sudomotor, microvascular, cardiovascular
  if(NHS == "Y" & LHS == "Y" & MHS == "Y" & HHS == "Y") {

    # Generate new lists dependent on the inputs
    list_sample_sizes_BL <- list()
    list_sample_sizes_LH <- list()
    list_sample_sizes_MH <- list()
    list_sample_sizes_HH <- list()

    # For each variable, bind together the within subject and between subject elements
    # Need to have separate lists for each level of heat strain because they can be separate lengths
    for (i in 1: length(variable_names)){

      list_sample_sizes_BL[[i]] <- cbind(
        df[[i]][[2]],
        df[[i]][[3]])

      list_sample_sizes_LH[[i]] <- cbind(
        df[[i]][[5]],
        df[[i]][[6]])

      list_sample_sizes_MH[[i]] <- cbind(
        df[[i]][[8]],
        df[[i]][[9]])

      list_sample_sizes_HH[[i]] <- cbind(
        df[[i]][[11]],
        df[[i]][[12]])


    }

    # Combine the lists
    full_list <- list(list_sample_sizes_BL,
                      list_sample_sizes_LH,
                      list_sample_sizes_MH,
                      list_sample_sizes_HH)

    ##------------------------------##

    ###
    # Now need to extract a dataframe with just the sample sizes and effect size inputs
    ###

    # Need to know how long each of the column lengths are for each level of heat strain and each outcome measure
    column_lengths <- c()

    # There will be the number of major indicees equal to the numberer of heat strain levels.
    # So for each heating stage (i.e., major index in the list), and each outcome measure, add the column lengths to the running vector
    for(stages in 1: length(full_list)) {

      for (measure in 1: length(full_list[[stages]])) {

        column_lengths <- c(column_lengths, nrow(full_list[[stages]][[measure]]))
      }

    }

    ## Repeat the names of the variables the number of times required (i.e., length of the columns) for all of the baseline

    # The column lengths are a list of equal to the number of variables * number of levels of heat strain
    # So it will repeat the name of the ith vector equal to the length of the ith vector of the column lengths

    names_df_BL <- c()
    for(i in 1: length(variable_names)){
      names_df_BL <- c(names_df_BL, rep(variable_names[i], times = column_lengths[i]))

    }

    no_heat_within <- c()
    no_heat_between <- c()

    # Repeat the names of the variables the number of times required (i.e., length of the columns) for all of the low heat
    names_df_LH <- c()
    for(i in 1: length(variable_names)){
      names_df_LH <- c(names_df_LH, rep(variable_names[i],
                                        times = column_lengths[length(variable_names)+i]))

    }

    low_heat_within <- c()
    low_heat_between <- c()

    # Repeat the names of the variables the number of times required (i.e., length of the columns) for all of the medium heat
    # Because the column_lengths vector is the combination of LH and MH, need to start the index at the level corresponding to the second half
    names_df_MH <- c()
    for(i in 1: length(variable_names)){
      names_df_MH <- c(names_df_MH, rep(variable_names[i],
                                        times = column_lengths[c(length(variable_names) +
                                                                   length(variable_names) + i)]))

    }

    med_heat_within <- c()
    med_heat_between <- c()

    # Repeat the names of the variables the number of times required (i.e., length of the columns) for all of the medium heat
    # Because the column_lengths vector is the combination of LH and MH, need to start the index at the level corresponding to the second half
    names_df_HH <- c()
    for(i in 1: length(variable_names)){
      names_df_HH <- c(names_df_HH, rep(variable_names[i],
                                        times = column_lengths[c(length(variable_names) +
                                                                   length(variable_names)+
                                                                   length(variable_names) + i)]))

    }

    high_heat_within <- c()
    high_heat_between <- c()

    ##------------------------------##

    ## Now need to add the rows of the data together dependent on the function input
    # Few combinations in manuscripts, so again will only provide those options
    for (i in 1: length(variable_names)){

      no_heat_within <- rbind(no_heat_within, full_list[[1]][[i]][,c(1:2)])
      no_heat_between <- rbind(no_heat_between, full_list[[1]][[i]][,c(3:4)])
      low_heat_within <- rbind(low_heat_within, full_list[[2]][[i]][,c(1:2)])
      low_heat_between <- rbind(low_heat_between, full_list[[2]][[i]][,c(3:4)])
      med_heat_within <- rbind(med_heat_within, full_list[[3]][[i]][,c(1:2)])
      med_heat_between <- rbind(med_heat_between, full_list[[3]][[i]][,c(3:4)])
      high_heat_within <- rbind(high_heat_within, full_list[[4]][[i]][,c(1:2)])
      high_heat_between <- rbind(high_heat_between, full_list[[4]][[i]][,c(3:4)])

    }

    # repeated the number of times equal to 1: number of variables (BL), number of variables + 1 to number of variables * 2 (LH), etc...
    data_frame_return <- data.frame(stages = c(
      rep(c("BL"), times = sum(column_lengths[1:length(variable_names)])),
      rep(c("LH"), times = sum(column_lengths[(length(variable_names) + 1):(length(variable_names) * 2)])),
      rep(c("MH"), times = sum(column_lengths[((length(variable_names) * 2) + 1): (length(variable_names) * 3)])),
      rep(c("HH"), times = sum(column_lengths[((length(variable_names) * 3) + 1): (length(variable_names) * 4)]))),
      outcome = c(names_df_BL, names_df_LH, names_df_MH, names_df_HH),
      within = rbind(no_heat_within, low_heat_within, med_heat_within, high_heat_within),
      between = rbind(no_heat_between, low_heat_between, med_heat_between, high_heat_between))

  }

  ##------------------------------##
  ##------------------------------##
  ## NEXT OPTION
  ##------------------------------##
  ##------------------------------##


  # If there are just 3 levels (i.e., no baseline) - sudomotor
  if(NHS == "N" & LHS == "Y" & MHS == "Y" & HHS == "Y") {
    # Generate new lists dependent on the inputs
    list_sample_sizes_LH <- list()
    list_sample_sizes_MH <- list()
    list_sample_sizes_HH <- list()

    # For each variable, bind together the within subject and between subject elements
    # Need to have separate lists for each level of heat strain because they can be separate lengths
    for (i in 1: length(variable_names)){

      list_sample_sizes_LH[[i]] <- cbind(
        df[[i]][[2]],
        df[[i]][[3]])

      list_sample_sizes_MH[[i]] <- cbind(
        df[[i]][[5]],
        df[[i]][[6]])

      list_sample_sizes_HH[[i]] <- cbind(
        df[[i]][[8]],
        df[[i]][[9]])


    }

    # Combine the lists
    full_list <- list(list_sample_sizes_LH,
                      list_sample_sizes_MH,
                      list_sample_sizes_HH)

    ##------------------------------##

    ###
    # Now need to extract a dataframe with just the sample sizes and effect size inputs
    ###

    # Need to know how long each of the column lengths are for each level of heat strain and each outcome measure
    column_lengths <- c()

    # There will be the number of major indicees equal to the numberer of heat strain levels.
    # So for each heating stage (i.e., major index in the list), and each outcome measure, add the column lengths to the running vector
    for(stages in 1: length(full_list)) {

      for (measure in 1: length(full_list[[stages]])) {

        column_lengths <- c(column_lengths, nrow(full_list[[stages]][[measure]]))
      }

    }

    ## Repeat the names of the variables the number of times required (i.e., length of the columns) for all of the baseline

    # The column lengths are a list of equal to the number of variables * number of levels of heat strain
    # So it will repeat the name of the ith vector equal to the length of the ith vector of the column lengths

    names_df_LH <- c()
    for(i in 1: length(variable_names)){
      names_df_LH <- c(names_df_LH, rep(variable_names[i], times = column_lengths[i]))

    }

    low_heat_within <- c()
    low_heat_between <- c()

    # Repeat the names of the variables the number of times required (i.e., length of the columns) for all of the low heat
    names_df_MH <- c()
    for(i in 1: length(variable_names)){
      names_df_MH <- c(names_df_MH, rep(variable_names[i],
                                        times = column_lengths[length(variable_names)+i]))

    }

    med_heat_within <- c()
    med_heat_between <- c()

    # Repeat the names of the variables the number of times required (i.e., length of the columns) for all of the medium heat
    # Because the column_lengths vector is the combination of LH and MH, need to start the index at the level corresponding to the second half
    names_df_HH <- c()
    for(i in 1: length(variable_names)){
      names_df_HH <- c(names_df_HH, rep(variable_names[i],
                                        times = column_lengths[c(length(variable_names) +
                                                                   length(variable_names) + i)]))

    }

    high_heat_within <- c()
    high_heat_between <- c()

    ##------------------------------##

    ## Now need to add the rows of the data together dependent on the function input
    # Few combinations in manuscripts, so again will only provide those options
    for (i in 1: length(variable_names)){
      low_heat_within <- rbind(low_heat_within, full_list[[1]][[i]][,c(1:2)])
      low_heat_between <- rbind(low_heat_between, full_list[[1]][[i]][,c(3:4)])
      med_heat_within <- rbind(med_heat_within, full_list[[2]][[i]][,c(1:2)])
      med_heat_between <- rbind(med_heat_between, full_list[[2]][[i]][,c(3:4)])
      high_heat_within <- rbind(high_heat_within, full_list[[3]][[i]][,c(1:2)])
      high_heat_between <- rbind(high_heat_between, full_list[[3]][[i]][,c(3:4)])
    }

    # repeated the number of times equal to 1: number of variables (BL), number of variables + 1 to number of variables * 2 (LH), etc...
    data_frame_return <- data.frame(stages = c(
      rep(c("LH"), times = sum(column_lengths[1:length(variable_names)])),
      rep(c("MH"), times = sum(column_lengths[(length(variable_names) + 1): (length(variable_names) * 2)])),
      rep(c("HH"), times = sum(column_lengths[((length(variable_names) * 2) + 1): (length(variable_names) * 3)]))),
      outcome = c(names_df_LH, names_df_MH, names_df_HH),
      within = rbind(low_heat_within, med_heat_within, high_heat_within),
      between = rbind(low_heat_between, med_heat_between, high_heat_between))

  }

  ##------------------------------##
  ##------------------------------##
  ## NEXT OPTION
  ##------------------------------##
  ##------------------------------##


  # If there are just 3 levels (i.e., no high heat strain) - microvascular
  if(NHS == "Y" & LHS == "Y" & MHS == "Y" & HHS == "N") {
    # Generate new lists dependent on the inputs
    list_sample_sizes_BL <- list()
    list_sample_sizes_LH <- list()
    list_sample_sizes_MH <- list()

    # For each variable, bind together the within subject and between subject elements
    # Need to have separate lists for each level of heat strain because they can be separate lengths
    for (i in 1: length(variable_names)){

      list_sample_sizes_BL[[i]] <- cbind(
        df[[i]][[2]],
        df[[i]][[3]])

      list_sample_sizes_LH[[i]] <- cbind(
        df[[i]][[5]],
        df[[i]][[6]])

      list_sample_sizes_MH[[i]] <- cbind(
        df[[i]][[8]],
        df[[i]][[9]])


    }

    # Combine the lists
    full_list <- list(list_sample_sizes_BL,
                      list_sample_sizes_LH,
                      list_sample_sizes_MH)


    ##------------------------------##

    ###
    # Now need to extract a dataframe with just the sample sizes and effect size inputs
    ###

    # Need to know how long each of the column lengths are for each level of heat strain and each outcome measure
    column_lengths <- c()

    # There will be the number of major indicees equal to the numberer of heat strain levels.
    # So for each heating stage (i.e., major index in the list), and each outcome measure, add the column lengths to the running vector
    for(stages in 1: length(full_list)) {

      for (measure in 1: length(full_list[[stages]])) {

        column_lengths <- c(column_lengths, nrow(full_list[[stages]][[measure]]))
      }

    }

    ## Repeat the names of the variables the number of times required (i.e., length of the columns) for all of the baseline

    # The column lengths are a list of equal to the number of variables * number of levels of heat strain
    # So it will repeat the name of the ith vector equal to the length of the ith vector of the column lengths

    names_df_BL <- c()
    for(i in 1: length(variable_names)){
      names_df_BL <- c(names_df_BL, rep(variable_names[i], times = column_lengths[i]))

    }

    no_heat_within <- c()
    no_heat_between <- c()

    # Repeat the names of the variables the number of times required (i.e., length of the columns) for all of the low heat
    names_df_LH <- c()
    for(i in 1: length(variable_names)){
      names_df_LH <- c(names_df_LH, rep(variable_names[i],
                                        times = column_lengths[length(variable_names)+i]))

    }

    low_heat_within <- c()
    low_heat_between <- c()

    # Repeat the names of the variables the number of times required (i.e., length of the columns) for all of the medium heat
    # Because the column_lengths vector is the combination of LH and MH, need to start the index at the level corresponding to the second half
    names_df_MH <- c()
    for(i in 1: length(variable_names)){
      names_df_MH <- c(names_df_MH, rep(variable_names[i],
                                        times = column_lengths[c(length(variable_names) +
                                                                   length(variable_names) + i)]))

    }

    med_heat_within <- c()
    med_heat_between <- c()


    ##------------------------------##

    ## Now need to add the rows of the data together dependent on the function input
    # Few combinations in manuscripts, so again will only provide those options
    for (i in 1: length(variable_names)){

      no_heat_within <- rbind(no_heat_within, full_list[[1]][[i]][,c(1:2)])
      no_heat_between <- rbind(no_heat_between, full_list[[1]][[i]][,c(3:4)])
      low_heat_within <- rbind(low_heat_within, full_list[[2]][[i]][,c(1:2)])
      low_heat_between <- rbind(low_heat_between, full_list[[2]][[i]][,c(3:4)])
      med_heat_within <- rbind(med_heat_within, full_list[[3]][[i]][,c(1:2)])
      med_heat_between <- rbind(med_heat_between, full_list[[3]][[i]][,c(3:4)])

    }

    # repeated the number of times equal to 1: number of variables (BL), number of variables + 1 to number of variables * 2 (LH), etc...
    data_frame_return <- data.frame(stages = c(
      rep(c("BL"), times = sum(column_lengths[1:length(variable_names)])),
      rep(c("LH"), times = sum(column_lengths[(length(variable_names) + 1):(length(variable_names) * 2)])),
      rep(c("MH"), times = sum(column_lengths[((length(variable_names) * 2) + 1): (length(variable_names) * 3)]))),
      outcome = c(names_df_BL, names_df_LH, names_df_MH),
      within = rbind(no_heat_within, low_heat_within, med_heat_within),
      between = rbind(no_heat_between, low_heat_between, med_heat_between))


  }


  ##------------------------------##
  ##------------------------------##
  ## NEXT OPTION
  ##------------------------------##
  ##------------------------------##

  # If there are just 2 levels (i.e., no low or moderate) - cardiovascular squats
  if(NHS == "Y" & LHS == "N" & MHS == "N" & HHS == "Y") {
    # Generate new lists dependent on the inputs
    list_sample_sizes_BL <- list()
    list_sample_sizes_HH <- list()

    # For each variable, bind together the within subject and between subject elements
    # Need to have separate lists for each level of heat strain because they can be separate lengths
    for (i in 1: length(variable_names)){

      list_sample_sizes_BL[[i]] <- cbind(
        df[[i]][[2]],
        df[[i]][[3]])

      list_sample_sizes_HH[[i]] <- cbind(
        df[[i]][[5]],
        df[[i]][[6]])


    }

    # Combine the lists
    full_list <- list(list_sample_sizes_BL,
                      list_sample_sizes_HH)


    ##------------------------------##

    ###
    # Now need to extract a dataframe with just the sample sizes and effect size inputs
    ###

    # Need to know how long each of the column lengths are for each level of heat strain and each outcome measure
    column_lengths <- c()

    # There will be the number of major indicees equal to the numberer of heat strain levels.
    # So for each heating stage (i.e., major index in the list), and each outcome measure, add the column lengths to the running vector
    for(stages in 1: length(full_list)) {

      for (measure in 1: length(full_list[[stages]])) {

        column_lengths <- c(column_lengths, nrow(full_list[[stages]][[measure]]))
      }

    }

    ## Repeat the names of the variables the number of times required (i.e., length of the columns) for all of the baseline

    # The column lengths are a list of equal to the number of variables * number of levels of heat strain
    # So it will repeat the name of the ith vector equal to the length of the ith vector of the column lengths

    names_df_BL <- c()
    for(i in 1: length(variable_names)){
      names_df_BL <- c(names_df_BL, rep(variable_names[i], times = column_lengths[i]))

    }

    no_heat_within <- c()
    no_heat_between <- c()

    # Repeat the names of the variables the number of times required (i.e., length of the columns) for all of the low heat
    names_df_HH <- c()
    for(i in 1: length(variable_names)){
      names_df_HH <- c(names_df_HH, rep(variable_names[i],
                                        times = column_lengths[length(variable_names)+i]))

    }

    high_heat_within <- c()
    high_heat_between <- c()

    ##------------------------------##

    ## Now need to add the rows of the data together dependent on the function input
    # Few combinations in manuscripts, so again will only provide those options
    for (i in 1: length(variable_names)){

      no_heat_within <- rbind(no_heat_within, full_list[[1]][[i]][,c(1:2)])
      no_heat_between <- rbind(no_heat_between, full_list[[1]][[i]][,c(3:4)])
      high_heat_within <- rbind(high_heat_within, full_list[[2]][[i]][,c(1:2)])
      high_heat_between <- rbind(high_heat_between, full_list[[2]][[i]][,c(3:4)])

    }

    # repeated the number of times equal to 1: number of variables (BL), number of variables + 1 to number of variables * 2 (LH), etc...
    data_frame_return <- data.frame(stages = c(
      rep(c("BL"), times = sum(column_lengths[1:length(variable_names)])),
      rep(c("HH"), times = sum(column_lengths[(length(variable_names) + 1):(length(variable_names) * 2)]))),
      outcome = c(names_df_BL, names_df_HH),
      within = rbind(no_heat_within, high_heat_within),
      between = rbind(no_heat_between, high_heat_between))

  }


  ##------------------------------##
  ##------------------------------##
  ## NEXT OPTION
  ##------------------------------##
  ##------------------------------##


  # If there is just 1 level - sudomotor
  if(NHS == "Y" & LHS == "N" & MHS == "N" & HHS == "N") {
    # Generate new lists dependent on the inputs
    list_sample_sizes_BL <- list()

    # For each variable, bind together the within subject and between subject elements
    # Need to have separate lists for each level of heat strain because they can be separate lengths
    for (i in 1: length(variable_names)){

      list_sample_sizes_BL[[i]] <- cbind(
        df[[i]][[2]],
        df[[i]][[3]])


    }

    # Combine the lists
    full_list <- list(list_sample_sizes_BL)


    ##------------------------------##

    ###
    # Now need to extract a dataframe with just the sample sizes and effect size inputs
    ###

    # Need to know how long each of the column lengths are for each level of heat strain and each outcome measure
    column_lengths <- c()

    # There will be the number of major indicees equal to the numberer of heat strain levels.
    # So for each heating stage (i.e., major index in the list), and each outcome measure, add the column lengths to the running vector
    for(stages in 1: length(full_list)) {

      for (measure in 1: length(full_list[[stages]])) {

        column_lengths <- c(column_lengths, nrow(full_list[[stages]][[measure]]))
      }

    }

    ## Repeat the names of the variables the number of times required (i.e., length of the columns) for all of the baseline

    # The column lengths are a list of equal to the number of variables * number of levels of heat strain
    # So it will repeat the name of the ith vector equal to the length of the ith vector of the column lengths

    names_df_BL <- c()
    for(i in 1: length(variable_names)){
      names_df_BL <- c(names_df_BL, rep(variable_names[i], times = column_lengths[i]))

    }

    no_heat_within <- c()
    no_heat_between <- c()

    ##------------------------------##

    ## Now need to add the rows of the data together dependent on the function input
    # Few combinations in manuscripts, so again will only provide those options
    for (i in 1: length(variable_names)){

      no_heat_within <- rbind(no_heat_within, full_list[[1]][[i]][,c(1:2)])
      no_heat_between <- rbind(no_heat_between, full_list[[1]][[i]][,c(3:4)])

    }

    # repeated the number of times equal to 1: number of variables (BL), number of variables + 1 to number of variables * 2 (LH), etc...
    data_frame_return <- data.frame(stages = c(
      rep(c("BL"), times = sum(column_lengths[1:length(variable_names)]))),
      outcome = c(names_df_BL),
      within = rbind(no_heat_within),
      between = rbind(no_heat_between))

  }


  ##------------------------------##
  ##------------------------------##

  names(data_frame_return) <- c("stages", "outcome", "within_sample_size", "effect_size1", "between_sample_size", "effect_size2")

  return(data_frame_return)

}
