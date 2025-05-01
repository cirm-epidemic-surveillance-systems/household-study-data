#' 
#' linelist_i <- structure(list(id = "HH_0001_02", house_id = "HH_0001", status = 1,age = 60, sex = 2, onset_date = 72.05, incub_period = 3.25,infec_date = 68.8, t = 68.8, inf_mean=3,inf_var=1,cur_inf = structure(2L, levels = c("0","1"), class = "factor"), cur_symp = structure(1L, levels = c("0", "1"), class = "factor"), onset_inf = structure(2L, levels = c("0","1"), class = "factor"), onset_symp = structure(1L, levels = c("0", "1"), class = "factor"), id_num = 478L), row.names = c(NA, -1L), class = c("tbl_df", "tbl", "data.frame"))
#' design=list(test_pos_height=1,specificity_test=1,sensitivity_symptoms=0.99,specificity_symptoms=0.95)
#' probability_recruit(linelist_i, symptom_based=TRUE, test_based=TRUE, design=list(test_pos_height=1))
#' 

probability_recruit <- function(linelist_i, symptom_based=TRUE, test_based=TRUE, design=list(test_pos_height=1,specificity_test=1,sensitivity_symptoms=1,specificity_symptoms=1,probability_enrollment=1)){
  ## Check if symptomatic
  is_symp <- linelist_i$cur_symp == 1
  
  ## Check TSI and if positive given TSI
  tsi <- linelist_i$t - linelist_i$infec_date
  is_inf <- tsi > 0
  sensitivity_test <- probability_positive(tsi, linelist_i$inf_mean, linelist_i$inf_var, design)
  
  ## Check if recruited
  if(symptom_based){
    prob_recruit <- is_symp*design$sensitivity_symptoms + (1-is_symp)*(1-design$specificity_symptoms)
    if(test_based){
      prob_recruit <- prob_recruit * (is_inf*sensitivity_test + (1-is_inf)*(1-design$specificity_test))
    } 
  } else {
    prob_recruit <- (is_inf*sensitivity_test + (1-is_inf)*(1-design$specificity_test))
  }
  return(prob_recruit*probability_enrollment)
}


probability_positive <- function(tsi, inf_mean, inf_var, design=list(test_pos_height=1)){
  pars <- jahR::convert_dist_pars(dist="gamma",mean_par=inf_mean,var_par=inf_var)
  ## Probability of being positive given TSI
  prob_pos <- design$test_pos_height*dgamma(tsi, shape=pars$pars["shape"], scale=pars$pars["scale"])
  return(prob_pos)
}

resolve_coprimaries <- function(household){
  
}

