#' 
#' linelist_i <- structure(list(id = "HH_0001_02", house_id = "HH_0001", status = 1,age = 60, sex = 2, onset_date = 72.05, incub_period = 3.25,infec_date = 68.8, t = 68.8, inf_mean=3,inf_var=1,cur_inf = structure(2L, levels = c("0","1"), class = "factor"), cur_symp = structure(1L, levels = c("0", "1"), class = "factor"), onset_inf = structure(2L, levels = c("0","1"), class = "factor"), onset_symp = structure(1L, levels = c("0", "1"), class = "factor"), id_num = 478L), row.names = c(NA, -1L), class = c("tbl_df", "tbl", "data.frame"))
#' design=list(test_pos_height=1,specificity_test=1,sensitivity_symptoms=0.99,specificity_symptoms=0.95)
#' probability_recruit(linelist_i, symptom_based=TRUE, test_based=TRUE, design=list(test_pos_height=1))
#' 
probability_recruit <- function(linelist_i, symptom_based=TRUE, test_based=TRUE, design=list(test_pos_height=1,specificity_test=1,sensitivity_symptoms=1,specificity_symptoms=1,probability_enrollment=1)){
  ## Check if symptomatic
  is_symp <- linelist_i$cur_symp == 1
  ## Check TSI and if positive given TSI
  tsi <- linelist_i$t - linelist_i$infect_date_start
  is_inf <- ifelse(!is.na(tsi) & tsi > 0, TRUE, FALSE)
  sensitivity_test <- probability_positive(tsi, linelist_i$inf_mean, linelist_i$inf_var, design)
  sensitivity_test <- ifelse(is.na(sensitivity_test), 0, sensitivity_test)
  ## Check if recruited
  if(symptom_based){
    prob_recruit <- is_symp*design$sensitivity_symptoms + (1-is_symp)*(1-design$specificity_symptoms)
    if(test_based){
      prob_recruit <- prob_recruit * (is_inf*sensitivity_test + (1-is_inf)*(1-design$specificity_test))
    } 
  } else {
    prob_recruit <- (is_inf*sensitivity_test + (1-is_inf)*(1-design$specificity_test))
  }
  return(prob_recruit*design$probability_enrollment)
}

probability_recruit_dt <- function(linelist, symptom_based=TRUE, test_based=TRUE, design=list(test_pos_height=1,specificity_test=1,sensitivity_symptoms=1,specificity_symptoms=1,probability_enrollment=1)){
  if(symptom_based == TRUE & test_based == TRUE){
    linelist[, prob_recruit := (cur_symp*design$sensitivity_symptoms + (1-cur_symp)*(1-design$specificity_symptoms)) * 
                (cur_inf*prob_test_pos + (1-cur_inf)*(1-design$specificity_test)) * 
                design$probability_enrollment]
  } else if (symptom_based == FALSE & test_based == FALSE){
    linelist[, prob_recruit := design$probability_enrollment]
  } else if (symptom_based == FALSE & test_based == TRUE){
    linelist[, prob_recruit := (cur_inf*prob_test_pos + (1-cur_inf)*(1-design$specificity_test)) * design$probability_enrollment]
  } else if (symptom_based == TRUE & test_based == FALSE){
    linelist[, prob_recruit := (cur_symp*design$sensitivity_symptoms + (1-cur_symp)*(1-design$specificity_symptoms)) * design$probability_enrollment]
  } else {
    linelist[, prob_recruit := design$probability_enrollment]
  } 
  return(linelist)
}



probability_positive <- function(tsi, shape, scale, design=list(test_pos_height=1)){
  #browser()
  #pars <- gamma_mean_var_to_shape_scale(mean_par=inf_mean,var_par=inf_var)
  ## Probability of being positive given TSI
  prob_pos <- design$test_pos_height*dgamma(tsi, shape=shape, scale=scale)
  return(prob_pos)
}

expand_linelist <- function(linelist, dt, tmax,round_dates=FALSE){
  timeframe <- seq(0, tmax, by=dt)
  
  ## Round dates to nearest dt
  if(round_dates){
  linelist$onset_date <- round(linelist$onset_date / dt)*dt
  linelist$infec_date <- round(linelist$infec_date / dt)*dt
  linelist$incub_period <- round(linelist$incub_period / dt)*dt
  }
  expanded_linelist <- expand_grid(linelist, t=timeframe)
  
  expanded_linelist <- expanded_linelist %>%
    mutate(cur_inf = (ifelse(!is.na(infec_date) & t >= infec_date, 1, 0)),
           cur_symp = (ifelse(!is.na(onset_date) & t >= onset_date, 1, 0)),
           onset_inf = (ifelse(!is.na(infec_date) & t == infec_date, 1, 0)),
           onset_symp = (ifelse(!is.na(onset_date) & t == onset_date, 1, 0))
    )
  return(expanded_linelist)
}

assign_numeric_indiv_key <- function(linelist){
  indiv_key <- linelist %>% select(id, infec_date) %>% 
    arrange(infec_date) %>% distinct() %>% 
    mutate(id_num = 1:n())
  linelist <- linelist %>% left_join(indiv_key)
  return(linelist)
}

assign_numeric_hh_key <- function(linelist){
  indiv_key <- linelist %>% select(house_id) %>% 
    distinct() %>% 
    mutate(house_id_num = 1:n())
  linelist <- linelist %>% left_join(indiv_key)
  return(linelist)
}

resolve_coprimaries <- function(household){
  
}

gamma_mean_var_to_shape_scale <- function(mean_par, var_par){
  scale <- var_par/mean_par
  shape <- mean_par/scale
  return(c(shape=shape,scale=scale))
}
