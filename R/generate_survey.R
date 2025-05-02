library(tidyverse)
library(data.table)
setwd("~/Documents/GitHub/household-study-data/")
source("R/recruitment_functions.R")
# Read in the linelist
linelist <- read_csv("simulation/outputs/simulated_infections.csv")
head(linelist)

## Simulation design decisions
dt <- 1 ## Time step of simulation (1 = 1 day)
## Decide on testing cadence and delay from recruitment to first test
test_interval = 3 ## Frequency of testing in days 
test_delay = 2 ## Delay from first detection to enrollment in days
test_max = 30 ## Maximum number of testing days following enrollment AND testing beginning

## Set timeframe for simulation
tmax <- max(linelist$onset_date,na.rm=TRUE)
timeframe <- seq(0, tmax,by=dt)

## Test characteristics
test_based <- TRUE ## Must test positive to trigger recruitment
symptom_based <- TRUE ## Must be symptomatic to trigger recruitment
test_pos_var_scale <- 5
max_prob_pos <- 2

inf_mean <- mean(linelist$ip_mean,na.rm=TRUE)
inf_var <- mean(linelist$ip_var,na.rm=TRUE)
scale=inf_var/inf_mean
shape=inf_mean/scale
inf_mode <- if (shape > 1) (shape - 1) * scale else 0

shape_new <- find_alpha_new(inf_mode, inf_var*test_pos_var_scale)
scale_new <- inf_mode / (shape_new - 1)

study_design <- list(test_pos_height=max_prob_pos/max(dgamma(seq(0,30,by=0.1),shape=shape_new,scale=scale_new)), ## Scaling constant to give maximum probability of being positive (here 0.9)
                     specificity_test=1, ## Prob false positive test
                     sensitivity_symptoms=1, ## Prob of recording someone as symptomatic, given symptomatic
                     specificity_symptoms=1, ## Prob of reporting someone as symptomatic, given asymptomatic
                     probability_enrollment=1 ## Probability of enrollment, given all other criteria met
                     )

## Look at probability of testing positive over time-since-infection
tsi <- seq(0,30,by=0.1)
prob_pos_dat <- data.frame(x=tsi,y=probability_positive(tsi,shape=shape_new,scale=scale_new,
                                                         design=study_design))
study_design1 <- study_design
study_design1$test_pos_height <- 1
prob_inf_dat <- data.frame(x=tsi,y=probability_positive(tsi,shape=shape,scale=scale,
                                                        design=study_design1))
p_pos_tsi <- ggplot(data=prob_pos_dat) + 
  geom_line(aes(x=x,y=y,col="Probability positive")) +
  geom_line(data=prob_inf_dat,aes(x=x,y=y,col="Probability infectious")) +
  theme_minimal() +
  ylab("Probability of being positive") +
  xlab("Time since infection (days)") +
  scale_color_manual(name="",values=c("Probability infectious"="red","Probability positive"="blue")) +
  theme(legend.position="bottom")
ggsave("figures/probability_positive_tsi.png",plot=p_pos_tsi,width=7,height=4,dpi=300)

## Expand linelist to show individual status at all timepoints
linelist <- linelist %>%
  ungroup() %>%
  rowwise() %>%
  mutate(scale=inf_var/inf_mean,
         shape=inf_mean/scale,
         inf_mode=if_else(shape > 1, (shape - 1) * scale, 0),
         shape_new=find_alpha_new(inf_mode,inf_var*10),
         scale_new=inf_mode / (shape_new - 1)) %>%
  rename(shape_old=shape,scale_old=scale,shape=shape_new,scale=scale_new)

expanded_linelist <- expand_linelist(linelist, dt, tmax)

## Assign new unique, numeric IDs to individuals and households
expanded_linelist <- assign_numeric_indiv_key(expanded_linelist)
expanded_linelist <- assign_numeric_hh_key(expanded_linelist)
linelist <- assign_numeric_indiv_key(linelist)

## Rename some variables and give "infection end" and "onset end"
expanded_linelist <- expanded_linelist %>% 
  rename(infect_date_start=infec_date) %>%
  ## Placeholder for now
  mutate(infect_date_end=infect_date_start+7,
         onset_date_start=onset_date,
         onset_date_end=onset_date_start+7)
data_hh <- as.data.table(expanded_linelist)


## Calculate TSI and gamma parameters for infectivity profile for all individuals
data_hh[,tsi := t - infect_date_start]

## Calculate probability of testing positive for all TIS
data_hh <- data_hh %>% mutate(prob_test_pos = probability_positive(tsi,shape,scale, design=study_design))
data_hh[is.na(prob_test_pos), prob_test_pos:=0]

## Calculate probability of recruitment for all times
data_hh <- probability_recruit_dt(data_hh,test_based = test_based,symptom_based = symptom_based)

## Plots of individual statuses over time
p1 <- ggplot(expanded_linelist) + 
  geom_tile(aes(x=t,y=id_num, fill=as.factor(cur_inf))) +
  geom_point(data=linelist,aes(x=onset_date,y=id_num),col="red",shape=8) +
  theme_minimal() +
  ylab("Individual ID") +
  xlab((paste0("Time (per ",dt, " days)"))) +
  scale_fill_viridis_d()+
  theme(legend.position="bottom")

p_pos <- ggplot(data_hh) + 
  geom_tile(aes(x=t,y=id_num, fill=prob_test_pos)) +
  #geom_point(data=linelist,aes(x=onset_date,y=id_num),col="red",shape=8) +
  theme_minimal() +
  ylab("Individual ID") +
  xlab((paste0("Time (per ",dt, " days)"))) +
  scale_fill_viridis_c()+
  theme(legend.position="bottom")

## Look at probability positive over time and probability of getting recruited
p_recruit <- ggplot(data_hh) + 
  geom_tile(aes(x=t,y=id_num, fill=prob_recruit)) +
  #geom_point(data=linelist,aes(x=onset_date,y=id_num),col="red",shape=8) +
  theme_minimal() +
  ylab("Individual ID") +
  xlab((paste0("Time (per ",dt, " days)"))) +
  scale_fill_viridis_c() +
  theme(legend.position="bottom")
p_all <- p1 | p_pos | p_recruit
ggsave("figures/linelist_expansion.png", p_all, width=10, height=10, dpi=300)

########################################
# 1. Simulate recruitment process
########################################
## For each time period, check if an individual is recruited, and if so, record
## them as an index case for that household
data_hh[, index_hh:='0']
for(t in seq_along(timeframe)){
  t1 = timeframe[t]
  print(paste0("Timestep: ", t1))
  # identify which households have not yet to be recruited
  test_id = unique(data_hh[index_hh==0 & t==t1, id_num])
  
  ## Calculate recruitment probabilities of each unrecruited individual
  recruit_probs <- data_hh[id_num %in% test_id & t == t1,prob_recruit]
  ## Simulate if individual is recruited
  recruit_id = test_id[which(runif(length(test_id)) < recruit_probs)]

  ## If we have recruited any individuals
  if(length(recruit_id) > 0){
    
    ## Flag the recruited individuals as index case for their household (own entry)
    data_hh[id_num %in% recruit_id, index_hh:=as.character(id_num)]
    ## Record the current time as the detection date
    data_hh[id_num %in% recruit_id, date_detected:=t1]
    
    ## For each recruited household, record the ID of all recruited individuals in that household, accounting for co-index cases    test_hh = unique(data_hh[id_num %in% recruit_id, house_id_num])
    index_id = lapply(1:length(test_hh), function(x){
      rep(paste(unique(data_hh[house_id_num == test_hh[x] & index_hh!='0',]$index_hh), collapse=','),
          data_hh[house_id_num == test_hh[x] & index_hh=='0',.N])
      
    })
    
    index_id = unlist(index_id)
    data_hh[house_id_num %in% test_hh & index_hh=='0', index_hh:=index_id]
    
    
    }
}

########################################
# Simulate household testing
########################################
test_intervals <- seq(1+test_delay, test_max+test_delay, by=test_interval)

## For each household, find the date of first detection, and then filter to only 
## keep dates where a test was performed, OR you are the index case and you were detected that day
data_hh <- data_hh %>% 
  group_by(house_id_num) %>% 
  mutate(date_detected = min(date_detected,na.rm=TRUE)) %>%
  filter(t >= date_detected) %>%
  ungroup() %>% 
  group_by(id_num) %>%
  filter(t %in% (test_intervals + date_detected[1]) | ## Check if time is in testing cadence
           ## Or it's the detection date, and the current ID is an index case
           (t == date_detected[1] & id_num %in% as.numeric(str_split(index_hh,",")[[1]])) )%>%
  ungroup()

## Simulate binary test outcomes
data_hh <- data_hh %>%
  mutate(test_outcome = if_else(t == date_detected,test_based,rbinom(n(),1,prob_test_pos)))
#data_hh$test_outcome <- runif(nrow(data_hh)) < data_hh$prob_test_pos
data_hh$time_since_enrollment <- data_hh$t - data_hh$date_detected

p_hh_tests <- ggplot(data_hh %>% filter(house_id %in% 1:25)) + geom_tile(aes(x=time_since_enrollment,y=id,fill=as.factor(test_outcome))) +
  facet_wrap(~house_id,scales="free_y") + theme_minimal()

p_all_tests <- ggplot(data_hh) + geom_tile(aes(x=time_since_enrollment,y=id,fill=as.factor(test_outcome))) + theme_minimal()

ggsave("figures/p_hh_tests.png",p_hh_tests,height=7,width=8,dpi=300)
ggsave("figures/p_all_tests.png",p_all_tests,height=10,width=6,dpi=300)

## Add date noising
data_hh_report <- data_hh %>% select(id,id_num,age,house_id, house_id_num, t, onset_date, test_outcome)
data_hh_report <- data_hh_report %>% mutate(onset_date = floor(onset_date),
                                            t=floor(t)) %>%
  mutate(is_symp = !is.na(onset_date))

## Remove onsets after testing period
data_hh_report <- data_hh_report %>% 
  group_by(house_id) %>%
  mutate(end_study_date=max(t)) %>%
  ungroup() %>%
  mutate(onset_date = if_else(onset_date > end_study_date,NA,onset_date))

## Save dataset
write_csv(data_hh_report, "data/simulated_perfect_survey.csv")
