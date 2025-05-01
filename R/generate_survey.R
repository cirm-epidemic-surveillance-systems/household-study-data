library(tidyverse)
library(data.table)
setwd("~/Documents/GitHub/household-study-data/")
source("R/recruitment_functions.R")
# Read in the linelist
linelist <- read_csv("data/dummy-output.csv")
head(linelist)

dt <- 1
tmax <- max(linelist$onset_date,na.rm=TRUE)
timeframe <- seq(0, tmax,by=dt)
expanded_linelist <- expand_linelist(linelist, dt, tmax)
expanded_linelist <- assign_numeric_indiv_key(expanded_linelist)
expanded_linelist <- assign_numeric_hh_key(expanded_linelist)
linelist <- assign_numeric_indiv_key(linelist)
# ggsave("figures/linelist_expansion.png", p1, width=5, height=10, dpi=300)


expanded_linelist <- expanded_linelist %>% 
  rename(infect_date_start=infec_date) %>%
  ## Placeholder for now
  mutate(infect_date_end=infect_date_start+7,
         onset_date_start=onset_date,
         onset_date_end=onset_date_start+7)
data_hh <- as.data.table(expanded_linelist)

study_design <- list(test_pos_height=7,specificity_test=1,sensitivity_symptoms=1,specificity_symptoms=1,probability_enrollment=1)

inf_mean <- 5
inf_var <- 10
scale=inf_var/inf_mean
shape=inf_mean/scale

ggplot(data.frame(x=seq(0,30,by=0.1),y=probability_positive(seq(0,30,by=0.1),shape=shape,scale=scale,design=study_design))) + 
  geom_line(aes(x=x,y=y)) +
  theme_minimal() +
  ylab("Probability of being positive") +
  xlab("Time since infection (days)") +
  scale_fill_viridis_c() +
  theme(legend.position="bottom")

data_hh[,inf_mean:=5]
data_hh[,inf_var:=10]


data_hh[,tsi := t - infect_date_start]
data_hh <- data_hh %>% mutate(scale=inf_var/inf_mean,
                              shape=inf_mean/scale)

data_hh <- data_hh %>% mutate(prob_test_pos = probability_positive(tsi,shape,scale, design=study_design))
#data_hh[,prob_test_pos := probability_positive(tsi, inf_mean, inf_var, design=study_design)]
data_hh[is.na(prob_test_pos), prob_test_pos:=0]
data_hh <- probability_recruit_dt(data_hh,test_based = TRUE,symptom_based = TRUE)

## Look at probability positive over time and probability of getting recruited

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
  geom_point(data=linelist,aes(x=onset_date,y=id_num),col="red",shape=8) +
  theme_minimal() +
  ylab("Individual ID") +
  xlab((paste0("Time (per ",dt, " days)"))) +
  scale_fill_viridis_c() +
  theme(legend.position="bottom")
p1 | p_pos | p_recruit

########################################
# simulate recruitment process
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
    
    # identify recruited household
    test_hh = unique(data_hh[id_num %in% recruit_id, house_id_num])
    index_id = lapply(1:length(test_hh), function(x){
      rep(paste(unique(data_hh[house_id_num == test_hh[x] & index_hh!='0',]$index_hh), collapse=','),
          data_hh[house_id_num == test_hh[x] & index_hh=='0',.N])
      
    })
    
    index_id = unlist(index_id)
    data_hh[house_id_num %in% test_hh & index_hh=='0', index_hh:=index_id]
    
    
    ## For each recruited household, record the ID of all recruited individuals in that household, accounting for co-index cases
    
    }
  }
data_hh_store <- data_hh
########################################
# Simulate household testing
########################################

## Decide on testing cadence and delay from recruitment to first test
test_interval = 3
test_delay = 2
test_max = 30

test_intervals <- seq(1+test_delay, test_max, by=test_interval)



data_hh <- data_hh %>% 
  group_by(house_id_num) %>% 
  mutate(date_detected = min(date_detected,na.rm=TRUE)) %>%
  filter(t >= date_detected) %>%
  ungroup() %>% 
  group_by(id_num) %>%
  filter(t %in% (test_intervals + date_detected) | (t == date_detected & id_num %in% as.numeric(str_split(index_hh,",")[[1]])) )%>%
  ungroup()

data_hh$test_outcome <- runif(nrow(data_hh)) < data_hh$prob_test_pos

ggplot(data_hh) + geom_tile(aes(x=t,y=id_num,fill=test_outcome))
