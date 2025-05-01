library(tidyverse)
setwd("~/Documents/GitHub/household-study-data/")

# Read in the linelist
linelist <- read_csv("data/dummy-output.csv")
head(linelist)

dt <- 0.05
tmax <- max(linelist$onset_date,na.rm=TRUE)
timeframe <- seq(0, tmax, by=dt)

## Round dates to nearest dt
linelist$onset_date <- round(linelist$onset_date / dt)*dt
linelist$infec_date <- round(linelist$infec_date / dt)*dt
linelist$incub_period <- round(linelist$incub_period / dt)*dt

expanded_linelist <- expand_grid(linelist, t=timeframe)

indiv_key <- expanded_linelist %>% select(id, infec_date) %>% arrange(infec_date) %>% distinct() %>% mutate(id_num = 1:n())

expanded_linelist <- expanded_linelist %>%
  mutate(cur_inf = as.factor(ifelse(!is.na(infec_date) & t >= infec_date, 1, 0)),
         cur_symp = as.factor(ifelse(!is.na(onset_date) & t >= onset_date, 1, 0)),
         onset_inf = as.factor(ifelse(!is.na(infec_date) & t == infec_date, 1, 0)),
         onset_symp = as.factor(ifelse(!is.na(onset_date) & t == onset_date, 1, 0))
         )
expanded_linelist <- expanded_linelist %>% left_join(indiv_key)
linelist <- linelist %>% left_join(indiv_key)

p1 <- ggplot(expanded_linelist) + 
  geom_tile(aes(x=t,y=id_num, fill=cur_inf)) +
  geom_point(data=linelist,aes(x=onset_date,y=id_num),col="red",shape=8) +
  theme_minimal() +
  ylab("Individual ID") +
  xlab((paste0("Time (per ",dt, " days)"))) +
  scale_fill_viridis_d()

ggsave("figures/linelist_expansion.png", p1, width=5, height=10, dpi=300)
