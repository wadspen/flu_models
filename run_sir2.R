library(cdcfluview)
library(ggplot2)
library(dplyr)
# library(rstan)
library(cmdstanr)
library(rstanarm)
library(boot)
library(lubridate)
library(forecast)
library(tidyr)
source('~/Documents/flu_forecast/flu_models/sir2_ili_mod.R')
stan_mod <- stan_model(model_code=sir2_ili_mod)

# source('~/Documents/flu_forecast/flu_forecast_competition/get_data_functions.R')


ILINet_state <- ilinet(region = 'state')
ILINet_us <- ilinet(region = 'national') %>% mutate(region = 'US')
ILINet <- rbind(ILINet_state,ILINet_us)


ILINet <- ILINet %>% 
  mutate(year=year(week_start), week2=week(week_start),location_name=region)




ILINet <- ILINet %>% 
  filter(week %in% c(40:53,1:30)) %>% 
  mutate(season=factor(ifelse(week %in% 40:53,year,year-1))) %>% 
  mutate(season_week = ifelse(week %in% 40:53, week-39,
                              ifelse(year %in% c(2021,2023),week+12,week+13)))


ILINet_samp <- ILINet %>% 
  filter(season %in% 2015:2018, region == 'Georgia') %>% 
  group_by(season) %>% 
  mutate(n_weeks = length(season_week))

write.csv(ILINet_samp, '~/Documents/flu_forecast/flu_models/ILINet_samp.csv')

dat = list(N = nrow(ILINet_samp),
           n_seasons = length(unique(ILINet_samp$season)),
           seasons = as.numeric(factor(ILINet_samp$season)),
           ili = ILINet_samp$unweighted_ili/100,
           x = ILINet_samp$season_week,
          # n_weeks = ILINet_samp$n_weeks,
           n_weeks = 1:43,
           ps = max(as.numeric(factor(unique(ILINet_samp$season)))),
           S0 = .9,
           t0 = 0)
# init <- list(list(mu=21))
game_on = sampling(stan_mod, dat, chains=1,iter=5000)


