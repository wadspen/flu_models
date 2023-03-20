

library(cdcfluview)
library(ggplot2)
library(dplyr)
library(rstan)
library(rstanarm)
library(boot)
library(lubridate)
library(forecast)
library(tidyr)
source('~/Documents/flu_forecast/flu_models/sir_ili_mod.R')
stan_mod <- stan_model(model_code=sir_ili_mod)

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
  mutate(n_weeks = max(season_week))








dat = list(N = nrow(ILINet_samp),
           n_seasons = length(unique(ILINet_samp$season)),
           seasons = as.numeric(factor(ILINet_samp$season)),
           ili = ILINet_samp$unweighted_ili/100,
           x = ILINet_samp$season_week,
           ps = max(as.numeric(factor(unique(ILINet_samp$season)))),
           last_week = 18)
# init <- list(list(mu=21))
rmle = sampling(stan_mod, dat, chains=1,iter=5000)
traceplot(rmle)
# saveRDS(rmle,file = '~/Documents/flu_forecast/flu_models/sirtest.rds')
post <- rstan::extract(rmle)
postpr <- post$ili_pred
ints <- as.data.frame(predictive_interval(postpr,prob = .5))
colnames(ints) <- c('lower','upper')
ints$season_week = 1:22


ILINet_samp %>% 
  filter(season == 2018) %>% 
  ggplot() +
  geom_point(aes(x=season_week, y=unweighted_ili/100)) +
  geom_line(data=ints,aes(x=season_week,y=upper))

test <- ILINet_samp2 %>% 
  group_by(season) %>% 
  mutate(n_weeks = max(season_week))
  group_by(season) %>% 
  mutate()

rm <- readRDS('~/Documents/flu_forecast/flu_models/sirtest.rds')
post <- rstan::extract(rm)
postpr <- post$ili_pred
ints <- predictive_interval(postpr)



beta <- 95.12
rho <- 94.12
gamma <- beta*rho
I0 <- 0.05

S0 <- .9
n <- 44

r <- RK4SIR(n, beta, gamma, S0, I0)
apply(r[,-1],MARGIN=1,FUN=sum)
plot(r$I~c(0:200))
















