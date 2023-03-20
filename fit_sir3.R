library(outbreaks)
library(tidyverse)
library(ggplot2)
# library(rstan)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(rstanarm)
library(cdcfluview)
library(lubridate)
# source('~/Documents/flu_forecast/flu_models/sir4_ili_mod.R')


niter <- 2000


########################
#######ILI data#########
########################

ILINet_state <- ilinet(region = 'state')
ILINet_us <- ilinet(region = 'national') %>% mutate(region = 'US')
ILINet <- rbind(ILINet_state,ILINet_us)


ILINet <- ILINet %>% 
  mutate(year=year(week_start), week2=week(week_start),location_name=region)

years53 <- ILINet %>% 
  group_by(year) %>% 
  summarise(
    m = max(week2)
  ) %>% 
  filter(m > 52) %>% 
  select(year) %>% 
  as.vector()
years53 <- years53$year


ILINet <- ILINet %>% 
  filter(week2 %in% c(40:53,1:30)) %>% 
  mutate(season=factor(ifelse(week2 %in% 40:53,year,year-1))) %>% 
  mutate(season_week = ifelse(week2 %in% 40:53, week2-39,
                              ifelse(season %in% years53, 
                                     week2 + 14, week2 + 13)))
  


ILINet_samp <- ILINet %>% 
  filter(season %in% c(2012:2019,2022), region == 'US') %>% 
  group_by(season) %>%
  mutate(n_weeks = max(season_week))

ILINet_samp <- ILINet_samp[
  with(ILINet_samp, order(season, season_week)),
]

# test_df <- ILINet_samp[,c('year','week','week2','week_start','season','season_week','test_week')]
# table(test_df$season_week)
# which(test_df$season_week == 14)
# View(test_df)
# test <- ILINet_samp %>% 
#   group_by(season,season_week) %>% 
#   summarise(
#     n = n()
#   )

ILINet_samp %>% 
  group_by(season) %>% 
  summarise(
    ma = max(season_week)
  )

ILINet_samp %>% 
  ggplot() +
  geom_point(aes(x = season_week, y = unweighted_ili,colour=season))


week_inds <- ILINet_samp %>% 
  group_by(season) %>% 
  summarise(
    max_week = max(season_week)
  )
is_last <- function(x) !duplicated(x, fromLast = TRUE)
week_inds$index <- which(is_last(ILINet_samp$season))
dat_inds <- week_inds %>% 
  mutate(for_ind = index - max_week + 1)




# ILINet_samp <- ILINet_samp %>% 
#   filter(season_week <= 43)

data_sir <- list(M = nrow(ILINet_samp),
                 n_seasons = length(unique(ILINet_samp$season)),
                 n_weeks = 44,
                 seasons = as.numeric(factor(unique(ILINet_samp$season))),
                 seg_ind_start = dat_inds$for_ind,
                 seg_ind_length = dat_inds$max_week,
                 seg_ind_max = dat_inds$index,
                 weeks = ILINet_samp$season_week,
                 S0 = .9, t0 = 0, 
                 ts = ILINet_samp$season_week, 
                 N = 1, 
                 ili = ILINet_samp$unweighted_ili/100,
                 ps = max(as.numeric(factor(unique(ILINet_samp$season)))))

# number of MCMC steps
niter <- 2000

stan_mod <- cmdstanr::cmdstan_model(stan_file = 
              '~/Documents/flu_forecast/flu_models/sir6_ili_mod.stan')


samples <- stan_mod$sample(data = data_sir,
                    chains = 1, 
                    seed = 0)

# mcmc_trace(draws, pars = 'pred_ili[23]')


# samples$summary()
draws <- samples$draws(format = 'df')
mcmc_trace(draws, pars = 'beta[8]')
# draws$pred_ili
# traceplot(samples)
postpr <- draws[,1586:1608]
pi <- as.data.frame(predictive_interval(as.matrix(postpr), prob=.95))
colnames(pi) <- c('lower','upper')
pi$season_week <- 1:23


ILINet_samp %>% 
  filter(season %in% 2022) %>% 
  ggplot() +
  geom_point(aes(x=season_week, y=unweighted_ili/100,colour=season)) +
  geom_line(data=pi,aes(x=season_week,y=lower),colour='orange') +
  geom_line(data=pi,aes(x=season_week,y=upper),colour='orange')

post <- rstan::extract(samples,pars='pred_cases')
postpr <- post$pred_cases

pred_ints <- as.data.frame(predictive_interval(postpr,prob=.95))
colnames(pred_ints) <- c('lower','upper')
pred_ints$season_week <- sort(ILINet_samp$season_week)

ILINet_samp %>% 
  ggplot() +
  geom_point(aes(x=season_week, y=unweighted_ili/100)) +
  geom_line(data = pred_ints, aes(x=season_week, y=upper), colour='orange') +
  geom_line(data = pred_ints, aes(x=season_week, y= lower), colour='orange')




smr_pred <- cbind(as.data.frame(summary(
  samples, pars = "pred_cases", probs = c(0.05, 0.5, 0.95))$summary), t, cases)
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names

ggplot(smr_pred, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = 'orange', alpha = 0.35) +
  geom_line(mapping = aes(x = t, y = X50.), color = 'orange') + 
  geom_point(mapping = aes(y = cases)) +
  labs(x = "Day", y = "Number of students in bed")
traceplot(samples)




sort(ILINet_samp$season_week)

week_inds <- ILINet_samp %>% 
  group_by(season) %>% 
  summarise(
    max_week = max(season_week)
  )
is_last <- function(x) !duplicated(x, fromLast = TRUE)
week_inds$index <- which(is_last(ILINet_samp$season))
dat_inds <- week_inds %>% 
  mutate(for_ind = index - max_week + 1)

data_sir <- list(n_weeks = nrow(ILINet_samp), 
                 S0 = .9, t0 = t0, 
                 ts = sort(unique(ILINet_samp$season_week)), 
                 N = 1, 
                 cases = ILINet_samp$unweighted_ili/100)



which(is_last(ILINet_samp$season))


# Chain 1 0.252544 
# Chain 1 1.22789 
# Chain 1 2.20871 

I <- .252544
mu <- -1.22789 

I*exp(mu)/(1 + (exp(mu) - 1)*I)
exp(mu)/(exp(mu) + (-1 -I))


