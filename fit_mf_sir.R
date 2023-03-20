library(cmdstanr)
library(cdcfluview)
library(lubridate)
library(ggplot2)
library(posterior)
library(bayesplot)
library(rstanarm)
library(stringr)
library(dplyr)
# mf_mod <- cmdstan_model(stan_file = 
#               '~/Documents/flu_forecast/flu_models/mf_sir_disc.stan')

raw_flu <- read.csv(
  'https://raw.githubusercontent.com/cdcepi/Flusight-forecast-data/master/data-truth/truth-Incident%20Hospitalizations.csv')

raw_flu <- raw_flu %>% 
  mutate(year = year(date),week = week(date)-1)

ILINet_state <- ilinet(region = 'state')
ILINet_us <- ilinet(region = 'national') %>% mutate(region = 'US')
ILINet <- rbind(ILINet_state,ILINet_us)

ILINet <- ILINet %>% 
  mutate(year=year(week_start), week2=week(week_start),location_name=region)

both_flu <- raw_flu %>% 
  left_join(ILINet, by = c('year','week','location_name'))


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

both_flu <- both_flu %>% 
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

table(ILINet_samp$season_week)
# ILINet_samp[,c('year','season','season_week','week','week2')] %>% View()

both_flu_samp <- both_flu %>% 
  filter(season == 2022, region == 'US') %>% 
  group_by(season) %>%
  mutate(n_weeks = max(season_week))

both_flu %>% 
  filter(season %in% 2022, location_name == 'US') %>% 
  ggplot(aes(x = season_week, y = unweighted_ili)) +
  # ggplot(aes(x=unweighted_ili/100, y=value)) + 
  geom_point() + 
  geom_smooth(method = 'lm') #+ 
  facet_wrap(~region_type)
  
both_flu_test <- both_flu %>% 
  filter(season == 2022, location_name == 'US', region_type == 'States')

lmmod <- lm(value~unweighted_ili, data=both_flu_samp)

both_flu_samp %>% 
  ggplot() +
  geom_point(aes(x = season_week, y = value), colour='blue') + 
  geom_point(aes(x = season_week, y = -173.18 + 131.32*unweighted_ili))

both_flu_samp %>% 
  ggplot() +
  geom_point(aes(x = season_week, y = unweighted_ili))

week_inds <- ILINet_samp %>% 
  group_by(season) %>% 
  summarise(
    max_week = max(season_week)
  )
is_last <- function(x) !duplicated(x, fromLast = TRUE)
week_inds$index <- which(is_last(ILINet_samp$season))
dat_inds <- week_inds %>% 
  mutate(for_ind = index - max_week + 1)


data_sir <- list(M = nrow(ILINet_samp),
                 n_seasons = length(unique(ILINet_samp$season)),
                 n_weeks = 44,
                 cur_yr_n_weeks = max(both_flu_samp$season_week),
                 seasons = as.numeric(factor(unique(ILINet_samp$season))),
                 seg_ind_start = dat_inds$for_ind,
                 seg_ind_length = dat_inds$max_week,
                 seg_ind_max = dat_inds$index,
                 weeks = ILINet_samp$season_week,
                 S0 = .9, t0 = 0, 
                 ts = ILINet_samp$season_week, 
                 N = 1, 
                 ili = ILINet_samp$unweighted_ili/100,
                 hosp = both_flu_samp$value,
                 ps = max(as.numeric(factor(unique(ILINet_samp$season)))))


mf_mod <- cmdstan_model(stan_file = 
                          '~/Documents/flu_forecast/flu_models/mf_sir2_disc.stan')

samples <- mf_mod$sample(data = data_sir,
                           chains = 1, 
                           seed = 0)

draws <- samples$draws(format = 'df')
mcmc_trace(draws, pars = 'alpha1')
postpr <- draws[,which(str_detect(colnames(draws),pattern = 'pred_hosp'))]
pi <- as.data.frame(predictive_interval(as.matrix(postpr), prob=.15))
colnames(pi) <- c('lower','upper')
pi$season_week <- 1:23

# pi$lower <- -10786.7 + 4492*(100*pi$lower)
# pi$upper <- -10786.7 + 4492*(100*pi$upper)


both_flu_samp %>% 
  ggplot() +
  # geom_line(aes(x=season_week, y=value/exp(unweighted_ili/100)))
  geom_point(aes(x = season_week, y = value),colour='red') +
  # geom_point(aes(x = season_week, y = unweighted_ili/100),colour='red') +
  geom_line(data = pi, aes(x = season_week, y = lower),colour = 'orange') +
  geom_line(data = pi, aes(x = season_week, y = upper),colour = 'orange')
  # geom_line(aes(x = season_week, y = 2000*unweighted_ili))





seg_ind_start = dat_inds$for_ind
seg_ind_length = dat_inds$max_week
# seg_ind_max = dat_inds$index
ts = ILINet_samp$season_week
for (i in 1:length(seg_ind_start)) {
  print(i)
  print(seg_ind_start[i])
  print(seg_ind_length[i])
  print(ts[seq(from = seg_ind_start[i],length.out = seg_ind_length[i])])
}








