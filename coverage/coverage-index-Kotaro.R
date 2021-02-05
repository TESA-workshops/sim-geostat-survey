library(dplyr)
library(ggplot2)
library(SimSurvey)
library(sdmTMB)
future::plan(future::multisession) # for the furrr package; or `plan(sequential)`

sim_and_fit <- function(iter, set_den, sd, phi_year, year_toremove, fraction_x, fraction_y) {
  set.seed(iter * 283028)
  sim <- sim_abundance(ages = seq(1, 10), years = seq(1, 10)) %>%
    sim_distribution(
      grid = make_grid(res = c(10, 10), depth_range = c(10, 500), strat_splits = 1, strat_breaks = seq(0, 1000, by = 100)),
      ays_covar = sim_ays_covar(phi_age = 0.8, phi_year = phi_year, range=600, sd=sd),
      depth_par = sim_parabola(mu = 200, sigma = 30)
    )

  ## Adding a few plotting functions to help us understand what we are doing
  # To look at the grid and how the strata look like (only one preset for depth arangemet)
  # a <- make_grid(res = c(10, 10), depth_range = c(10, 500), strat_splits = 1, strat_breaks = seq(10, 1000, by = 100))
  # plot_grid(a)
  # To plot fish distribution
  #plot_distribution(sim, ages=1:3, year=1:3)

  # For the survey simulation. To examine the effect of set_den (which in essence controls the number of set per strate)
  # but min_sets sets the minimum number so that it is not equal to zero by error
  #a  <- sim_survey(sim, n_sims = 1, set_den=set_den, min_sets = 2)$setdet
  #a %>% ggplot(aes(x, y)) + geom_point() + facet_wrap(~year)


  survey <- sim_survey(sim, n_sims = 1, set_den=set_den) %>% run_strat()
  xy <- as_tibble(survey$grid_xy)
  dat <- as_tibble(survey$setdet) %>%
    select(x, y, set, year, N = n, tow_area)
  dat <- left_join(dat, xy, by = c("x", "y"))
  dat$offset <- log(dat$tow_area)

  grid_dat <- as_tibble(select(sim$grid_xy, x, y, depth)) %>% distinct()
  grid_dat <- purrr::map_dfr(sort(unique(dat$year)), ~ bind_cols(grid_dat, year = .))
  grid_dat$offset <- mean(dat$offset)

  dat_filtered <- dplyr::filter(dat, !(x %in% c(min(dat$x), min(dat$x)+diff(range(dat$x))*fraction_x)  &
                                         y %in% c(min(dat$y), min(dat$y)+diff(range(dat$y))*fraction_y) &
                                         year %in% sample(1:10, year_toremove)))
  # ggplot(dat_filtered, aes(x, y)) + geom_point() + facet_wrap(~year)

  mesh <- sdmTMB::make_mesh(dat_filtered, xy_cols = c("x", "y"), cutoff = 20)
  #mesh <- sdmTMB::make_mesh(dat, xy_cols = c("x", "y"), cutoff = 20)
  fit <- sdmTMB(N ~ 0 + as.factor(year) + offset,
                data = dat_filtered,
                family = nbinom2(link = "log"), spde = mesh,
                ar1_fields = TRUE,
                include_spatial = TRUE, time = "year"
  )

  pred <- predict(fit, newdata = grid_dat, return_tmb_object = TRUE, area = 100)
  index <- get_index(pred)

  true_abund <- tibble(year = unique(dat$year), N = as.numeric(colSums(survey$I))) %>%
    mutate(type = "True")
  strat_abund <- tibble::as_tibble(survey$total_strat) %>%
    mutate(N = total, type = "Design-based") %>%
    select(year, N, type)
  mutate(index, type = "Model-based", N = est) %>%
    bind_rows(strat_abund) %>%
    bind_rows(true_abund) %>%
    mutate(iter = iter,
           set_den=set_den,
           sd = sd,
           phi_year=phi_year,
           year_toremove = year_toremove,
           fraction_x=fraction_x,
           fraction_y=fraction_y,
           scenario=paste(iter,set_den,sd,phi_year,year_toremove,fraction_x,fraction_y,sep="_"))
}


### Scenario of interest that we came up with
  # 1. Vary effort level
  # set_den=1/10000, 1/1000, 1/500

  # 2. strata/area to remove
  # contrling the type: remove depth layer(s), remove 1 horizontal layer (i.e. across depth), hole in the middle (or percentage)
  # controling the size: Several small vs one big

  # 3. only few years (1 vs 3) vs whole years

  # 4. type of population structure:
  #wide distribution :  ays_covar = sim_ays_covar(phi_age = 0.8, phi_year = 0.1, range=600, sd=1)
  #more pactchy: ays_covar = sim_ays_covar(phi_age = 0.8, phi_year = 0.1, range=600, sd=5)

  # 5. Temporal correlation in species distribution
  # phi_year= 0.1 = uncorrelated between year, 0.9 = correlated between years

scen <- expand.grid(iter=1:30, set_den=c(1/1000, 1/500), sd = c(1,5), phi_year =c(0.1, 0.9),
                    year_toremove = c(0, 3), fraction_x = c(0, 0.1), fraction_y = c(0, 0.1))


# sequential:
# result <- purrr::pmap_dfr(scen, sim_and_fit)
#
#


# parallel:
result <- furrr::future_pmap_dfr(scen, sim_and_fit,
                                 .options = furrr::furrr_options(seed = TRUE))




not_converged <- dplyr::filter(result, (bad_eig | max_gradient > 0.001) & type == "Model-based")
stopifnot(nrow(not_converged) == 0L)

result_scaled <- result %>%
  group_by(type, iter, scenario) %>%
  mutate(geo_mean = exp(mean(log(N), na.rm = TRUE)),
         lwr_scaled = lwr / geo_mean, N_scaled = N / geo_mean, upr_scaled = upr / geo_mean,
         type = factor(type, levels = c("True", "Design-based", "Model-based"))) %>%
  ungroup() %>%
  arrange(type)

result_scaled_model <- result_scaled %>% filter(type == "Model-based")
result_scaled_design <- result_scaled %>% filter(type == "Design-based")
result_scaled_true <- result_scaled %>% filter(type == "True")

result_scaled_model$RE <- (result_scaled_model$N_scaled - result_scaled_true$N_scaled)/result_scaled_true$N_scaled
result_scaled_design$RE <- (result_scaled_design$N_scaled - result_scaled_true$N_scaled)/result_scaled_true$N_scaled

results_scaled_new <-  result_scaled_model %>%
  bind_rows(result_scaled_design) %>%
  group_by(scenario, type, year) %>%
  summarize(RE_median = median(RE),
            RE_lwr = quantile(RE, 0.05, na.rm=T),
            RE_upr = quantile(RE, 0.95, na.rm=T))

ggplot(results_scaled_new, aes(x=year, y=RE_median, col=type, fill=type, linetype=type)) + geom_ribbon(aes(ymin=RE_lwr, ymax=RE_upr),alpha=0.3) +
  facet_wrap(~scenario, scales = "free_y") + geom_line() +geom_hline(yintercept=0)+
  scale_color_manual(values = c("Model-based" = "red", "Design-based" = "steelblue")) +
  scale_fill_manual(values = c("Model-based" = "red", "Design-based" = "steelblue")) +
  scale_linetype_manual(values = c("Model-based" = 1, "Design-based" = 2))


#### running conditinal random forest to examine simulation results

library(party)
new_dat = result_scaled_model %>% bind_rows(result_scaled_design)
Analysis <- cforest(data= new_dat, formula = RE ~ set_den + sd + phi_year+year_toremove+fraction_x+fraction_y)
barplot(varimp(Analysis))



#
# result_scaled %>%
#   ggplot(aes(year, N_scaled, group = type)) +
#   geom_line(aes(colour = type, size = type)) +
#   geom_ribbon(aes(ymin = lwr_scaled, ymax = upr_scaled, fill = type), alpha = 0.3) +
#   labs(x = "Year", y = "Relative abundance", colour = "Type", fill = "Type", size = "Type") +
#   scale_color_manual(values = c("Model-based" = "grey30", "Design-based" = "steelblue", "True" = "red")) +
#   scale_fill_manual(values = c("Model-based" = "grey30", "Design-based" = "steelblue", "True" = "red")) +
#   scale_size_manual(values = c("Model-based" = 0.5, "Design-based" = 0.5, "True" = 1)) +
#   facet_wrap(~scenario, scales = "free_y") +
#   theme_minimal() + theme(legend.position = "bottom")
#
# summary_stats <- result_scaled %>%
#   group_by(year, type, scenario) %>%
#   summarise(
#     est_lwr = lwr_scaled[type == "Model-based"],
#     est_upr = upr_scaled[type == "Model-based"],
#     est = N_scaled[type == "Model-based"],
#     true = N_scaled[type == "True"], .groups = "drop") %>%
#   mutate(covered = est_lwr < true & est_upr > true)
# mean(summary_stats$covered)
#
# ggplot(summary_stats, aes(log(true), log(est))) + geom_point() +
#   coord_fixed() + geom_abline(intercept = 0, slope = 1)
#
