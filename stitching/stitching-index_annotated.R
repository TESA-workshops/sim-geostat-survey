library(dplyr)
library(ggplot2)
library(SimSurvey)
library(sdmTMB)
future::plan(future::multisession) # for the furrr package; or `plan(sequential)`
set.seed(1) # Explore simulation code line-by-line
sim_and_fit <- function(iter) {
  set.seed(iter * 283028)
  sim <- sim_abundance(ages = seq(1, 10), years = seq(1, 10)) %>%
    sim_distribution(
      grid = make_grid(res = c(10, 10), depth_range = c(10, 500)),
      ays_covar = sim_ays_covar(phi_age = 0.8, phi_year = 0.1),
      depth_par = sim_parabola(mu = 200, sigma = 30)
    )


  # this is the first survey
  survey <- sim_survey(sim, n_sims = 1) %>% run_strat()
  xy <- as_tibble(survey$grid_xy) # Names of strata
  dat <- as_tibble(survey$setdet) %>% # Numbers caught and tow area
    select(x, y, set, year, N = n, tow_area)
  dat <- left_join(dat, xy, by = c("x", "y"))
  dat$offset <- log(dat$tow_area)

  # catchability by age:
  # logistic_fun <- sim_logistic(k = 2, x0 = 3, plot = TRUE)
  # logistic_fun(x = 1:10)
  # logistic_fun2 <- sim_logistic(k = 2, x0 = 6, plot = TRUE)
  # logistic_fun2(x = 1:10)
  # survey2 <- sim_survey(sim, n_sims = 1, q = sim_logistic(x0 = 6)) %>% run_strat()

  # or, more simply, use a scalar q here:
  # This is a second survey data set over the same operating model with q = 66 % of that in survey 1
  survey2 <- sim_survey(sim, n_sims = 1) %>% run_strat()
  xy2 <- as_tibble(survey2$grid_xy)
  dat2 <- as_tibble(survey2$setdet) %>%
    select(x, y, set, year, N = n, tow_area)
  dat2 <- left_join(dat2, xy2, by = c("x", "y"))
  dat2$offset <- log(dat2$tow_area)
  dat2$N <- round(dat2$N * 0.66) # 66% of q from survey "A"

  # Compare to a third survey where the entire survey area is sampled
  # Value of information analysis: is it worth it to sample entire area annually vs. alternating sides
  # and stitching? With same effort (# of sets) in either case

  # Grid of survey area
  grid_dat <- as_tibble(select(sim$grid_xy, x, y, depth)) %>% distinct()
  grid_dat <- purrr::map_dfr(sort(unique(dat$year)), ~ bind_cols(grid_dat, year = .))
  grid_dat$offset <- mean(dat$offset)

  dat_filtered1 <- dplyr::filter(dat, x > 0) %>% mutate(survey = "A") # Survey #1 (A) samples right side of grid
  dat_filtered2 <- dplyr::filter(dat2, x < 0) %>% mutate(survey = "B") # Survey #2 (B) samples left side of grid
  dat_combined <- bind_rows(dat_filtered1, dat_filtered2)
  # ggplot(dat_combined, aes(x, y, colour = survey)) + geom_point() + facet_wrap(~year)

  mesh <- sdmTMB::make_mesh(dat_combined, xy_cols = c("x", "y"), cutoff = 20)
  # plot(mesh)

  fit <- sdmTMB(N ~ 0 + as.factor(year) + offset + as.factor(survey),
    data = dat_combined,
    family = nbinom2(link = "log"), spde = mesh,
    include_spatial = TRUE, time = "year", silent = TRUE
  )
  grid_dat$survey <- "A"
  pred <- predict(fit, newdata = grid_dat, return_tmb_object = TRUE, area = 100)
  index <- get_index(pred)

  # q_A = q_B
  true_A <- sim_survey(sim, n_sims = 1)
  true_A$setdet <- true_A$setdet %>% dplyr::filter(x > 0)
  true_A <- run_strat(true_A)
  #matplot(t(true_A$I), typ = 'l')
  plot(mean ~ year, true_A$total_strat, typ = "o")
  #true_A$I %>% colSums() %>% plot(type = 'o')

  true_B <- sim_survey(sim, n_sims = 1)
  true_B$setdet <- true_B$setdet %>% dplyr::filter(x < 0)
  true_B <- run_strat(true_B)
  #matplot(true_B$I, typ = 'l')
  lines(0.66 * mean ~ year, true_B$total_strat, typ = "o", col = "red")
  #true_A$I %>% colSums() %>% plot(type = 'o')
  lines(est/mean(est) ~ year, index, type = "o", col = "blue")
  legend("topleft", c("A", "B", "stitched"), col = c("black", "red", "blue"), lwd = 2)

  plot(survey$I %>% colSums() %>% `/`(mean(.)), ylab = "Global survey", typ = 'o', ylim = c(0, 2))
  lines(est/mean(est) ~ year, index, type = "o", col = "blue")
  legend("topleft", c("true", "stitched"), col = c("black", "blue"), lwd = 2)

  # extract q ratio estimate:
  q_hat <- tidy(fit) %>% dplyr::filter(term == "as.factor(survey)B")

  true_abund <- tibble(year = unique(dat$year), N = as.numeric(colSums(survey$I))) %>%
    mutate(type = "True")
  strat_abund <- tibble::as_tibble(survey$total_strat) %>%
    mutate(N = total, type = "Design-based") %>%
    select(year, N, type)
  mutate(index, type = "Model-based", N = est) %>%
    bind_rows(strat_abund) %>%
    bind_rows(true_abund) %>%
    mutate(iter = iter) %>%
    mutate(q_est = q_hat$estimate, q_se = q_hat$std.error)
}

# sequential:
# result <- purrr::map_dfr(seq_len(8), sim_and_fit)

# parallel:
result <- furrr::future_map_dfr(seq_len(4), sim_and_fit,
  .options = furrr::furrr_options(seed = TRUE))

not_converged <- dplyr::filter(result, (bad_eig | max_gradient > 0.001) & type == "Model-based")
stopifnot(nrow(not_converged) == 0L)

result_scaled <- result %>%
  group_by(type, iter) %>%
  mutate(geo_mean = exp(mean(log(N), na.rm = TRUE)),
    lwr_scaled = lwr / geo_mean, N_scaled = N / geo_mean, upr_scaled = upr / geo_mean,
    type = factor(type, levels = c("True", "Design-based", "Model-based"))) %>%
  ungroup() %>%
  arrange(type)

result_scaled %>%
  ggplot(aes(year, N_scaled, group = type)) +
  geom_line(aes(colour = type, size = type)) +
  geom_ribbon(aes(ymin = lwr_scaled, ymax = upr_scaled, fill = type), alpha = 0.3) +
  labs(x = "Year", y = "Relative abundance", colour = "Type", fill = "Type", size = "Type") +
  scale_color_manual(values = c("Model-based" = "grey30", "Design-based" = "steelblue", "True" = "red")) +
  scale_fill_manual(values = c("Model-based" = "grey30", "Design-based" = "steelblue", "True" = "red")) +
  scale_size_manual(values = c("Model-based" = 0.5, "Design-based" = 0.5, "True" = 1)) +
  facet_wrap(~iter, scales = "free_y") +
  theme_minimal() + theme(legend.position = "bottom")

summary_stats <- result_scaled %>%
  group_by(year, iter) %>%
  summarise(
    est_lwr = lwr_scaled[type == "Model-based"],
    est_upr = upr_scaled[type == "Model-based"],
    est = N_scaled[type == "Model-based"],
    true = N_scaled[type == "True"], .groups = "drop") %>%
  mutate(covered = est_lwr < true & est_upr > true)
mean(summary_stats$covered)

ggplot(summary_stats, aes(log(true), log(est))) + geom_point() +
  coord_fixed() + geom_abline(intercept = 0, slope = 1)

result %>% select(iter, q_est, q_se) %>%
  distinct() %>%
  mutate(est = exp(q_est), lwr = exp(q_est - 1.96 * q_se), upr = exp(q_est + 1.96 * q_se)) %>%
  ggplot(aes(est, iter, xmin = lwr, xmax = upr)) + geom_pointrange() +
  geom_vline(xintercept = 0.66, lty = 2)

