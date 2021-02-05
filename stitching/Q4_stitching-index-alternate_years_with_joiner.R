
## Question 4:
#- Same vessel samples different areas in alternating years, one year in common
#- Compare to global survey with half the set density (same effort)
#- Different q between areas

library(dplyr)
library(ggplot2)
library(SimSurvey)
library(sdmTMB)
future::plan(future::multisession) # for the furrr package; or `plan(sequential)`
#set.seed(1) # Explore simulation code line-by-line
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
  # Global survey with half the set density compared to the other two "half-surveys"
  survey_global <- sim_survey(sim, n_sims = 1, set_den = 1/1000) %>% run_strat()
  xy_global <- as_tibble(survey_global$grid_xy) # Names of strata
  dat_global <- as_tibble(survey_global$setdet) %>% # Numbers caught and tow area
    select(x, y, set, year, N = n, tow_area)
  dat_global <- left_join(dat_global, xy, by = c("x", "y"))
  dat_global$offset <- log(dat_global$tow_area)

  # Grid of survey area
  grid_dat <- as_tibble(select(sim$grid_xy, x, y, depth)) %>% distinct()
  grid_dat <- purrr::map_dfr(sort(unique(dat$year)), ~ bind_cols(grid_dat, year = .))
  grid_dat$offset <- mean(dat$offset)

  dat_filtered1 <- dplyr::filter(dat, x > 0, year %in% seq(2, 10, 2)) %>%
    mutate(survey = "A") # Survey #1 (A) samples right side of grid in even years and terminal year
  dat_filtered2 <- dplyr::filter(dat2, x < 0, year %in% c(seq(1, 9, 2), 10)) %>%
    mutate(survey = "B") # Survey #2 (B) samples left side of grid in odd years and terminal year
  dat_combined <- bind_rows(dat_filtered1, dat_filtered2)
  # ggplot(dat_combined, aes(x, y, colour = survey)) + geom_point() + facet_wrap(~year) + theme_bw()
  # ggsave("stitching/stitched_design.png", height = 4, width = 6)

  mesh <- sdmTMB::make_mesh(dat_combined, xy_cols = c("x", "y"), cutoff = 20)
  # plot(mesh)

  # Stitched survey
  fit <- sdmTMB(N ~ 0 + as.factor(year) + offset + as.factor(survey),
    data = dat_combined,
    family = nbinom2(link = "log"), spde = mesh,
    include_spatial = TRUE, time = "year", silent = TRUE
  )
  grid_dat$survey <- "A"
  pred <- predict(fit, newdata = grid_dat, return_tmb_object = TRUE, area = 100)
  index <- get_index(pred)

  # Global survey
  fit_global <- sdmTMB(N ~ 0 + as.factor(year) + offset,
                       data = dat_global,
                       family = nbinom2(link = "log"), spde = mesh,
                       include_spatial = TRUE, time = "year", silent = TRUE
  )
  pred_global <- predict(fit_global, newdata = grid_dat, return_tmb_object = TRUE, area = 100)
  index_global <- get_index(pred_global)

  # True underlying abundance in each area
  xy_sim <- as_tibble(sim$grid_xy)
  df_sim <- as_tibble(sim$sp_N)
  df_sim <- left_join(df_sim, xy_sim, by = "cell")

  N_true <- df_sim %>%
    group_by(x, y, year, cell) %>%
    summarise(N = sum(N), .groups = "drop")

  N_A <- N_true %>% dplyr::filter(x > 0) %>% group_by(year) %>% summarise(N = sum(N)) %>%
    mutate(type = "True (A)", type2 = "subarea")
  N_B <- N_true %>% dplyr::filter(x < 0) %>% group_by(year) %>% summarise(N = sum(N)) %>%
    mutate(type = "True (B)", type2 = "subarea")
  #plot(N/mean(N) ~ year, N_A, typ = 'o', ylim = c(0, 4))
  #lines(N/mean(N) ~ year, N_B, typ = 'o', col = "red")
  #lines(est/mean(est) ~ year, index, type = "o", col = "blue")
  #lines(est/mean(est) ~ year, index_global, type = "o", col = "green")
  #legend("topleft", c("A", "B", "stitched", "global"), col = c("black", "red", "blue", "green"), lwd = 2)

  #plot(survey$I %>% colSums() %>% `/`(mean(.)), ylab = "Global survey", typ = 'o', ylim = c(0, 3))
  #lines(est/mean(est) ~ year, index, type = "o", col = "blue")
  #lines(est/mean(est) ~ year, index_global, type = "o", col = "green")
  #legend("topleft", c("true", "stitched", "green"), col = c("black", "blue", "green"), lwd = 2)

  # extract q ratio estimate:
  q_hat <- tidy(fit) %>% dplyr::filter(term == "as.factor(survey)B")

  true_abund <- N_true %>% group_by(year) %>% summarise(N = sum(N)) %>%
    mutate(type = "True (global)", type2 = "global")

  #strat_abund <- tibble::as_tibble(survey$total_strat) %>%
  #  mutate(N = total, type = "Design-based") %>%
  #  select(year, N, type)
  mutate(index, type = "Model (stitch A/B)", N = est, type2 = "global") %>%
    bind_rows(mutate(index_global, type = "Model (global)", N = est, type2 = "global")) %>%
    #bind_rows(strat_abund) %>%
    bind_rows(true_abund) %>% bind_rows(N_A) %>% bind_rows(N_B) %>%
    mutate(iter = iter) #%>%
    #mutate(q_est = q_hat$estimate, q_se = q_hat$std.error)
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
         lwr_scaled = lwr / geo_mean, N_scaled = N / geo_mean, upr_scaled = upr / geo_mean) %>%
  ungroup() %>%
  arrange(type)

result_scaled %>%
  ggplot(aes(year, N_scaled, group = type)) +
  geom_line(aes(colour = type)) + geom_point(aes(colour = type, size = grepl("True", type) %>% as.character())) +
  geom_linerange(aes(ymin = lwr_scaled, ymax = upr_scaled, colour = type), size = 1.5) +
  labs(x = "Year", y = "Relative abundance", colour = "Type", fill = "Type", size = "Type") +
  scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
  facet_grid(type2~iter, scales = "free_y") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("stitching/stitch_iter_with_errorbar.png", width = 7, height = 4)


result_scaled %>%
  ggplot(aes(year, N_scaled, group = type)) +
  geom_line(aes(colour = type)) + geom_point(aes(colour = type, size = grepl("True", type) %>% as.character())) +
  labs(x = "Year", y = "Relative abundance", colour = "Type", fill = "Type", size = "Type") +
  scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
  facet_grid(type2~iter, scales = "free_y") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("stitching/stitch_iter.png", width = 7, height = 4)

result %>%
  ggplot(aes(year, N, group = type)) +
  geom_line(aes(colour = type)) + geom_point(aes(colour = type)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = type), alpha = 0.3) +
  labs(x = "Year", y = "Relative abundance", colour = "Type", fill = "Type", size = "Type") +
  facet_grid(type2~iter, scales = "free_y") +
  theme_bw() + theme(legend.position = "bottom")
