library(dplyr)
library(ggplot2)
library(SimSurvey)
library(sdmTMB)
future::plan(future::multisession) # for the furrr package; or `plan(sequential)`

# iter <- 1

sim_and_fit <- function(iter) {
  set.seed(iter * 283028)
  sim <- sim_abundance(ages = seq(1, 10), years = seq(1, 10)) %>%
    sim_distribution(
      grid = make_grid(res = c(10, 10),
        # shelf_width = 20,
        # strat_breaks = seq(0, 1000, by = 200),
        depth_range = c(10, 500)),
      ays_covar = sim_ays_covar(phi_age = 0.8, phi_year = 0.1),
      depth_par = sim_parabola(mu = 200, sigma = 30)
    )

  # View(sim$grid_xy)
  # sim$grid_xy %>%
  #   ggplot(aes(x, y, fill = depth)) +
  #   geom_raster() +
  #   scale_fill_viridis_c(#limits=c(150, 300),
  #     direction = -1)

  survey <- sim_survey(sim, n_sims = 1) %>% run_strat()
  xy <- as_tibble(survey$grid_xy)
  dat <- as_tibble(survey$setdet) %>%
    select(x, y, set, year, N = n, tow_area)
  dat <- left_join(dat, xy, by = c("x", "y"))
  dat$offset <- log(dat$tow_area)

  grid_dat <- as_tibble(select(sim$grid_xy, x, y, depth)) %>% distinct()
  grid_dat <- purrr::map_dfr(sort(unique(dat$year)), ~ bind_cols(grid_dat, year = .))
  grid_dat$offset <- mean(dat$offset)



  mesh <- sdmTMB::make_mesh(dat, xy_cols = c("x", "y"), cutoff = 20)
  # plot(mesh)
  fit <- sdmTMB(N ~ 0 + as.factor(year) + offset,
    data = dat,
    family = nbinom2(link = "log"), spde = mesh,
    include_spatial = TRUE, time = "year"
  )
  fit_depth <- sdmTMB(N ~ 0 + as.factor(year) + s(depth, k = 3) + offset,
    data = dat,
    family = nbinom2(link = "log"), spde = mesh,
    include_spatial = TRUE, time = "year"
  )

  # nd <- data.frame(depth = seq(min(dat$depth), max(dat$depth), length.out = 80))
  # # You'll need at least one time element:
  # nd$year <- dat$year[1]
  # nd$offset <- mean(dat$offset)
  # p <- predict(fit, newdata = nd, se_fit = TRUE, re_form = NA)
  # ggplot(p, aes(depth, exp(est),
  #   ymin = exp(est - 1.96 * est_se), ymax = exp(est + 1.96 * est_se))) +
  #   geom_line() + geom_ribbon(alpha = 0.4)
  # ggplot(p, aes(depth, est,
  #   ymin = est - 1.96 * est_se, ymax = est + 1.96 * est_se)) +
  #   geom_line() + geom_ribbon(alpha = 0.4)

  pred <- predict(fit, newdata = grid_dat, return_tmb_object = TRUE, area = 100)
  pred_depth <- predict(fit_depth, newdata = grid_dat, return_tmb_object = TRUE, area = 100)
  index <- get_index(pred)
  index_depth <- get_index(pred_depth)

  true_abund <- tibble(year = unique(dat$year), N = as.numeric(colSums(survey$I))) %>%
    mutate(type = "True")
  strat_abund <- tibble::as_tibble(survey$total_strat) %>%
    mutate(N = total, type = "Design-based") %>%
    select(year, N, type)
  mutate(index, type = "Model-based", N = est) %>%
    bind_rows(mutate(index_depth, type = "Model-based-depth", N = est)) %>%
    bind_rows(strat_abund) %>%
    bind_rows(true_abund) %>%
    mutate(iter = iter)
}

# sequential:
# result <- purrr::map_dfr(seq_len(8), sim_and_fit)

# parallel:
result <- furrr::future_map_dfr(seq_len(8), sim_and_fit,
  .options = furrr::furrr_options(seed = TRUE))

not_converged <- dplyr::filter(result, (bad_eig | max_gradient > 0.001) & type == "Estimated")
stopifnot(nrow(not_converged) == 0L)

result_scaled <- result %>%
  group_by(type, iter) %>%
  mutate(geo_mean = exp(mean(log(N), na.rm = TRUE)),
    lwr_scaled = lwr / geo_mean, N_scaled = N / geo_mean, upr_scaled = upr / geo_mean,
    type = factor(type, levels = c("True", "Design-based", "Model-based", "Model-based-depth"))) %>%
  ungroup() %>%
  arrange(type)

result_scaled %>%
  group_by(type) %>%
  mutate(true_index = ifelse("True" %in% type, "yes", "no")) %>%
  ggplot(aes(year, N_scaled, group = type)) +
  geom_line(aes(colour = type, lty = true_index)) +
  geom_ribbon(aes(ymin = lwr_scaled, ymax = upr_scaled, fill = type), alpha = 0.3) +
  labs(x = "Year", y = "Relative abundance", colour = "Type", fill = "Type") +
  facet_wrap(~iter, scales = "free_y") +
  theme_minimal() + theme(legend.position = "bottom")

# is true value within the CI of model-based?
summary_stats1 <- result_scaled %>%
  group_by(year, iter) %>%
  summarise(
    est_lwr = lwr_scaled[type == "Model-based"],
    est_upr = upr_scaled[type == "Model-based"],
    est = N_scaled[type == "Model-based"],
    true = N_scaled[type == "True"], .groups = "drop") %>%
  mutate(covered = est_lwr < true & est_upr > true)
mean(summary_stats1$covered)

# is true value within the CI of model-based with depth?
summary_stats2 <- result_scaled %>%
  group_by(year, iter) %>%
  summarise(
    est_lwr = lwr_scaled[type == "Model-based-depth"],
    est_upr = upr_scaled[type == "Model-based-depth"],
    est = N_scaled[type == "Model-based-depth"],
    true = N_scaled[type == "True"], .groups = "drop") %>%
  mutate(covered = est_lwr < true & est_upr > true)
mean(summary_stats2$covered)

# compare true with predictions as scatter plot
ggplot(summary_stats1, aes(log(true), log(est))) + geom_point() +
  coord_fixed() + geom_abline(intercept = 0, slope = 1)

ggplot(summary_stats2, aes(log(true), log(est))) + geom_point() +
  coord_fixed() + geom_abline(intercept = 0, slope = 1)

