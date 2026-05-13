#' Fluctuations in direct human presence, not predictable weekly cycles, influence avoidance behaviour in ravens
#' Script by: Varalika Jain

#' This code supplements '1 data_prep.R'

#### KDE DOWNSAMPLING EXERCISE ----
library(segmented)

##### (1) Create sample size set ----
#' Testing out minimum number of fixes needed in a day
set.seed(123)
xx <- final_5yr_filtered %>% group_by(id_date, date, mmdd, year) %>% tally() %>%
  filter(n>30)

df_day <- final_5yr_filtered %>% filter(id_date %in% xx$id_date)

#' sample sizes to test
sample_sizes <- c(6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30)
n_reps <- 20  # number of random repeats for each sample size

names(df_day)

##### (2) Create track ----
# Convert day's points to an amt track object
track_day <- make_track(
  df_day,
  .id = id_date,
  .x = location_long,
  .y = location_lat,
  .t = timestamps,
  crs = "+proj=longlat +datum=WGS84"
)

# Split full track object by individual
track_list <- split(track_day, track_day$.id)
set.seed(180)
sampled_ids <- sample(names(track_list), 200)
track_list <- track_list[sampled_ids]

##### (3) Reference & downsampled KDEs ----
# get 95% isopleth polygon from an amt track (full iso)
sf_use_s2(FALSE)
run_kde_downsampling <- function(track_obj, sample_sizes, n_reps) {
  id <- track_obj$.id[1]
  message("Processing individual: ", id)
  
  # reference KDE -
  get_kde_iso50 <- function(track_obj) {
    # Compute KDE (still in lon/lat)
    hr <- hr_kde(track_obj, levels = 0.95)
    iso <- hr_isopleths(hr)
    return(iso)
  }
  
  # downsampling KDE
  downsample_kde_iso <- function(track_obj, n, rep_id) {
    if (nrow(track_obj) < n) return(NULL)
    sampled <- track_obj %>% slice_sample(n = n) %>% arrange(t_)
    iso <- get_kde_iso50(sampled)
    
    full_iso <- hr_isopleths(hr_kde(track_obj, levels = 0.95))
    full_area <- hr_area(hr_kde(track_obj, levels = 0.95), levels = 0.95)
    full_area <- full_area$area
    
    inter <- tryCatch(st_intersection(full_iso, iso), error = function(e) NULL)
    inter_area <- if (!is.null(inter) && length(inter) > 0) 
      sum(as.numeric(st_area(st_make_valid(inter)))) else 0
    
    iso_area <- hr_area(hr_kde(sampled, levels = 0.95))$area
    
    # summary
    tibble(
      new_ID = id,
      sample_size = n,
      rep = rep_id,
      iso_area = iso_area,
      inter_area = inter_area,
      prop_overlap_full = ifelse(full_area > 0, inter_area / full_area, NA)
    )
  }
  
  # run for all sample sizes + reps
  all_metrics <- purrr::map_dfr(sample_sizes, function(n) {
    purrr::map_dfr(1:n_reps, function(rep_id) {
      downsample_kde_iso(track_obj, n, rep_id)
    })
  })
  
  return(all_metrics)
}

#' running into some individuals with errors, so skipping those that are 
#' problematic
safe_run <- purrr::possibly(
  run_kde_downsampling,
  otherwise = NULL,   # return NULL on error
  quiet = FALSE
)

#' applying to whole list
all_results <- purrr::map_dfr(
  track_list,
  safe_run,
  sample_sizes = sample_sizes,
  n_reps = n_reps
)

all_results %>% group_by(new_ID) %>% tally() #88 individuals

##### (4) Breakpoint analysis ----
#' breakpoint analysis 
fit_segmented_breakpoint <- function(df) {
  # fit basic linear model
  fit <- lm(iso_area ~ sample_size, data = df)
  
  # segmented regression with starting breakpoint (psi = 15)
  seg_fit <- segmented(fit, seg.Z = ~sample_size, psi = 15)
  
  # extract breakpoint and SE
  bp <- summary(seg_fit)$psi
  
  # summary
  tibble(
    new_ID = unique(df$new_ID),
    breakpoint = bp[,"Est."],
    se = bp[,"St.Err"],
    response = "iso_area"
  )
}

# split data by individual
data_list <- split(all_results, all_results$new_ID)

#' running into some individuals with errors, so skipping those that are 
#' problematic
safe_run_2 <- purrr::possibly(
  fit_segmented_breakpoint,
  otherwise = NULL,   # return NULL on error
  quiet = FALSE
)

# run segmented regression for iso_area
breakpoints_iso <- map_dfr(data_list, safe_run_2)

fit_segmented_breakpoint <- function(df) {
  # Fit basic linear model
  fit <- lm(prop_overlap_full ~ sample_size, data = df)
  
  # Segmented regression with starting breakpoint (psi = 15)
  seg_fit <- segmented(fit, seg.Z = ~sample_size, psi = 15)
  
  # Extract breakpoint and SE
  bp <- summary(seg_fit)$psi
  
  tibble(
    new_ID = unique(df$new_ID),
    breakpoint = bp[,"Est."],
    se = bp[,"St.Err"],
    response = "prop_overlap_full"
  )
}

#' running into some individuals with errors, so skipping those that are 
#' problematic
safe_run_2 <- purrr::possibly(
  fit_segmented_breakpoint,
  otherwise = NULL,   # return NULL on error
  quiet = FALSE
)

# Run segmented regression for prop_overlap_full
breakpoints_overlap <- map_dfr(data_list, safe_run_2)

# View results
breakpoints_iso
breakpoints_overlap

breakpoints_iso %>% ggplot(aes(x = breakpoint))+geom_histogram()
breakpoints_overlap %>% ggplot(aes(x = breakpoint))+geom_histogram()

summary(breakpoints_iso$breakpoint)
summary(breakpoints_overlap$breakpoint)

#' example from 1 individual
df = data_list[["Rollo2021-06-20"]]
fit <- lm(iso_area ~ sample_size, data = df)
# Segmented regression with starting breakpoint (psi = 15)
seg_fit <- segmented(fit, seg.Z = ~sample_size, psi = 15)

pdf("SIplotx6.pdf",
    width = 6, height = 6)
# Plot points
plot(df$sample_size, df$iso_area,
     pch = 19, col = "gray40",
     xlab = "Sample size",
     ylab = "Area")
# Add segmented regression fit
plot(seg_fit, add = TRUE, col = "red", lwd = 2)
# Optional: add breakpoint line
bp <- summary(seg_fit)$psi[ , "Est."]
abline(v = bp, col = "blue", lty = 2, lwd = 2)
dev.off()

df = data_list[["Rollo2021-06-20"]]
fit <- lm(prop_overlap_full ~ sample_size, data = df)
# Segmented regression with starting breakpoint (psi = 15)
seg_fit <- segmented(fit, seg.Z = ~sample_size, psi = 15)

pdf("SIplotx7.pdf",
    width = 6, height = 6)
# Plot points
plot(df$sample_size, df$prop_overlap_full,
     pch = 19, col = "gray40",
     xlab = "Sample size",
     ylab = "Proportion overlap")
# Add segmented regression fit
plot(seg_fit, add = TRUE, col = "red", lwd = 2)
# Optional: add breakpoint line
bp <- summary(seg_fit)$psi[ , "Est."]
abline(v = bp, col = "blue", lty = 2, lwd = 2)
dev.off()


##### (5) Plotting for supplementary information ----
all_metrics <- all_results %>% 
  filter(new_ID == 'Rollo2021-06-20')

# First compute summary stats
summary_metrics <- all_results %>% 
  filter(new_ID == 'Rollo2021-06-20') %>%
  group_by(sample_size) %>%
  summarise(
    mean_area = mean(iso_area, na.rm = TRUE),
    lower_ci = quantile(iso_area, 0.025, na.rm = TRUE),
    upper_ci = quantile(iso_area, 0.975, na.rm = TRUE),
    mean_overlap = mean(prop_overlap_full, na.rm = TRUE),
    lower_ci_overlap = quantile(prop_overlap_full, 0.025, na.rm = TRUE),
    upper_ci_overlap = quantile(prop_overlap_full, 0.975, na.rm = TRUE)
  )

p1<- ggplot() +
  geom_jitter(
    data = all_metrics,
    aes(x = sample_size, y = iso_area/1000000, color = as.factor(sample_size)),
    width = 2, height = 0, alpha = 0.6
  ) +
  geom_ribbon(
    data = summary_metrics,
    aes(x = sample_size, ymin = lower_ci/1000000, ymax = upper_ci/1000000),
    inherit.aes = FALSE, fill = "grey80", alpha = 0.3
  ) +
  geom_line(
    data = summary_metrics,
    aes(x = sample_size, y = mean_area/1000000),
    color = "black", linewidth = 0.5
  ) +
  scale_x_continuous(breaks = unique(all_metrics$sample_size)) +
  labs(
    x = "Number of points",
    y = "KDE area (km²)",
    title = "a)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")+
  geom_vline(xintercept = 15, linetype="dotted", 
             color = "red", size=1)


p2 <- ggplot() +
  geom_jitter(
    data = all_metrics,
    aes(x = sample_size, y = prop_overlap_full, color = as.factor(sample_size)),
    width = 2, height = 0, alpha = 0.6
  ) +
  geom_ribbon(
    data = summary_metrics,
    aes(x = sample_size, ymin = lower_ci_overlap, ymax = upper_ci_overlap),
    inherit.aes = FALSE, fill = "grey80", alpha = 0.3
  ) +
  geom_line(
    data = summary_metrics,
    aes(x = sample_size, y = mean_overlap),
    color = "black", linewidth = 0.5
  ) +
  scale_x_continuous(breaks = unique(all_metrics$sample_size)) +
  labs(
    x = "Number of points",
    y = "Proportion overlap",
    title = "b)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")+
  geom_vline(xintercept = 15, linetype="dotted", 
             color = "red", size=1) 


library(ggpubr)
quartz(height = 6, width = 12)
ggarrange(p1, p2,ncol=2, nrow=1)
dev.copy2pdf(file="SIplotx1.pdf")

#' map    
n_reps <- 5  # for visualisation

downsample_kde <- function(track_obj, n, rep_id) {
  #if (nrow(track_obj) < n) return(NULL)
  
  sampled_track <- track_obj %>%
    slice_sample(n = n) %>%
    arrange(t_)
  
  kde_res <- amt::hr_kde(sampled_track, levels = 0.95)
  
  # Extract the 50% isopleth polygon as sf
  kde_df <- hr_isopleths(kde_res) %>%
    mutate(sample_size = n, rep = rep_id)
  
  return(kde_df)
}

# Run downsampling loop
kde_results <- purrr::map_dfr(sample_sizes, function(n) {
  purrr::map_dfr(1:n_reps, function(rep_id) {
    downsample_kde(track_list[['Rollo2021-06-20']], n, rep_id)
  })
})

# Plot with full-data overlay
quartz(height = 8, width = 18)
ggplot() +
  geom_sf(data = kde_results, aes(fill = factor(sample_size)), color = "grey", alpha = 0.5) +
  facet_grid(rep ~ sample_size) +
  theme_minimal() +
  labs(
    fill = "Sample size"
  )+
  geom_point(data = track_list[['Rollo2021-06-20']], aes(x = x_, y = y_), size = 0.1, alpha = 0.4)+
  labs(x = NULL, y = NULL, title = "c)") +
  scale_x_continuous(labels = ~ .x) +
  scale_y_continuous(labels = ~ .x) +
  theme(legend.position = "none")
dev.copy2pdf(file="SIplotx2.pdf")

#### MCP DOWNSAMPLING EXERCISE ----
##### (1) Reference & downsampled KDEs ----
run_mcp_downsampling <- function(track_obj, sample_sizes, n_reps) {
  id <- track_obj$.id[1]
  message("Processing individual: ", id)
  
  # reference MCP
  get_mcp_iso50 <- function(track_obj) {
    # Compute KDE (still in lon/lat)
    hr <- hr_mcp(track_obj, levels = 0.95)
    iso <- hr_isopleths(hr)
    return(iso)
  }
  
  # downsampling MCP function
  downsample_kde_iso <- function(track_obj, n, rep_id) {
    if (nrow(track_obj) < n) return(NULL)
    sampled <- track_obj %>% slice_sample(n = n) %>% arrange(t_)
    iso <- get_mcp_iso50(sampled)
    
    full_iso <- hr_isopleths(hr_mcp(track_obj, levels = 0.95))
    full_area <- hr_area(hr_mcp(track_obj, levels = 0.95), levels = 0.95)
    full_area <- full_area$area
    
    inter <- tryCatch(st_intersection(full_iso, iso), error = function(e) NULL)
    inter_area <- if (!is.null(inter) && length(inter) > 0) 
      sum(as.numeric(st_area(st_make_valid(inter)))) else 0
    
    iso_area <- hr_area(hr_kde(sampled, levels = 0.95))$area
    
    # Summary
    tibble(
      new_ID = id,
      sample_size = n,
      rep = rep_id,
      iso_area = iso_area,
      inter_area = inter_area,
      prop_overlap_full = ifelse(full_area > 0, inter_area / full_area, NA)
    )
  }
  
  # run for all sample sizes + reps
  all_metrics <- purrr::map_dfr(sample_sizes, function(n) {
    purrr::map_dfr(1:n_reps, function(rep_id) {
      downsample_kde_iso(track_obj, n, rep_id)
    })
  })
  
  return(all_metrics)
}
n_reps <- 20  # for visualisation

safe_run_3 <- purrr::possibly(
  run_mcp_downsampling,
  otherwise = NULL,   # return NULL on error
  quiet = FALSE
)

#' applying to first item of list
mcp_results <- purrr::map_dfr(
  track_list['Rollo2021-06-20'],
  safe_run_3,
  sample_sizes = sample_sizes,
  n_reps = n_reps
)

# First compute summary stats
summary_metrics_mcp <- mcp_results %>% 
  group_by(sample_size) %>%
  summarise(
    mean_area = mean(iso_area, na.rm = TRUE),
    lower_ci = quantile(iso_area, 0.025, na.rm = TRUE),
    upper_ci = quantile(iso_area, 0.975, na.rm = TRUE),
    mean_overlap = mean(prop_overlap_full, na.rm = TRUE),
    lower_ci_overlap = quantile(prop_overlap_full, 0.025, na.rm = TRUE),
    upper_ci_overlap = quantile(prop_overlap_full, 0.975, na.rm = TRUE)
  )

##### (2) Plotting for supplementary information ----
p1<- ggplot() +
  geom_jitter(
    data = mcp_results,
    aes(x = sample_size, y = iso_area/1000000, color = as.factor(sample_size)),
    width = 2, height = 0, alpha = 0.6
  ) +
  geom_ribbon(
    data = summary_metrics_mcp,
    aes(x = sample_size, ymin = lower_ci/1000000, ymax = upper_ci/1000000),
    inherit.aes = FALSE, fill = "grey80", alpha = 0.3
  ) +
  geom_line(
    data = summary_metrics_mcp,
    aes(x = sample_size, y = mean_area/1000000),
    color = "black", linewidth = 0.5
  ) +
  scale_x_continuous(breaks = unique(mcp_results$sample_size)) +
  labs(
    x = "Number of points",
    y = "MCP area (km²)",
    title = "a)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")+
  geom_vline(xintercept = 15, linetype="dotted", 
             color = "red", size=1)


p2 <- ggplot() +
  geom_jitter(
    data = mcp_results,
    aes(x = sample_size, y = prop_overlap_full, color = as.factor(sample_size)),
    width = 2, height = 0, alpha = 0.6
  ) +
  geom_ribbon(
    data = summary_metrics_mcp,
    aes(x = sample_size, ymin = lower_ci_overlap, ymax = upper_ci_overlap),
    inherit.aes = FALSE, fill = "grey80", alpha = 0.3
  ) +
  geom_line(
    data = summary_metrics_mcp,
    aes(x = sample_size, y = mean_overlap),
    color = "black", linewidth = 0.5
  ) +
  scale_x_continuous(breaks = unique(mcp_results$sample_size)) +
  labs(
    x = "Number of points",
    y = "Proportion overlap",
    title = "b)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")+
  geom_vline(xintercept = 15, linetype="dotted", 
             color = "red", size=1) 

quartz(height = 6, width = 12)
ggarrange(p1, p2,ncol=2, nrow=1)
dev.copy2pdf(file="SIplotx3.pdf")

#' map    
n_reps <- 5  # for visualisation

downsample_kde <- function(track_obj, n, rep_id) {
  #if (nrow(track_obj) < n) return(NULL)
  
  sampled_track <- track_obj %>%
    slice_sample(n = n) %>%
    arrange(t_)
  
  kde_res <- amt::hr_mcp(sampled_track, levels = 0.95)
  
  # Extract the 50% isopleth polygon as sf
  kde_df <- hr_isopleths(kde_res) %>%
    mutate(sample_size = n, rep = rep_id)
  
  return(kde_df)
}

# Run downsampling loop
kde_results <- purrr::map_dfr(sample_sizes, function(n) {
  purrr::map_dfr(1:n_reps, function(rep_id) {
    downsample_kde(track_list[['Rollo2021-06-20']], n, rep_id)
  })
})


# Plot with full-data overlay
quartz(height = 8, width = 18)

ggplot() +
  geom_sf(data = kde_results, aes(fill = factor(sample_size)), color = "grey", alpha = 0.5) +
  facet_grid(rep ~ sample_size) +
  theme_minimal() +
  labs(
    fill = "Sample size"
  )+
  geom_point(data = track_list[['Rollo2021-06-20']], aes(x = x_, y = y_), size = 0.1, alpha = 0.4)+
  labs(x = NULL, y = NULL, title = "c)") +
  scale_x_continuous(labels = ~ .x) +
  scale_y_continuous(labels = ~ .x) +
  theme(legend.position = "none")

dev.copy2pdf(file="SIplotx4.pdf")


##### (3) Breakpoint analysis ----
# get 95% isopleth polygon from an amt track (full iso)
sf_use_s2(FALSE)

n_reps <- 20
track_list_mcp <- track_list[names(data_list)]

run_mcp_downsampling <- function(track_obj, sample_sizes, n_reps) {
  id <- track_obj$.id[1]
  message("Processing individual: ", id)
  
  # reference MCP
  get_mcp_iso50 <- function(track_obj) {
    # Compute KDE (still in lon/lat)
    hr <- hr_mcp(track_obj, levels = 0.95)
    iso <- hr_isopleths(hr)
    return(iso)
  }
  
  # downsampling MCP function
  downsample_mcp_iso <- function(track_obj, n, rep_id) {
    if (nrow(track_obj) < n) return(NULL)
    sampled <- track_obj %>% slice_sample(n = n) %>% arrange(t_)
    iso <- get_mcp_iso50(sampled)
    
    full_iso <- hr_isopleths(hr_mcp(track_obj, levels = 0.95))
    full_area <- hr_area(hr_mcp(track_obj, levels = 0.95), levels = 0.95)
    full_area <- full_area$area
    
    inter <- tryCatch(st_intersection(full_iso, iso), error = function(e) NULL)
    inter_area <- if (!is.null(inter) && length(inter) > 0) 
      sum(as.numeric(st_area(st_make_valid(inter)))) else 0
    
    iso_area <- hr_area(hr_kde(sampled, levels = 0.95))$area
    
    # summary
    tibble(
      new_ID = id,
      sample_size = n,
      rep = rep_id,
      iso_area = iso_area,
      inter_area = inter_area,
      prop_overlap_full = ifelse(full_area > 0, inter_area / full_area, NA)
    )
  }
  
  # run for all sample sizes + reps
  all_metrics_mcp <- purrr::map_dfr(sample_sizes, function(n) {
    purrr::map_dfr(1:n_reps, function(rep_id) {
      downsample_mcp_iso(track_obj, n, rep_id)
    })
  })
  
  return(all_metrics_mcp)
}

#' running into some individuals with errors, so skipping those that are 
#' problematic
safe_run <- purrr::possibly(
  run_mcp_downsampling,
  otherwise = NULL,   # return NULL on error
  quiet = FALSE
)

#' applying to whole list
all_results_mcp <- purrr::map_dfr(
  track_list_mcp,
  safe_run,
  sample_sizes = sample_sizes,
  n_reps = n_reps
)

all_results_mcp %>% group_by(new_ID) %>% tally() 
beep(2)

# Split data by individual
data_list_mcp <- split(all_results_mcp, all_results_mcp$new_ID)

#' breakpoint analysis 
fit_segmented_breakpoint <- function(df) {
  # Fit basic linear model
  fit <- lm(iso_area ~ sample_size, data = df)
  
  # Segmented regression with starting breakpoint (psi = 15)
  seg_fit <- segmented(fit, seg.Z = ~sample_size, psi = 15)
  
  # Extract breakpoint and SE
  bp <- summary(seg_fit)$psi
  
  tibble(
    new_ID = unique(df$new_ID),
    breakpoint = bp[,"Est."],
    se = bp[,"St.Err"],
    response = "iso_area"
  )
}

#' running into some individuals with errors, so skipping those that are 
#' problematic
safe_run_2 <- purrr::possibly(
  fit_segmented_breakpoint,
  otherwise = NULL,   # return NULL on error
  quiet = FALSE
)

# Run segmented regression for iso_area
breakpoints_iso_mcp <- map_dfr(data_list_mcp, safe_run_2)

fit_segmented_breakpoint <- function(df) {
  # Fit basic linear model
  fit <- lm(prop_overlap_full ~ sample_size, data = df)
  
  # Segmented regression with starting breakpoint (psi = 15)
  seg_fit <- segmented(fit, seg.Z = ~sample_size, psi = 15)
  
  # Extract breakpoint and SE
  bp <- summary(seg_fit)$psi
  
  tibble(
    new_ID = unique(df$new_ID),
    breakpoint = bp[,"Est."],
    se = bp[,"St.Err"],
    response = "prop_overlap_full"
  )
}

#' running into some individuals with errors, so skipping those that are 
#' problematic
safe_run_2 <- purrr::possibly(
  fit_segmented_breakpoint,
  otherwise = NULL,   # return NULL on error
  quiet = FALSE
)

# Run segmented regression for prop_overlap_full
breakpoints_overlap_mcp <- map_dfr(data_list_mcp, safe_run_2)

# View results
breakpoints_iso_mcp
breakpoints_overlap_mcp

breakpoints_iso_mcp %>% ggplot(aes(x = breakpoint))+geom_histogram()
breakpoints_overlap_mcp %>% ggplot(aes(x = breakpoint))+geom_histogram()

summary(breakpoints_iso_mcp$breakpoint)
summary(breakpoints_overlap_mcp$breakpoint)

##### (4) Plotting ----
#' example from 1 individual
df = data_list_mcp[["Rollo2021-06-20"]]
fit <- lm(iso_area ~ sample_size, data = df)
# Segmented regression with starting breakpoint (psi = 15)
seg_fit <- segmented(fit, seg.Z = ~sample_size, psi = 15)

pdf("SIplotx8.pdf",
    width = 6, height = 6)
# Plot points
plot(df$sample_size, df$iso_area,
     pch = 19, col = "gray40",
     xlab = "Sample size",
     ylab = "Area")
# Add segmented regression fit
plot(seg_fit, add = TRUE, col = "red", lwd = 2)
# Optional: add breakpoint line
bp <- summary(seg_fit)$psi[ , "Est."]
abline(v = bp, col = "blue", lty = 2, lwd = 2)
dev.off()

df = data_list_mcp[["Rollo2021-06-20"]]
fit <- lm(prop_overlap_full ~ sample_size, data = df)
# Segmented regression with starting breakpoint (psi = 15)
seg_fit <- segmented(fit, seg.Z = ~sample_size, psi = 15)

pdf("SIplotx9.pdf",
    width = 6, height = 6)
# Plot points
plot(df$sample_size, df$prop_overlap_full,
     pch = 19, col = "gray40",
     xlab = "Sample size",
     ylab = "Proportion overlap")
# Add segmented regression fit
plot(seg_fit, add = TRUE, col = "red", lwd = 2)
# Optional: add breakpoint line
bp <- summary(seg_fit)$psi[ , "Est."]
abline(v = bp, col = "blue", lty = 2, lwd = 2)
dev.off()

