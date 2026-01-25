#' Fluctuations in direct human presence, not predictable weekly cycles, influence avoidance behaviour in ravens
#' Script by: Varalika Jain

#' The final processed datasets are available on the Zenodo data
#' repository associated with this study

#load("final_data_prep.RData")
#rm(list = setdiff(ls(), "ravens_5yr"))

####LOAD LIBRARIES----
library(rgdal)
library(sp)
library(terra)
library(move)
library(move2)
#library(lubridate)
library(dplyr)
#library(suncalc)
#library(sf)
#library(plyr)
library(tidyr)
library(ggplot2)
#library(glmmTMB)
#library(splines)
#library(performance)
#library(sjPlot)
#library(effects)
#library(ggeffects)
#library(gtsummary)
library(amt)
#library(emmeans)
library(scales)
library(purrr)
library(ggpubr)
library(recurse)

####DOWNLOAD DATA FROM MOVEBANK (MOVE OBJECT)----
#' Permissions required to download data
#' Input movebank login details
login <- movebankLogin() 

##### (1) 5 year period----
#' Source data from GPS tagged ravens from 1st Jan 2018 to 1st Jan 2023 
#' (duration of available ticket sales data)
ravens_5yr <- getMovebankData(study="Common Ravens in the Eastern Alps", 
                              login=login,
                              timestamp_start = "20180101000000000", 
                              timestamp_end = "20230101000000000",
                              removeDuplicatedTimestamps=TRUE)

#' Assigning the correct time zone
timestamps(ravens_5yr) <- lubridate::with_tz(timestamps(ravens_5yr), tz="Europe/Vienna")

#' Check the timezone
head(timestamps(ravens_5yr))

#' Check individuals
levels(ravens_5yr@trackId)

##### (2) 6 week period----
#' Source data from GPS tagged ravens from 2nd May 2023 to 3rd Jul 2023 
#' (duration of available hourly ticket sales data)
ravens_6wk <- getMovebankData(study="Common Ravens in the Eastern Alps", 
                              login=login,
                              timestamp_start = "20230522000000000", 
                              timestamp_end = "20230703000000000",
                              removeDuplicatedTimestamps=TRUE)

#' Assigning the correct time zone
timestamps(ravens_6wk) <- lubridate::with_tz(timestamps(ravens_6wk), tz="Europe/Vienna")

#' Check the timezone
head(timestamps(ravens_6wk))

#' Check individuals
levels(ravens_6wk@trackId)

####CONVERT MOVE OBJECT TO DATAFRAME----
ravens_5yr_df <- as(ravens_5yr, "data.frame")
ravens_6wk_df <- as(ravens_6wk, "data.frame")

####MODIFY DATAFRAME FOR REQUIRED INFO----
#####(1) Function for creating and filtering columns----
#' Check to see if timestamps are still in local time
head(ravens_5yr_df$timestamps)
head(ravens_6wk_df$timestamps)

#' Function
prep_df_func = function(data){
  data %>%
    mutate(date = as.Date(timestamps),
           time = format(timestamps, format = "%H:%M:%S"),
           weekday = factor(weekdays(date), levels = c("Monday", "Tuesday", "Wednesday", 
                                                       "Thursday", "Friday", "Saturday", "Sunday")),
           mmdd = format(as.Date(date), "%m-%d"),
           year = format(timestamps, format = "%Y"),
           local_identifier = as.factor(local_identifier)) %>%
    dplyr::select(c("local_identifier", "location_lat", 
                    "location_long", "timestamps", "date",
                    "time", "weekday", "mmdd", "year"))
}

##### (2) Apply prep function----
ravens_5yr_df <- prep_df_func(ravens_5yr_df)
ravens_6wk_df <- prep_df_func(ravens_6wk_df)
head(ravens_5yr_df)
head(ravens_6wk_df)

####SELECT DAYTIME POINTS ONLY----
#####(1) Function for daytime points---- 
#' To remove bias of night time points being higher in colder months
#' Assume AFS use during daytime

#' Function
daytime_func = function(data){
  #' Create a dataframe with timestamps and location for 'suncalc' package
  sun <- data.frame(date = as.Date(data$timestamps, tz = "Europe/Vienna"), 
                    lat = data$location_lat, lon = data$location_long)
  
  #' Caluclate sunrise and sunset times
  sunrise <- suncalc::getSunlightTimes(data = sun, keep="sunrise", tz = "Europe/Vienna") 
  sunset <- suncalc::getSunlightTimes(data = sun, keep="sunsetStart", tz = "Europe/Vienna") 
  
  #' Creating a new column for the move object, where fixes that have timestamps
  #' before sunrise and after sunset (i.e., night) should be marked as 1, else 0
  data$night_day <- ifelse(data$timestamps < sunrise$sunrise | 
                             data$timestamps > sunset$sunsetStart,
                           night <- 1, night <- 0)
  
  data <- data %>% filter(night_day == 0)
  return(data)
}

#####(2) Apply daytime function---- 
ravens_5yr_day <- daytime_func(ravens_5yr_df)
ravens_6wk_day <- daytime_func(ravens_6wk_df)

# Check
range(ravens_5yr_day$time)
range(ravens_6wk_day$time)

####STANDARDIZE SAMPLING FOR ALL DATA----
#####(1) Minimum 15 minute sampling rate; function----
#' In 2018, some birds had much higher sampling rates than in other years 
#' This was because some experiments were going on 
#' We select minimum of 15 minute sampling rate as that is what the 
#' loggers are set to otherwise (this is to standardize ALL the data)
#' first extract periods where sampling was very high

#' Calculate timestamp lags 
ravens_5yr_day <- ravens_5yr_day %>% group_by(local_identifier, date) %>%
  mutate(diff = as.numeric(timestamps - lag(timestamps), "mins"))
ravens_6wk_day <- ravens_6wk_day %>% group_by(local_identifier, date) %>%
  mutate(diff = as.numeric(timestamps - lag(timestamps), "mins"))

#' create id-date column
ravens_5yr_day$id_date <- paste0(ravens_5yr_day$local_identifier, ravens_5yr_day$date)
ravens_6wk_day$id_date <- paste0(ravens_6wk_day$local_identifier, ravens_6wk_day$date)

sampling_5yr <- ravens_5yr_day %>% group_by(id_date, date, mmdd) %>% tidyr::drop_na() %>%
  summarise(mean_interval = mean(diff))

sampling_6wk <- ravens_6wk_day %>% group_by(id_date, date, mmdd) %>% tidyr::drop_na() %>%
  summarise(mean_interval = mean(diff))
range(sampling_6wk$mean_interval) #15 min minimum sampling rate

quartz(height = 6, width = 10)
sampling_5yr %>% filter(mean_interval < 60) %>% 
  ggplot(aes(x = as.POSIXct(date), y = mean_interval)) + geom_point(alpha = 0.5, size = 0.2)+
  scale_x_datetime(labels = date_format("%m-%y"),
                   date_breaks = "month")+
  geom_hline(yintercept = 15, linetype = "dashed", colour = "red")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(y = "Sampling interval", x = "Date")
dev.copy2pdf(file="~/Documents/Ravens/weekend/240312 Final/revisions/Supplementary Information/SIplot1.pdf")

#' first set sampling to a minimum of 15 minutes
#' convert to move2 object
ravens_5yr_move2 =  mt_as_move2(ravens_5yr_day,
                                coords = c("location_long", "location_lat"), time_column = "timestamps",
                                track_id_column = "id_date") |> sf::st_set_crs(4326L)
ravens_5yr_15min = mt_filter_per_interval(ravens_5yr_move2, criterion = "first", unit = "15 mins")
#' note mt_filter per interval:  The time lag between the selected events does not 
#' necessarily correspond to the defined interval. For example, if the defined time interval 
#' is "1 hour" with the criterion "first", the function will select the event that is closest 
#' to every full hour, so if the first event of a track is at 10:45 and the second at 11:05, 
#' both of them will be selected, as they fall into different hour windows, but the time lag 
#' between them is 20 minutes.  

#' so we need to check lags 
ravens_5yr_15min$lags <- mt_time_lags(ravens_5yr_15min, units = "min")
hist(ravens_5yr_15min$lags)
#' still some points at very low sampling 
#' remove those points and recaculate lags (with 2 minute tolerance)
ravens_5yr_15min <- ravens_5yr_15min %>% filter(as.numeric(lags) >= 13)
#' recalculate lags
ravens_5yr_15min$lags <- mt_time_lags(ravens_5yr_15min, units = "min")

##### (2) Testing for autocorrelation (5yr data) ----
#running with 15 minute data 
ravens_5yr_15min$speeds <- mt_speed(ravens_5yr_15min)
ravens_5yr_15min$angles <- mt_turnangle(ravens_5yr_15min)
ravens_5yr_15min$dist_m <- mt_distance(ravens_5yr_15min)

# extracting individuals with consistent 15 min sampling in a day
ravens_5yr_15only <- ravens_5yr_15min %>% ungroup() %>%
  mutate(location_lat=st_coordinates(.)[,2], location_long=st_coordinates(.)[,1]) %>%
  st_drop_geometry()
ravens_5yr_15only <- as.data.frame(ravens_5yr_15only)
names(ravens_5yr_15only)

#' need to get rid of diff and angles columns because of NAs in the first row for each day
xx <- ravens_5yr_15only[,-c(9,13)] 
xx <- xx %>% 
  group_by(id_date) %>% 
  tidyr::drop_na() %>%
  summarise(min = min(lags), max = max(lags))
#' extracting individuals with 15 minutes of regular sampling within the day
xx <- xx %>% filter(as.numeric(min) < 17 & as.numeric(max) < 17)
ravens_5yr_15only <- ravens_5yr_15only %>% filter(id_date %in% xx$id_date)

#' double check lags again
xx <- ravens_5yr_15only %>% group_by(local_identifier, date) %>%
  mutate(diff = as.numeric(timestamps - lag(timestamps), "mins"))

ravens_5yr_15only %>% group_by(id_date) %>% tally()

# function to compute decorrelation time (first non-sig autocorr)
get_decor_time <- function(x, dt = 15, lag.max = 48) {
  # take metric of interest (speed / dist / andgles) and 
  # remove NA at start/end safely
  x <- as.numeric(x)
  x <- x[!is.na(x)] 
  if(length(x) < 10) return(NA_real_)  # not enough points
  # run acf function
  a <- acf(x, na.action = na.contiguous, lag.max = lag.max, plot = FALSE)$acf
  # skip lag0
  a <- a[-1]
  # 95% CI threshold for zero ACF
  crit <- qnorm(0.975) / sqrt(length(x))
  # find first lag where ACF ~ 0
  idx <- which(abs(a) < crit)[1]
  if(is.na(idx)) return(NA_real_)
  return(idx * dt)
}

# apply per id-date
acf_results_5yr <- ravens_5yr_15only %>%
  group_by(id_date) %>%
  summarise(
    lag_speed = get_decor_time(speeds),
    lag_dist  = get_decor_time(dist_m),
    lag_angle = get_decor_time(angles),
    .groups = "drop"
  )

#' two results with NA
acf_results_5yr <- acf_results_5yr %>% tidyr::drop_na()

#' investigate results
quartz(height = 4, width = 4)
acf_results_5yr %>% 
  tidyr::pivot_longer(-id_date, names_to = "metric", values_to = "lag_min") %>%
  ggplot(aes(metric, lag_min)) + geom_violin()+
  ylim(15, 150)+theme_minimal()+
  geom_hline(yintercept = 30, linetype = "dashed", color = "red", size = 0.5)+
  labs(x = "Metric", y = "Lag in minutes")
dev.copy2pdf(file="~/Documents/Ravens/weekend/240312 Final/revisions/Supplementary Information/SIplot1a.pdf")

summary(as.numeric(acf_results_5yr$lag_speed))
summary(as.numeric(acf_results_5yr$lag_angle))
summary(as.numeric(acf_results_5yr$lag_dist))

##### (3) ACF plots ---- 
get_acf_curve <- function(x, dt = 15, lag.max = 10) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if(length(x) < (lag.max + 1)) return(NULL)  # skip too-short series
  
  a <- acf(x, na.action = na.contiguous, lag.max = lag.max, plot = FALSE)$acf
  a <- a[-1]  # remove lag zero
  # pad if needed
  if(length(a) < lag.max) {
    a <- c(a, rep(NA, lag.max - length(a)))
  }
  
  tibble(
    lag = 1:lag.max,
    acf = a,
    dt = dt
  )
}

acf_curve_df1 <- ravens_5yr_15only %>%
  group_by(id_date) %>%
  summarise(
    acf_data = list(get_acf_curve(dist_m, dt = 15, lag.max = 10)),
    .groups = "drop"
  ) %>%
  filter(!sapply(acf_data, is.null)) %>%  # remove skipped groups
  unnest(acf_data)

acf_curve_df2 <- ravens_5yr_15only %>%
  group_by(id_date) %>%
  summarise(
    acf_data = list(get_acf_curve(angles, dt = 15, lag.max = 10)),
    .groups = "drop"
  ) %>%
  filter(!sapply(acf_data, is.null)) %>%  # remove skipped groups
  unnest(acf_data)

acf_curve_df3 <- ravens_5yr_15only %>%
  group_by(id_date) %>%
  summarise(
    acf_data = list(get_acf_curve(speeds, dt = 15, lag.max = 10)),
    .groups = "drop"
  ) %>%
  filter(!sapply(acf_data, is.null)) %>%  # remove skipped groups
  unnest(acf_data)

acf_summary1 <- acf_curve_df1 %>%
  group_by(lag) %>%
  summarise(
    mean = mean(acf, na.rm = TRUE),
    sd   = sd(acf, na.rm = TRUE),
    se   = sd / sqrt(n()),
    .groups = "drop"
  )

acf_summary2 <- acf_curve_df2 %>%
  group_by(lag) %>%
  summarise(
    mean = mean(acf, na.rm = TRUE),
    sd   = sd(acf, na.rm = TRUE),
    se   = sd / sqrt(n()),
    .groups = "drop"
  )

acf_summary3 <- acf_curve_df3 %>%
  group_by(lag) %>%
  summarise(
    mean = mean(acf, na.rm = TRUE),
    sd   = sd(acf, na.rm = TRUE),
    se   = sd / sqrt(n()),
    .groups = "drop"
  )

acf1 <- ggplot(acf_summary1, aes(x = lag, y = mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.6, size = 1.2) +
  geom_point(shape = 22, fill = "white", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red")+
  labs(x = "Lag at a 15-minute interval", y = "Step length autocorrelation",
       title = "a)") +
  theme_bw()+
  ylim(-0.5, 0.5)+
  scale_x_continuous(breaks = seq(0, 10, 1))

acf2 <- ggplot(acf_summary2, aes(x = lag, y = mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.6, size = 1.2) +
  geom_point(shape = 22, fill = "white", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red")+
  labs(x = "Lag at a 15-minute interval", y = "Turning angle autocorrelation",
       title = "b)") +
  theme_bw()+
  ylim(-0.5, 0.5)+
  scale_x_continuous(breaks = seq(0, 10, 1))

acf3 <- ggplot(acf_summary3, aes(x = lag, y = mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.6, size = 1.2) +
  geom_point(shape = 22, fill = "white", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red")+
  labs(x = "Lag at a 15-minute interval", y = "Speed autocorrelation",
       title = "c)") +
  theme_bw()+
  ylim(-0.5, 0.5)+
  scale_x_continuous(breaks = seq(0, 10, 1))

quartz(height = 6, width = 14)
ggarrange(acf1, acf2, acf3, ncol=3, nrow=1)
dev.copy2pdf(file="~/Documents/Ravens/weekend/240312 Final/revisions/Supplementary Information/SIplot2.pdf")

##### (4) ACF example for one individual ----
xx <- ravens_5yr_15only %>% group_by(local_identifier, date) %>%
  mutate(diff = as.numeric(timestamps - lag(timestamps), "mins"))

acf_example <- xx %>% filter(id_date=="Alfon2019-05-18") 
acf_example <- acf_example[c(-60), ]

pdf("~/Documents/Ravens/weekend/240312 Final/revisions/Supplementary Information/SIplot3.pdf",
    width = 12, height = 4)  # wide format for 3 columns
par(mfrow = c(1, 3),
    mar = c(4, 4, 2, 1))

xx <- acf_example %>% pull(as.numeric(dist_m))
acf(xx, na.action=na.contiguous, lag.max=48,
    xlab = "Lag interval",
    ylab = "Step length autocorrelation") # plot up to 48 lags (12 hrs)
mtext("a)", side = 3, adj = 0, line = 1)

xx <- acf_example %>% pull(as.numeric(angles))
acf(xx, na.action=na.contiguous, lag.max=48,
    xlab = "Lag interval",
    ylab = "Turning angle autocorrelation") # plot up to 48 lags (12 hrs)
mtext("b)", side = 3, adj = 0, line = 1)

xx <- acf_example %>% pull(as.numeric(speeds))
acf(xx, na.action=na.contiguous, lag.max=48,
    xlab = "Lag interval",
    ylab = "Speed autocorrelation") # plot up to 48 lags (12 hrs)
mtext("c)", side = 3, adj = 0, line = 1)
dev.off()

##### (5) Thinning all data to 30-minute intervals ----
#' back to the complete daytime dataset 
ravens_5yr_move2 =  mt_as_move2(ravens_5yr_day,
                                coords = c("location_long", "location_lat"), time_column = "timestamps",
                                track_id_column = "id_date") |> sf::st_set_crs(4326L)
ravens_6wk_move2 =  mt_as_move2(ravens_6wk_day,
                                coords = c("location_long", "location_lat"), time_column = "timestamps",
                                track_id_column = "id_date") |> sf::st_set_crs(4326L)

#' filtering at 30 minutes (same way as earlier)
ravens_5yr_30min = mt_filter_per_interval(ravens_5yr_move2, criterion = "first", unit = "30 mins")
ravens_5yr_30min$lags <- mt_time_lags(ravens_5yr_30min, units = "min")
ravens_5yr_30min <- ravens_5yr_30min %>% filter(as.numeric(lags) >= 25)
#' recalculate lags
ravens_5yr_30min$lags <- mt_time_lags(ravens_5yr_30min, units = "min")
#' add other metrics back
ravens_5yr_30min$speeds <- mt_speed(ravens_5yr_30min)
ravens_5yr_30min$angles <- mt_turnangle(ravens_5yr_30min)
ravens_5yr_30min$dist_m <- mt_distance(ravens_5yr_30min)

ravens_6wk_30min = mt_filter_per_interval(ravens_6wk_move2, criterion = "first", unit = "30 mins")
ravens_6wk_30min$lags <- mt_time_lags(ravens_6wk_30min, units = "min")
ravens_6wk_30min <- ravens_6wk_30min %>% filter(as.numeric(lags) >= 25)
#' recalculate lags
ravens_6wk_30min$lags <- mt_time_lags(ravens_6wk_30min, units = "min")
#' add other metrics back
ravens_6wk_30min$speeds <- mt_speed(ravens_6wk_30min)
ravens_6wk_30min$angles <- mt_turnangle(ravens_6wk_30min)
ravens_6wk_30min$dist_m <- mt_distance(ravens_6wk_30min)

##### (6) Convert back to dataframe ----
ravens_5yr_30min <- ravens_5yr_30min %>% ungroup() %>%
  mutate(location_lat=st_coordinates(.)[,2],location_long=st_coordinates(.)[,1]) %>%
  st_drop_geometry()
ravens_5yr_30min <- as.data.frame(ravens_5yr_30min)

ravens_6wk_30min <- ravens_6wk_30min %>% ungroup() %>%
  mutate(location_lat=st_coordinates(.)[,2],location_long=st_coordinates(.)[,1]) %>%
  st_drop_geometry()
ravens_6wk_30min <- as.data.frame(ravens_6wk_30min)

####BUFFER SIZE ESTIMATION ----
#' An individual was considered to be using the AFS if its GPS fix intersected 
#' with the buffer around the enclosure
#' doing the initial buffer size estimation with the 5yr dataset

##### (1) Convert raven data into UTM projection----
head(ravens_5yr_30min)

ravens_coords <- ravens_5yr_30min[c("location_long", "location_lat", "timestamps", "local_identifier")]
coordinates(ravens_coords) <- c("location_long", "location_lat")
proj4string(ravens_coords) <- CRS("+proj=longlat +datum=WGS84")

#' The units for the radius are those of the (x,y) coordinates 
#' (e.g., meters in the case of a UTM projection)
ravens_coords_utm <- spTransform(ravens_coords, CRS("+proj=utm +zone=33 ellps=WGS84"))
ravens_coords_utm <- as.data.frame(ravens_coords_utm)
ravens_coords_utm <- ravens_coords_utm[c(3,4,1,2)] #long, lat, timestamp, id
head(ravens_coords_utm)

##### (2) Convert into an sf object----
ravens_sf <- st_as_sf(x = ravens_coords_utm, 
                      coords = c("coords.x1", "coords.x2"),
                      crs = "+proj=utm +zone=33 +datum=WGS84")

##### (3) Read in the AFS location data----
#' AFS locations were identified using clusters of GPS fixes (from tracked ravens)
#' near obvious anthropogenic structures
afs <- read.csv('./input_df/AFSs.csv') #locs - locations

#' only interested in the game park
afs <- afs[1,]
head(afs)

##### (4) Select for columns: name, long, lat and type of AFS----
afs_latlong <- afs[c(1,3,2,6)]
head(afs_latlong)

##### (5) Convert location data into UTM projection----
afs_coords <- afs_latlong[c("FID", "Longitude", "Latitude", "Category")]
coordinates(afs_coords) <- c("Longitude", "Latitude")
proj4string(afs_coords) <- CRS("+proj=longlat +datum=WGS84")
afs_coords_utm <- spTransform(afs_coords, CRS("+proj=utm +zone=33 +datum=WGS84"))

##### (6) Convert location coordinates into dataframe----
afs_coords_df <- as.data.frame(afs_coords_utm@coords)
head(afs_coords_df)

##### (7) Using recurse to determine buffer/ radius size----
#' Calculate the number of revisits at AFSs at different radii
#' Same method as in Jain et al (2022) 
#' https://link.springer.com/article/10.1186/s40462-022-00335-4

afs_revisits = list()
for (i in 1:150) {
  afs_revisits$revisits[i] = getRecursionsAtLocations(x = ravens_coords_utm, 
                                                      locations = afs_coords_df,
                                                      radius = i,
                                                      threshold = 0) #threshold in hours
}
afs_revisits <- as.data.frame(afs_revisits$revisits)

##### (8) Plotting the revisits according to radius size ---- 
radii <- c(1:150) 
titles <- paste(afs_coords$FID)

pdf("~/Documents/Ravens/weekend/240312 Final/revisions/Supplementary Information/SIplot4.pdf",
    width = 6, height = 6)  # wide format for 3 columns
plot(radii,afs_revisits[1,],
     ylab = "Total number of revisits", xlab = "Radii")
abline(v = 80, col = "red", lty = 2, lwd = 2)   # vertical red dashed line at x = 80
dev.off()

#' This finds the smallest radius where the curve reaches a 
#' percentage of its max.
y <- afs_revisits[1,]
x <- radii
target <- 0.95 * max(y)  # 95% of plateau
idx <- which(y >= target)[1]
r_asymptote <- x[idx]
r_asymptote #80 is okay as a buffer

####GPS FIX & BUFFER INTERSECTION ----
##### (1) Function for projection---- 
project_func = function(data){
  #' duplicate data for manipulation into spatial dataframe
  data_sp <- data
  
  #' Reprojecting the data
  coordinates(data_sp) <-  c("location_long", "location_lat")
  proj4string(data_sp) <- CRS("+proj=longlat +datum=WGS84")
  data_proj <- spTransform(data_sp, CRS("EPSG:3416"))
  
  #'Convert to sf object
  data_sf <- sf::st_as_sf(x = data_proj, 
                          coords = c("location_long", "location_lat"),
                          crs = "EPSG:3416")
  return(data_sf)
} 

##### (2) Apply proj function to raven data----
ravens_5yr_proj <- project_func(ravens_5yr_30min)
ravens_6wk_proj <- project_func(ravens_6wk_30min)

#####(3) Read in the AFS location data again----
afs <- read.csv('./input_df/AFSs.csv') #locs - locations

#' only interested in the game park
afs <- afs[1,]
head(afs)

#####(4) Select for columns: name, long, lat and type of AFS----
afs_latlong <- afs[,c("FID", "Longitude", "Latitude", "Category")]
#' rename columns
colnames(afs_latlong)[colnames(afs_latlong) == "Longitude"] ="location_long"
colnames(afs_latlong)[colnames(afs_latlong) == "Latitude"] ="location_lat"

head(afs_latlong)

#####(5) Apply proj function to raven data----
afs_proj <- project_func(afs_latlong)

#####(6) Create intersection buffer----
#' In UTM, so dist is in meters
buf <- sf::st_buffer(afs_proj, dist = 80)

#####(7) Create intersection buffer----
#' 5 year data 
#' For each GPS fix, calculate whether or not it intersects with any of the bufs
#' #' If it does not intersect, in the intersection column, paste " ", else "AFS"
intersect_5yr_dat <- ravens_5yr_proj %>% mutate(
  intersection = as.integer(sf::st_intersects(geometry, buf))
  , location = if_else(is.na(intersection), " ", paste0("AFS"))) 

#' 6 weeks data 
intersect_6wk_dat <- ravens_6wk_proj %>% mutate(
  intersection = as.integer(sf::st_intersects(geometry, buf))
  , location = if_else(is.na(intersection), " ", paste0("AFS"))) 

#####(8) Convert to dataframe ----
intersect_5yr_df <- intersect_5yr_dat %>% 
  mutate(location_lat=st_coordinates(.)[,2],
         location_long=st_coordinates(.)[,1]) %>%
  st_drop_geometry()

intersect_6wk_df <- intersect_6wk_dat %>% 
  mutate(location_lat=st_coordinates(.)[,2],
         location_long=st_coordinates(.)[,1]) %>%
  st_drop_geometry()

intersect_5yr_df <- as.data.frame(intersect_5yr_df)
intersect_6wk_df <- as.data.frame(intersect_6wk_df)

####CALCULATE INDV PT DISTANCES TO WB SITE----
#####(1) Convert to sp object----
#' Difference in distance in meters stored in a vector
dist_5yr <- ravens_5yr_30min[,c("local_identifier", "timestamps", 
                                "location_lat", "location_long")]
coordinates(dist_5yr) <-  c("location_long", "location_lat")
proj4string(dist_5yr) <- CRS("+proj=longlat +datum=WGS84")

dist_6wk <- ravens_6wk_30min[,c("local_identifier", "timestamps", 
                                "location_lat", "location_long")]
coordinates(dist_6wk) <-  c("location_long", "location_lat")
proj4string(dist_6wk) <- CRS("+proj=longlat +datum=WGS84")

afs_deg <- afs_latlong
coordinates(afs_deg) <-  c("location_long", "location_lat")
proj4string(afs_deg) <- CRS("+proj=longlat +datum=WGS84")

#####(2) Apply distVincentyEllipsoid----
#' Difference in distance in meters stored in a vector
pt_dist_5yr <- geosphere::distVincentyEllipsoid(dist_5yr, afs_deg,
                                                a=6378137, b=6356752.3142, f=1/298.257223563)
pt_dist_6wk <- geosphere::distVincentyEllipsoid(dist_6wk, afs_deg,
                                                a=6378137, b=6356752.3142, f=1/298.257223563)
#####(3) Add to df----
intersect_5yr_df$wb_dist <- pt_dist_5yr
intersect_6wk_df$wb_dist <- pt_dist_6wk

#### SAFETY BACK-UP ----
save.image("revision_data_prep_part_1.RData")
rm(ravens_5yr) #also saved separately as an RData
# "~/Documents/Ravens/weekend/240312 Final/ravens_5yr.RData"
rm(acf_curve_df1, acf_curve_df2, acf_curve_df3,
   acf_example, acf_results_5yr, acf_summary1,
   acf_summary2, acf_summary3, acf1, acf2, acf3)
rm(ravens_5yr_15min, ravens_5yr_15only, ravens_5yr_day,
   ravens_6wk_day)
gc()

####MERGE DATAFRAMES WITH DEMOGRAPHIC DATA----
#####(1) Add location long and location lat columns to intersect data----
str(ravens_5yr_30min)
str(intersect_5yr_df)

#' checking for differences between the dfs
vars = c("local_identifier", "timestamps")
setdiff(ravens_5yr_30min[,vars], intersect_5yr_df[,vars])
setdiff(ravens_6wk_30min[,vars], intersect_6wk_df[,vars])

#' no difference so can add lat long columns to intersect_df
intersect_5yr_df$location_lat <- ravens_5yr_30min$location_lat
intersect_5yr_df$location_long <- ravens_5yr_30min$location_long

intersect_6wk_df$location_lat <- ravens_6wk_30min$location_lat
intersect_6wk_df$location_long <- ravens_6wk_30min$location_long

#####(2) Read in the demography datafile----
demog <- read.csv('./input_df/GPS Tagging information 2021.csv') 
head(demog)

demog<- demog[!duplicated(demog$Name), ]

demog <- demog[,c("Name", "ID.ring", "Sex", "Origin", "Year.hatch", "Tracking.status", "GPS.year")]

#' id name as factor
class(demog$Name)
demog$id_name <- as.factor(demog$Name)

#####(3) Checking to see if all the individual names match----
levels(intersect_5yr_df$local_identifier)
levels(demog$id_name)
setdiff(demog$id_name, intersect_5yr_df$local_identifier)

#' Check which birds are to be excluded, and which ones are not included in tracking period
demog %>% filter(grepl("excluded", Tracking.status))
demog %>% filter(grepl("GPS23", GPS.year))
demog %>% filter(grepl("GPS17", GPS.year))

#' Change the names in the demog file as it's smaller & faster to do
demog$id_name <- plyr::revalue(demog$id_name, c("Lisa 2" = "Lisa"))

#####(4) Remove excluded indv----
excl <- demog %>% filter(!grepl("excluded", Tracking.status))
demog <- demog %>% filter(id_name %in% excl$id_name)

#####(5) Merging the dataframe with demography data----
#' 5 year dataframe
names(demog)
names(intersect_5yr_df)
demog_5yr_intersect <- merge.data.frame(demog, intersect_5yr_df, by.x = c("id_name"),
                                        by.y = c('local_identifier'),
                                        all.x = FALSE, sort = TRUE)

demog_5yr_intersect <- droplevels(demog_5yr_intersect)
levels(demog_5yr_intersect$id_name)
names(demog_5yr_intersect)

#' 6 week dataframe
demog_6wk_intersect <- merge.data.frame(demog, intersect_6wk_df, by.x = c("id_name"),
                                        by.y = c('local_identifier'),
                                        all.x = FALSE, sort = TRUE)

demog_6wk_intersect <- droplevels(demog_6wk_intersect)
levels(demog_6wk_intersect$id_name)
#' rumble removed from dataset so becomes smaller. 
#' rumble has no sex so it's okay to leave out.

#####(6) Add age and other columns function----
add_cols_func = function(data){
  data <- data %>%
    mutate(
      date = as.Date(date),
      month = format(date, "%m"),
      doy = lubridate::yday(date),
      hour = format(timestamps, "%H"),
      Year.hatch = plyr::revalue(Year.hatch, 
                           c("<2007" = "2006", "<2009" = "2008","<2014"= "2013", "<2015"= "2014",
                             "<2020"= "2019", "2015/16"= "2016")),
      bday = as.Date(paste0(Year.hatch,"-03-21")),
      age = lubridate::time_length(difftime(date, bday), "years"),
      age_class = ifelse(age <1, "juvenile", 
                         ifelse(age>= 1 & age < 3, "sub", "adult")))
}

#####(7) Apply function----
demog_5yr_intersect <- add_cols_func(demog_5yr_intersect)
demog_6wk_intersect <- add_cols_func(demog_6wk_intersect)

#' check age class
plot(as.factor(demog_5yr_intersect$age_class), demog_5yr_intersect$age)
plot(as.factor(demog_6wk_intersect$age_class), demog_6wk_intersect$age)

#####(8) Export as csv----
#' Export
write.csv(demog_5yr_intersect, "demog_5yr_intersect.csv")
write.csv(demog_6wk_intersect, "demog_6wk_intersect.csv")

####TOTAL NUMBER OF INDIVIDUALS; TOTAL NUMBER OF FIXES; IN-OUT COLS----
#####(1) 5YR: Daily totals----
#' Find total number of individuals tagged in a day
indv_tot <- demog_5yr_intersect %>% droplevels() %>%
  group_by(id_name, date) %>% tally() %>%
  group_by(date) %>% tally()
#' rename column
names(indv_tot)[names(indv_tot) == "n"] <- "indv_tot"
head(indv_tot)

#' Total number of fixes per indv and day
fix_tot <- demog_5yr_intersect %>%
  group_by(id_name, date) %>% tally()
#' rename column
names(fix_tot)[names(fix_tot) == "n"] <- "fix_tot"
head(fix_tot)

#' merge the two
daily_tot <- merge(fix_tot, indv_tot, by = c("date"), all.x = TRUE) 

#' Merge with 5yr df
final_5yr <- merge(demog_5yr_intersect, daily_tot, by = c("id_name","date"), all.x = TRUE) 

#####(2) 6WK: Hourly totals----
#' Add hour column
demog_6wk_intersect <- demog_6wk_intersect %>% mutate(hour = format(as.POSIXct(timestamps), "%H"))

#' Sampling rate of movement data is too coarse for hourly resolution
#' All ticket sales from 9AM to 1PM = AM 
#' All ticket sales from 1PM to 5PM = PM
#' convert to bi-hourly
demog_6wk_intersect$bi_hourly <- plyr::revalue(as.character(demog_6wk_intersect$hour), 
                                              c("07" = "7",
                                                "08" = "7",
                                                "09" = "9",
                                                "10" = "9", 
                                                "11" = "11",
                                                "12" = "11",
                                                "13" = "13",
                                                "14" = "13",
                                                "15" = "15",
                                                "16" = "15",
                                                "17" = "17",
                                                "18" = "17"))

#' Find total number of individuals tagged per hour per day
indv_tot_bihrly <- demog_6wk_intersect %>%
  group_by(id_name, date, bi_hourly) %>% tally() %>%
  group_by(date, bi_hourly) %>% tally()
#' rename column
names(indv_tot_bihrly)[names(indv_tot_bihrly) == "n"] <- "indv_tot_bihrly"
head(indv_tot_bihrly)

#' Total number of fixes per indv and hour
fix_tot_bihrly <- demog_6wk_intersect %>%
  group_by(id_name, date, bi_hourly) %>% tally()
#' rename column
names(fix_tot_bihrly)[names(fix_tot_bihrly) == "n"] <- "fix_tot_bihrly"
head(fix_tot_bihrly)

#' merge the two
tot_bihrly <- merge(fix_tot_bihrly, indv_tot_bihrly, by = c("date", "bi_hourly"), all.x = TRUE) 

#' Merge with 6wk df
final_6wk <- merge(demog_6wk_intersect, tot_bihrly, by = c("id_name","date", "bi_hourly"), all.x = TRUE) 

#####(3) GPS fixes in versus out----
#' out and in columns
final_5yr$location <- ifelse(is.na(final_5yr$intersection), "out", "in")
final_6wk$location <- ifelse(is.na(final_6wk$intersection), "out", "in")

#####(4) Export as csv----
write.csv(final_5yr, "final_5yr.csv")
write.csv(final_6wk, "final_6wk.csv")

####ASSIGN WEEKDAY-WEEKEND & SEASON----
#####(1) Weekend function----
#' Weekends = saturdays and sundays 
#' Weekdays = all other days

weekend_func = function(data){
  data <- data %>%
    mutate(type_of_day = case_when(weekday == 'Sunday' | weekday == 'Saturday' ~ 'weekend',
                                   weekday == 'Monday' | weekday == 'Tuesday' | 
                                     weekday == "Thursday" | weekday == "Friday" |
                                     weekday == "Wednesday" ~ 'weekday'))
}

#####(2) Apply function----
final_5yr <- weekend_func(final_5yr)
final_6wk <- weekend_func(final_6wk)

#####(3) Season function----
assign_season <- function(date) {
  # ensure date is Date class
  date <- as.Date(date)
  year <- as.integer(format(date, "%Y"))
  
  # define astronomical cutoffs for each year
  spring_start  <- as.Date(paste0(year, "-03-21"))
  summer_start  <- as.Date(paste0(year, "-06-21"))
  fall_start    <- as.Date(paste0(year, "-09-23"))
  winter_start  <- as.Date(paste0(year, "-12-21"))
  
  # winter extends into next year:
  next_winter_end <- as.Date(paste0(year + 1, "-03-19"))
  
  ifelse(date >= spring_start & date < summer_start, "Spring",
         ifelse(date >= summer_start & date < fall_start,   "Summer",
                ifelse(date >= fall_start   & date < winter_start, "Autumn",
                       "Winter")))
}

#####(4) Apply function----
final_5yr$season <- assign_season(final_5yr$date)
final_6wk$season <- assign_season(final_6wk$date)

####REMOVE BIRDS FAR FROM SITE----
hist(final_5yr$wb_dist)
names(final_5yr)

##### (1) Filter birds that were at the site, how far did they travel ----
max_dist_park = final_5yr %>% group_by(id_date, mmdd, year,date) %>% 
  summarise(min = min(wb_dist),
            max = max(wb_dist)) %>%
  filter(min <= 80)

range(max_dist_park$max)
hist(max_dist_park$max)
summary(max_dist_park$max)
quantile(x=max_dist_park$max, prob=0.99) #23km
xx <- max_dist_park %>% filter(max > 80) #25km
summary(xx$max)
quantile(x=xx$max, prob=0.99) #at 99th percentile of distances above 80m, 25km 

quantile(x=max_dist_park$max, prob=0.999) #23km


##### (2) plotting ---- 
dist_test <- ravens_5yr_30min %>% filter(id_date %in% max_dist_park$id_date)

# Convert to sf points
dist_sf <- dist_test %>%
  st_as_sf(coords = c("location_long", "location_lat"), crs = 4326)

# Group by trajectory and convert to LINESTRING
tracks_sf <- dist_sf %>%
  group_by(id_date) %>%
  summarise(geometry = st_combine(geometry)) %>% 
  st_cast("LINESTRING")

st_write(tracks_sf, "raven_max_dist_tracks.kml", driver = "KML")

#' Visually, I think a 50km buffer is meaningful - it encompasses sites outside 
#' of the park, while covering this network of alternative AFSs where informaiton
#' exchange can happen

##### (3) extracting individuals not visiting the park in a day ----
xx = final_5yr %>% filter(intersection == 1) %>%
  group_by(id_date) %>% tally()

no_visit_5yr <- final_5yr %>% filter(!id_date %in% xx$id_date) %>%
  group_by(id_date, mmdd, year,date) %>% 
  summarise(min = min(wb_dist),
            max = max(wb_dist))

xx = final_6wk %>% filter(intersection == 1) %>%
  group_by(id_date) %>% tally()

no_visit_6wk <- final_6wk %>% filter(!id_date %in% xx$id_date) %>%
  group_by(id_date, mmdd, year,date) %>% 
  summarise(min = min(wb_dist),
            max = max(wb_dist))

##### (4) remove them if all of their movement points exceed 50km away from park ----
no_visit_5yr <- no_visit_5yr %>% filter(min > 50000)
final_5yr_filtered <- final_5yr %>% filter(!id_date %in% no_visit_5yr$id_date)

no_visit_6wk <- no_visit_6wk %>% filter(min > 50000)
final_6wk_filtered <- final_6wk %>% filter(!id_date %in% no_visit_6wk$id_date)

##### (5) Sanity check ----
xx = final_5yr_filtered %>% group_by(id_date, intersection) %>% 
  summarise(min = min(wb_dist),
            max = max(wb_dist))

xx = final_6wk_filtered %>% group_by(id_date, intersection) %>% 
  summarise(min = min(wb_dist),
            max = max(wb_dist))

####CALCULATE METRIC 1: PROBABILITY OF BIRDS VISITING SITE----
#####(1) 5YR: Daily no. of GPS fixes in versus out of site----
names(final_5yr_filtered)

prob_daily_5yr <- final_5yr_filtered %>%
  group_by(id_name, date, Sex, Origin, weekday,
           type_of_day, month, year, location, doy, age_class,
           fix_tot) %>% 
  tally(name = "number_pts") %>%
  tidyr::drop_na() %>%
  ungroup() %>%
  dplyr::select(id_name, date, Sex, Origin, weekday, type_of_day,
         month, year, doy, age_class, fix_tot, location, number_pts) %>%
  tidyr::pivot_wider(
    names_from = location,
    values_from = number_pts,
    values_fill = 0   # <-- ensures individuals with only inside or outside remain
  ) %>%
  rename(
    'outside' = 'out',
    'inside' = 'in'
  ) %>%
  mutate(
    prob = inside / fix_tot
  )

#####(2) 6WK: Bi-hourly no. of GPS fixes in versus out of site----
names(final_6wk_filtered)
prob_bihourly_6wk <- final_6wk_filtered %>%
  group_by(id_name, date, bi_hourly, Sex, Origin, weekday,
           type_of_day,
           year, location, month, doy, age_class,
           fix_tot_bihrly) %>% 
  tally(name = "number_pts") %>%
  tidyr::drop_na() %>%
  ungroup() %>%
  dplyr::select(id_name, date, bi_hourly, Sex, Origin, weekday,
                type_of_day,
                year, location, month, doy, age_class,
                fix_tot_bihrly, number_pts) %>%
  tidyr::pivot_wider(
    names_from = location,
    values_from = number_pts,
    values_fill = 0   # <-- ensures individuals with only inside or outside remain
  ) %>%
  rename(
    'outside' = 'out',
    'inside' = 'in'
  ) %>%
  mutate(
    prob = inside / fix_tot_bihrly
  )


#####(3) Export as csv----
write.csv(prob_daily_5yr, "prob_daily_5yr.csv")
write.csv(prob_bihourly_6wk, "prob_bihourly_6wk.csv")

####CALCULATE METRIC 2: KDE DAILY ----
#' Need to decide a minimum number of points in a day for KDE estimation
#' see script 'WE data prep KDE'
prob_daily_5yr 

yy <- final_5yr_filtered %>% group_by(id_date, date, mmdd, year) %>% tally() %>%
  filter(n>14)
xx <- prob_daily_5yr %>% mutate(id_date = paste0(id_name,date)) %>%
  filter(fix_tot > 14) 

setdiff(as.character(yy$id_date), as.character(xx$id_date))

final_kde <- final_5yr_filtered %>% filter(id_date %in% yy$id_date)
final_kde <- final_kde %>% arrange(id_date, timestamps)

#####----(1) Get projection per day ----#####
# Convert your day's points to an amt track object
track_day <- make_track(
  final_kde,
  .id = id_date,
  .x = location_long,
  .y = location_lat,
  .t = timestamps,
  crs = "+proj=longlat +datum=WGS84"
)

# Split full track object by individual
track_list <- split(track_day, track_day$.id)
levels(as.factor(final_kde$id_date))

sf::sf_use_s2(FALSE)
kde_5yr <- map(track_list, ~ hr_kde(.x, levels = 0.95))

xx = track_day %>% filter(.id =="Napoleon2022-03-19")

xx = hr_kde(xx, levels = 0.95)
plot(xx)

#####----(2) Extract area----####
kde_5yr_area_df <- lapply(kde_5yr, function(x) {
  tryCatch(
    hr_area(x),
    error = function(e) NA
  )
})

kde_5yr_area_full <- do.call(rbind, kde_5yr_area_df)

kde_5yr_area_df <- tibble(
  id = names(kde_5yr_area_df),
  area_m2 = kde_5yr_area_full$area
  )

kde_5yr_area_df <- as.data.frame(kde_5yr_area_df)

quartz(height = 4, width = 5)
xx = merge(max_dist_park, kde_5yr_area_df, by.x = "id_date", by.y = "id")
ggplot(xx, aes(max/1000, area_m2/1000000))+geom_point(alpha = 0.5, size = 2)+
  theme_minimal(base_size = 14)+
  labs(y= expression("KDE area" ~(km^2)),
       x = "Maximum displacement from park (km)")+
  stat_smooth()
dev.copy2pdf(file="~/Documents/Ravens/weekend/240312 Final/revisions/SIplotx14.pdf")


####----COMBINE MOVEMENT METRICS----####
head(prob_daily_5yr)
head(prob_bihourly_6wk)
head(kde_5yr_area_df)

prob_daily_5yr$id_date <- paste0(prob_daily_5yr$id_name, prob_daily_5yr$date)
prob_bihourly_6wk$id_date <- paste0(prob_bihourly_6wk$id_name, prob_bihourly_6wk$date)
prob_bihourly_6wk$id.bi_hour <- paste0(prob_bihourly_6wk$id_date,"-", prob_bihourly_6wk$bi_hourly)
prob_daily_5yr$id_date
kde_5yr_area_df$id
prob_bihourly_6wk$id.bi_hour

prob_daily_5yr <- merge(prob_daily_5yr, kde_5yr_area_df,
                        by.x = "id_date",
                        by.y = "id",
                        all.x = TRUE)

####FINAL 5YR WEEKLY DATASET----
#####(1) Read in ticket sales data----
ticket_sales<- read.csv("./input_df/Ticket sales.csv")
head(ticket_sales)
ticket_sales$date <- as.Date(ticket_sales$date, format = "%d/%m/%Y")

#' Select only columns needed
ticket_sales <- ticket_sales[,c(1:2)]
range(ticket_sales$date)
range(prob_daily_5yr$date)
prob_daily_5yr$date <- as.Date(prob_daily_5yr$date)

#####(2) Add stringency index data----
strin_idx <- read.csv("./input_df/owid-covid-data.csv", stringsAsFactors = T)
str(strin_idx)
strin_idx <- strin_idx %>% filter(location == "Austria") %>%
  dplyr::select(c(date, stringency_index))
#' after 2022-12-30 assume SI over
strin_idx[is.na(strin_idx)] <- 0

#' convert to date
strin_idx$date <- as.Date(strin_idx$date)

#####(3) Merge ticket sales, stringency index and dataset----
extra_cols <- merge(ticket_sales, strin_idx, all.x = TRUE, by = "date")
apply(is.na(extra_cols), FUN = sum, MARGIN = 2)
extra_cols[is.na(extra_cols)] <- 0

prob_final_5yr <- merge(prob_daily_5yr, extra_cols, by = "date", all.x = TRUE)
prob_final_5yr <- ungroup(prob_final_5yr)

####FINAL 6WK BIHOURLY DATASET----
#####(1) Read in hourly sales data----
hourly_sales <- read.csv("./input_df/add_hourly_sales.csv")

hourly_sales <- hourly_sales %>% 
  mutate(
    date = as.Date(date, format = "%d/%m/%Y")
  )


#' Hourly sales bi-hour
hourly_sales$bi_hourly <- plyr::revalue(as.character(hourly_sales$hour_start), 
                                          c("7" = "7",
                                            "8" = "7",
                                            "9" = "9",
                                            "10" = "9", 
                                            "11" = "11",
                                            "12" = "11",
                                            "13" = "13",
                                            "14" = "13",
                                            "15" = "15",
                                            "16" = "15",
                                            "17" = "17",
                                            "18" = "17"))

#' Calculate AM-PM visitor data 
visitors_bi <- hourly_sales %>% group_by(date, bi_hourly, visitors) %>%
  tally() %>%
  group_by(date, bi_hourly) %>%
  dplyr::summarise(visitors_original = sum(visitors)) %>%
  mutate(bi_hourly = as.numeric(bi_hourly)) %>% arrange(date, bi_hourly) %>%
  group_by(date) %>%
  mutate(visitors_shifted = lag(visitors_original, default = 0)) %>%
  mutate(visitors_cumulative = visitors_shifted + visitors_original)


#####(2) Merge hourly sales with prob_daily_6wk----
prob_final_6wk <- merge(prob_bihourly_6wk, visitors_bi,
                        by.x = c("date", "bi_hourly"),
                        by.y = c("date", "bi_hourly"),
                        all.x = TRUE)


#' hourly range of ticket sales is only 
#' from 7am to 6pm. need to filter out other vals
levels(as.factor(prob_final_6wk$bi_hourly))
prob_final_6wk <- prob_final_6wk %>% filter(bi_hourly != "05" &
                                              bi_hourly != "06" &
                                              bi_hourly != "19" &
                                              bi_hourly != "20" &
                                              bi_hourly != "21") %>%droplevels()

levels(as.factor(prob_final_6wk$bi_hourly))
prob_final_6wk$bi_hourly <- factor(prob_final_6wk$bi_hourly,
                                 levels = c("7", "9", "11", "13", "15", "17"),
                                 ordered = TRUE)
#levels(as.factor(hourly_sales$bi_hourly))

#' remove all rows with less than 3 total fixes in both dataframes
#prob_final_5yr <- prob_final_5yr %>% filter(fix_tot > 3)
#prob_final_6wk <- prob_final_6wk %>% filter(fix_tot_bihrly >3)

####EXPORT FINAL DATASETS----
write.csv(prob_final_5yr, "prob_final_5yr.csv")
write.csv(prob_final_6wk, "prob_final_6wk.csv")

####SAVE IMAGE----
#save.image("REVISION_final_data_prep.RData")
#load('REVISION_final_data_prep.RData')
