## 
##
## January 2026
##
## 
library(tidyverse)
library(dplyr)
library(ggplot2)
library(diptest)
library(gridExtra)
library(BioSIM)
library(fitdistrplus)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(reldist)

###########################################
### Read flight data from Marc's files
###########################################
# Identify record source based on ID last digit
#   0: Original files 1-7 from Marc
#   1: Maritimes light trap catches
#   2: Last batch of files recovered by Marc
# Records were selected based on a set of criteria
# Criteria #1: at least 30 moths captured in total
# Criteria #2: at least one survey record of absence before emergence
# According to Marc: Among the site*year that Criteria #2 would eliminate, several are very close to 0
# or have 0 but it wasn't entered. Based on his information, Marc made the following list of site*year
# to eliminate
# Marc.Eliminate <- c("4599667572","4668604372","4544634673","4614659173","4614659175","4487673677",
#               "4665683977","4708689378","4742673478","4589699779")
# Criteria #3: more than 10 survey dates
# Criteria #4: No more than 30% of the total sample size was trapped in one single day 
# After Discussion with Marc in March 2017, it was decided to increase the limit of max sumday to 50% of
# population caught but the duration of the flight has to be > 6 days to limit the risk of having
# migration flights 

Obs <- read.csv2("Observations.csv")
# Relocate a trap that appears in the ocean
Obs <- Obs |>
  mutate(lat = if_else(lat == 46.38, 46.40, lat))

##################################
# Figure 1: Map of trap locations
##################################
Traps.location <- Obs %>%
  dplyr::select(year,lat,long) %>%
  distinct()

# Get provinces and states
canada <- ne_states(country = "canada", returnclass = "sf")
usa <- ne_states(country = "united states of america", returnclass = "sf")

# Filter for Eastern regions
eastern_canada <- canada |> filter(longitude > -85 & longitude < -55)
eastern_usa <- usa |> filter(longitude > -85 & longitude < -55)

# Get centroids for labels
canada_labels <- eastern_canada |> 
  st_centroid() |> 
  st_coordinates() |> 
  as.data.frame() |> 
  bind_cols(name = eastern_canada$name)

usa_labels <- eastern_usa |> 
  st_centroid() |> 
  st_coordinates() |> 
  as.data.frame() |> 
  bind_cols(name = eastern_usa$name)

# Create map
ggplot() +
  geom_sf(data = eastern_canada, fill = "lightgray", color = "white") +
  geom_sf(data = eastern_usa, fill = "lightgray", color = "white") +
  geom_text(data = canada_labels, aes(x = X, y = Y, label = name), size = 3.5) +
  geom_text(data = usa_labels, aes(x = X, y = Y, label = name), size = 3.5) +
  geom_point(data = Traps.location, 
             aes(x = long, y = lat),
             color = "red", size = 2, alpha = 0.7) +
  coord_sf(xlim = c(-72, -62), ylim = c(44, 48)) +
  facet_wrap(~year) +
  labs(
    title = "Trap Locations by Year",
    x = "Longitude",
    y = "Latitude"
  ) +
  annotation_scale(location = "br", width_hint = 0.2) +
  theme_bw()

########################################################
# Figure Annexe 1: Distribution of capture in each trap
########################################################

pdf(file="ObservationPlots.pdf")
# Plot flight dynamics
Plot.obs <- function(i) {
  Obs.dat <- Obs[Obs$id == i,]
  tmp <- rep(Obs.dat$date,Obs.dat$sumday)
  Obs.dat$median <- wtd.quantile(tmp,q=0.5)
  g <- ggplot(Obs.dat,aes(x=date,y=sumday)) + 
    geom_col(fill="steelblue") +
    geom_bar(aes(x=median,y=max(sumday)),position=position_dodge(),stat="identity",colour="red",width=.1,fill="red") +
    labs(
      title =paste("Lat:",substr(Obs.dat$id,1,2),".",substr(Obs.dat$id,3,4),sep=""," Long:-",substr(Obs.dat$id,5,6),".",
                  substr(Obs.dat$id,7,8)," Year:19",substr(Obs.dat$id,9,10)),
      x = "Date (DOY)",
      y = "Moths captured"
    ) +
    theme(plot.title=element_text(size=7))

  return(g)
}

ids <- unique(Obs$id)
for (l in 1:17) {
  for (k in 1:9) {
    i <- (l-1)*9 + k
    g <- Plot.obs(ids[i])    
    assign(paste("p",k,sep=""),g)
  }
  grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,nrow=3)
}
dev.off()
rm(g,list=ls(pattern="p[:123456789:]"))

##############################################################
# Figure 2: Distribution of 5%, 95%, median and duration
##############################################################
# Calculate DOY of 5% flight for observations
Obs.summary <- Obs %>% 
  group_by(id) %>% 
  summarize(
    first_5pct.obs = date[which(sumday.cumsum > 0.05)[1]],
    last_95pct.obs = date[max(which(sumday.cumsum < 0.95))],
    median.obs = date[max(which(sumday.cumsum < 0.50))],
    duration.obs = last_95pct.obs - first_5pct.obs,
    .groups = "drop"
  ) %>%
  mutate(year = paste0("19",substr(id,9,10)),
         lat = as.integer(substr(id,1,4))/100,
         lon = as.integer(substr(id,5,8))/100*-1)

p1 <- ggplot(Obs.summary, aes(x = first_5pct.obs)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", color = "black") +
  geom_density(color = "red", linewidth = 1) +
  labs(
    title = "Distribution of First 5% Flight Date",
    x = "Day of Year",
    y = "Density"
  )
p2 <- ggplot(Obs.summary, aes(x = last_95pct.obs)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", color = "black") +
  geom_density(color = "red", linewidth = 1) +
  labs(
    title = "Distribution of 95% Flight Date",
    x = "Day of Year",
    y = "Density"
  )
p3 <- ggplot(Obs.summary, aes(x = median.obs)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", color = "black") +
  geom_density(color = "red", linewidth = 1) +
  labs(
    title = "Distribution of median Flight Date",
    x = "Day of Year",
    y = "Density"
  )
p4 <- ggplot(Obs.summary, aes(x = duration.obs)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", color = "black") +
  geom_density(color = "red", linewidth = 1) +
  labs(
    title = "Distribution of Flight Duration",
    x = "Days",
    y = "Density"
  )
grid.arrange(p1,p2,p3,p4,nrow=2)

######################################################################
# Exploring the distribution of flight dates: normality and bimodality
######################################################################

# Test normality and bimodality for each trap*year
normality_bimodality_tests <- Obs %>%
  summarize(
    n = n(),
    shapiro_statistic = shapiro.test(sumday)$statistic,
    shapiro_p_value = shapiro.test(sumday)$p.value,
    is_normal = shapiro.test(sumday)$p.value > 0.05,
    dip_statistic = dip.test(sumday)$statistic,
    dip_p_value = dip.test(sumday)$p.value,
    is_multimodal = dip.test(sumday)$p.value < 0.05,

    .by = id
  )

normality_bimodality_tests

# Add transformed versions of sumday
Obs_transformed <- Obs |>
  mutate(
    sumday_log = log(sumday + 1),        # log transformation (add 1 to handle zeros)
    sumday_sqrt = sqrt(sumday),           # square root transformation
    sumday_cuberoot = sumday^(1/3),      # cube root transformation
    sumday_reciprocal = 1/(sumday + 1)   # reciprocal transformation
  )

# Test normality of transformations for each group
transformation_tests <- Obs_transformed |>
  summarize(
    n = n(),
    # Original
    original_p = shapiro.test(sumday)$p.value > 0.05,
    # Log transform
    log_p = shapiro.test(sumday_log)$p.value > 0.05,
    # Square root
    sqrt_p = shapiro.test(sumday_sqrt)$p.value > 0.05,
    # Cube root
    cuberoot_p = shapiro.test(sumday_cuberoot)$p.value > 0.05,
    # Reciprocal
    reciprocal_p = shapiro.test(sumday_reciprocal)$p.value > 0.05,
    .by = id
  )

transformation_tests

# log(sumday+1) seems to be the transformation that gives the most 
# normalized series
transformation_tests |>
  summarize(across(where(is.logical), sum))

# Add log_sumday
Obs <- Obs %>% mutate(log_sumday = log(sumday+1))

# Test normality and bimodality on log(sumday+1)
normality_bimodality_tests <- Obs |>
  summarize(
    n = n(),
    shapiro_statistic = shapiro.test(log_sumday)$statistic,
    shapiro_p_value = shapiro.test(log_sumday)$p.value,
    is_normal = shapiro.test(log_sumday)$p.value > 0.05,
    dip_statistic = dip.test(log_sumday)$statistic,
    dip_p_value = dip.test(log_sumday)$p.value,
    is_multimodal = dip.test(log_sumday)$p.value < 0.05,

    .by = id
  )

normality_bimodality_tests

####################################################################################
# Run BioSIM
####################################################################################
# ids <- unique(Obs$id)
# Biosim.pred <- do.call("rbind",lapply(ids,function(i){
#   year <- as.integer(paste0("19",substr(i,9,10)))
#   lat <- as.integer(substr(i,1,4))/100
#   lon <- as.integer(substr(i,5,8))/100*-1
#   generateWeather(modelNames = "Spruce_Budworm_Biology",
#                   fromYr = year,
#                   toYr = year,
#                   id = i,
#                   latDeg = lat,
#                   longDeg = lon,
#                 )$Spruce_Budworm_Biology
# }))
# write.csv2(Biosim.pred,file="Biosim.pred.csv")

Biosim.pred <- read.csv2(file="Biosim.pred.csv")

####################################################################################
# Make BioSIM predictions 
# 
####################################################################################
Biosim.simple <- Biosim.pred %>%
  mutate(DOY = as.integer(format(as.Date(paste(Year, Month, Day, sep = "-")), "%j")),
         Adult.Flight = MaleFlight+FemaleFlight) %>%
  group_by(KeyID) %>%
  mutate(Adult.Flight.Percent = Adult.Flight/(sum(Adult.Flight))) %>%
  ungroup() %>%
  dplyr::select(KeyID,Year,DOY,Adult.Flight,Adult.Flight.Percent)

Obs.Adults.sum <- Obs %>% 
  summarise(SumAdult=sum(sumday),.by=id)

Biosim.pred.2 <- left_join(Obs.Adults.sum,Biosim.simple,by=c("id"="KeyID")) %>%
  mutate(Biosim.Flight=Adult.Flight.Percent*SumAdult) %>%
  dplyr::select(id,Year,DOY,Biosim.Flight)

Obs.with.Biosim <- full_join(Obs,Biosim.pred.2,by=c("id"="id","date"="DOY")) %>%
  dplyr::select(year,lat,long,id,date,sumday,Biosim.Flight) %>%
  group_by(id) %>%
  arrange(date,.by_group = TRUE) %>%
  mutate(sumday=replace_na(sumday,0)) %>%
  dplyr::filter(date>=160 & date <=250)

# Plot comparison between Observations and Biosim
pdf(file="Obs_vs_Biosim.pdf",paper="USr")
# Plot flight dynamics
ids <- unique(Obs$id)
test <- lapply(ids, function(i) {
  Obs.dat <- Obs.with.Biosim[Obs.with.Biosim$id == i,]
  ggplot(Obs.dat,aes(x=date,y=sumday)) + 
    geom_col(fill="steelblue") +
    geom_line(aes(y=Biosim.Flight),color="red") +
    xlim(160,250) +
    labs(
      title =paste("Lat:",substr(Obs.dat$id,1,2),".",substr(Obs.dat$id,3,4),sep=""," Long:-",substr(Obs.dat$id,5,6),".",
                  substr(Obs.dat$id,7,8)," Year:19",substr(Obs.dat$id,9,10)),
      x = "Date (DOY)",
      y = "Moths"
    ) +
    theme(plot.title=element_text(size=7))
})
dev.off()

# Calculate the difference between Obs and Biosim 
Diff <- Obs.with.Biosim %>%
  ungroup() %>%
  mutate(Diff = abs(sumday - Biosim.Flight)) %>%
  summarise(
    Approval.rate = 1-(sum(Diff)/(2*sum(sumday))),
    .by=id
  ) %>%
  mutate(year = as.factor(paste0("19",substr(id,9,10))),
  lat = as.integer(substr(id,1,4))/100 ,
  lon = as.integer(substr(id,5,8))/100*-1
  )

# Histogram of differences between observations and Biosim
ggplot(Diff,aes(x=Approval.rate)) + geom_histogram()

# Differences vs lat, long and Year
p1 <- ggplot(Diff,aes(x=lat,y=Approval.rate)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue")

p2 <- ggplot(Diff,aes(x=lon,y=Approval.rate)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue")

p3 <- ggplot(Diff, aes(x = year, y = Approval.rate)) +
  geom_violin(fill = "steelblue", alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(
    title = "Approval Rate Distribution by Year",
    x = "Year",
    y = "Approval Rate"
  )

grid.arrange(p1,p2,p3,nrow=2)

# Map of the worst predictions
Diff.worst <- Diff %>% dplyr::filter(Approval.rate<0.1)
ggplot() +
  geom_sf(data = eastern_canada, fill = "lightgray", color = "white") +
  geom_sf(data = eastern_usa, fill = "lightgray", color = "white") +
  geom_text(data = canada_labels, aes(x = X, y = Y, label = name), size = 3.5) +
  geom_text(data = usa_labels, aes(x = X, y = Y, label = name), size = 3.5) +
  geom_point(data = Diff.worst, 
             aes(x = lon, y = lat),
             color = "red", size = 2, alpha = 0.7) +
  coord_sf(xlim = c(-72, -62), ylim = c(44, 48)) +
  facet_wrap(~year) +
  labs(
    title = "Trap Locations by Year",
    x = "Longitude",
    y = "Latitude"
  ) +
  annotation_scale(location = "br", width_hint = 0.2) +
  theme_bw()

# Combine Observations and Biosim predictions
Obs.with.Biosim <- Obs.with.Biosim %>% 
  dplyr::select(id,date,sumday,Biosim.Flight) %>%
  mutate(year = as.factor(paste0("19",substr(id,9,10))),
  lat = as.integer(substr(id,1,4))/100 ,
  lon = as.integer(substr(id,5,8))/100*-1) %>%
  ungroup()

write.csv2(a,file="ObsVsPred.csv",row.names=FALSE)

Flight.cumsum <- Obs.with.Biosim %>%
  group_by(id) %>%
  mutate(
    Obs.cumsum = cumsum(sumday)/sum(sumday),
    Biosim.cumsum = cumsum(Biosim.Flight)/sum(Biosim.Flight)
  ) %>%
ungroup()

# Plot comparison between Observations and Biosim
# Plot flight dynamics
ids <- unique(Obs$id)
test.2 <- lapply(ids, function(i) {
  Obs.dat <- Flight.cumsum[Flight.cumsum$id == i,]
  ggplot(Obs.dat,aes(x=date)) + 
    geom_line(aes(y=Obs.cumsum),color="blue") +
    geom_line(aes(y=Biosim.cumsum),color="red") +
    xlim(170,225) +
    labs(
      title =paste("Lat:",substr(Obs.dat$id,1,2),".",substr(Obs.dat$id,3,4),sep=""," Long:-",substr(Obs.dat$id,5,6),".",
                  substr(Obs.dat$id,7,8)," Year:19",substr(Obs.dat$id,9,10)),
      x = "Date (DOY)",
      y = "Moths"
    ) +
    theme(plot.title=element_text(size=7))
})

xx <- marrangeGrob(test.2,ncol=5,nrow=5)
ggsave(file="Obs_vs_Biosim_cumsum.pdf",width=11,height=7,unit="in",xx)

###############################
# Fitting logistic regressions
###############################
# Fit logistic curves for each id
logistic_fits <- Flight.cumsum |>
  group_by(id) |>
  summarize(
    # Fit logistic to Obs.cumsum
    obs_fit = list(tryCatch(
      nls(Obs.cumsum ~ SSlogis(date, Asym, xmid, scal)),
      error = function(e) NULL
    )),
    # Fit logistic to Biosim.cumsum
    biosim_fit = list(tryCatch(
      nls(Biosim.cumsum ~ SSlogis(date, Asym, xmid, scal)),
      error = function(e) NULL
    )),
    .groups = "drop"
  ) |>
  mutate(
    # Extract parameters for Obs
    obs_params = map(obs_fit, ~if(!is.null(.x)) tidy(.x) else NULL),
    # Extract parameters for Biosim
    biosim_params = map(biosim_fit, ~if(!is.null(.x)) tidy(.x) else NULL)
  )

# Extract the midpoint (xmid) parameter for each curve
midpoints <- logistic_fits |>
  mutate(
    obs_midpoint = map_dbl(obs_params, ~if(!is.null(.x)) .x$estimate[.x$term == "xmid"] else NA_real_),
    biosim_midpoint = map_dbl(biosim_params, ~if(!is.null(.x)) .x$estimate[.x$term == "xmid"] else NA_real_)
  ) |>
  dplyr::select(id, obs_midpoint, biosim_midpoint)

midpoints

logistic_plots <- lapply(ids, function(i) {
  # Get data for this id
  plot_data <- Flight.cumsum[Flight.cumsum$id == i, ]
  
  # Get fitted parameters
  params <- logistic_fits[logistic_fits$id == i, ]
  
  # Create base plot with both observations and biosim as points
  p <- ggplot(plot_data, aes(x = date)) + 
    geom_point(aes(y = Obs.cumsum), color = "blue", alpha = 0.6, size = 1.5) +
    geom_point(aes(y = Biosim.cumsum), color = "green", alpha = 0.6, size = 1.5) +
    xlim(170, 225) +
    ylim(0, 1) +
    labs(
      title = paste("Lat:", substr(i, 1, 2), ".", substr(i, 3, 4), 
                    " Long:-", substr(i, 5, 6), ".", substr(i, 7, 8), 
                    " Year:19", substr(i, 9, 10)),
      x = "Date (DOY)",
      y = "Cumulative Proportion"
    ) +
    theme_bw() +
    theme(plot.title = element_text(size = 7))
  
  # Add obs logistic curve if fit converged
  if (!is.null(params$obs_fit[[1]])) {
    pred_dates <- seq(170, 225, length.out = 200)
    pred_values <- predict(params$obs_fit[[1]], newdata = data.frame(date = pred_dates))
    pred_df <- data.frame(date = pred_dates, pred = pred_values)
    p <- p + geom_line(data = pred_df, aes(x = date, y = pred), color = "darkblue", linewidth = 1)
  }
  
  # Add Biosim logistic curve if fit converged
  if (!is.null(params$biosim_fit[[1]])) {
    pred_dates <- seq(170, 225, length.out = 200)
    pred_values <- predict(params$biosim_fit[[1]], newdata = data.frame(date = pred_dates))
    pred_df <- data.frame(date = pred_dates, pred = pred_values)
    p <- p + geom_line(data = pred_df, aes(x = date, y = pred), color = "red", linewidth = 1)
  }
  
  p
})

# Arrange in grid and save
zz <- marrangeGrob(logistic_plots, ncol = 5, nrow = 5)
ggsave("Logistic_fits_both.pdf", zz, width = 11, height = 8.5, units = "in")

# Add midpoint difference, lat, long and year
midpoints <- midpoints %>% 
  mutate(ObsMinusBio=obs_midpoint-biosim_midpoint,
    year = as.factor(paste0("19",substr(id,9,10))),
    lat = as.integer(substr(id,1,4))/100 ,
    lon = as.integer(substr(id,5,8))/100*-1)

ggplot(midpoints,aes(x=ObsMinusBio)) + geom_histogram()

midpoints.worst <- midpoints %>%
  dplyr::filter(ObsMinusBio < -7 | ObsMinusBio > 7)

# Map of the worst predictions
ggplot() +
  geom_sf(data = eastern_canada, fill = "lightgray", color = "white") +
  geom_sf(data = eastern_usa, fill = "lightgray", color = "white") +
  geom_text(data = canada_labels, aes(x = X, y = Y, label = name), size = 3.5) +
  geom_text(data = usa_labels, aes(x = X, y = Y, label = name), size = 3.5) +
  geom_point(data = midpoints.worst, 
             aes(x = lon, y = lat),
             color = "red", size = 2, alpha = 0.7) +
  coord_sf(xlim = c(-72, -62), ylim = c(44, 48)) +
  facet_wrap(~year) +
  labs(
    title = "Trap Locations by Year",
    x = "Longitude",
    y = "Latitude"
  ) +
  annotation_scale(location = "br", width_hint = 0.2) +
  theme_bw()

summary(lm(ObsMinusBio~lon+lat,data=midpoints))

####################################################################################
# Part 1: Compare 5%, 25%, 50%, 75%, 95% and duration of observations vs Biosim
####################################################################################

# Calculate DOY of 5% flight for Biosim
Flight5.bio <- Biosim.pred %>% 
  group_by(KeyID) %>%
  mutate(
    TotalFlight = MaleFlight + FemaleFlight,
    TotalFlightCumsum = cumsum(TotalFlight)/sum(TotalFlight),
    DOY = as.integer(format(as.Date(paste(Year, Month, Day, sep = "-")), "%j"))
  ) %>%
  summarize(
    x5pct.bio = DOY[which(TotalFlightCumsum > 0.05)[1]],
    x25pct.bio = DOY[which(TotalFlightCumsum > 0.25)[1]],
    x50pct.bio = DOY[which(TotalFlightCumsum > 0.50)[1]],
    x75pct.bio = DOY[which(TotalFlightCumsum > 0.75)[1]],
    x95pct.bio = DOY[max(which(TotalFlightCumsum < 0.95))],
    duration.bio = x95pct.bio - x5pct.bio,
    .groups = "drop"
  )
 
# Calculate DOY of 5% flight for observations
Obs5.bio <- Obs %>% 
  group_by(id) %>% 
  summarize(
    x5pct.obs = date[which(sumday.cumsum > 0.05)[1]],
    x25pct.obs = date[which(sumday.cumsum > 0.05)[1]],
    x50pct.obs = date[which(sumday.cumsum > 0.05)[1]],
    x75pct.obs = date[which(sumday.cumsum > 0.05)[1]],
    x95pct.obs = date[max(which(sumday.cumsum < 0.95))],
    duration.obs = x95pct.obs - x5pct.obs,
    .groups = "drop"
  )
  
# Combine results
Compare.pred <- right_join(Obs5.bio,Flight5.bio,by=c("id"="KeyID")) %>%
  mutate(ObsBioDiff.first = x5pct.obs-x5pct.bio,
         ObsBioDiff.last = x95pct.obs-x95pct.bio,
         ObsBioDiff.duration = duration.obs-duration.bio) %>%
  mutate(year = as.integer(paste0("19",substr(as.character(id),9,10))),
         lat = as.integer(substr(id,1,4))/100,
         lon = as.integer(substr(id,5,8))/100*-1 )

p1 <- ggplot(Compare.pred, aes(x = ObsBioDiff.first)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", color = "black") +
  stat_function(
    fun = dnorm,
    args = list(mean = mean(Compare.pred$ObsBioDiff.first), sd = sd(Compare.pred$ObsBioDiff.first)),
    color = "red",
    linewidth = 1
  ) +
  xlim(-25,25) +
  geom_vline(aes(xintercept = mean(ObsBioDiff.first)), color = "darkred", linewidth = 1, linetype = "dashed") +
  labs(
    title = "Differences Observed vs Biosim 5%",
    x = "Difference (days)",
    y = "Density"
  )

p2 <- ggplot(Compare.pred, aes(x = ObsBioDiff.last)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", color = "black") +
  stat_function(
    fun = dnorm,
    args = list(mean = mean(Compare.pred$ObsBioDiff.last), sd = sd(Compare.pred$ObsBioDiff.last)),
    color = "red",
    linewidth = 1
  ) +
  xlim(-25,25) +
  geom_vline(aes(xintercept = mean(ObsBioDiff.last)), color = "darkred", linewidth = 1, linetype = "dashed") +
  labs(
    title = "Differences Observed vs Biosim 95%",
    x = "Difference (days)",
    y = "Density"
  )

p3 <- ggplot(Compare.pred, aes(x = duration.obs)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", color = "black") +
  stat_function(
    fun = dnorm,
    args = list(mean = mean(Compare.pred$duration.obs), sd = sd(Compare.pred$duration.obs)),
    color = "red",
    linewidth = 1
  ) +
  xlim(0,30) +
  geom_vline(aes(xintercept = mean(duration.obs)), color = "darkred", linewidth = 1, linetype = "dashed") +
  labs(
    title = "Duration between 5% and 95% in observed",
    x = "Duration (days)",
    y = "Density"
  )

p4 <- ggplot(Compare.pred, aes(x = duration.bio)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", color = "black") +
  stat_function(
    fun = dnorm,
    args = list(mean = mean(Compare.pred$duration.bio), sd = sd(Compare.pred$duration.bio)),
    color = "red",
    linewidth = 1
  ) +
  xlim(0,30) +
  geom_vline(aes(xintercept = mean(duration.bio)), color = "darkred", linewidth = 1, linetype = "dashed") +
  labs(
    title = "Duration between 5% and 95% in Biosim",
    x = "Duration (days)",
    y = "Density"
  )

grid.arrange(p1,p3,p2,p4)

p5 <- ggplot(Compare.pred,aes(x=x5pct.obs,y=x5pct.bio)) + 
  geom_point(aes(color=as.factor(year))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "5% flight date Observations vs Biosim",
    x = "Observations",
    y = "Biosim") +
  xlim(170,230) + ylim(160,240)


p25 <- ggplot(Compare.pred,aes(x=x25pct.obs,y=x25pct.bio)) + 
  geom_point(aes(color=as.factor(year))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "25% flight date Observations vs Biosim",
    x = "Observations",
    y = "Biosim") +
  xlim(170,230) + ylim(160,240)


p50 <- ggplot(Compare.pred,aes(x=x50pct.obs,y=x50pct.bio)) + 
  geom_point(aes(color=as.factor(year))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "50% flight date Observations vs Biosim",
    x = "Observations",
    y = "Biosim") +
  xlim(170,230) + ylim(160,240)


p75 <- ggplot(Compare.pred,aes(x=x75pct.obs,y=x75pct.bio)) + 
  geom_point(aes(color=as.factor(year))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "75% flight date Observations vs Biosim",
    x = "Observations",
    y = "Biosim") +
  xlim(170,230) + ylim(160,240)

p95 <- ggplot(Compare.pred,aes(x=x95pct.obs,y=x95pct.bio)) + 
  geom_point(aes(color=as.factor(year))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "95% flight date Observations vs Biosim",
    x = "Observations",
    y = "Biosim") +
  xlim(170,230) + ylim(160,240)

pduration <- ggplot(Compare.pred,aes(x=duration.obs,y=duration.bio)) + 
  geom_point(aes(color=as.factor(year))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Flight duration Observations vs Biosim",
    x = "Observations",
    y = "Biosim")


grid.arrange(p5,p25,p50,p75,p95)

ggplot(Compare.pred,aes(x=last_95pct.obs,y=last_95pct.bio)) + 
  geom_point(aes(color=as.factor(year))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "95% flight date Observations vs Biosim",
    x = "Observations",
    y = "Biosim")

ggplot(Compare.pred,aes(x=lat,y=ObsBioDiff.first)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") 

ggplot(Compare.pred,aes(x=lon,y=ObsBioDiff.first)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") 

#################################################################################
# The comparisons between observations and Biosim are difficult to interpret
# Let's go back to visualizing observations and Biosim predictions for each trap
#################################################################################

#### 5%
shapiro.test(Compare.pred$ObsBioDiff.first)
fit.first <- fitdist(Compare.pred$ObsBioDiff.first,"norm")
fit.first.mean <- fit.first$estimate["mean"]
fit.first.sd <- fit.first$estimate["sd"]
qnorm(0.05,mean=fit.first.mean,sd=fit.first.sd)
qnorm(0.95,mean=fit.first.mean,sd=fit.first.sd)
qnorm(0.25,mean=fit.first.mean,sd=fit.first.sd)
qnorm(0.75,mean=fit.first.mean,sd=fit.first.sd)

#### 95%
shapiro.test(Compare.pred$ObsBioDiff.last)
fit.last <- fitdist(Compare.pred$ObsBioDiff.last,"norm")
fit.last.mean <- fit.last$estimate["mean"]
fit.last.sd <- fit.last$estimate["sd"]
qnorm(0.05,mean=fit.last.mean,sd=fit.last.sd)
qnorm(0.95,mean=fit.last.mean,sd=fit.last.sd)
qnorm(0.25,mean=fit.last.mean,sd=fit.last.sd)
qnorm(0.75,mean=fit.last.mean,sd=fit.last.sd)


####################################################################################
# Part 2: Download weather for all trap*year
####################################################################################
library(openmeteo)

# Obs.T <- do.call("rbind",lapply(ids,function(i) {
#   year <- paste0("19",substr(i,9,10))
#   lat <- as.integer(substr(i,1,4))/100
#   lon <- as.integer(substr(i,5,8))/100*-1
#   weather_history(location = c(lat,lon),
#     start = paste0(year,"-01-01"),
#     end = paste0(year,"-12-31"),
#     hourly="temperature_2m") %>% mutate(id=i)
# }))
# write.csv2(Obs.T,file="2m_temperature_at_sites.csv")

Obs.T <- read.csv2("2m_temperature_at_sites.csv")

# Calulate DD over 5C
Obs.DD <- Obs.T %>% 
  mutate(DOY=yday(datetime)) %>%
  dplyr::filter(DOY < 182) %>%
  summarise(
    T=mean(hourly_temperature_2m,na.rm=TRUE),
    .by=c(id,DOY)
  ) %>% 
  summarise(
    DD=sum(ifelse(T>5,T,0)),
    .by=c(id)
  )
  
midpoints.DD <- left_join(midpoints,Obs.DD,by="id")

ggplot(midpoints.DD,aes(x=DD,y=ObsMinusBio)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue")

summary(lm(ObsMinusBio~DD,data=midpoints.DD))
  summarise(
    TN=min(hourly_temperature_2m,na.rm=TRUE),
    T=mean(hourly_temperature_2m,na.rm=TRUE),
    TX=max(hourly_temperature_2m,na.rm=TRUE),
    .by=c(id,DOY)
  )

##
## Option 2: Use Biosim to get the weather data
##
ids <- unique(Obs$id)
Biosim.T <- do.call("rbind",lapply(ids,function(i){
   year <- as.integer(paste0("19",substr(i,9,10)))
   lat <- as.integer(substr(i,1,4))/100
   lon <- as.integer(substr(i,5,8))/100*-1
   generateWeather(modelNames = "TminTairTmax_Daily",
                   fromYr = year,
                   toYr = year,
                   id = i,
                   latDeg = lat,
                   longDeg = lon
                 )$TminTairTmax_Daily
 }))
write.csv2(Biosim.T,file="Biosim.T.csv")

Biosim.DD <- Biosim.T %>% 
  mutate(
    DOY = as.integer(format(as.Date(paste(Year, Month, Day, sep = "-")), "%j")),
    id = KeyID
  ) %>%
  dplyr::filter(DOY < 182) %>%
  summarise(
    DD.Biosim=sum(ifelse(Tair>5,Tair,0)),
    .by=c(id)
  )

midpoints.DD <- left_join(midpoints.DD,Biosim.DD,by="id")

ggplot(midpoints.DD,aes(x=DD.Biosim,y=ObsMinusBio)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue")





library(SBWPheno)
Obs.T.summary
test <- Obs.T.summary[Obs.T.summary$id == ids[1],]
help(Regniere.devage)

relrate <- Regniere.devvar(n=100,gc="field")
sex <- rep("m",100)
stage <- rep("L2o",100)
temp <- test[,c(3,5)]
res <- Regniere.devage(temp=temp,relrate=relrate,sex=sex,period=24,start.stage=stage,gc="field")

generateWeather(modelNames = "DegreeDay_Daily",
                   fromYr = 1980,
                   toYr = 1980,
                   id = 999,
                   latDeg = 43.87,
                   longDeg = -70.02,
                   LowerThreshold=5
                 )
