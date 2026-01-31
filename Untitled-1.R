source("./test_mgvc.R")

temp_2m <- result[,"id"]
temp_2m$temp <- as.matrix(result[, grepl("temp", names(result))])
temp_2m$lag <- as.matrix(result[, grepl("lag", names(result))])

points.temp <- merge(logistic.points,temp_2m,by=c("id")) %>%
  group_by(year) %>%
  mutate(id2=1:n())

tplot <- data.frame(points.temp$temp)
colnames(tplot) <- paste0(points.temp$lag[1,])
tplot$year <- points.temp$year
tplot$id2 <- points.temp$id2

tplot <- tplot %>%
  pivot_longer(-c(year, id2), names_to="lag", values_to="temperature") %>%
  left_join(points.temp[, c("id2", "obs_midpoint")], by=c("id2")) %>%
  mutate(lag = as.numeric(lag)) %>%
  mutate(ycoord = lag) %>%
  group_by(id2) %>%
  arrange(lag) %>%
  mutate(temperature = cumsum(ifelse(temperature < 5 | is.na(temperature),0,temperature)))


ggplot(points.temp, aes(x=id2, y=obs_midpoint)) +
  geom_tile(aes(x=id2, y=ycoord, fill=temperature), data=tplot) +
  geom_segment(aes(x=id2, yend=max(tplot$lag)+1, y=obs_midpoint), alpha=0.3) +
  geom_point(col="red") + 
  facet_wrap(~year) +
  ylim(150,210) +
  scale_fill_viridis_c() +
  scale_colour_gradient(low=muted("red"), high="orange") +
  labs(y="Day of year",
       fill="Cumulative\ntemperature (C)",
       colour="Jan-Feb average\ntemperature (C)") +
  theme_minimal() +
  theme(legend.position = "bottom")

