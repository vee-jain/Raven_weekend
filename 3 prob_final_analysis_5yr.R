#' Fluctuations in direct human presence, not predictable weekly cycles, influence avoidance behaviour in ravens
#' Script by: Varalika Jain

#' The final processed datasets are available on the Zenodo data
#' repository associated with this study
#' This code relies on custom functions currently not provided with the code - 
#' please contact me for more information


#### (1) Prepare working space----
rm(list = ls())
dev.off()

library(dplyr)
library(ggplot2)
library(purrr)
library(forcats)

#### (2) Load data ----
xdata = read.csv(file="prob_final_5yr.csv", stringsAsFactors=T)

#### (3) Plots ----
thresholds <- c(3, 4, 6, 15)

# compute counts per age class at each threshold
plot.data <- map_df(thresholds, ~{
  xdata %>%
    filter(fix_tot >= .x) %>%
    group_by(id_date, age_class) %>%
    tally() %>%
    group_by(age_class) %>%
    tally(name = "n_days") %>%
    mutate(threshold = paste0("fix_tot ≥ ", .x))
}) %>%
  mutate(
    threshold = fct_relevel(threshold,
                            "fix_tot ≥ 3",
                            "fix_tot ≥ 4",
                            "fix_tot ≥ 6",
                            "fix_tot ≥ 15")
  )

# compute percent lost relative to baseline (fix_tot ≥ 3)
baseline <- plot.data %>%
  filter(threshold == "fix_tot ≥ 3") %>%
  select(age_class, baseline_n = n_days)

plot.data <- plot.data %>%
  left_join(baseline, by = "age_class") %>%
  mutate(
    pct_lost = 100 * (1 - n_days / baseline_n),
    pct_lab  = paste0(round(pct_lost, 1), "%")
  )

quartz(height = 5, width = 6)

ggplot(plot.data, aes(x = age_class, y = n_days, fill = threshold)) +
  geom_col(position = position_identity(), alpha = 0.5) +
  geom_text(
    aes(label = pct_lab),
    position = position_identity(),
    vjust = +1.3,
    size = 3.4
  ) +
  labs(
    x = "Age class",
    y = "Number of individual-date points",
    fill = "Threshold"
  ) +
  theme_minimal(base_size = 14)
dev.copy2pdf(file="SIplotx5.pdf")

xdata$weekday <- factor(xdata$weekday,
                        levels = c("Monday", 
                                   "Tuesday", 
                                   "Wednesday", 
                                   "Thursday", 
                                   "Friday", 
                                   "Saturday", 
                                   "Sunday"))

cb_cols <- c("weekday" = "#7F7F7F",   # grey
             "weekend" = "#E69F00")   # orange

quartz(height = 5, width = 7)
ggplot(xdata, aes(x = weekday, y = ticket_sales, fill = type_of_day)) +
  geom_boxplot() +
  scale_fill_manual(values = cb_cols, name = "Type of day") +
  theme_minimal() +
  theme(
    legend.position = "right",
  )+
  labs(x ="Day of week",
       y = "Ticket sales",
       title = "b)")
dev.copy2pdf(file="SIplotx10.pdf")

quartz(height = 5, width = 3)
ggplot(xdata, aes(x = type_of_day, y = ticket_sales, fill = type_of_day)) +
  geom_boxplot() +
  scale_fill_manual(values = cb_cols) +
  theme_minimal() +
  theme(
    legend.position = "none",
  )+
  labs(x =" Type of day",
       y = "Ticket sales",
       title = "a)")
dev.copy2pdf(file="SIplotx11.pdf")

#### (4) check NA and str ----
#' are there any NAs
any(is.na(xdata))
str(xdata)

# drop area column
xdata <- xdata %>% dplyr::select(-area_m2)
any(is.na(xdata))

#### (5) dataframe touch-ups ----
#' select a minimum of 6 fixes in a day (covers whole day in winter)
#' with hourly sampling
xdata <- xdata %>% filter(fix_tot >= 6)

#' season as radians
xdata$date = as.Date(xdata$date)
xdata$date.rad=2*pi*as.numeric(xdata$date)/365.25
xdata$doy.rad=2*pi*as.numeric(xdata$doy)/365.25

#### (6) look at predictors ----
hist(xdata$ticket_sales)
hist(xdata$stringency_index)

#' range of ticket sales wrt to type of day
plot(xdata$type_of_day, xdata$ticket_sales) 

#' range of ticket sales wrt to origin
plot(xdata$Origin, xdata$ticket_sales) 

#' check for each indv. only 1 sex
xx = table(xdata$id_name, xdata$Sex)
#' check if each row has a value >0
range(apply(X=xx>0, # if larger than 0 = TRUE, else FALSE
            MARGIN = 1, #for each row (don't want 1 indv with both M & F)
            FUN= sum))
apply(X=xx>0, # if larger than 0 = TRUE, else FALSE
      MARGIN = 2, #for each col
      FUN= sum)

#' check for each indv. only 1 origin
xx = table(xdata$id_name, xdata$Origin)
#' check if each row has a value >0
range(apply(X=xx>0, # if larger than 0 = TRUE, else FALSE
            MARGIN = 1, #for each row (don't want 1 indv with both origin and captive)
            FUN= sum))
apply(X=xx>0, # if larger than 0 = TRUE, else FALSE
      MARGIN = 2, #for each col
      FUN= sum)

#' check for each indv. can have multiple age class
xx = table(xdata$id_name, xdata$age_class)
#' check if each row has a value >0
range(apply(X=xx>0, # if larger than 0 = TRUE, else FALSE
            MARGIN = 1, 
            FUN= sum))
apply(X=xx>0, # if larger than 0 = TRUE, else FALSE
      MARGIN = 2, #for each col
      FUN= sum)

#### (7) look at the response ----
hist(xdata$prob)
range(xdata$prob)

#### (8) grouping factors----
#' look at grouping in id_name and date 
#' should be only 1 daily measure per indv per date
xx = aggregate(x=1:nrow(xdata), by = xdata[, c("id_name", "date")], FUN= length)
head(xx)
plot(table(xx$x))
table(xx$x) # for 35828 combinations of indv and date, 1 prob measure
xdata$date = as.Date(xdata$date)

#### (9) checking random slopes ----
source("diagnostic_fcns.r")
xx.re.tab = fe.re.tab(fe.model = "prob~ type_of_day*ticket_sales+
                      Origin*ticket_sales+
                      stringency_index*type_of_day+
                      stringency_index*Origin+
                      Sex + 
                      age_class + 
                      date.rad+
                      doy.rad+
                      doy",
                      re = "(1|id_name)+(1|date)",
                      data = xdata, other.vars=c("inside", "outside", "fix_tot"))

xx=xx.re.tab$detailed[13]
head(xx)
xx.re.tab$summary[17] 

#### (10) setting up model ----
str(xx.re.tab$data)
t.data=xx.re.tab$data

#' z-transform the covariates
range(t.data$ticket_sales)
t.data$z.ticket_sales=as.vector(scale(t.data$ticket_sales))
range(t.data$stringency_index)
t.data$z.stringency_index=as.vector(scale(t.data$stringency_index))

#' response
hist(t.data$prob)
range(t.data$prob)
t.data$tr.prob = (t.data$prob*(length(t.data$prob)-1)+0.5)/length(t.data$prob)
#' squeezes everything closer to 0.5 
#' use when you have zero or 1 in the response, 
#' the beta distr cannot cope with it,
#' will get infinite likelihoods 
pdf("SIplot5.pdf",
    width = 6, height = 6)
plot(t.data$prob, t.data$tr.prob, xlab = "Original probability", ylab = "Transformed probability")
abline(a =0, b = 1)
dev.off()

range(t.data$tr.prob)
hist(t.data$tr.prob)
t.data$weights=nrow(t.data)*t.data$fix_tot/sum(t.data$fix_tot) 
plot(t.data$fix_tot, t.data$weights)
#' weight term allows to give more weight to observations that are more 
#' accurate (i.e., based on a larger sample)

#### (11) full with all correlations ----
library(glmmTMB)
library(beepr)
start=Sys.time()
full.wac = glmmTMB(tr.prob ~ 
                     type_of_day*z.ticket_sales+
                     Origin*z.ticket_sales+
                     z.stringency_index*type_of_day+
                     z.stringency_index*Origin+
                     Sex+ 
                     age_class+ 
                     sin(doy.rad) + cos(doy.rad)+ 
                     sin(2*doy.rad) + cos(2*doy.rad)+
                     (1+ type_of_day*z.ticket_sales+
                        z.stringency_index*type_of_day|id_name)+
                     (1+ Origin+
                        Sex+
                        age_class|date),
                   data = t.data,
                   weights = weights,
                   family=beta_family)
end=Sys.time()
end-start
beep(2)
#' did not converge
#' 26 minutes

#### (12) full with no correlations ----
start=Sys.time()
full = glmmTMB(tr.prob ~ 
                 type_of_day*z.ticket_sales+
                 Origin*z.ticket_sales+
                 z.stringency_index*type_of_day+
                 z.stringency_index*Origin+
                 Sex+ 
                 age_class+ 
                 sin(doy.rad) + cos(doy.rad)+ 
                 sin(2*doy.rad) + cos(2*doy.rad)+
                 (1+ type_of_day*z.ticket_sales+
                    z.stringency_index*type_of_day||id_name)+
                 (1+ Origin+
                    Sex+
                    age_class||date),
               data = t.data,
               weights = weights,
               family=beta_family)
end=Sys.time()
end-start
beep(2)
# 2.4 mins

summary(full)$varcor
summary(full)
confint(full)

round(summary(full)$coefficients$cond, 3)

round(confint(full)[1:16,], 3)

#### (13) model diagnostics ----
overdisp.test(full) #the model is not overdispersed
#'      chisq    df         P dispersion.parameter
#' 1 35815.04 35775 0.4395414             1.001119

ranef.diagn.plot(full)
#' BLUPs are normally distributed

plot(fitted(full), t.data$tr.prob, pch=19, cex=0.5, col=grey(level=0.2, alpha=0.2))

#note the exclusion of the interaction for vif: 
full.vif = glmmTMB(tr.prob ~ 
                     type_of_day+
                     z.ticket_sales+
                     z.stringency_index+
                     Origin+
                     Sex+ 
                     age_class+ 
                     sin(doy.rad) + cos(doy.rad)+ 
                     sin(2*doy.rad) + cos(2*doy.rad)+
                     (1+ type_of_day*z.ticket_sales+
                        z.stringency_index*type_of_day||id_name)+
                     (1+ Origin+
                        Sex+
                        age_class||date),
             data = t.data,
             weights = weights,
             family=beta_family)
performance::check_collinearity(full.vif)
beep(2)

#### (14) null model----
start=Sys.time()
null = glmmTMB(tr.prob ~ 
                 Origin+
                 Sex+ 
                 age_class+ 
                 sin(doy.rad) + cos(doy.rad)+ 
                 sin(2*doy.rad) + cos(2*doy.rad)+
                 (1+ type_of_day*z.ticket_sales+
                    z.stringency_index*type_of_day||id_name)+
                 (1+ Origin+
                    Sex+
                    age_class||date),
               data = t.data,
               weights = weights,
               family=beta_family)
end=Sys.time()
end-start
beep(2)
#1.5 mins 
round(as.data.frame(anova(full, null, test="Chisq")), 3)
#' full-null model comparison is significant
#'      Df       AIC       BIC   logLik  deviance    Chisq Chi Df Pr(>Chisq)
#' null 21 -230425.8 -230247.6 115233.9 -230467.8       NA     NA         NA
#' full 28 -230428.6 -230191.0 115242.3 -230484.6 16.76514      7 0.01897504

#### (15) summary ----
round(summary(full)$coefficients$cond, 3)
summary(full)$varcor
summary(full)
confint(full)

#' AIC
summary(full)$AICtab["AIC"]
#' LogLik
logLik(full)
#' no. levels per grouping factor 
summary(full)$ngrps
#' number random effects
length(summary(full)$varcor)
#' total no. of estimated effects
length(fixef(full)) +
  length(summary(full)$varcor) + 1

#' observations per estimated term
length(residuals(full))/
  (length(fixef(full)) + length(summary(full)$varcor) + 1)
#'  5118.286 observations per estimated term

plot(effects::allEffects(full))

#### (16) fixed effects inference ----
#' effect of stringency index 
#' drop 1 from the model 
tests=as.data.frame(drop1(full, test="Chisq"))
round(tests, 3)
beep(2)
#23 mins

#' need the full effect of season
start=Sys.time()
full.doy = glmmTMB(tr.prob ~ 
                     type_of_day*z.ticket_sales+
                     Origin*z.ticket_sales+
                     z.stringency_index*type_of_day+
                     z.stringency_index*Origin+
                     Sex+ 
                     age_class+
                     (1+ type_of_day*z.ticket_sales+
                        z.stringency_index*type_of_day||id_name)+
                     (1+ Origin+
                        Sex+
                        age_class||date),
                   data = t.data,
                   weights = weights,
                   family=beta_family)
end=Sys.time()
end-start
beep(2)
round(as.data.frame(anova(full, full.doy, test="Chisq")), 3)

#' we see the type of day and ticket sales interaction is not significant 
#' neither is stringency index nor type of day. 
#' neither is ticket sales and origin
#' neither origin nor stringency index
#' so we can basically remove all the interactions from the model and interpret
#' this would be identical to vif model

#### (17) final model ----
full.red = full.vif

plot(effects::allEffects(full.red))
ranef.diagn.plot(full.red) #BLUPS are normally distributed
overdisp.test(full.red)
#'      chisq    df         P dispersion.parameter
#'      1 35818.83 35779 0.4398574             1.001113

tests.red=as.data.frame(drop1(full.red, test="Chisq"))
round(tests.red, 3)
beep(2)

full.red.doy = glmmTMB(tr.prob ~ 
                     type_of_day+
                     z.ticket_sales+
                     z.stringency_index+
                     Origin+
                     Sex+ 
                     age_class+ 
                     (1+ type_of_day*z.ticket_sales+
                        z.stringency_index*type_of_day||id_name)+
                     (1+ Origin+
                        Sex+
                        age_class||date),
                   data = t.data,
                   weights = weights,
                   family=beta_family)
beep(2)

round(as.data.frame(anova(full.red, full.red.doy, test="Chisq")), 3)

round(summary(full.red)$coefficients$cond, 3)
round(confint(full.red), 3)
round(summary(full)$varcor,3)

plot(effects::allEffects(full.red))

#save.image("REVISION FINAL prob_5yr_new.RData")

####---- (18) plotting ----
library(ggplot2)
plot1 = predict_response(full.red, "type_of_day")

p1 = ggplot(plot1, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, col = group)) +
  geom_point(size = 3)+ geom_errorbar(width = 0.1, linewidth = 1)+
  scale_colour_manual(values = "black")+
  theme_minimal() +
  labs(title = "a)",
    x = "Type of day",
    y = "Predicted foraging probability") +
  theme_minimal(base_size = 18)+
  scale_y_continuous(breaks=seq(0, 0.5, 0.1), 
                     limits = c(0, 0.3))+theme(
                       legend.title = element_blank(),
                       legend.position = "bottom",
                       legend.text = element_blank(),
                       legend. = element_blank())+
  guides(colour = guide_legend(override.aes = list(colour = "white", fill = "white")))


plot2 = predict_response(full.red, "z.ticket_sales [all]")
plot2$ticket_sales <- plot2$x * 360.4542+433.9258

p2 <- ggplot(plot2, aes(x = ticket_sales, y = predicted, col = group))+
         geom_line(size=1) +
    scale_colour_manual(values = "black")+
         geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                     alpha = 0.2, fill = "black") +
         labs(
           title = "b)",
           x = "Ticket sales",
           y = "Predicted foraging probability") +
         theme_minimal(base_size = 18)+
  scale_y_continuous(breaks=seq(0, 0.5, 0.1), 
                     limits = c(0, 0.3))+theme(
                       legend.title = element_blank(),
                       legend.position = "bottom",
                       legend.text = element_blank(),
                       legend. = element_blank())+
  guides(colour = guide_legend(override.aes = list(colour = "white", fill = "white")))


plot3 = predict_response(full.red, "z.stringency_index [all]")
plot3$stringency_index <- plot3$x * 29.29588+28.31897

p3 <- ggplot(plot3, aes(x = stringency_index, y = predicted, col = group))+
  geom_line(size=1) +
  scale_colour_manual(values = "black")+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "black") +
  labs(
    title = "c)",
    x = "Stringency index",
    y = "Predicted foraging probability") +
  theme_minimal(base_size = 18)+
  scale_y_continuous(breaks=seq(0, 0.5, 0.1), 
                     limits = c(0, 0.3))+theme(
                       legend.title = element_blank(),
                       legend.position = "bottom",
                       legend.text = element_blank(),
                       legend. = element_blank())+
  guides(colour = guide_legend(override.aes = list(colour = "white", fill = "white")))

plot4 = predict_response(full.red, "doy.rad [all]")

ggplot(plot4, aes(x = x, y = predicted))+
  geom_line(size=1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "black") +
  labs(
    title = "b)",
    x = "Stringency index",
    y = "Predicted foraging probability") +
  theme_minimal(base_size = 18)+
  scale_y_continuous(breaks=seq(0, 0.5, 0.1), 
                     limits = c(0, 0.4))


library(ggpubr)
combined_plot <- ggarrange(
  p1, p2, p3,
  ncol = 3,
  nrow = 1,
  labels =NULL 
)

quartz(7,18)
combined_plot
dev.copy2pdf(file="plot1.pdf")

#### (19) exporting ---- 
#save.image("prob_5yr_new.RData")
