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

#### (3) Age-sample size plot ----
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

#### (4) check NA and str ----
#' are there any NAs
any(is.na(xdata))
str(xdata)

# drop NAs in area column
xdata <- xdata %>% tidyr::drop_na(area_m2)
any(is.na(xdata))

#### (5) dataframe touch-ups ----
#' check that minimum fixes in a day are at least 15 
range(xdata$fix_tot)

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
xx.re.tab = fe.re.tab(fe.model = "area_m2~ type_of_day*ticket_sales+
                      Origin*ticket_sales+
                      stringency_index*type_of_day+
                      stringency_index*Origin+
                      Sex + 
                      age_class + 
                      date.rad+
                      doy.rad",
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
#' 365.9858, 407.033
range(t.data$stringency_index)
t.data$z.stringency_index=as.vector(scale(t.data$stringency_index))
#' 30.40507, 28.35908

#' response
hist(t.data$area_m2)
range(t.data$area_m2)
t.data$tr.area_m2 = log(t.data$area_m2)

hist(t.data$tr.area_m2)
range(t.data$tr.area_m2)
t.data$weights=t.data$fix_tot
#' weight term allows to give more weight to observations that are more 
#' accurate (i.e., based on a larger sample)

#### (11) full with all correlations ----
library(glmmTMB)
library(beepr)
start=Sys.time()
full.wac = glmmTMB(tr.area_m2 ~ 
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
                   REML=F,
                   weights = weights,
                   family=gaussian)
end=Sys.time()
end-start
beep(2)
#' converged
#' 1 minute

#### (13) model diagnostics ----
diagnostics.plot(full.wac)
#' some structure in the data
ranef.diagn.plot(full.wac)
#' BLUPs are normally distributed

#note the exclusion of the interaction for vif: 
full.vif = glmmTMB(tr.area_m2 ~ 
                     type_of_day+
                     z.ticket_sales+
                     z.stringency_index+
                     Origin+
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
                   REML=F,
                   weights = weights,
                   family=gaussian)
performance::check_collinearity(full.vif)

#### (14) null model----
null = glmmTMB(tr.area_m2 ~ 
                 Origin+
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
               REML=F,
               weights = weights,
               family=gaussian)

round(as.data.frame(anova(full.wac, null, test="Chisq")), 3)
#' full-null model comparison is significant
#'         Df     AIC     BIC   logLik deviance  Chisq Chi Df Pr(>Chisq)
#'null     46 2834504 2834887 -1417206  2834412     NA     NA         NA
#'full.wac 53 2834481 2834923 -1417188  2834375 36.503      7          0

#### (15) summary ----
round(summary(full.wac)$coefficients$cond, 3)
summary(full.wac)$varcor
round(confint(full.wac),3)

#' AIC
summary(full.wac)$AICtab["AIC"]
#' LogLik
logLik(full.wac)
#' no. levels per grouping factor 
summary(full.wac)$ngrps
#' number random effects
length(summary(full)$varcor)
#' total no. of estimated effects
length(fixef(full.wac)) +
  length(summary(full.wac)$varcor) + 1

#' observations per estimated term
length(residuals(full.wac))/
  (length(fixef(full.wac)) + length(summary(full.wac)$varcor) + 1)
#'  4397.286 observations per estimated term

plot(effects::allEffects(full.wac))

#### (16) fixed effects inference ----
#' drop 1 from the model 
tests=as.data.frame(drop1(full.wac, test="Chisq"))
round(tests, 3)
beep(2)

#' need the full effect of season
start=Sys.time()
full.doy = glmmTMB(tr.area_m2 ~ 
                     type_of_day*z.ticket_sales+
                     Origin*z.ticket_sales+
                     z.stringency_index*type_of_day+
                     z.stringency_index*Origin+
                     Sex+ 
                     age_class+ 
                     (1+ type_of_day*z.ticket_sales+
                        z.stringency_index*type_of_day|id_name)+
                     (1+ Origin+
                        Sex+
                        age_class|date),
                   data = t.data,
                   REML=F,
                   weights = weights,
                   family=gaussian)
end=Sys.time()
end-start
beep(2)
round(as.data.frame(anova(full.wac, full.doy, test="Chisq")), 3)

#' we see the type of day and ticket sales interaction is not significant 
#' neither is ticket sales and origin
#' so we can basically remove those two interactions from the model and interpret

#### (17) final model ----
full.red = glmmTMB(tr.area_m2 ~ 
                    type_of_day+
                    z.ticket_sales+
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
                  REML=F,
                  weights = weights,
                  family=gaussian)
beep(2)
  
#' model diagnostics  
plot(effects::allEffects(full.red))
ranef.diagn.plot(full.red) #BLUPS are normally distributed

tests.red=as.data.frame(drop1(full.red, test="Chisq"))
round(tests.red, 3)
beep(2)

full.red.doy = glmmTMB(tr.area_m2 ~ 
                         type_of_day+
                         z.ticket_sales+
                         z.stringency_index*type_of_day+
                         z.stringency_index*Origin+
                         Sex+ 
                         age_class+ 
                         (1+ type_of_day*z.ticket_sales+
                            z.stringency_index*type_of_day|id_name)+
                         (1+ Origin+
                            Sex+
                            age_class|date),
                       data = t.data,
                       REML=F,
                       weights = weights,
                       family=gaussian)
beep(2)

round(as.data.frame(anova(full.red, full.red.doy, test="Chisq")), 3)

round(summary(full.red)$coefficients$cond, 3)
round(confint(full.red), 3)
summary(full.red)$varcor

#save.image("REVISION FINAL kde_5yr_new.RData")

####---- (18) plotting ----
library(ggeffects)

plot1 = predict_response(full.red, c("z.stringency_index [all]", "type_of_day"))

scale(t.data$stringency_index)
plot1$stringency_index = plot1$x * 28.35908+ 30.40507


p1 <- ggplot(plot1, aes(x = stringency_index, y = predicted, col = group))+
  geom_line(size=1, aes(linetype = group)) +
  scale_color_manual(values = c("#555555", "black"))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "black") +
  labs(
    title = "a)",
    x = "Stringency index",
    y = "Log (space-use)") +
theme_minimal(base_size = 18) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  )+
  scale_y_continuous(breaks=seq(6, 16, 2), 
                     limits = c(6, 16))


plot2 = predict_response(full.red, c("z.stringency_index [all]", "Origin"))
plot2$stringency_index = plot2$x * 28.35908+ 30.40507

p2 <- ggplot(plot2, aes(x = stringency_index, y = predicted, col = group))+
  geom_line(size=1, aes(linetype = group)) +
  scale_color_manual(values = c("#555555", "black"))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "black") +
  labs(
    title = "c)",
    x = "Stringency index",
    y = "Log (space-use)") +
  theme_minimal(base_size = 18)+
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  )+
  scale_y_continuous(breaks=seq(6, 16, 2), 
                     limits = c(6, 16))


plot3 = predict_response(full.red, c("z.ticket_sales [all]"))
plot3$z.ticket_sales = plot3$x * 365.9858+407.033

p3 <- ggplot(plot3, aes(x = z.ticket_sales, y = predicted, col = group))+
  geom_line(size=1) + 
  scale_colour_manual(values = "black")+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "black") +
  labs(
    title = "b)",
    x = "Ticket sales",
    y = "Log (space-use)") +
  theme_minimal(base_size = 18)+
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_blank(),
    legend. = element_blank()
  )+guides(colour = guide_legend(override.aes = list(colour = "white", fill = "white")))+
  scale_y_continuous(breaks=seq(6, 16, 2), 
                     limits = c(6, 16))


plot4 = predict_response(full.red, "doy.rad [all]")
ggplot(plot4, aes(x = x, y = predicted))+
  geom_line(size=1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "black")# +
#  labs(
#    title = "b)",
#    x = "Stringency index",
#    y = "Predicted foraging probability") +
  #theme_minimal(base_size = 14)#+
 # scale_y_continuous(breaks=seq(0, 0.5, 0.1), 
            #         limits = c(0, 0.4))


library(ggpubr)
combined_plot <- ggarrange(
  p1, p3, p2,
  ncol = 3,
  nrow = 1,
  labels =NULL 
)

quartz(7,18)
combined_plot
dev.copy2pdf(file="plot2.pdf")
  
#### (19) exporting ---- 
#save.image("mcp_5yr_new.RData")


