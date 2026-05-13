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
library(stringr)

#### (2) Load data ----
xdata = read.csv(file="prob_final_6wk.csv", stringsAsFactors=T)
str(xdata)

#### (3) Visitor data plot ----
df_plot <- xdata %>%
  mutate(
    bi_hour_label = sprintf("%02d–%02d", bi_hourly, bi_hourly + 2),
    bi_hour_label = factor(bi_hour_label, 
                           levels = sort(unique(sprintf("%02d–%02d", bi_hourly, bi_hourly + 2))))
  )

df_long <- df_plot %>%
  select(type_of_day, bi_hour_label, visitors_original, visitors_cumulative) %>%
  tidyr::pivot_longer(
    cols = c(visitors_original, visitors_cumulative),
    names_to = "metric",
    values_to = "visitors"
  ) %>%
  mutate(metric = recode(metric,
                         visitors_original = "Short visit",
                         visitors_cumulative = "Long visit"),
         metric = factor(metric, levels = c("Short visit", "Long visit")))

cb_cols <- c("weekday" = "#7F7F7F",   # grey
             "weekend" = "#E69F00")   # orange

xx = df_long %>% filter(metric == "Short visit")
quartz(height = 5, width = 4.5)
ggplot(xx,
       aes(x = bi_hour_label, y = visitors,
           fill = type_of_day)) +
  geom_boxplot() +
  scale_fill_manual(values = cb_cols, name = "Type of day") +
  labs(
    x = "Bi-hourly interval",
    y = "Number of visitors",
    title = "c)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
dev.copy2pdf(file="SIplotx12.pdf")

xx = df_long %>% filter(metric == "Long visit")

quartz(height = 5, width = 5.5)
ggplot(xx,
       aes(x = bi_hour_label, y = visitors,
           fill = type_of_day)) +
  geom_boxplot() +
  scale_fill_manual(values = cb_cols, name = "Type of day") +
  labs(
    x = "Bi-hourly interval",
    y = "Number of visitors",
    title = "d)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.copy2pdf(file="SIplotx13.pdf")


#### (4) check NA and str ----
#' are there any NAs
any(is.na(xdata))
str(xdata)

#' get bihour from 1-6
xdata$bi_hourly <- (as.integer(as.factor(xdata$bi_hourly)))
class(xdata$bi_hourly)

#### (5) dataframe touch-ups ----
#' fixes range from 1-4
range(xdata$fix_tot_bihrly)
xdata <- xdata %>% filter(fix_tot_bihrly >= 4) #lose about 1200 points
levels(droplevels(as.factor(xdata$id_name)))

#### (6) look at predictors ----
hist(xdata$bi_hourly)
hist((xdata$visitors_original))
hist((xdata$visitors_cumulative))

#' ID per weekend - each indv has weekday and weekend obs
xx = table(xdata$id_name, xdata$type_of_day)
head(xx)
table(apply(X = xx>0, MARGIN = 2, FUN = sum))

#' range of visitors wrt to type of day
plot(xdata$type_of_day, xdata$visitors_original) 
plot(xdata$type_of_day, xdata$visitors_cumulative) 

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
xx = aggregate(x=1:nrow(xdata), by = xdata[, c("id.bi_hour", "date")], FUN= length)
head(xx)
plot(table(xx$x))
table(xx$x) # for 4475 combinations of indv and date-bihour, 1 prob measure

#' obs per ID
table(xdata$id_name) # how many datapts per indv
plot(table(table(xdata$id_name)))

####(9) checking random slopes ----
source("diagnostic_fcns.r")
xx.re.tab.1 = fe.re.tab(fe.model = "prob~ type_of_day*visitors_original+
                      Origin*visitors_original+
                      bi_hourly+
                      Sex + 
                      age_class",
                      re = "(1|id_name)+(1|date)",
                      data = xdata, other.vars=c("inside", "outside", "fix_tot_bihrly"))

xx=xx.re.tab.1$detailed[4]
xx
xx.re.tab.1$summary[6] 


xx.re.tab.2 = fe.re.tab(fe.model = "prob~ type_of_day*visitors_cumulative+
                      Origin*visitors_cumulative+
                      bi_hourly+
                      Sex + 
                      age_class",
                      re = "(1|id_name)+(1|date)",
                      data = xdata, other.vars=c("inside", "outside", "fix_tot_bihrly"))

#### (10) setting up model ----
str(xx.re.tab.1$data)
str(xx.re.tab.2$data)

t.data.original =xx.re.tab.1$data
t.data.cumulative =xx.re.tab.2$data

#' z-transform the covariates
range(t.data.original$visitors_original)
t.data.original$z.visitors=as.vector(scale(t.data.original$visitors_original))
range(t.data.original$z.visitors)

range(t.data.cumulative$visitors_cumulative)
t.data.cumulative$z.visitors=as.vector(scale(t.data.cumulative$visitors_cumulative))
range(t.data.cumulative$z.visitors)

#' response
t.data.original$tr.prob = (t.data.original$prob*(length(t.data.original$prob)-1)+0.5)/length(t.data.original$prob)
#' squeezes everything closer to 0.5 
#' use when you have zero or 1 in the response, the beta distr cannot cope with it,
#' will get infinite likelihoods 
plot(t.data.original$prob, t.data.original$tr.prob)
abline(a =0, b = 1)
t.data.original$weights=nrow(t.data.original)*t.data.original$fix_tot_bihrly/sum(t.data.original$fix_tot_bihrly)
#' weight term allows to give more weight to observations that are more 
#' accurate (i.e., based on a larger sample)

t.data.cumulative$tr.prob = (t.data.cumulative$prob*(length(t.data.cumulative$prob)-1)+0.5)/length(t.data.cumulative$prob)
t.data.cumulative$weights=nrow(t.data.cumulative)*t.data.cumulative$fix_tot_bihrly/sum(t.data.cumulative$fix_tot_bihrly)

class(t.data.original$bi_hourly)
class(t.data.cumulative$bi_hourly)

#### (11) full with all correlations ----
library(gamm4)
library(beepr)

start=Sys.time()
full.wac.original = gamm4(cbind(inside, outside) ~ 
                            type_of_day*z.visitors+
                            Origin*z.visitors+
                            s(bi_hourly, k = 6)+
                            Sex+ 
                            age_class, 
                          data = t.data.original,
                          family=binomial(),
                          random = ~ 
                            (1|id_name)+(0+ type_of_day*z.visitors|id_name)+
                            (1|date)+(0+ z.visitors+Origin + Sex +age_class|date),
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
end=Sys.time()
end-start
beep(2)

start=Sys.time()
full.wac.cumulative = gamm4(cbind(inside, outside) ~ 
                              type_of_day*z.visitors+
                              Origin*z.visitors+
                              s(as.numeric(bi_hourly), k = 6)+
                              Sex+ 
                              age_class, 
                            data = t.data.cumulative,
                            family=binomial(),
                            random = ~ (1|id_name)+(0+ type_of_day*z.visitors|id_name)+
                              (1|date)+(0+ z.visitors+ Origin + Sex + age_class|date),
                            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
end=Sys.time()
end-start
beep(2)

#### (12) model diagnostics ----
overdisp.test(full.wac.original$mer)
overdisp.test(full.wac.cumulative$mer)
ranef.diagn.plot(full.wac.original$mer)
ranef.diagn.plot(full.wac.cumulative$mer)
#' BLUPs are normally distributed

#note the exclusion of the interaction for vif: 
full.original.vif = gamm4(cbind(inside, outside) ~ 
                            type_of_day+z.visitors+
                            Origin+
                            bi_hourly+
                            Sex+ 
                            age_class, 
                          data = t.data.original,
                          family=binomial(),
                          random = ~ 
                            (1|id_name)+(0+ type_of_day*z.visitors|id_name)+
                            (1|date)+(0+ z.visitors+Origin + Sex +age_class|date),
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
performance::check_collinearity(full.original.vif)

full.cumulative.vif = gamm4(cbind(inside, outside) ~ 
                              type_of_day+z.visitors+
                              Origin+
                              bi_hourly+
                              Sex+ 
                              age_class, 
                            data = t.data.cumulative,
                            family=binomial(),
                            random = ~ 
                              (1|id_name)+(0+ type_of_day*z.visitors|id_name)+
                              (1|date)+(0+ z.visitors+Origin + Sex +age_class|date),
                            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
performance::check_collinearity(full.cumulative.vif)

#### (13) null model----
null.original = gamm4(cbind(inside, outside) ~ 
                        Origin+
                        s(as.numeric(bi_hourly), k = 6)+
                        Sex+ 
                        age_class, 
                      data = t.data.original,
                      family=binomial(),
                      random = ~ (1|id_name)+(0+ type_of_day*z.visitors|id_name)+
                        (1|date)+(0+ z.visitors+
                                    Origin + 
                                    Sex +
                                    age_class|date),
                      control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

round(as.data.frame(anova(full.wac.original$mer, null.original$mer, test="Chisq")), 3)
#' full-null model comparison is insignificant

null.cumulative = gamm4(cbind(inside, outside) ~ 
                          Origin+
                          s(as.numeric(bi_hourly), k = 6)+
                          Sex+ 
                          age_class, 
                        data = t.data.cumulative,
                        family=binomial(),
                        random = ~ (1|id_name)+(0+ type_of_day*z.visitors|id_name)+
                          (1|date)+(0+ z.visitors+
                                      Origin + 
                                      Sex +
                                      age_class|date),
                        control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
round(as.data.frame(anova(full.wac.cumulative$mer, null.cumulative$mer, test="Chisq")), 3)
#' full-null model comparison is insignificant

#' observations per estimated term
length(residuals(full.wac.cumulative$mer))/
  (length(fixef(full.wac.cumulative$mer)) + length(summary(full.wac.cumulative$mer)$varcor) + 1)
#'  275 observations per estimated term

length(residuals(full.wac.original$mer))/
  (length(fixef(full.wac.original$mer)) + length(summary(full.wac.original$mer)$varcor) + 1)

#save.image('prob_final_analysis_6wk.RData')
