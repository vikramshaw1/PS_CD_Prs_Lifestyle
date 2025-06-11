# 1. Load packages ----
library(data.table)
library(dplyr)
library(survival)
library(survminer)
library(cmprsk)
library(MatchIt)
library(tidyverse)

# 2. Load data ----

# 2.1 Load downloaded data from UKB DNA nexus (this information must be downloaded and analyzed by an individual with an active UKB application who has obtained DNA Nexus access through the UKB Access Management (AMS) system)

# 2.2 Processing data
data <- data[data$p131742 != "Code has event date matching participant's date of birth", ] #502,121
data <- data[data$p131742 != "Code has event date after participant's date of birth and falls in the same calendar year as date of birth", ] #502,116
data$p131742 <- as.Date(data$p131742)
data$p53_i0 <- as.Date(data$p53_i0)
data$p40000_i0 <- as.Date(data$p40000_i0)
data <- data[data$p131626 != "Code has event date matching participant's date of birth", ] #502,115
data$p131626 <- as.Date(data$p131626)

count_rows <- sum(is.na(data$p131626) | is.na(data$p131742)) #219 with both

#Exclude if both occurred before assessment 
data <- data[!(data$p131742 < data$p53_i0 & data$p131626 < data$p53_i0) | is.na(data$p131742) | is.na(data$p131626), ] #502,031, removed 84 where both occurred before

data$psoriasis <- ifelse(!is.na(data$p131742), 1, 0)
data$crohn <- ifelse(!is.na(data$p131626), 1, 0)
data$psor_first <- ifelse(data$p131742 < data$p131626, 1, 0)
data$cd_first <- ifelse(data$p131742 > data$p131626, 1, 0)

#Keep only European background
data <- data[data$p21000_i0 == "Any other white background" | data$p21000_i0 == "British" | data$p21000_i0 == "White" | data$p21000_i0 == "Irish", ] #472,264

count_rows <- sum(data$psoriasis == 1 & data$crohn == 1, na.rm = TRUE) #132
data_both <- data[data$psoriasis == 1 & data$crohn == 1, ]

#Times (study end date 2024-01-01)
end_date <- as.Date("2024-01-01", format = "%Y-%m-%d")

#Create the time_to_both column with the date comorbid disease diagnosed
data$p131742 <- as.character(data$p131742)
data$p131626 <- as.character(data$p131626)

data$time_to_both <- ifelse(
  data$psoriasis == 1 & data$crohn == 1,
  ifelse(
    as.Date(data$p131742, format = "%Y-%m-%d") >= as.Date(data$p131626, format = "%Y-%m-%d"),
    data$p131742,
    data$p131626
  ),
  NA
)

#Time to development of comorbid disease
data$tsurv <- with(data, ifelse(
  psoriasis == 1 & crohn == 1, 
  as.numeric(difftime(time_to_both, p53_i0, units = "days")) / 30.44, #Months between `assessment date` and development of both
  ifelse(
    !is.na(p40000_i0), 
    as.numeric(difftime(p40000_i0, p53_i0, units = "days")) / 30.44, #Months between `assessment date` and `date of death`
    as.numeric(difftime(end_date, p53_i0, units = "days")) / 30.44 #Months between `assessment date` and `end_date`
  )
))

#Calculating the PRS quartiles
#Remove rows with NA
data <- data[!is.na(data$p26269), ] #457,828
data$p26269 <- as.numeric(data$p26269)
data <- data %>%
  mutate(psor_prs_half = ntile(p26269, 2))

data <- data[!is.na(data$p26229), ] #457,828
data$p26229 <- as.numeric(data$p26229)
data <- data %>%
  mutate(cd_prs_half = ntile(p26229, 2))

#Physical activity
data <- data[data$p884_i0 != "", ] #457,828
data <- data[data$p884_i0 != "Do not know", ] #437,192
data <- data[data$p884_i0 != "Prefer not to answer", ] #435,751

#Smoking
data <- data[data$p20160_i0 != "", ] #434,438

#Drinking
data <- data[data$p20117_i0 != "", ] #434,438
data <- data[data$p20117_i0 != "Prefer not to answer", ] #434,170

#Sleep
data <- data[data$p1160_i0 != "Prefer not to answer", ] #434,083
data <- data[data$p1160_i0 != "Do not know", ] #432,565

#TDI, age (removing NAs)
data <- subset(data, !is.na(p22189)) #432,053 #TDI
data <- subset(data, !is.na(p21022)) #432,053 #age
data <- subset(data, !is.na(p21001_i0)) #430,758 #bmi

# 2.3 Creating modifiable lifestyle variable
data <- data %>%
  mutate(
    lifestyle_raw = 
      as.integer(p1160_i0 >= 7) + 
      as.integer(p884_i0 %in% c(5, 6, 7)) + 
      as.integer(p20117_i0 == "Never") + 
      as.integer(p20160_i0 == "No") + 
      as.integer(p21001_i0 < 30),
    lifestyle_factor = case_when(
      lifestyle_raw %in% c(0, 1) ~ "Unhealthy",
      lifestyle_raw %in% c(2, 3) ~ "Intermediate",
      lifestyle_raw %in% c(4, 5) ~ "Healthy"
    ),
    lifestyle_factor = factor(
      lifestyle_factor, 
      levels = c("Unhealthy", "Intermediate", "Healthy")  
    )
  )

data$comorbid <- ifelse(data$psoriasis == 1 & data$crohn == 1, 1, 0)

# 3. Survival and density plots (Figure 2) ----
#Coding variable types
data$p31 <- as.factor(data$p31)
data$p26201_a0 <- as.numeric(data$p26201_a0)
data$p26201_a1 <- as.numeric(data$p26201_a1)
data$p26201_a2 <- as.numeric(data$p26201_a2)
data$p26201_a3 <- as.numeric(data$p26201_a3)
data$p21022 <- as.numeric(data$p21022)
data$p22189 <- as.numeric(data$p22189)
data$p1160_i0 <- as.numeric(data$p1160_i0)
data$p21001_i0 <- as.numeric(data$p21001_i0)
data$p884_i0 <- as.numeric(data$p884_i0)
data$p20117_i0 <- as.factor(data$p20117_i0)
data$p20160_i0 <- as.factor(data$p20160_i0)
data$psor_prs_half <- as.factor(data$psor_prs_half)
data$cd_prs_half <- as.factor(data$cd_prs_half)
data$p26229 <- as.numeric(data$p26229)
data$lifestyle_factor <- as.factor(data$lifestyle_factor)

#Combining CD and PsO PRS into single variable with 4 levels
data$genetics <- paste(data$psor_prs_half, data$cd_prs_half, sep="_")
data$genetics <- as.factor(data$genetics)

#Survival analysis (univariate) for PRS groups
surv_obj <- Surv(time = data$tsurv, event = data$comorbid)
km_fit <- survfit(surv_obj ~ genetics, data = data)

ggsurvplot(
  km_fit,
  data = data,
  pval = TRUE,
  pval.method = TRUE, 
  conf.int = TRUE,                 
  risk.table = TRUE,               
  legend.title = "",  
  legend.labs = c("PsO-Lo & CD-Lo", "PsO-Lo & CD-Hi", "PsO-Hi & CD-Lo", "PsO-Hi & CD-Hi"), 
  xlab = "Time (Months)",          
  ylab = "Disease-Free Probability", 
  title = "",
  palette = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"),  
  ylim = c(0.999, 1),
  pval.coord = c(175, 0.999),
  pval.method.coord = c(175, 0.99905),
  ggtheme = theme_survminer(base_size = 14) + 
    theme(
      plot.title = element_text(size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 16)
    )
)

#Survival analysis (univariate) for lifestyle factor groups
surv_obj <- Surv(time = data$tsurv, event = data$comorbid)
km_fit <- survfit(surv_obj ~ lifestyle_factor, data = data)

ggsurvplot(
  km_fit,
  data = data,
  pval = TRUE,
  pval.method = TRUE, 
  conf.int = TRUE,                 
  risk.table = TRUE,               
  legend.title = "",  
  legend.labs = c("Unhealthy", "Intermediate", "Healthy"), 
  xlab = "Time (Months)",          
  ylab = "Disease-Free Probability", 
  title = "Disease-Free Survival by Lifestyle Factor",
  palette = c("#66c2a5", "#fc8d62", "#8da0cb"),  
  ylim = c(0.999, 1),
  pval.coord = c(175, 0.999),
  pval.method.coord = c(175, 0.99905),
  ggtheme = theme_survminer(base_size = 14) + 
    theme(
      plot.title = element_text(size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 16)
    )
)

#Density plots for PRS analysis
data <- data %>%
  rename(`PsO PRS` = p26269, `CD PRS` = p26229)

data_long <- data %>%
  pivot_longer(cols = c(`PsO PRS`, `CD PRS`), names_to = "PRS_Type", values_to = "Score")
data_long$Status <- ifelse(data_long$comorbid == 1, "Case", "Control")
table(interaction(data_long$PRS_Type, data_long$Status)) 
legend_labels <- c("Comorbid - CD PRS", "Comorbid - PsO PRS", "Control - CD PRS", "Control - PsO PRS")

ggplot(data_long, aes(x = Score, fill = interaction(PRS_Type, Status))) +
  geom_density(alpha = 0.4, color = "black") +
  facet_wrap(~PRS_Type, scales = "free") +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = legend_labels) +  
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"), 
    panel.border = element_blank()  
  ) +
  labs(x = "PRS Score", y = "Density", title = "")



# 4. Multivariable Cox proportional hazards regression analysis (Table 2) ----

data <- data %>%
  mutate(
    lifestyle_factor = factor(
      lifestyle_factor, 
      levels = c("Healthy", "Intermediate", "Unhealthy")  
    )
  )

cox_model <- coxph(Surv(tsurv, comorbid) ~ genetics + p31 + p26201_a0 + p26201_a1 + p26201_a2 + p26201_a3 + p21022 + p22189 + lifestyle_factor, data = data)
summary(cox_model)

#Ptrend
data$genetics <- factor(data$genetics,
                                  levels = c("1_1",
                                             "1_2",
                                             "2_1",
                                             "2_2"),
                                  ordered = TRUE)
cox_trend <- coxph(Surv(tsurv, comorbid) ~ genetics, data = data)
summary(cox_trend)

#Ptrend lifestyle factor
data$lifestyle_factor <- factor(data$lifestyle_factor,
                        levels = c("Unhealthy",
                                   "Intermediate",
                                   "Healthy"),
                        ordered = TRUE)
cox_trend <- coxph(Surv(tsurv, comorbid) ~ lifestyle_factor, data = data)
summary(cox_trend)


# 5. Multivariable Cox proportional hazards regression analysis, stratified by disease of first onset (Table 3) ----
#Note: Clear environment, re-run blocks 2 and 3 before running this code
# 5.1 Psoriasis onset first analysis
data <- data %>%
  mutate(
    lifestyle_factor = factor(
      lifestyle_factor, 
      levels = c("Healthy", "Intermediate", "Unhealthy")  
    )
  )

psor_first <- data %>%
  filter(psor_first == "1" | is.na(psor_first))
cox_model_psor <- coxph(Surv(tsurv, comorbid) ~ genetics + p31 + p26201_a0 + p26201_a1 + p26201_a2 + p26201_a3 + p21022 + p22189 + lifestyle_factor, data = psor_first)
summary(cox_model_psor)

#Ptrend genetics
psor_first$genetics <- factor(psor_first$genetics,
                        levels = c("1_1",
                                   "1_2",
                                   "2_1",
                                   "2_2"),
                        ordered = TRUE)
cox_trend <- coxph(Surv(tsurv, comorbid) ~ genetics, data = psor_first)
summary(cox_trend)

#Ptrend lifestyle factor
psor_first$lifestyle_factor <- factor(psor_first$lifestyle_factor,
                                levels = c("Unhealthy",
                                           "Intermediate",
                                           "Healthy"),
                                ordered = TRUE)
cox_trend <- coxph(Surv(tsurv, comorbid) ~ lifestyle_factor, data = psor_first)
summary(cox_trend)

# 5.2 CD onset first analysis
cd_first <- data %>%
  filter(psor_first == "0" | is.na(psor_first))
cox_model_cd <- coxph(Surv(tsurv, comorbid) ~ genetics + p31 + p26201_a0 + p26201_a1 + p26201_a2 + p26201_a3 + p21022 + p22189 + lifestyle_factor, data = cd_first)
summary(cox_model_cd)

#Ptrend genetics
cd_first$genetics <- factor(cd_first$genetics,
                              levels = c("1_1",
                                         "1_2",
                                         "2_1",
                                         "2_2"),
                              ordered = TRUE)
cox_trend <- coxph(Surv(tsurv, comorbid) ~ genetics, data = cd_first)
summary(cox_trend)

#Ptrend lifestyle factor
cd_first$lifestyle_factor <- factor(cd_first$lifestyle_factor,
                                      levels = c("Unhealthy",
                                                 "Intermediate",
                                                 "Healthy"),
                                      ordered = TRUE)
cox_trend <- coxph(Surv(tsurv, comorbid) ~ lifestyle_factor, data = cd_first)
summary(cox_trend)

# 6. Comparing highest risk group to all others (Supplementary Table 1) ----
#Note: Clear environment, re-run blocks 2 and 3 before running this code
data$genetics_binary <- ifelse(data$genetics != "2_2",0,1)

data$genetics_binary <- factor(data$genetics_binary,
                               levels = c(0,
                                          1),
                               ordered = TRUE)

data$lifestyle_factor <- factor(data$lifestyle_factor,
                                levels = c("Healthy",
                                           "Intermediate",
                                           "Unhealthy"),
                                ordered = TRUE)

cox_model <- coxph(Surv(tsurv, comorbid) ~ genetics_binary + p31 + p26201_a0 + p26201_a1 + p26201_a2 + p26201_a3 + p21022 + p22189 + lifestyle_factor, data = data)
summary(cox_model)
