##############################################
# 02_survival_analysis_MRS_Cox_KM.R
##############################################

library(survival)
library(survminer)
library(tidyverse)

dat <- readRDS("outputs/dat_TCGA_MRS_clinical.rds")

# Make sure variables are correctly formatted
dat_cox <- dat %>%
  mutate(
    stage_simple = factor(stage_simple, levels = c("Early", "Late")),
    HPV          = factor(HPV,          levels = c("negative", "positive")),
    alcohol      = factor(Alcohol.History, levels = c("No", "Yes")),
    smoker       = factor(Smoking.History, levels = c("No", "Yes")),
    sub_site_grouped = factor(
      sub_site_grouped,
      levels = c("Oral cavity", "Oropharynx", "Larynx", "Hypopharynx")
    )
  )

# Cox model
fit_cox <- coxph(
  Surv(OS.time, OS.event) ~ risk_score + Age + stage_simple +
    HPV + alcohol + smoker + sub_site_grouped,
  data = dat_cox
)

summary(fit_cox)

# Export coefficients for supplement
cox_coef <- broom::tidy(fit_cox, exponentiate = TRUE, conf.int = TRUE)
write.csv(cox_coef, "outputs/cox_model_coefficients_MRS.csv", row.names = FALSE)

# Median split for KM
dat_cox <- dat_cox %>%
  mutate(risk_group = ifelse(risk_score > median(risk_score, na.rm = TRUE),
                             "High", "Low"))

fit_km <- survfit(Surv(OS.time, OS.event) ~ risk_group, data = dat_cox)

# KM plot
p_km <- ggsurvplot(
  fit_km,
  data      = dat_cox,
  risk.table = TRUE,
  pval      = TRUE,
  legend.labs = c("Low risk", "High risk"),
  xlab      = "Time (months)",
  ylab      = "Overall survival probability"
)

ggsave("outputs/Fig_KM_MRS_OS.pdf", p_km$plot, width = 5, height = 4)
