library(papaja)
library(dplyr)
library(gtools)
library(tidyverse)
library(pander)
library(psych)
library(ggplot2)
library(ggh4x)
library(lme4)
library(clubSandwich)
library(lmerTest)
library(brms)
library(broom.mixed)
library(gtools)
library(haven)
library(gridExtra)
library(CR2)

## readin dataset
mlbook_red <- read.table("data/mlbook2_r.dat", header=TRUE)
mlbook_red <- mlbook_red %>%
  group_by(schoolnr) %>%
  mutate(sch_iq = mean(IQ_verb),
         sch_iqc = IQ_verb - sch_iq,
         sch_ses = mean(ses),
         sch_sesc = ses - sch_ses)

# randomly sample a smaller datset 
set.seed(177)
ran_id <- sample(mlbook_red$schoolnr, 15)
mlbook_sam <- mlbook_red %>%
  filter(schoolnr %in% ran_id)
fit_em <- lmer(langPOST ~ 1 + (1| schoolnr), data = mlbook_sam)
icc <-  as.data.frame(VarCorr(fit_em))$vcov[1]/sum(as.data.frame(VarCorr(fit_em))$vcov) 

## fit different models
fit_ols <- lm(langPOST ~ sch_iqc + sch_sesc + sch_iq + sch_ses, data = mlbook_sam)
fit_ri <- lmer(langPOST ~ sch_iqc + sch_sesc + sch_iq + sch_ses + (1| schoolnr), data = mlbook_sam)
fit_rs <- lmer(langPOST ~ sch_iqc + sch_sesc + sch_iq + sch_ses + (IQ_verb| schoolnr), data = mlbook_sam)
## test homoscedasticity assumption
schiq_bp <- ncvMLM(lmer(langPOST ~ sch_iq + (1| schoolnr), data = mlbook_sam))
stuiq_bp <- ncvMLM(lmer(langPOST ~ sch_iqc + (1| schoolnr), data = mlbook_sam))
stuses_bp <- ncvMLM(lmer(langPOST ~ sch_sesc + (1| schoolnr), data = mlbook_sam))

## test fixed effects
KR_res <- summary(fit_rs, ddf = "Kenward-Roger")
sandt_ols <- coef_test(fit_ols, vcov = "CR2", cluster = mlbook_sam$schoolnr, test = "Satterthwaite")
sandt_ri <- coef_test(fit_ri, vcov = "CR2", cluster = mlbook_sam$schoolnr, test = "Satterthwaite")
sandt_rs <- coef_test(fit_rs, vcov = "CR2", cluster = mlbook_sam$schoolnr, test = "Satterthwaite")

## visualizing the variance pattern of IQ and SES
resid_lv2 <- ranef(fit_rs)$schoolnr
resid_lv2_var <- attr(resid_lv2, "postVar")
resid_lv2_sd <- sqrt(cbind(
  resid_lv2_var[1, 1, ],
  resid_lv2_var[2, 2, ]
))
std_resid_lv2 <- resid_lv2 %>%
  rownames_to_column("schoolnr") %>%
  transmute(schoolnr = as.double(schoolnr),
            u0 = `(Intercept)` / resid_lv2_sd[, 1],
            u1 = IQ_verb / resid_lv2_sd[, 2])
# Merge lv-2 data with lv-2 standardized residuals
mlbook_lv2 <- mlbook_sam %>%
  group_by(schoolnr) %>%
  summarise(sch_iq = sch_iq[1]) %>%
  left_join(
    std_resid_lv2
  )
p1 <- augment(fit_rs) %>%
  mutate(.std_resid = resid(fit_rs, scaled = TRUE)) %>%
  ggplot(aes(x = sch_iqc, y = .std_resid)) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_smooth(se = FALSE) +
  labs(x = "IQ scores", y = "Level-1 Residuals") +
  scale_color_brewer(palette = "Set1")+
  labs(title = '(a)') +
  theme_apa(box = TRUE)+
  theme(legend.position = c(0.2, 0.8))
p2 <- augment(fit_rs) %>%
  mutate(.std_resid = resid(fit_rs, scaled = TRUE)) %>%
  ggplot(aes(x = sch_sesc, y = .std_resid)) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_smooth(se = FALSE)+
  labs(x = "SES", y = "Level-1 Residuals") +
  scale_color_brewer(palette = "Set1")+
  ggtitle('(b)') +
  theme_apa(box = TRUE)+
  theme(legend.position = c(0.2, 0.8))
p3 <- ggplot(mlbook_lv2, aes(x = sch_iq, y = u0)) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_smooth(se = FALSE)+
  labs(x = "IQ scores", y = "Random Intercept u0") +
  scale_color_brewer(palette = "Set1")+
  ggtitle('(c)') +
  theme_apa(box = TRUE)+
  theme(legend.position = c(0.2, 0.8))
p4 <- ggplot(mlbook_lv2, aes(x = sch_iq, y = u1)) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_smooth(se = FALSE) +
  labs(x = "IQ scores", y = "Random Slope u1") +
  scale_color_brewer(palette = "Set1") +
  ggtitle('(d)') +
  theme_apa(box = TRUE)+
  theme(legend.position = c(0.2, 0.8))
grid.arrange(p1, p2, p3, p4, nrow = 2)




