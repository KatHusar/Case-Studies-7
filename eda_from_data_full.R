if (!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, magrittr, ggpubr, latex2exp, xtable, tidybayes, 
               survival, mice, ggsurvfit, gtsummary, survminer)
dir.create("Figure", showWarnings = F)
theme_set(theme_classic(base_size = 12))
# Define colors
purple <- "#756bb1"
orange <- "#ffb50f"

data_full <- readRDS("data_full.RDS")
data.subject <- data_full %>% 
  mutate(DEATH2 = ifelse(death_adj == 0, "Censored", "Dead"),
         DAYS2LKA = survtime, DEATH = death_adj)

PlotKMCurve <- function(group_var, var_name = NULL) {
  if (is.null(var_name)) var_name <- group_var
  survfit2(Surv(DAYS2LKA + 0.1, DEATH) ~ get(group_var), data = data.subject) %>% 
    ggsurvfit() +
    add_confidence_interval() +
    labs(
      x = "Days",
      y = "Survival probability",
      fill = element_blank(), color = element_blank()
    )
}

fig1 <- PlotKMCurve("S3", "Third Heart Sound") +
  labs(subtitle = "Third Heart Sound") +
  scale_color_manual(values = c(purple, orange), labels = c("No (96%)", "Yes (4%)")) +
  scale_fill_manual(values = c(purple, orange), labels = c("No (96%)", "Yes (4%)"))
fig2 <- PlotKMCurve("past_CABG", "Past CABG") +
  labs(subtitle = "Closest Coronary Artery Bypass Surgery") +
  scale_color_manual(values = c("skyblue2", purple, orange), labels = c("Never (76%)", "More than 2 years ago (19%)",  "Within 2 years (5%)")) +
  scale_fill_manual(values = c("skyblue2", purple, orange), labels = c("Never (76%)", "More than 2 years ago (19%)",  "Within 2 years (5%)"))
fig3 <- PlotKMCurve("CATHAPPR", "CATHAPPR") +
  labs(subtitle = "Type of Cardiac Catheterization") +
  scale_color_manual(values = c(orange, "skyblue2", purple, "darkgrey"),
                     labels = c("Unknown (2%)", "Right (2%)",
                                "Left (85%)",
                                "Right and Left (11%)")) +
  scale_fill_manual(values = c(orange, "skyblue2", purple, "black"),
                    labels = c("Unknown (2%)", "Right (2%)",
                               "Left (85%)",
                               "Right and Left (11%)"))

pdf("eda1.pdf", height = 3.5, width = 6)
fig1
dev.off()
pdf("eda2.pdf", height = 3.5, width = 6)
fig2
dev.off()
pdf("eda3.pdf", height = 3.5, width = 6)
fig3
dev.off()

















