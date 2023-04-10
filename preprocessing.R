library(tidyverse)
library(dplyr)
library(mice)
library(survival)

# Pre-processing

data = read.csv('dukecathr.csv')

ID.vec = unique(data$RSUBJID)

death.adj = vector(length = nrow(data))
time = vector(length = nrow(data))
prev.cat = vector(length = nrow(data))

count = 0
cutoff = 2*365
for (i in 1:length(ID.vec)){
  if (i %% 300 == 0)
  {
    cat('progress: ', (i/length(ID.vec))*100, '% \n')
  }
  id = ID.vec[i]
  cath.dat = data[data$RSUBJID == id,]
  tot.cath = nrow(cath.dat)
  
  #cat('id: ', id, ' tot.cath: ', tot.cath, '\n')
  for (j in 1:tot.cath)
  {
    count = count + 1
    if (j == 1)
    {
      prev.cat[count] = 0
    }
    if (j > 1)
    {
      t.diff = cath.dat$DAYS2LKA[j-1] - cath.dat$DAYS2LKA[j]
      if (t.diff >= cutoff)
      {
        prev.cat[count] = 1
      }
      if (t.diff < cutoff)
      {
        prev.cat[count] = 2
      }
    }
    if (j == tot.cath)
    {
      time[count] = cath.dat$DAYS2LKA[j]
      death.adj[count] = cath.dat$DEATH[j]
    }
    if (j < tot.cath)
    {
      time[count] = cath.dat$DAYS2LKA[j] - cath.dat$DAYS2LKA[j+1]
      death.adj[count] = 0
    }
  }
}

data$prev_cat = prev.cat
data$death_adj = death.adj
data$survtime = time

# recode: DPCABG, DPMI, DPPCI, 
# collapse: CHFSEV

# weight and height encode as BMI and maybe interact with age 
# recde age to the median of each category


# potentially exclude: CATHAPPR (only take left heart)

# TARGET: DAYS2LKA

# Make binary and set as did not happen if there is a catherization happening before it: DSCABG, DSMI, DSPCI,DSSTROKE,



data$past_CABG = 1*(!is.na(data$DPCABG) & data$DPCABG > 2*365) + 2*(!is.na(data$DPCABG) & data$DPCABG <= 2*365)
data$past_MI = 1*(!is.na(data$DPMI) & data$DPMI > 365) + 2*(!is.na(data$DPMI) & data$DPMI <= 365)
data$past_PCI = 1*(!is.na(data$DPPCI) & data$DPPCI > 365) + 2*(!is.na(data$DPPCI) & data$DPPCI <= 365)


data$CHF_severity = 0*(data$CHFSEV == 0) + 1*(data$CHFSEV == 1 | data$CHFSEV == 2) + 2*(data$CHFSEV == 3 | data$CHFSEV == 4)
data$CHF_severity[is.na(data$CHFSEV)] = NA

data$BMI = data$WEIGHT_R / ((data$HEIGHT_R/100)^2)


data$age = 21*(data$AGE_G == 1) + 27*(data$AGE_G == 2) + 32*(data$AGE_G == 3) + 37*(data$AGE_G == 4) + 
  42*(data$AGE_G == 5) + 47*(data$AGE_G == 6) + 52*(data$AGE_G == 7) + 57*(data$AGE_G == 8) + 
  62*(data$AGE_G == 9) + 67*(data$AGE_G == 10) + 72*(data$AGE_G == 11) + 77*(data$AGE_G == 12) + 
  85*(data$AGE_G == 13)


# Change year into 3 categories: before coronary stents, before drug-eluting stents, and after both
data = data|> mutate(year = case_when(
  YRCATH_G %in% c(1,2) ~ 1, # coronary stents were approved by the FDA in 1994
  YRCATH_G %in% c(3,4) ~ 2, # first drug-eluting stents were approved in 2003
  TRUE ~ 3
))

# DSCABG, DSMI, DSPCI,DSSTROKE

data$subsequent_CABG = 1*(!is.na(data$DSCABG) & (data$DSCABG <= data$survtime))
data$subsequent_MI = 1*(!is.na(data$DSMI) & (data$DSMI <= data$survtime))
data$subsequent_PCI = 1*(!is.na(data$DSPCI) & (data$DSPCI <= data$survtime))
data$subsequent_stroke = 1*(!is.na(data$DSSTROKE) & (data$DSSTROKE <= data$survtime))

####################################

# Definitely exclude: FUPROTCL, RDAYSFROMINDEX, HXCHF (we already have CHF severity), NUMPRMI ?

data_proc = data %>% select(c('survtime', 'death_adj', 'age','year', 'GENDER', 'RACE_G', 'ACS', 'CHF_severity', 'past_CABG', 'past_MI', 'past_PCI', 
                              'HXANGINA', 'HXCEREB', 'HXCOPD', 'HXDIAB', 'HXHTN', 'HXHYL', 'HXMI', 'HXSMOKE',
                              'NUMPRMI', 'DIASBP_R', 'PULSE_R', 'SYSBP_R', 'CBRUITS', 'BMI', 'S3', 
                              'CREATININE_R', 'HDL_R', 'LDL_R', 'TOTCHOL_R',
                              'CATHAPPR', 'DIAGCATH', 'INTVCATH', 'CORDOM',
                              'GRAFTST', 'LADST', 'LCXST', 'LMST', 'LVEF_R',
                              'NUMDZV', 'PRXLADST', 'RCAST'))

data_complete = data_proc %>% na.omit()

head(data_complete)

nrow(data_complete)

data_mice = data_proc %>%
  mutate(GENDER = factor(GENDER), RACE_G = factor(RACE_G), ACS = factor(ACS), year = factor(year),
         CHF_severity = factor(CHF_severity), past_CABG = factor(past_CABG), past_MI = factor(past_MI), 
         past_PCI = factor(past_PCI), HXANGINA = factor(HXANGINA), HXCEREB = factor(HXCEREB), HXCOPD= factor(HXCOPD),
         HXDIAB = factor(HXDIAB), HXHTN = factor(HXHTN), HXHYL = factor(HXHYL), HXMI = factor(HXMI),
         HXSMOKE = factor(HXSMOKE), CBRUITS = factor(CBRUITS), S3 = factor(S3), CATHAPPR = factor(CATHAPPR),
         DIAGCATH = factor(DIAGCATH), INTVCATH = factor(INTVCATH), CORDO = factor(CORDOM))


mice_obj <- mice(data = data_mice, m = 1)
data_full <- complete(mice_obj)



survmodel = survreg(Surv(survtime + 0.1, death_adj)~1 + age + GENDER + RACE_G + year +  ACS + 
                      CHF_severity + past_CABG + past_MI + past_PCI +  
                      HXANGINA + HXCEREB + HXCOPD + HXDIAB + HXHTN +
                      HXHYL  + HXSMOKE + NUMPRMI + DIASBP_R + PULSE_R + SYSBP_R + #+ as.factor(HXMI)
                      + CBRUITS + BMI + S3 + CREATININE_R + HDL_R + LDL_R + TOTCHOL_R + CATHAPPR +
                      DIAGCATH + INTVCATH + CORDOM + GRAFTST + LADST + LCXST + LMST + LVEF_R +
                      NUMDZV + PRXLADST + RCAST, data = data_full, dist='weibull')
summary(survmodel)



coeffs = survmodel$coefficients

low = confint(survmodel, level=.95)[-1,1]
high = confint(survmodel, level=.95)[-1,2]

df_results = data.frame(value = coeffs[-1], Q2.5= low, Q97.5 = high)
df_results$variable = rownames(df_results)

df_results$variable <- factor(df_results$variable, 
                              levels = df_results$variable[order(df_results$value)])
df_signif = df_results[df_results$Q2.5 > 0 | df_results$Q97.5 < 0,]
ggplot(df_signif, aes(x =variable, y = value)) + 
  geom_errorbar(aes(ymax = Q97.5, ymin = Q2.5),width=0.2) +
  geom_point(position = position_dodge(0.9)) +
  geom_hline(yintercept=0, color = "red")+
  ylab("Estimate") +
  xlab("Variable") + 
  labs(title = "95% CI for coefficients")+
  theme_bw()+
  theme( legend.position = "none", axis.title = element_text(size = 20), 
         plot.title = element_text(size = 22),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
 

coxmodel <- coxph(Surv(survtime + 0.1, death_adj) ~ age + GENDER + RACE_G + year +  ACS + 
                    CHF_severity + past_CABG + past_MI + past_PCI +  
                    HXANGINA + HXCEREB + HXCOPD + HXDIAB + HXHTN +
                    HXHYL  + HXSMOKE + NUMPRMI + DIASBP_R + PULSE_R + SYSBP_R + #+ as.factor(HXMI)
                    + CBRUITS + BMI + S3 + CREATININE_R + HDL_R + LDL_R + TOTCHOL_R + CATHAPPR +
                    DIAGCATH + INTVCATH + CORDOM + GRAFTST + LADST + LCXST + LMST + LVEF_R +
                    NUMDZV + PRXLADST + RCAST, data = data_full)
summary(coxmodel)

coeffs2 = coxmodel$coefficients


low2 = confint(coxmodel, level=.95)[,1]
high2 = confint(coxmodel, level=.95)[,2]

df_results_cox = data.frame(value = coeffs2, Q2.5= low2, Q97.5 = high2)
df_results_cox$variable = rownames(df_results_cox)

df_results_cox$variable <- factor(df_results_cox$variable, 
                              levels = df_results_cox$variable[order(df_results_cox$value)])
df_signif_cox = df_results_cox[df_results_cox$Q2.5 > 0 | df_results_cox$Q97.5 < 0,]
ggplot(df_signif_cox, aes(x =variable, y = value)) + 
  geom_errorbar(aes(ymax = Q97.5, ymin = Q2.5),width=0.2) +
  geom_point(position = position_dodge(0.9)) +
  geom_hline(yintercept=0, color = "red")+
  ylab("Estimate") +
  xlab("Variable") + 
  labs(title = "95% CI for coefficients")+
  theme_bw()+
  theme( legend.position = "none", axis.title = element_text(size = 20), 
         plot.title = element_text(size = 22),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


saveRDS(data_full, "data_full.RDS")







### diagnoatics


source("http://myweb.uiowa.edu/pbreheny/7210/f18/notes/fun.R")

sfit <- survfit(coxmodel)
H0 <- -log(sfit$surv)
H <- approxfun(c(0, sfit$time), c(0, H0), method='constant')
e1 <- H(coxmodel$y[,1])*exp(coxmodel$linear.predictors)
e2 <- coxmodel$y[,2]-residuals(coxmodel)
head(e1)
head(e2)

# Slide 5: Diagnostic plot
efit <- survfit(Surv(e1, coxmodel$y[,2])~1)
lim <- c(0,15)
plot(efit, fun='cumhaz', mark.time=FALSE, bty='n', conf.int=FALSE, lwd=1, las=1,
     xlab='Residual', ylab='Cumulative hazard', xlim=lim, ylim=lim)
ciband(efit, fun=function(x) -log(x))
lines(lim, lim, col='red', lwd=1)
