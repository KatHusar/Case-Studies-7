# Setup
library(survival)
source("http://myweb.uiowa.edu/pbreheny/7210/f18/notes/fun.R")
PBC <- subset(pbc, !is.na(trt))
PBC <- PBC[order(PBC$time),]
PBC$trt <- PBC$trt == 1
GVHD <- read.delim("https://s3.amazonaws.com/pbreheny-data-sets/gvhd.txt")

# Residuals: Cox-Snell ----------------------------------------------------

# Slide 4: Construction
fit <- coxph(Surv(time/365.25, status!=0) ~ trt + stage + hepato + bili, PBC)
sfit <- survfit(fit)
H0 <- -log(sfit$surv)
H <- approxfun(c(0, sfit$time), c(0, H0), method='constant')
e1 <- H(fit$y[,1])*exp(fit$linear.predictors)
e2 <- fit$y[,2]-residuals(fit)
head(e1)
head(e2)

# Slide 5: Diagnostic plot
efit <- survfit(Surv(e1, fit$y[,2])~1)
lim <- c(0,2)
plot(efit, fun='cumhaz', mark.time=FALSE, bty='n', conf.int=FALSE, lwd=2, las=1,
     xlab='Residual', ylab='Cumulative hazard', xlim=lim, ylim=lim)
ciband(efit, fun=function(x) -log(x))
lines(lim, lim, col='red', lwd=2)

# Slide 9: Alternative display
plot(log(efit$time), (-log(efit$surv) - efit$time)/efit$std.err, type='l', col='red', lwd=2, bty='n', las=1, ylim=c(-3,3),
     xlab='Residual', ylab=expression((hat(Lambda)(e)-e)/SE), xaxt='n')
at <- seq(-4, 1, len=4)
axis(1, at=at, lab=round(exp(at), 2))
abline(h=c(-2,2), col='gray')

# Slide 6: Simulate PH violation
set.seed(14)
n <- 500
x <- rbinom(n, 1, 0.5)
y <- 2*x + rnorm(n)
fit <- coxph(Surv(exp(y)) ~ x)
e <- fit$y[,2]-residuals(fit)
efit <- survfit(Surv(e, fit$y[,2])~1)
lim <- c(0, 6)
plot(efit, fun='cumhaz', mark.time=FALSE, bty='n', conf.int=FALSE, lwd=2, las=1,
     xlab='Residual', ylab='Cumulative hazard', xlim=lim, ylim=lim)
ciband(efit, fun=function(x) -log(x))
lines(lim, lim, col='red', lwd=2)

# Slide 7: Zoom in
lim <- c(0, 0.25)
plot(efit, fun='cumhaz', mark.time=FALSE, bty='n', conf.int=FALSE, lwd=2, las=1,
     xlab='Residual', ylab='Cumulative hazard', xlim=lim, ylim=lim)
ciband(efit, fun=function(x) -log(x))
lines(lim, lim, col='red', lwd=2)

# Slide 8: Alternative display
plot(log(efit$time), (-log(efit$surv) - efit$time)/efit$std.err, type='l', col='red', lwd=2, bty='n', las=1, ylim=c(-3,3),
     xlab='Residual', ylab=expression((hat(Lambda)(e)-e)/SE), xaxt='n')
at <- seq(-6, 2, 2)
axis(1, at=at, lab=round(exp(at), 2))
abline(h=c(-2,2), col='gray')

# Slide 10: GVHD
fit <- coxph(Surv(Time, Status) ~ Group, GVHD)
e <- fit$y[,2]-residuals(fit)
efit <- survfit(Surv(e, fit$y[,2])~1)
lim <- c(0, max(e))
plot(efit, fun='cumhaz', mark.time=FALSE, bty='n', conf.int=FALSE, lwd=2, las=1,
     xlab='Residual', ylab='Cumulative hazard', xlim=lim, ylim=lim)
ciband(efit, fun=function(x) -log(x))
lines(lim, lim, col='red', lwd=2)
plot(log(efit$time), (-log(efit$surv) - efit$time)/efit$std.err, type='l', col='red', lwd=2, bty='n', las=1, ylim=c(-3,3))

# Slide 11: GVHD, stratified
efit <- survfit(Surv(e, fit$y[,2])~GVHD$Group)
lim <- c(0, max(e))
plot(efit, fun='cumhaz', mark.time=FALSE, bty='n', conf.int=FALSE, lwd=2, col=pal(2), las=1,
     xlab='Residual', ylab='Cumulative hazard', xlim=lim, ylim=lim)
lines(lim, lim, col='gray', lwd=2)

# Residuals: Martingale ---------------------------------------------------

# Slide 15: Martingales for PBC
fit <- coxph(Surv(time/365.25, status!=0) ~ trt + stage + hepato + bili, PBC)
r <- residuals(fit)
plot(PBC$time/365.25, residuals(fit), pch=19, las=1, bty='n', col=pal(2)[fit$y[,2]+1],
     xlab='Time', ylab='Martingale residual')
toplegend(legend=c("Censored", "Died / Liver failure"), pch=19, col=pal(2))
PBC[which.min(residuals(fit)),]
ecdf(PBC$bili)(PBC$bili[which.min(residuals(fit))])
ecdf(fit$linear.predictors)(fit$linear.predictors[which.min(residuals(fit))])

# Residuals: Deviance -----------------------------------------------------

# Slide 19: Construction
fit <- coxph(Surv(time/365.25, status!=0) ~ trt + stage + hepato + bili, PBC)
m <- residuals(fit)
d <- fit$y[,2]
r <- sign(m)*sqrt(-2*(m + d*log(d-m)))
head(r)
head(residuals(fit, type='deviance'))

# Slide 20: Base plot
plot(PBC$time/365.25, r, pch=19, bty='n', las=1, col=pal(2)[fit$y[,2]+1],
     xlab='Time', ylab='Deviance residual')
toplegend(legend=c("Censored", "Died / Liver failure"), pch=19, col=pal(2))

# Slide 22: Outliers
cbind(fit$y, model.matrix(fit), r)[order(r)[1:5],]
cbind(fit$y, model.matrix(fit), r, rank(-fit$linear.predictors))[order(-r)[1:5],]

# Slide 24: vs albumin
plot(PBC$albumin, r, pch=19, bty='n', las=1, col='gray60',
     xlab='Albumin', ylab='Deviance residual')
lines(lowess(PBC$albumin, r), col='red', lwd=3)

# Functional forms --------------------------------------------------------

# Slide 26-27: Albumin
library(visreg)
f <- function(x) {pmin(x, 3.5)}
fit2 <- coxph(Surv(time/365.25, status!=0) ~ trt + stage + hepato + bili + f(albumin), PBC)
f2 <- function(x) {cbind(x, (x-3.5)*(x-3.5 > 0))}
fit3 <- coxph(Surv(time/365.25, status!=0) ~ trt + stage + hepato + bili + f2(albumin), PBC)
visreg(fit2, 'albumin', xlab='Albumin', ylab="Linear predictor")
visreg(fit3, 'albumin', xlab='Albumin', ylab="Linear predictor")
visreg(fit2, 'albumin', xlab='Albumin', ylab="Hazard ratio", trans=exp, ylim=c(0, 8))

# Slide 28: Comparing AIC for albumin
fit0 <- coxph(Surv(time, status!=0) ~ trt + stage + hepato + bili, PBC)
fit1 <- coxph(Surv(time, status!=0) ~ trt + stage + hepato + bili + albumin, PBC)
fit2 <- coxph(Surv(time, status!=0) ~ trt + stage + hepato + bili + f(albumin), PBC)
fit3 <- coxph(Surv(time, status!=0) ~ trt + stage + hepato + bili + f2(albumin), PBC)
AIC(fit0)
AIC(fit1)
AIC(fit2)
AIC(fit3)

# Slide 29: Bilirubin
fit0 <- coxph(Surv(time, status!=0) ~ trt + stage + hepato + f(albumin), PBC)
fit1 <- coxph(Surv(time, status!=0) ~ trt + stage + hepato + bili + f(albumin), PBC)
r0 <- residuals(fit0, type='deviance')
r1 <- residuals(fit1, type='deviance')
plot(PBC$bili, r0, pch=19, bty='n', las=1, col='gray60',
     xlab='Bilirubin', ylab='Deviance residual')
mtext('without bili')
lines(lowess(PBC$bili, r0), col='red', lwd=3)
plot(PBC$bili, r1, pch=19, bty='n', las=1, col='gray60',
     xlab='Bilirubin', ylab='Deviance residual')
mtext('with bili')
lines(lowess(PBC$bili, r1), col='red', lwd=3)

# Slide 30-31: log(Bilirubin)
fit <- coxph(Surv(time, status!=0) ~ trt + stage + hepato + log(bili) + f(albumin), PBC)
visreg(fit, 'bili', xlab='Bilirubin', ylab="Linear predictor")
2^coef(fit)['log(bili)']

# Slide 33: Spline(Bilirubin)
fit <- coxph(Surv(time, status!=0) ~ trt + stage + hepato + pspline(bili) + f(albumin), PBC)
visreg(fit, 'bili', xlab='Bilirubin', ylab="Linear predictor")

# Slide 34: Comparing AIC for bilirubin
fit0 <- coxph(Surv(time, status!=0) ~ trt + stage + hepato + f(albumin), PBC)
fit1 <- coxph(Surv(time, status!=0) ~ trt + stage + hepato + bili + f(albumin), PBC)
fit2 <- coxph(Surv(time, status!=0) ~ trt + stage + hepato + log(bili) + f(albumin), PBC)
fit3 <- coxph(Surv(time, status!=0) ~ trt + stage + hepato + pspline(bili) + f(albumin), PBC)
AIC(fit0)
AIC(fit1)
AIC(fit2)
AIC(fit3)

# Slide 35: Stage
PBC$Stage <- factor(PBC$stage)
fit <- coxph(Surv(time, status!=0) ~ trt + Stage + hepato + log(bili) + f(albumin), PBC)
visreg(fit, 'Stage', ylab="Linear predictor")

# Influence (delta-beta) --------------------------------------------------

# Slide 38: Setup/syntax
fit <- coxph(Surv(time/365.25, status!=0) ~ trt + stage + hepato + log(bili) + f(albumin), PBC)
D <- residuals(fit, type='dfbeta')

# Slide 39: Treatment
plot(D[,1], type='h', col=pal(2)[PBC$trt], bty='n', las=1, lwd=1.5,
     xlab='Index (ordered by time on study)', ylab=expression(hat(Delta)[ij]))
toplegend(legend=1:2, lwd=1.5, col=pal(2))
mtext("Treatment:", line=3)
cbind(fit$y, model.matrix(fit))[order(-D[,1])[1],]
PBC[order(-D[,1])[1],]

# Slide 40: Stage
plot(D[,2], type='h', col=pal(4)[PBC$Stage], bty='n', las=1, lwd=1.5,
     xlab='Index (ordered by time on study)', ylab=expression(hat(Delta)[ij]))
toplegend(legend=levels(PBC$Stage), lwd=1.5, col=pal(4))
mtext("Stage:", line=3)
PBC[order(D[,2])[1],]

# Slide 41: Bili, log scale
xx <- log(PBC$bili)
xx <- (xx - min(xx)) / (diff(range(xx)))
RGB <- colorRamp(c(pal(2)[2], 'gray60', pal(2)[1]))(xx)
col <- apply(RGB, 1, function(x) rgb(x[1], x[2], x[3], max=255))
plot(D[,4], type='h', col=col, bty='n', las=1, lwd=1.5,
     xlab='Index (ordered by time on study)', ylab=expression(hat(Delta)[ij]))
toplegend(legend=c("High bili", "Low bili"), lwd=1.5, col=pal(2))
mtext("Log scale", 3, line=4)
PBC[order(D[,4])[1],]

# Slide 41: Bili (original scale)
fit2 <- coxph(Surv(time/365.25, status!=0) ~ trt + stage + hepato + bili + f(albumin), PBC)
D2 <- residuals(fit2, type='dfbeta')
xx <- PBC$bili
xx <- (xx - min(xx)) / (diff(range(xx)))
RGB <- colorRamp(c(pal(2)[2], 'gray60', pal(2)[1]))(xx)
col <- apply(RGB, 1, function(x) rgb(x[1], x[2], x[3], max=255))
plot(D2[,4], type='h', col=col, bty='n', las=1, lwd=1.5,
     xlab='Index (ordered by time on study)', ylab=expression(hat(Delta)[ij]))
toplegend(legend=c("High bili", "Low bili"), lwd=1.5, col=pal(2))
mtext("Original scale", 3, line=4)


#  https://myweb.uiowa.edu/pbreheny/7210/f19/notes/11-12.pdf