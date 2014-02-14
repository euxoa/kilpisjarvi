# Stan modelling and plotting with a time series with potential
# outliers and autocorrelations. 
# The data set here are yearly average summer temperatures from Kilpisjarvi, Finland

d <- read.table("kilpisjarvi-summer-temp.txt", sep=";", header=T)
d <- with(d, data.frame(year=year, temp=temp.summer))

library(ggplot2)
library(reshape2)
par(ask=F)
ggplot(d, aes(x=year, y=temp)) + geom_point() + geom_smooth(method="lm") + theme_bw(18)
ggsave("figure/kilpisjarvi-raw.png") 

# Linear model with R tools
m <- lm(temp ~ year, data=d)
summary(m)
# 95% confidence intervals for decadal trend
# Laplace approximation, gaussian posterior assumed. 
# (This is how all conventional stats programs do it.)
local({p <- .05; z <- qnorm(p/2)
       10*(c(z, -z)*sqrt(diag(vcov(m))["year"]) + coef(m)["year"]) })

library(rstan)
# Data for stan. D is the forthcoming MA model lag.
sdat <- list(N = nrow(d), y = d$temp, D=10)

trend.stats <- function (fit) {
  trends <- extract(fit, "trend")[[1]]
  list (p = 2*sum(trends<0)/length(trends),
        cfi = 10*quantile(trends, c(0.025, 1-0.025))) }

# Linear regression
# Normally 1000—5000 MCMC iterations are enough for a well-mixing chain, 
# but here we use longer chains for accurate tail p-values. 
# (There are also analytical ways to get accurate tail p's; I don't know them.)
m.slm <- stan_model("kj-lm.stan")
system.time(fit <- sampling(m.slm, data=sdat, refresh=500, iter=1000))
trend.stats(fit)

# With t residuals
m.t <- stan_model("kj-t.stan")
fit <- sampling(m.t, data=sdat, chains=1, init=0) # Fast, inaccurate p, overall ok
fit <- sampling(m.t, data=sdat, iter=100000) # Takes some time
trend.stats(fit)

# First MA try. Does not converge too well.
m.c <- stan_model("kj-correlated.stan")

# MA with pooled coefs.
m.c2 <- stan_model("kj-correlated-pooled-theta.stan")

# Check that the model is ok. Chains should be roughly similar.
fit <- sampling(m.c2, data=sdat, iter=3000, init="random", chains=4, seed=5, nondiag_mass=T)
# Long run for an accurate p.
fit <- sampling(m.c2, data=sdat, warmup=30000, iter=100000, init=0, chains=1, nondiag_mass=T)
summary(fit)
trend.stats(fit)

# Whole model, with t residuals and MA.
m.k <- stan_model("kj.stan")

# Quick look.
system.time(fit1q <- sampling(m.k, data=sdat, iter=3000, init=0, chains=1, nondiag_mass=T))
# Does it really converge? Probably yes, but with random inits, some chains get stuck.
# Adding the priors not in the model does not really help.
system.time(fit14 <- sampling(m.k, data=sdat, iter=3000, init="random", chains=4, nondiag_mass=T))
# Initialize at zero. stepsize_jitter allows the sampler to escape.
system.time(fit1 <- sampling(m.k, data=sdat, iter=100000, init=0, chains=1, nondiag_mass=T, 
                             control=list(stepsize_jitter=0.3)))

trend.stats(fit1)
# Plot posteriors. 
ggplot(as.data.frame(extract(fit1, "trend")), aes(x=trend)) + 
  geom_histogram(binwidth=.0002) + theme_bw(18) +  ylab("") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("figure/trend.png")

ggplot(as.data.frame(extract(fit1, "df")), aes(x=df+3)) + 
  geom_histogram(binwidth=1) + theme_bw(18) +  ylab("") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("figure/df.png")

thetas <- as.data.frame(extract(fit1, "theta"))
ggplot(melt(thetas), aes(x=tanh(value))) + 
  geom_vline(xintercept=0, color="red") + geom_histogram(binwidth=.003) +
  facet_wrap(~variable) + 
  theme_bw(18) + xlab("") + ylab("") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("theta-posteriors.png")

# Check residual. There might be something in the middle. 
res <- data.frame(err=d$temp-get_posterior_mean(fit1, "eta")[,1], t=(1:nrow(d)))
ggplot(res, aes(y=err, x=t)) + geom_point() #+ geom_smooth(method="gam")
ggsave("figure/residuals.png")


# Lag-decaying MA, quadratic trend.
m.k2 <- stan_model("kj2.stan")
fit2 <- sampling(m.k2, data=sdat, iter=20000, init=0, chains=1, 
                 nondiag_mass=T, control=list(stepsize_jitter=0.3))

hist(extract(fit2, "k2")[[1]], n=100) # Quadratic trend 
sum(extract(fit2, "k2")[[1]]<0)/10000 # Quadratic trend >0 ?

# Just out of curiosity...
# With the quadratic model, when was the underlying trend likely positive?
d$t <- (1:nrow(d))
coeffs <- as.data.frame(extract(fit2, c("trend", "k2")))
plot(d$year, apply(apply(coeffs, 1, 
                         function (c) 2*c["k2"]*(d$t-nrow(d)/2) + c["trend"])<0, 1, sum)/nrow(coeffs), type="l")


# See how the posterior *residual* model looks like compared to normal distribution.
x <- seq(-5, 5, .1)
df=extract(fit, "df")[[1]]
ggplot(data.frame(x=x, 
                  p=dnorm(x, 0, 1),
                  p2=apply(sapply(df,  function (df) dt(x, df)), 1, mean)), 
       aes(x=x, y=p)) + geom_line() + geom_line(aes(y=p2), color="red")

ggplot(data.frame(df=df), aes(x=df)) + geom_histogram(binwidth=1) + theme_bw(18)
ggsave("figure/residual-posterior.png")
