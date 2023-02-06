# Libraries Needed:
library(INLA)
library(sf)
# Set your working directory to where HW2data.Rdata and VR.graph are
# located
setwd("./Data")
load("HW2data.Rdata")

# Reading in the expected and observed lung cancer vectors
expected <- Exp.mv3[, "Lung"]
cases <- Obs.mv3[, "Lung"]
# Adding them as columns on VR.cart
VR.cart$expected <- expected
VR.cart$cases <- cases
# Saving it as a dataframe object
data <- VR.cart@data
# Adding SMR column
data$SMR <- data$cases/data$expected
data$geometry <- st_as_sfc(VR.cart)
dmap <- st_as_sf(data)

# Plotting SMR and Estimates:
pal = function(n) brewer.pal(n, "Oranges")
plot(dmap["cases"], nbreaks = 8, pal = pal, breaks = "equal", main = "Observed Lung Cancer Cases",
sub = "Yi for Valencia")

pal = function(n) brewer.pal(n, "Blues")
plot(dmap["expected"], nbreaks = 8, pal = pal, breaks = "equal", main = "Expected Lung Cancer Cases",
sub = "Ei for Valencia")
summary(dmap$expected)

pal = function(n) brewer.pal(n, "Purples")
plot(dmap["SMR"], nbreaks = 8, pal = pal, breaks = "equal", main = "SMRs for Lung Cancer Cases")
summary(dmap$SMR)

ggplot(data.frame(se=sqrt(data$SMR/data$expected),SMR=data$SMR),aes(x=se,y=SMR)) + geom_point() + labs(y="SMR",x="Standard Error", title='Plot of SMR vs Estimated Standard Errors') +
  theme(plot.title = element_text(hjust = 0.5))

# Fit Poisson-Lognormal model in INLA:
model.fit0 <- inla(cases ~ 1 + f(MUNI_ID, model = "iid"), data = dmap, family = "poisson", E = expected) 

# Results for the \beta_0: 
model.fit0$summary.fixed[, 1:5]
# Results for the precision: 
model.fit0$summary.hyper[, 1:5]
# Results for the sigma_e:
sigma <- 1/sqrt(model.fit0$summary.hyper[, 4])
lower <- 1/sqrt(model.fit0$summary.hyper[, 3])
upper <- 1/sqrt(model.fit0$summary.hyper[, 5])
print(c(sigma, lower, upper))
# Extracting the posterior medians using our fitted model
dmap$fit0fitted <- model.fit0$summary.fitted.values$`0.5quant`
# Mapping our results
pal = function(n) brewer.pal(n,"Purples")
plot(dmap["fit0fitted"], pal = pal, nbreaks = 8, breaks = "equal", main="Map of Posterior Medians")

# Find the Standard Errors
se <- sqrt(dmap$SMR/dmap$expected)
# Plotting Posterior Median vs SMR
ggplot(data.frame(pmedian=model.fit0$summary.fitted.values$`0.5quant`,SMR=dmap$SMR),
       aes(y=pmedian,x=SMR, size = se)) + geom_point(alpha = 0.5) + labs(y="Posterior Median",x="SMR") + 
  labs(y="Posterior Median",x="SMR", title='Posterior Median Relative Risk against SMR') +
  theme(plot.title = element_text(hjust = 0.5)) + geom_abline(intercept=0,slope=1,color="red") + xlim(0,2.5) + ylim(0,2.5)

# Find the Standard Errors
se <- sqrt(dmap$SMR/dmap$expected)
# Plotting standard deviations vs standard errors
ggplot(data.frame(psd=model.fit0$summary.fitted.values$`sd`,se=se),
       aes(y=psd,x=se)) + geom_point() + labs(y="Standard Deviation",x="Standard Error", title='Posterior Standard Deviations against Standard Errors of SMRs') +
  theme(plot.title = element_text(hjust = 0.5))

# INLA with BYM2:

# Assigning a region column with a numeric identifier
dmap$region <- 1:nrow(data)
# Running INLA with given specifications
formula <- cases ~ 1 + f(region, model="bym2", graph="VR.graph", scale.model=T, constr=T,
                         hyper=list(phi=list(prior="pc", param=c(0.5, 0.5), initial=1), prec=list(prior="pc.prec", param=c(0.3,0.01),initial=5)))

model.fit1 <- inla(formula, data=dmap, family="poisson",E=expected, control.predictor=list(compute=TRUE),
                   control.compute=list(return.marginals.predictor=TRUE, config = TRUE))
# Beta and Sigma:
model.fit1$summary.fixed[,1:5]
model.fit1$summary.hyper[,1:5]

1/sqrt(model.fit1$summary.hyper[1,4])
model.fit1$summary.hyper[2,4]

# Mapping Relative Risks
dmap$model1fitted <- model.fit1$summary.fitted.values$`0.5quant`
plot(dmap["model1fitted"], pal = pal, nbreaks=8, breaks = "equal", main="Map of Posterior Medians (Spatial BYM2 Model)")

# Extracting the posterior medians using our fitted model
dmap$fit0fitted <- model.fit0$summary.fitted.values$`0.5quant`


