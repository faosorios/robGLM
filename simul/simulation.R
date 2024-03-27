# read R scripts and R packages
source("simul_Poi.R")
source("../code/robLoglin.R")
library(robustbase)

# 0% contamination
z025.00 <- simul.Poi(Nsize = 1000, nobs =  25, eps = 0.0)
z100.00 <- simul.Poi(Nsize = 1000, nobs = 100, eps = 0.0)
z400.00 <- simul.Poi(Nsize = 1000, nobs = 400, eps = 0.0)

# 5% contamination
z025.05.2 <- simul.Poi(Nsize = 1000, nobs =  25, eps = 0.05, nu = 2)
z100.05.2 <- simul.Poi(Nsize = 1000, nobs = 100, eps = 0.05, nu = 2)
z400.05.2 <- simul.Poi(Nsize = 1000, nobs = 400, eps = 0.05, nu = 2)
z025.05.5 <- simul.Poi(Nsize = 1000, nobs =  25, eps = 0.05, nu = 5)
z100.05.5 <- simul.Poi(Nsize = 1000, nobs = 100, eps = 0.05, nu = 5)
z400.05.5 <- simul.Poi(Nsize = 1000, nobs = 400, eps = 0.05, nu = 5)

# 10% contamination
z025.10.2 <- simul.Poi(Nsize = 1000, nobs =  25, eps = 0.10, nu = 2)
z100.10.2 <- simul.Poi(Nsize = 1000, nobs = 100, eps = 0.10, nu = 2)
z400.10.2 <- simul.Poi(Nsize = 1000, nobs = 400, eps = 0.10, nu = 2)
z025.10.5 <- simul.Poi(Nsize = 1000, nobs =  25, eps = 0.10, nu = 5)
z100.10.5 <- simul.Poi(Nsize = 1000, nobs = 100, eps = 0.10, nu = 5)
z400.10.5 <- simul.Poi(Nsize = 1000, nobs = 400, eps = 0.10, nu = 5)

# 25% contamination
z025.25.2 <- simul.Poi(Nsize = 1000, nobs =  25, eps = 0.25, nu = 2)
z100.25.2 <- simul.Poi(Nsize = 1000, nobs = 100, eps = 0.25, nu = 2)
z400.25.2 <- simul.Poi(Nsize = 1000, nobs = 400, eps = 0.25, nu = 2)
z025.25.5 <- simul.Poi(Nsize = 1000, nobs =  25, eps = 0.25, nu = 5)
z100.25.5 <- simul.Poi(Nsize = 1000, nobs = 100, eps = 0.25, nu = 5)
z400.25.5 <- simul.Poi(Nsize = 1000, nobs = 400, eps = 0.25, nu = 5)

# summary
print(rbind(z025.00$bias, z100.00$bias, z400.00$bias), digits = 4)
print(rbind(z025.00$iqr, z100.00$iqr, z400.00$iqr), digits = 4)

print(rbind(z025.05.2$bias, z100.05.2$bias, z400.05.2$bias), digits = 4)
print(rbind(z025.05.5$bias, z100.05.5$bias, z400.05.5$bias), digits = 4)
print(rbind(z025.10.2$bias, z100.10.2$bias, z400.10.2$bias), digits = 4)
print(rbind(z025.10.5$bias, z100.10.5$bias, z400.10.5$bias), digits = 4)
print(rbind(z025.25.2$bias, z100.25.2$bias, z400.25.2$bias), digits = 4)
print(rbind(z025.25.5$bias, z100.25.5$bias, z400.25.5$bias), digits = 4)

print(rbind(z025.05.2$iqr, z100.05.2$iqr, z400.05.2$iqr), digits = 4)
print(rbind(z025.05.5$iqr, z100.05.5$iqr, z400.05.5$iqr), digits = 4)
print(rbind(z025.10.2$iqr, z100.10.2$iqr, z400.10.2$iqr), digits = 4)
print(rbind(z025.10.5$iqr, z100.10.5$iqr, z400.10.5$iqr), digits = 4)
print(rbind(z025.25.2$iqr, z100.25.2$iqr, z400.25.2$iqr), digits = 4)
print(rbind(z025.25.5$iqr, z100.25.5$iqr, z400.25.5$iqr), digits = 4)
