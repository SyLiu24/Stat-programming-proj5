# Shuying Liu
# s2436365

# The project aims to obtain an excess deaths time series from 2020 onwards and
# analyze the series.

# Excess deaths are the number of deaths over some time, relative to
# what would have been expected given some data from previous period.
# We will first implement a demographic model for England and Wales to generate
# the expected deaths.
# Considering death rates have a seasonal pattern and the ageing process, 
# we build a model with 101 age classes - N_i, i=1,...,101 to represent the
# population in age 0 to age 100, and update N_i by week, to get the expected 
# number of deaths per week. For each age class we have the instantaneous per 
# capita death rate per year m_i. In the absence of seasonality, the expected
# proportion of N_i dying in a week is q_i=1-e^{-m_i/52}. To allow for seasonal
# variation d_j in mortality rates, the proportion dying is d_jq_i in week j.
# So for week j, applying the dying and ageing process, the model becomes
# D_i = d_jq_iN_i, N^*_i = N_i-D_i, N^+_i = 51N^*_i/52+N^*_{i-1}/52, i=1,...,101
# If N_i is the population in age class i at the start of week j, then D_i is 
# the number of deaths for class i, and N^+_i is the population in class i at
# the start of week j+1. Assuming a constant birth rate, we set N^*_0 to the
# initial value of age 0 class, i.e. N^*_0 is always set the value of N_1 at
# the start of simulation. And changing to D_i=0.9885d_jq_iN_i reduces
# the crudeness of the model. Summing up the D_i over age classes we get the
# expected deaths each week. And since m_i is rather different between sexes,
# model will be run separately for each sex and the resulting weekly expected
# deaths summed.
# We iterate this model with the populations from 2020 and the death rates from
# 2017-2019 to predict the expected deaths from 2020 if death rates had stayed
# the same as previous 3 years. Then compare this to the actual deaths per week
# from 2020 to obtain the excess deaths time series. The annual death rates and
# populations of of each age group data are available from the Office for 
# National Statistics (ONS).
# We then then model this series x_i, excess deaths for week i, using a simple
# Bayesian model
# mu_1=alpha, mu_{i+1}=(x_i-alpha)rho+alpha, x_i~t_k(mu_i,tau), i=1,...,n.
# t_k(mu_i,tau) denotes a scaled t distribution with mean mu_i ,precision tau
# and k degrees of freedom. Priors tau~exp(1), rho~U(0,0.9),
# alpha~N(0,tau=0.0001) and k~U(2,100).
# We implement the model in JAGS and draw samples from the posterior densities
# of mu vector, rho and k given x_i data to highlight the parts of the series
# that are most unusual. Note that we will not use x_i for weeks 51, 52, 53, 105
# and 106 for inference as these are Christmas/New Year data with various 
# recording problems.



setwd("D:/Edinburgh/Statistical Programming/proj5")

# Set up source data
# lt1720uk.dat contains columns:
#   age - the age classes 0:0-1 years, 1: 1-2 years, etc
#   fpop17 and mpop17 (fpop20 and mpop20) - the female and male populations in
#   each age class at the start of 2017 (2020)
#   mf and mm - female and male annual death rates for each 1-year age band
lt <- read.table("lt1720uk.dat",header=TRUE)
# death1722uk.dat contains columns:
#   deaths - the number of deaths that week
#   week - the week since the start of 2017, 2020 starts in week 157
#   d - the mortality rate modifier, d_j, value for that week
death <- read.table("death1722uk.dat",header=TRUE)

predict_death <- function(fpop,mpop,mf,mm,d) {
  # Simulate expected deaths by week time series
  
  # Input:
  #   fpop and mpop - female and male starting populations by 1-year age band
  #   mf and mm - female and male annual death rates for each 1-year age band
  #   d - the mortality rate modifier by week
  
  # Returns the predicted total number of deaths each week for length(d) weeks
  
  n <- length(d) # number of weeks forward to predict
  N <- rbind(fpop,mpop) # populations in each age group
  q <- 1-exp(-rbind(mf,mm)/52) # expected proportion dying each week
  birth <- c(fpop[1],mpop[1]) # constant birth
  death_pre <- rep(0,n) # predicted number of deaths each week
  
  # Iterate the demographic model
  for (week in 1:n) {
    D <- 0.9885*d[week]*q*N # deaths in each age group over the week
    # Calculate the predicted deaths for the week
    death_pre[week] <- sum(D)
    
    # Update population
    Nd <- N - D
    N <- (51*Nd + cbind(birth,Nd[,-101]))/52 # dying and ageing process
  }
  death_pre
}


d1=death$deaths[1:156]
d2=predict_death(lt$fpop17, lt$mpop17, lt$mf, lt$mm, death$d[1:156])
sum(d1-d2)

w22 <- nrow(death) # week at the end of data

# Compute excess deaths from the start of 2020 to the end of the data
ob_death2022 <- death$deaths[157:w22]
pre_death2022 <- predict_death(lt$fpop20,lt$mpop20,
                               lt$mf,lt$mm,death$d[157:w22])
ed2022_w <- ob_death2022 - pre_death2022 # excess deaths each week
ed2022 <- sum(ed2022_w)

# Compute excess deaths for 2020
pre_death20 <- predict_death(lt$fpop20,lt$mpop20,
                             lt$mf,lt$mm,death$d[157:209])
ed20 <- sum(death$deaths[157:209] - pre_death20)

# Plot the observed and predicted deaths against week
plot(ob_death2022,ylim=c(0,max(ob_death2022,pre_death2022)),
     xlab="Week",ylab="Deaths",
     main=sprintf("Deaths in 2020-2022\nExcess deaths: %f (2020)\n%f (Overall)",
                  ed20, ed2022))
lines(pre_death2022,col="blue")
legend("topright","predicted deaths",lty=1,col="blue")

# Plot the cumulative excess deaths by week
plot(cumsum(ed2022_w),xlab="Week",ylab="Cumulative excess deaths",
     main="Cumulative excess deaths in 2020-2022")


library(rjags)

# Omit data with recording problems
x <- ed2022_w # excess deaths from 2020
x[c(51,52,53,105,106)] <- NA

# Implement the Bayesian model for excess deaths in JAGS
nadapt <- 1000 # burn-in
wN <- w22-157+1 # total number of weeks from 2020 to the end of data
mod <- jags.model("model.jags",data=list(x=x,N=wN),n.adapt=nadapt)
# Draw 10000 samples form the posterior of parameters mu vector, rho and k
sam.coda <- coda.samples(mod,c("mu","rho","k"),n.iter=10000)

# Investigate the parameters of model
# Produce the trace plot and histogram of posterior for rho
traceplot(sam.coda[[1]][,"rho"],
          ylab=expression(rho),
          main="Trace of rho")
hist(sam.coda[[1]][,"rho"],probability=TRUE,
     xlab=expression(rho),main="Histogram of rho") # posterior for rho
rho <- seq(0,.9,length=50) 
lines(rho,dunif(rho,0,.9),col=2) # prior for rho

# Extract mu samples from mod
cols <- sprintf("mu[%d]",1:wN)
mu_sam <- sam.coda[[1]][,cols]
# Compute the posterior expected values for mu
mu_exp <- colMeans(mu_sam)

# Produce a plot showing every 50th sampled mu, with the estimated
# expectation for mu and observed excess deaths x_i
mu_win <- window(mu_sam,start=nadapt+50,thin=50) # every 50th sampled mu
matplot(t(mu_win),type="l",lty=1,col="grey",
        ylim=range(ed2022_w,mu_win,mu_exp),
        xlab="Week", ylab="Excess deaths",
        main="Mean mu of excess deaths in 2020-2022")
lines(mu_exp,col=" blue",type='l')
points(ed2022_w,col=(1:wN %in% c(51,52,53,105,106))+1)
legend("topright",c("every 50th sampled mu",
                    "estimated expectation for mu",
                    "used observed excess deaths",
                    "unused observed excess deaths"),
       lty=c(1,1,NA,NA),col=c("grey","blue","black","red"),pch=c(NA,NA,1,1))

# Plot the residuals from this model, x_i−E(mu_i), against time
plot(ed2022_w-mu_exp,xlab="Week",ylab="Residuals")
abline(0,0,lty=2)