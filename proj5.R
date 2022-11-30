setwd("D:/Edinburgh/Statistical Programming/proj5")

lt <- read.table("lt1720uk.dat",header=TRUE)
death <- read.table("death1722uk.dat",header=TRUE)

predict_death <- function(fpop, mpop, mf, mm, d) {
  
  # Number of weeks forward to predict
  n <- length(d)
  
  # Expected proportion dying in a week
  q <- 1-exp(-rbind(mf,mm)/52)
  
  birth <- c(fpop[1],mpop[1])
  N <- rbind(fpop,mpop)
  death_pre <- rep(0,n) # predicted number of deaths each week
  for (week in 1:n) {
    D <- 0.9885*d[week]*q*N
    # Calculate the predicted deaths for the week
    death_pre[week] <- sum(D)
    
    # Update population
    Nd <- N - D
    N <- (51*Nd + cbind(birth,Nd[,-101]))/52
  }
  
  death_pre
}


d1=death$deaths[1:156]
d2=predict_death(lt$fpop17, lt$mpop17, lt$mf, lt$mm, death$d[1:156])
sum(d1-d2)

w22 <- nrow(death)
N=w22-157+1
weeks <- 1:N

ob_death2022 <- death$deaths[157:w22]
pre_death2022 <- predict_death(lt$fpop20, lt$mpop20,
                               lt$mf, lt$mm, death$d[157:w22])
ed2022_w <- ob_death2022 - pre_death2022
ed2022 <- sum(ed2022_w)

ob_death20 <- death$deaths[157:209]
pre_death20 <- predict_death(lt$fpop20, lt$mpop20,
                             lt$mf, lt$mm, death$d[157:209])
ed20 <- sum(ob_death20 - pre_death20)

plot(weeks, ob_death2022, ylim=c(0,max(ob_death2022,pre_death2022)),
     xlab="Week", ylab="Deaths",
     main=sprintf("Deaths in 2020-2022\nExcess deaths: %f (2020)\n%f (Overall)",
                  ed20, ed2022))
lines(weeks, pre_death2022, col="blue")
legend("topright","predicted deaths",lty=1, col="blue")

# ed2022_w <- ob_death2022 - pre_death2022
plot(weeks, cumsum(ed2022_w), xlab="Week", ylab="Cumulative excess deaths",
     main="Cumulative excess deaths in 2020-2022")

library(rjags)
x <- ed2022_w
x[c(51,52,53,105,106)] <- NA

N=w22-157+1
nadapt <- 1000

mod <- jags.model("model.jags", data=list(x=x,N=N), n.adapt=nadapt)
# sam <- jags.samples(mod, c("mu","rho","k"), n.iter=10000)
# str(sam)
# adapt(mod, n.iter=1000)
sam.coda <- coda.samples(mod, c("mu","rho","k"), n.iter=10000)
str(sam.coda)

traceplot(sam.coda[[1]][,"rho"], ylab=expression(rho),
          main="Trace of rho")
# densplot(sam.coda[[1]][,"rho"],type="h")
hist(sam.coda[[1]][,"rho"], probability=TRUE, xlab=expression(rho),
     main="Histogram of rho")
rho <- seq(0, .9, length=50)
lines(rho, dunif(rho,0,.9), col=2)

cols <- sprintf("mu[%d]",1:N)
mu_sam <- sam.coda[[1]][,cols]
mu_exp <- colMeans(mu_sam)
mu_win <- window(mu_sam, start=nadapt+50, thin=50)
weeks <- 1:N


# plot(weeks,ed2022_w,ylim=range(ed2022_w,mu_win,mu_exp),type="n")
# for (i in 1:200) {
#   lines(weeks,mu_sam[i*50,],col="grey")
# }
# lines(weeks,mu_exp,col=" blue",type='l')
# points(weeks,ed2022_w,col=(1:N %in% c(51,52,53,105,106))+1)

matplot(weeks, t(mu_win), type="l", lty=1, col="grey",
        ylim=range(ed2022_w, mu_win, mu_exp),
        xlab="Week", ylab="Excess deaths",
        main="Mean mu of excess deaths in 2020-2022")
# title(main="Mean mu of excess deaths in 2020-2022",
#       xlab="Week", ylab=expression(mu))
lines(weeks, mu_exp, col=" blue", type='l')
points(weeks, ed2022_w, col=(1:N %in% c(51,52,53,105,106))+1)
legend("topright",c("every 50th sampled mu",
                    "estimated expectation for mu",
                    "used observed excess deaths",
                    "unused observed excess deaths"),
       lty=c(1,1,NA,NA),col=c("grey","blue","black","red"), pch=c(NA,NA,1,1))
# legend(legend=c("used observed excess deaths","unused observed excess deaths"),
#        col=c("black","red"),
#        pch=1)


plot(weeks, ed2022_w-mu_exp, xlab="Week", ylab="Residuals")
abline(0, 0, lty=2)