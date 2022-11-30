setwd("D:/Edinburgh/Statistical Programming/proj5")

lt <- read.table("lt1720uk.dat",header=TRUE)
death <- read.table("death1722uk.dat",header=TRUE)

predict_death <- function(fpop,mpop,mf,mm,d) {
  
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
d2=predict_death(lt$fpop17,lt$mpop17,lt$mf,lt$mm,death$d[1:156])

w22 <- nrow(death)
ob_death2022 <- death$deaths[157:w22]
pre_death2022 <- predict_death(lt$fpop20,lt$mpop20,lt$mf,lt$mm,death$d[157:w22])
ed2022 <- sum(ob_death2022 - pre_death2022)

ob_death20 <- death$deaths[157:209]
pre_death20 <- predict_death(lt$fpop20,lt$mpop20,lt$mf,lt$mm,death$d[157:209])
ed20 <- sum(ob_death20 - pre_death20)
plot(157:w22,ob_death2022)
lines(157:w22,pre_death2022)

ed2022_w <- ob_death2022 - pre_death2022
plot(157:w22,cumsum(ed2022_w))

library(rjags)
x <- ed2022_w
x[c(51,52,53,105,106)] <- NA
N=w22-157+1
mod <- jags.model("model.jags",data=list(x=x,N=N))
# sam <- jags.samples(mod,c("mu","rho","k"),n.iter=10000)
# str(sam)
sam.coda <- coda.samples(mod,c("mu","rho","k"),n.iter=10000)
str(sam.coda)

traceplot(sam.coda[[1]][,"rho"])
hist(sam.coda[[1]][,"rho"],probability=TRUE)
rho <- seq(0,.9,length=50)
lines(rho,dunif(rho,0,.9),col=2)

cols <- sprintf("mu[%d]",1:N)
mu_sam <- sam.coda[[1]][,cols]
mu_exp <- colMeans(mu_sam)
mu_win <- window(mu_sam,thin=50)
weeks <- 1:N

plot(weeks,ed2022_w,type="n")
for (i in 1:200) {
  lines(weeks,mu_win[i,],col="grey")
}
lines(weeks,mu_exp,col=" blue",type='l')
points(weeks,ed2022_w,col=" black")

plot(weeks,ed2022_w-mu_exp)
abline(0,0,lty=2)