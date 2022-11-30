


mod2 <- jags.model("basic.jags",data=list(x=nhtemp,t=1:length(nhtemp),N=length(nhtemp)),n.adapt = 0,
                   inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed =1989))
sam2.coda <- coda.samples(mod2,c("alpha","beta","tau","df"),n.iter=20)
sam2.coda[,'tau']