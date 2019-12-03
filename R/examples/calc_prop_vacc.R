# Two dose measles vaccine: routine without SIA
calc_prop_vacc(V=c(0.9, 0.8),
               effectiveness=c(0.84, 0.941),
               independent=FALSE)

calc_prop_vacc(V=c(0.9, 0.8),
               effectiveness=c(0.84, 0.941),
               independent=TRUE)

# Three dose vaccine: routine without SIA
calc_prop_vacc(V=c(0.9, 0.8, 0.7),
               effectiveness=c(0.85, 0.9, 0.95))

# Estimate posterior distribution of proportion vaccinated given uncertainty around MCV1 and MCV2
sims <- foreach(i=1:5000, .combine='c') %do% {

  p <- calc_prop_vacc(V=c(rbeta(1,40,1), rbeta(1,4,2)),
                      effectiveness=c(0.85, 0.94),
                      independent = F)
}

q <- quantile(sims, c(0.025, 0.5, 0.975))

par(mfrow=c(1,1))
hist(sims, breaks=100, col='cyan', xlab='Proportion vaccinated')
abline(v=q[2], lwd=3)
abline(v=q[c(1,3)], lty=2, lwd=2)

