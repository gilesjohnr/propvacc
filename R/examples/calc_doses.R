# Two dose vaccine
calc_doses(0.9, 0.9)
calc_doses(1, 0.5)
calc_doses(0.5, 0.5)

# Three dose vaccine
calc_doses(0.9, 0.9, 0.8)
calc_doses(1, 0.5, 0)
calc_doses(1, 0, 0.5)

# Should be equivalent
all.equal(calc_doses(0.999, 0.5)[,2],
          calc_doses(0.999, 0.5, 0)[2:4,2],
          tolerance=0.02)

all.equal(calc_doses(0.9, 0.8)[,2],
          calc_doses(0.9, 0.8, 0)[2:4,2],
          tolerance=0.02)

# Boundary conditions
calc_doses(0, 0)
calc_doses(1, 0)
calc_doses(1, 0, 0)
calc_doses(1, 1)
calc_doses(1, 1, 0)
calc_doses(1, 1, 1)
calc_doses(1, 0, 1)
calc_doses(0, 1, 1)

# The assumptions may not hold well when v1 << v2 << v3.
# In this case its better to assume independence?
calc_doses(0.1, 0.9)
calc_doses(0.1, 0.5, 0.9)

# Calculate total proportion vaccinated with MCV1 and MCV2 with efficacy
p <- calc_doses(0.9, 0.8)
sum(0.93*p[2,2], 0.97*p[3,2])

# Estimate posterior distribution of proportion vaccinated given uncertainty around MCV1 and MCV2
n <- 1000
sims <- rep(NA, n)
for (i in 1:n) {

     p <- calc_doses(rbeta(1,40,1), rbeta(1,4,2))
     sims[i] <- 0.93*p[2,2] + 0.97*p[3,2]
}

q <- quantile(sims, c(0.025, 0.5, 0.975))

par(mfrow=c(1,1))
hist(sims, breaks=100, col='cyan', xlab='Proportion vaccinated')
abline(v=q[2], lwd=3)
abline(v=q[c(1,3)], lty=2, lwd=2)
