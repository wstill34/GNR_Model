trials <- 1000

#Teaching
xm <- rbinom(trials, 112680, 0.0002)
xl <- rbinom(trials, 112680, 0.000155)
xu <- rbinom(trials, 112680, 0.000245)

par(mfrow = c(4,1))
hist(xl, xlim = c(0, 50))
hist(xm, xlim = c(0, 50))
hist(xu, xlim = c(0, 50))

#Coding
h_vec <- rnorm(trials, mean = 0.0002, sd = 0.000045)
hosp <- rep(NA, times = 1000)

for (j in 1:trials) {
  hosp[j] <- rbinom(1, 112680, h_vec[j])
}


hist(hosp, xlim = c(0, 50))