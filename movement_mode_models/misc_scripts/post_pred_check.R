rwcauchy <- function(n, mu = 0, rho = 0) {
    u = runif(n)
    V = cos(2 * pi * u)
    c = 2 * rho/(1 + rho^2)
    t <- (sign(runif(n) - 0.5) * acos((V + c)/(1 + c * V)) + mu)%%(2 * pi)
    return(t)
}
n <- nrow(obs_1)

single_fit_sims <- rbind(single_fit[[1]], single_fit[[2]], single_fit[[3]])
n.sims <- nrow(single_fit_sims) # number of simulations
pps.single <- matrix(NA, n.sims, n)  # matrix to save steps
ppt.single <- matrix(NA, n.sims, n)  # matrix to save turns

v <- single_fit_sims[4]
lambda <- single_fit_sims[1]
mu <- single_fit_sims[2]
rho <- single_fit_sims[3]

for (i in 1:n.sims) {
    pps.single[i, ] <- rweibull(n, shape = v[i], scale = (1/lambda[i])^(1/v[i])) 
    ppt.single[i, ] <- rwcauchy(n, mu[i], rho[i])
}



ppac.s = matrix(NA, n.sims, 61)
for (i in 1:n.sims) {
    ppac.s[i, ] = acf(log(pps.single[i, ]), lag.max = 60, plot = F)$acf
}

oac = acf(log(l_single), lag.max = 60, plot = F)  #  observed acf

matplot(oac$lag,t(ppac.s),lty=1,pch=1,col='gray')
plot(oac$lag, oac$acf, type = "b", lwd = 2, col = 2, xlab = "lag", ylab = "autocorrelation")
miquantile <- function(x) quantile(x, prob = c(0.025, 0.975))
qs = apply(ppac.s, 2, miquantile)
lines(oac$lag, qs[1, ])
lines(oac$lag, qs[2, ])