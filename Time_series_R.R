#Simulate ARMA(1,1)
n<-1000
y <- arima.sim(n=n, model=list(ar=0.5, ma=-0.5), rand.gen=rnorm)#ar=p, ma=-q
plot(y)

iter<-1000
plot(0,0, xlim=c(-0.9, 0.9), ylim=c(-0.9,0.9), xlab="Theta", ylab = "Phi", main = "Likelihood")
for (i in 1:iter){
  y <- arima.sim(n=500, model=list(ar=-0.50, ma=0.5), rand.gen=rnorm)
  arma_ex<-artfima(y, glp = "ARIMA", arimaOrder = c(1, 0, 1),likAlg = c("exact"),  fixd = NULL, b0 = NULL)
  points(arma_ex$thetaHat, arma_ex$phiHat, col=2, pch=4, lwd=1)
}
abline(0,1, col=3)
print.artfima(arma_ex)

plot(0,0, xlim=c(-0.9, 0.9), ylim=c(-0.9,0.9), xlab="Theta", ylab = "Phi", main = "Whittle")
for (i in 1:iter){
  y <- arima.sim(n=500, model=list(ar=-0.50, ma=0.5), rand.gen=rnorm)
  arma_whittle<-artfima(y, glp = "ARIMA", arimaOrder = c(1, 0, 1),likAlg = c("Whittle"),  fixd = NULL, b0 = NULL)
  points(arma_whittle$thetaHat, arma_whittle$phiHat, col=2, pch=4, lwd=1)
}
abline(0,1, col=3)
print.artfima(arma_whittle)

#Contour plot for likelihood for time-domain
#Packages: ltsa: exactLoglikelihood, tacvfARMA
# Produce a contour plot
gxy<-function(beta){
  r <- tacvfARMA(phi=beta[1], theta=beta[2], maxLag=n-1)
  logLL<-exactLoglikelihood(r, y, innovationVarianceQ = TRUE)
  logLL$LL
}

x1grid <- seq(-0.1, 0.1, by=0.005)
x2grid <- seq(-0.1, 0.1, by=0.005)
dx1 <- length(x1grid)
dx2 <- length(x2grid)
dx  <- dx1*dx2
gx  <- matrix(rep(NA, dx), nrow=dx1)
for (i in 1:dx1)
  for (j in 1:dx2)
  {
    gx[i,j] <- gxy(c(x1grid[i],x2grid[j]))
  }
mgx <- matrix(gx-max(gx), nrow=dx1, ncol=dx2)#normalize
par(mfrow=c(1,1))
contour(x1grid, x2grid, mgx, nlevels=100)
abline(0,1, col=3)



