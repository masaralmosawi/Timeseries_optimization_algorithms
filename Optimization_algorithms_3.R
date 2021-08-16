### Optimisation algorithms                           
### Masar Al-Mosawi  

#Data
x<-c(0,0,0,0.1,0.1,0.3,0.3,0.9,0.9,0.9)
y<-c(0,0,1,0,1,1,1,0,1,1)

L<-function(x,y){
  L<-sum(sqrt(1+x^2))
  L
}
L(x,y)

s<-function(x,y){
  s<-1/length(x)*sum(((y+1)^2)+(y*x+1)^2)
  s
}
s(x,y)

f<-function(b,x,y){
  h1  <- 1/(1+exp(-b[1]-b[2]*x))
  dg1 <- (y - h1)
  dg2 <- (y - h1)*x
  num<-2/(length(x)^2)*(sum(dg1))^2+(sum(dg2)^2)
  den<-L(x,y)*length(x)*sum((dg1^2+dg2^2))
  num/den
}
# Produce a contour plot
x1grid <- seq(-6, 6, by=0.1)
x2grid <- seq(-6, 6, by=0.1)
dx1 <- length(x1grid)
dx2 <- length(x2grid)
dx  <- dx1*dx2
gx  <- matrix(rep(NA, dx), nrow=dx1)
for (i in 1:dx1)
  for (j in 1:dx2)
  {
    gx[i,j] <- f(c(x1grid[i],x2grid[j]),x,y)
  }
mgx <- matrix(gx, nrow=dx1, ncol=dx2)
par(mfrow=c(1,1))
contour(x1grid, x2grid, mgx, nlevels=34)

#Define the functions
g <- function(b,x,y)
{
  # to ensure that h1<1, we add a little 1e-9; otherwise log(1-h1) might crash
  h1 <- 1/(1+exp(-b[1]-b[2]*x)+1e-9)
  sum(log(h1)*y + log(1-h1)*(1-y))
}
gradient <- function(b,x,y)
{
  h1  <- 1/(1+exp(-b[1]-b[2]*x))
  dg1 <- sum((y - h1))/length(x)
  dg2 <- sum((y - h1)*x)/length(x)
  c(dg1,dg2)
}

# Produce a contour plot
x1grid <- seq(-20, 20, by=1)
x2grid <- seq(-40, 20, by=1)
dx1 <- length(x1grid)
dx2 <- length(x2grid)
dx  <- dx1*dx2
gx  <- matrix(rep(NA, dx), nrow=dx1)
for (i in 1:dx1)
  for (j in 1:dx2)
  {
    gx[i,j] <- g(c(x1grid[i],x2grid[j]),x,y)
  }
mgx <- matrix(gx, nrow=dx1, ncol=dx2)
par(mfrow=c(1,1))
contour(x1grid, x2grid, mgx, nlevels=34)

#Quasi-Newton algorithm with BFGS method with stepsize halving
Quasi_newton <- function(x0, eps=1e-8, alpha0=0.1)
{
  #start
  iter<-0
  alpha <- alpha0
  xt   <- x0
  points(xt[1], xt[2], col=2, pch=4, lwd=3)
  Mt <-alpha*diag(2)
  xt1 <- xt-alpha*solve(Mt)%*%gradient(xt,x,y)
  conv <- 999
  while(conv>eps)
  {
    yt <- gradient(xt1,x,y)-gradient(xt,x,y)
    zt <- xt1-xt
    Mt1 <- Mt-(Mt%*%zt%*%t(Mt%*%zt))/(t(zt)%*%Mt%*%zt)[1]+(yt%*%t(yt))/(t(yt)%*%zt)[1]
    #update
    Mt<-Mt1
    xt<-xt1
    xt1 <- xt-alpha*solve(Mt)%*%gradient(xt,x,y)
    #stepsize halving, minimize
    while (g(xt,x,y)>g(xt1,x,y))
    {
      alpha <- alpha/2
      xt <- xt1-alpha*solve(Mt)%*%gradient(xt1,x,y)
    }
    conv <- sum((xt-xt1)*(xt-xt1))
    iter<-iter+1
    points(xt[1], xt[2], col=4, pch=4, lwd=3)
  }
  cbind(iter,xt, conv)
}
Quasi_newton(c(-5,3))

#Steepest descent 
#Steepest ascent function:
iter<-0
steepestasc <- function(x0, eps=1e-8, alpha0=0.1)
{
  xt   <- x0
  conv <- 999
  points(xt[1], xt[2], col=2, pch=4, lwd=3)
  while(conv>eps)
  {
    alpha <- alpha0
    xt1   <- xt
    xt    <- xt1 + alpha*gradient(xt1,x,y)
    while (g(xt,x,y)<g(xt1,x,y))
    {
      alpha <- alpha/2
      xt    <- xt1 + alpha*gradient(xt1,x,y)
    }
    points(xt[1], xt[2], col=2, pch=4, lwd=1)
    conv <- sum((xt-xt1)*(xt-xt1))
    iter<-iter+1
  }
  points(xt[1], xt[2], col=4, pch=4, lwd=3)
  xt
  cbind(iter,xt, conv)
}
time<-proc.time()
steepestasc(c(-5,3))
proc.time()-time

