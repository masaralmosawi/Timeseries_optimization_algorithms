### Optimisation algorithms                                         
### Masar Al-Mosawi                                       

#data
x<-c(0.72, 5.50, -2.21, -0.35, -0.67, 0.16, 23.64, 1.00, 1.06, -495.81, -1.98, 37.72)
m<-length(x)
#Define the functions
cachy_log<-function(x,theta)
{
  s<-(-12*log(pi))
  for (i in 1:m)
  {
    s<-s-log(1+(x[i]-theta)^2)
  }
  s
}
cachy_log_deriv<-function(x,theta)
{
  s<-0
  for (i in 1:m)
  {
    s<-s+2*(x[i]-theta)/(1+(x[i]-theta)^2)
  }
  s
}
cachy_log_deriv2<-function(x, theta)
{
  s<-0
  for (i in 1:m)
  {
    s<-s+(-2+2*(x[i]-theta)^2)/(1+(x[i]-theta^2))^2
  }
  s
}

#test of function
cachy_log(c(1:12),6)
cachy_log_deriv(c(1:12),6)
cachy_log_deriv2(c(1:12),6)


#Plot the functions
par(mfrow=c(2,1))
a<-150
theta<-seq(-a,a,2*a/100)
y_log=cachy_log(x,theta)
y_log_deriv=cachy_log_deriv(x,theta)
plot(theta,y_log, type = 'l', main='Log likelihood')
plot(theta,y_log_deriv, type='l', main='deriviate of log likelihood')
abline(0,0)

#Newton-Raphson method
Newton_raphson <- function(x0, eps=1e-8)
{
  result<-c()
  xt   <- x0
  epsilon <- 999
  while(epsilon>eps)
  {
    xt1 <- xt
    xt <- xt1-cachy_log_deriv(x,xt1)/cachy_log_deriv2(x,xt1)
    epsilon <- sum((xt-xt1)*(xt-xt1))
    result<-c(result,xt)
  }
  result
}
Newton_raphson(0.5)

#Define the function
gxy <- function(x)
{
  (-x[1]^2+10*x[2]-2*x[2]^3+1/2*x[1]^2*x[2])
}

# Produce a contour plot
x1grid <- seq(-10, 10, by=0.05)
x2grid <- seq(-5, 5, by=0.05)
dx1 <- length(x1grid)
dx2 <- length(x2grid)
dx  <- dx1*dx2
gx  <- matrix(rep(NA, dx), nrow=dx1)
for (i in 1:dx1)
  for (j in 1:dx2)
  {
    gx[i,j] <- gxy(c(x1grid[i],x2grid[j]))
  }
mgx <- matrix(gx, nrow=dx1, ncol=dx2)
par(mfrow=c(1,1))
contour(x1grid, x2grid, mgx, nlevels=34)

#3D-plot
gxy_plot <- function(x,y)
{
  (-x^2+10*y-2*y^3+1/2*x^2*y)
}
xx<-seq(-10,10,0.1)
yy<-seq(-10,10,0.1)
par(mfrow=c(1,1))
persp(xx, yy, z=outer(xx,yy,gxy_plot), theta=100, zlab='BE/MV',phi=20, col="green", shade=0.1)

#eigenvalues for Hessian matrices
H1=matrix(c(sqrt(5/3)-2, 0, 0, -12*sqrt(5/3)), nrow=2, ncol=2, byrow=TRUE)
H2=matrix(c(-sqrt(5/3)-2, 0, 0, 12*sqrt(5/3)), nrow=2, ncol=2, byrow=TRUE)
H3=matrix(c(0, 2*sqrt(7), 2*sqrt(7), -12*2), nrow=2, ncol=2, byrow=TRUE)
H4=matrix(c(0, -2*sqrt(7), -2*sqrt(7), -12*2), nrow=2, ncol=2, byrow=TRUE)
eigen(H1)
eigen(H2)
eigen(H3)
eigen(H4)

#Define the functions
#constants
k=100
n=10*k
x<-rep(c(0,0,0,0.1,0.1,0.3,0.3,0.9,0.9,0.9), times=k)
y<-rep(c(0,0,1,0,1,1,1,0,1,1), times=k)

h<-function(b)
{
  s <- exp(b[1]+x*b[2])
  sum(y*log(s)-log(1+s))
}
hd_b1<-function(b)
{
  s <- exp(b[1]+x*b[2])
  sum(y-s/(1+s))
}
hd_b2<-function(b)
{
  s <- exp(b[1]+x*b[2])
  sum(y*x-x*s/(1+s))
}
hgradient <- function(b)
{
  c(hd_b1(b),hd_b2(b))
}

#test function
hd_b1(c(0.2,0.1))
hd_b2(c(0.2,0.1))
h(c(0.2,0.1))

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
    gx[i,j] <- h(c(x1grid[i],x2grid[j]))
  }
mgx <- matrix(gx, nrow=dx1, ncol=dx2)
par(mfrow=c(1,1))
contour(x1grid, x2grid, mgx, nlevels=34)

#Steepest ascent function:
steepestasc <- function(x0, eps=1e-8, alpha0=0.1)
{
  xt   <- x0
  conv <- 999
  points(xt[1], xt[2], col=2, pch=4, lwd=3)
  while(conv>eps)
  {
    alpha <- alpha0
    xt1   <- xt
    xt    <- xt1 + alpha*hgradient(xt1)
    while (h(xt)<h(xt1))
    {
      alpha <- alpha/2
      xt    <- xt1 + alpha*hgradient(xt1)
    }
    points(xt[1], xt[2], col=2, pch=4, lwd=1)
    conv <- sum((xt-xt1)*(xt-xt1))
  }
  points(xt[1], xt[2], col=4, pch=4, lwd=3)
  xt
}
time<-proc.time()
steepestasc(c(-5,3))
proc.time()-time

#GLM
#x<-c(0,0,0,0.1,0.1,0.3,0.3,0.9,0.9,0.9)
#y<-c(0,0,1,0,1,1,1,0,1,1)
glm_model<-glm(y~x, family=binomial)
summary(glm_model)