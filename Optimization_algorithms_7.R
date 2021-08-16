### Optimisation algorithms                                         
### Masar Al-Mosawi                                       

#Define the functions
g<-function(A,x)
{
  -1/2*x%*%A%*%x
}
g_deriv<-function(A,x)
{
  -A%*%x
}
#Define matrices
A1=matrix(c(8, 1, 1, 8), nrow=2, ncol=2, byrow=TRUE)
A2=matrix(c(1, -2, -2, 12), nrow=2, ncol=2, byrow=TRUE)
A3=matrix(c(1, 0, 0, 100), nrow=2, ncol=2, byrow=TRUE)

#test the function
g_deriv(A3,c(1,9))
g(A3,c(1,9))

# Produce a contour plot
#choose which matrix to contour plot
Ax=A3
x1grid <- seq(-1, 1, by=0.01)
x2grid <- seq(-1, 1, by=0.01)
dx1 <- length(x1grid)
dx2 <- length(x2grid)
dx  <- dx1*dx2
gx  <- matrix(rep(NA, dx), nrow=dx1)
for (i in 1:dx1)
  for (j in 1:dx2)
  {
    gx[i,j] <- g(Ax,c(x1grid[i],x2grid[j]))
  }
mgx <- matrix(gx, nrow=dx1, ncol=dx2)
par(mfrow=c(1,1))
contour(x1grid, x2grid, mgx, nlevels=34)

#Steepest ascent function:
#Algorithm is sensitive of values for large for alpha>1, for Beta=0, we get steepest ascent
steepestasc_polyak <- function(x0, eps=1e-6, alpha=0.033, beta=0.067)
{
  xt2<-c(0,0)
  iter<-0
  xt<-x0
  conv<-999
  points(xt[1], xt[2], col=2, pch=4, lwd=3)
  while(conv>eps)
  {
    xt1<-xt
    xt<-xt1+alpha*g_deriv(Ax,xt1)+beta*(xt1-xt2)
    xt2<-xt1
    points(xt[1], xt[2], col=2, pch=4, lwd=1)
    conv <- norm(xt-xt2,type = "2") #since xt optimal=0, otherwise norm(xt-xt_optimal, type="2")
    iter<-iter+1
    
  }
  points(xt[1], xt[2], col=4, pch=4, lwd=3)
  return(c(iter, conv))
}
time<-proc.time()
steepestasc_polyak(c(-4,2))
proc.time()-time

#eigenvalues for matrix Ax
eigen(A1)
eigen(A2)
eigen(A3)


#Define the functions and the code
#constants
n=1
x<-c(0,0,0,100,100,300,300,900,900,900)
y<-c(0,0,1,0,1,1,1,0,1,1)
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

# Produce a contour plot
x1grid <- seq(-0.02, 0.02, by=0.002)
x2grid <- seq(-0.02, 0.02, by=0.002)
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
steepestasc <- function(x0, eps=1e-16, alpha0=0.0005)
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
steepestasc(c(-5/1000,3/1000)) 
steepestasc(c(-2/1000,3/1000))
#change plot axis to see the search path. I choose different starting values, and starting values larger than these
#the same optimal value as the starting value, not good. These give -0.008008813  0.001260943 and -0.008010555  0.001260939
#Which is close to the optimal wirh GLM. Since the algorithm gives the same optimal value as the starting point
#This means that the point is not traveling, alpha is small enough but the convergence criteria is to large, i edited it 
#to be much smaller 

#Optimal value with GLM
glm_model<-glm(y~x, family=binomial)
summary(glm_model)

#Problem 2.3
h<-function(x,y,b)
{
  s <- exp(b[1]+x*b[2])
  sum(y*log(s)-log(1+s))
}
hd_b1<-function(x,y,b)
{
  s <- exp(b[1]+x*b[2])
  sum(y-s/(1+s))
}
hd_b2<-function(x,y,b)
{
  s <- exp(b[1]+x*b[2])
  sum(y*x-x*s/(1+s))
}
hgradient <- function(x,y,b)
{
  c(hd_b1(x,y,b),hd_b2(x,y,b))
}

#plot
axis<-1
plot(0,0, xlim=c(-axis,axis), ylim=c(-axis,axis))
#Stochastic steepest ascent function:
steepestasc_stoch <- function(x0, eps=1e-6, alpha0=0.001)
{
  xt   <- x0
  conv <- 999
  points(xt[1], xt[2], col=2, pch=4, lwd=3)
  for (i in 1:100000)
  {
    #Uniformly random observations
    n_unif<-sample(1:nrow(my_data),1)
    alpha <- alpha0
    xt1   <- xt
    xt    <- xt1 + alpha*hgradient(my_data[n_unif,"X3"],my_data[n_unif,"X1"],xt1)
    points(xt[1], xt[2], col=2, pch=4, lwd=1)
    conv <- sum((xt-xt1)*(xt-xt1))
    
  }
  points(xt[1], xt[2], col=4, pch=4, lwd=3)
  xt
}
steepestasc_stoch(c(0.2,0.5))

#Read the data into R. First column x (X3), second column y (X1), my_data[2,"X3] second value on column X3
my_data <- read.delim("logist.txt")