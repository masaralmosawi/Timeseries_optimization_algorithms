### Optimisation algorithms                          
### Masar Al-Mosawi   
### Packages: pso

#Define the functions
#Data
x<-cbind(c(3,6), c(-4,11), c(19,8), c(-1,4), c(1,28))
cauchy<-function(x,y,b)
{
  s<-1+(b[1]-x)^2+(b[2]-y)^2
  -sum(log(s)) #minus sign to make it a minimization problem with a maximization algorithm
}

h1<-function(x,y,b)
{
  s1<-b[1]-x
  s2<-1+(b[1]-x)^2+(b[2]-y)^2
  -sum(2*s1/s2) #minus sign to make it a minimization problem with a maximization algorithm
}

h2<-function(x,y,b)
{
  s1<-b[2]-y
  s2<-1+(b[1]-x)^2+(b[2]-y)^2
  -sum(2*s1/s2) #minus sign to make it a minimization problem with a maximization algorithm
}

gradient <- function(x,y,b)
{
  c(h1(x,y,b),h2(x,y,b))
}

#test of functions
test<-c(1,1)
cauchy(3,6,c(1,1))
gradient(3,6,test)

#Plot a contour
x1grid <- seq(-20, 30, by=1)
x2grid <- seq(-10, 40, by=1)
dx1 <- length(x1grid)
dx2 <- length(x2grid)
dx  <- dx1*dx2
gx  <- matrix(rep(NA, dx), nrow=dx1)
for (i in 1:dx1)
  for (j in 1:dx2)
  {
    gx[i,j] <- cauchy(x[1,],x[2,],c(x1grid[i],x2grid[j]))
  }
mgx <- matrix(gx, nrow=dx1, ncol=dx2)
par(mfrow=c(1,1))
contour(x1grid, x2grid, mgx, nlevels=34)

#Stochastic steepest descent with momentum: g(x*)=-g(x*)
steepestasc_stoch <- function(x0, alpha0=0.02, beta=1)
{
  xt2  <-c(0,0)
  xt   <- x0
  conv <- 999
  points(xt[1], xt[2], col=2, pch=4, lwd=3)
  for (i in 1:10000)
  {
    #Uniformly random observations
    n_unif<- sample(1:ncol(x),1)
    alpha <- alpha0
    xt1   <- xt
    xt    <- xt1 + alpha*gradient(x[1,n_unif],x[2,n_unif],xt1)+beta*(xt1-xt2)
    xt2   <- xt1
    while (cauchy(x[1,n_unif],x[2,n_unif],xt)<cauchy(x[1,n_unif],x[2,n_unif],xt1))
    {
      alpha <- alpha/2
      xt    <- xt1 + alpha*gradient(x[1,n_unif],x[2,n_unif],xt1)+beta*(xt1-xt2)
    }
    points(xt[1], xt[2], col=2, pch=4, lwd=1)
    conv <- norm(xt-xt2,type = "2")
    
  }
  points(xt[1], xt[2], col=4, pch=4, lwd=3)
  xt
}
steepestasc_stoch(c(-10,30))

#Re-define the function
f<-function(b)
{
  s<-1+(b[1]-x[1,])^2+(b[2]-x[2,])^2
  sum(log(s)) #minus sign to make it a minimization problem with a maximization algorithm
}

#Simulated annealing
temperature<-seq(0.05,5,by=0.05)
startv <- c(30,-10)
result<-matrix(nrow = length(temperature), ncol = 3)
for (i in 1:length(temperature)) {
  optim_sim_ann<-optim(startv, fn=f, method = "SANN", 
                       control = list(temp=temperature[i], tmax=10, maxit=10000))
  result[i,1]<-optim_sim_ann$par[1]
  result[i,2]<-optim_sim_ann$par[2]
  result[i,3]<-optim_sim_ann$value
}
cbind(result,temperature)

#optim_sim_ann<-optim(startv, fn=f, method = "SANN", control = list(temp=0.5, tmax=10, maxit=10000))
#optim_sim_ann$par
#optim_sim_ann$value
#result
#The temperatues are decreased according to the logarithmic cooling schedule; more specifically, the temperature is set
#to temp/log(((t-1)%/%tmax)*tmax+exp(1)), where t is the current iteration step and tmax are is the nr of function
#evaluations at each temperature
