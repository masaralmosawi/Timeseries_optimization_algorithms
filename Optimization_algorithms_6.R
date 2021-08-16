### Optimisation algorithms                                        
### Masar Al-Mosawi                                       

#Bimodel normal mixture density function
g <- function(x)
{
  s1  <- 0.6
  s2  <- 0.5
  mu1 <- 1.5
  mu2 <- 1.2
  (exp(-(x[1]^2+x[2]^2)/(2*s1))/s1 + exp(-((x[1]-mu1)^2+(x[2]-mu2)^2)/(2*s2))/s2) / (4*pi)
}

# Produce a contour plot; define first a grid where function is evaluated
x1grid <- seq(-1, 3, by=0.05)
x2grid <- seq(-1, 2.5, by=0.05)
dx1 <- length(x1grid)
dx2 <- length(x2grid)
dx  <- dx1*dx2
gx  <- matrix(rep(NA, dx), nrow=dx1)
for (i in 1:dx1)
  for (j in 1:dx2)
  {
    gx[i,j] <- g(c(x1grid[i],x2grid[j]))
  }
mgx <- matrix(gx, nrow=dx1, ncol=dx2)
contour(x1grid, x2grid, mgx, nlevels=40)

#Nelder-Mead algorithm
Nelder_mead<-function(x0,n)
{
  x_best<-c(x0[1],x0[2])
  x_bad<-c(x0[3],x0[4])
  x_worst<-c(x0[5],x0[6])
  #Initialize, 3 vertices
  #x_best<-c(-0.5,1)
  #x_bad<-c(0,1)
  #x_worst<-c(-0.5,0)
  #n=4 #nr of iterations
  alpha<-c(1,2,1/2,1/2)
  
  for (i in 1:n)
  {
    
    x01<-x_best
    x02<-x_bad
    x03<-x_worst
    
    #sort
    m<-matrix(c(x01,x02,x03), ncol=3)
    g_val<-c(g(x01),g(x02),g(x03))
    sort<-order(-g_val)
    x_best<-m[,sort[1]]
    x_bad<-m[,sort[2]]
    x_worst<-m[,sort[3]]
    
    #orient
    c<-1/(3-1)*(x01+x02+x03-x_worst)
    
    #Reflect
    x_r<-c+alpha[1]*(c-x_worst)
    if(g(x_best)>=g(x_r) && g(x_r)>g(x_bad))
    {
      x_worst<-x_r #stop
    } else if(g(x_r)>g(x_best))
    {
      #expansion
      x_e<-c+alpha[2]*(x_r-c)
      if(g(x_e)>g(x_r))
      {
        x_worst<-x_e #stop
      } else
      {
        x_worst<-x_r #stop 
      }
    } else
    {
      #contraction
      if(g(x_bad)>=g(x_r) && g(x_r)>=g(x_worst))
      {
        #outer contraction
        x_o<-c+alpha[3]*(x_r-c)
        if(g(x_o)>=g(x_r))
        {
          x_worst<-x_o #stop
        } else
        {
          #shrinking
          x_bad<-x_best+alpha[4]*(x_bad-x_best)
          x_worst<-x_best+alpha[4]*(x_worst-x_best) #stop
        }
        
      } else
      {
        #inner contraction
        x_i<-c+alpha[3]*(x_worst-c)
        if(g(x_i)>g(x_worst))
        {
          x_worst<-x_i #stop
        } else 
        {
          #shrinking
          x_bad<-x_best+alpha[4]*(x_bad-x_best)
          x_worst<-x_best+alpha[4]*(x_worst-x_best) #stop
        }
      }
    }
    points(x_best[1], x_best[2], col=2, pch=4, lwd=1)
    points(x_bad[1], x_bad[2], col=2, pch=4, lwd=1)
    points(x_worst[1], x_worst[2], col=2, pch=4, lwd=1)
    
    x_tri<-c(x_best[1], x_bad[1], x_worst[1])
    y_tri<-c(x_best[2], x_bad[2], x_worst[2])
    polygon(x_tri,y_tri)
  }
}
Nelder_mead(c(-0.5,1,0,1,-0.5,0),4)

#Define multivariate normal mixture density function, x in R^3
g2<-function(x)
{
  det_4<-det(cbind(c(0.98,0,0), c(0,0.98,0), c(0,0,0.98)))
  w<-1/(sqrt(2*pi)^3)*1/4
  f1<-exp(-1/2*((x[1]-0)^2+(x[2]-0)^2+(x[3]-0)^2))
  f2<-exp(-1/2*((x[1]-2)^2+(x[2]-2)^2+(x[3]-0)^2))
  f3<-exp(-1/2*((x[1]-2)^2+(x[2]-0)^2+(x[3]-2)^2))
  f4<-1/det_4*exp(-1/2*0.98*((x[1]-0)^2+(x[2]-2)^2+(x[3]-2)^2))
  w*sum(f1,f2,f3,f4)
}
#test function
g2(c(1,0,0))

#Optim Nelder-Mead
startvalue1<-c(1,0,0) #local
startvalue2<-c(1,0,2) #local
startvalue3<-c(-5,-5,2) #local
startvalue4<-c(5,5,-2) #global
NM_2<-optim(startvalue1, fn=g2, control = list(fnscale=-1))
round(NM_2$par,3)
NM_2$value

#Particle swarm optimization
pso2<-psoptim(par=rep(NA,3), fn=g2, lower=-2, upper=2, 
              control=list(fnscale=-1, maxit=1000, p=0.11, s=25)) 
round(pso2$par, 3)
pso2$value #matches with my startvalue4 on value and coordinates.

count_global<-0
for (i in 1:10)
{
  pso2_loop<-psoptim(par=rep(NA,3), fn=g2, lower=-2, upper=2, 
                     control=list(fnscale=-1, maxit=1000, p=0.30, s=25)) 
  if (pso2_loop$value < -0.017)
  {
    count_global<-count_global+1
  }
}
count_global

#Read the data into R. First column fertilizer (X1), second column yeild (X)
my_data <- read.delim("cressdata_masar3.txt")
fert<-my_data$X1
yeild<-my_data$X
A<-cbind(rep(1, length(fert)), fert, fert^2, fert^3, fert^4, fert^5)

#Define the function
g3<-function(x,l=0)
{
  A<-cbind(rep(1, length(fert)), fert, fert^2, fert^3, fert^4, fert^5)
  norm_value<-norm(A%*%x-yeild, type='2')^2
  absolut<-l*sum(abs(x))
  norm_value+absolut
}

#Particle swarm optimization
pso3<-psoptim(par=rep(NA,6), fn=g3, lower=-2000, upper=2000, 
              control=list(maxit=5000, p=0.11, s=15)) 
round(pso3$par,1)
pso3$value

#Optim with Nelder-Mead
startvalue_nm2<-c(10,-1,-1,0,1,-17)
NM<-optim(startvalue_nm2, fn=g3)
round(NM$par,1)
NM$value

#Plot the fifth-degree polynomial
xx<-seq(0, 1.2, 0.05)
yy_swarm<- pso3$par[1]+pso3$par[2]*xx+pso3$par[3]*xx^2+pso3$par[4]*xx^3+pso3$par[5]*xx^4+pso3$par[6]*xx^5
yy_neldermead<- NM$par[1]+NM$par[2]*xx+NM$par[3]*xx^2+NM$par[4]*xx^3+NM$par[5]*xx^4+NM$par[6]*xx^5


plot(fert, yeild, main='lamda=0')
lines(xx,yy_swarm, type = 'l', col="red")
lines(xx,yy_neldermead, type='l', col="blue")
legend(1,95, legend=c("PSO", "NM"), col=c("red", "blue"), lty=1:2, cex=0.8)

#Least square solution
sol<-solve(t(A)%*%A)%*%t(A)%*%yeild
g3(sol)