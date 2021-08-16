### Optimisation algorithms                                         
### Masar Al-Mosawi   
### Packages: pso

#Define the function, matrix and vectors
g<-function(s)
{
  x<-c(s[1], s[2], s[3], s[4])
  w<-c(s[5], s[6], s[7])
  A1<-c(1,x[1],x[1]^2,x[1]^3)%*%t(c(1,x[1],x[1]^2,x[1]^3))*w[1]
  A2<-c(1,x[2],x[2]^2,x[2]^3)%*%t(c(1,x[2],x[2]^2,x[2]^3))*w[2]
  A3<-c(1,x[3],x[3]^2,x[3]^3)%*%t(c(1,x[3],x[3]^2,x[3]^3))*w[3]
  A4<-c(1,x[4],x[4]^2,x[4]^3)%*%t(c(1,x[4],x[4]^2,x[4]^3))*(1-w[3]-w[2]-w[1])
  det(A1+A2+A3+A4)
  
}
#test of function
#<-c(1,2,3,4, -100, 1/2, 1/3)
#g(t)

U <- matrix(0, nrow=12, ncol=7)
U[1,1] <- U[3,2] <- U[5,3] <- U[7,4] <- U[9,5] <- U[10,6] <- U[11,7] <- 1
U[2,1] <- U[4,2] <- U[6,3] <- U[8,4] <- U[12,5] <- U[12,6] <-U[12,7] <- -1
d <- c(rep(c(0, -10), 4), 0, 0, 0, -1)
startv <- c(1, 2, 4, 7, 0.2, 0.1, 0.3)

# ConstrOptim: Nelder-Mead as inner optimisation method:
res <- constrOptim(startv, f=g, grad=NULL, ui=U, ci=d, 
                   control=list(fnscale=-1))
round(res$par, 3) 
g(res$par)
#Yes the values are in the range and for different 
#input values a got the same result, so yes it seems reasonale

g_tilde<-function(s, mu=0.001)
{
  b<-U[c(1:dim(U)[1]),]%*%s-d #log(b) crashes since b<0
  g(s)+mu*sum(log(b))
}

#test of function
#startv <- c(1, 2, 4, 7, 0.2, 0.1, 0.3)
#g_tilde(startv)

#Non contraint optimization
#Particle swarm optimization
startv <- c(3, 2, 4, 7, 0.2, 0.1, 0.3)
optim_tilde<-psoptim(par=rep(NA,7), fn=g_tilde, lower=c(0,0,0,0,0,0,0),
                     upper=c(10,10,10,10,0.3,0.3,0.3), 
                     control=list(fnscale=-1, maxit=1000, p=0.11, s=10)) 
g(optim_tilde$par)
round(optim_tilde$par,2)

#Nelder-mead algorithm, default, testing a different algorithm for the sake of it
startv_nelder_mead <- c(5, 2, 3, 6, 0.1, 0.4, 0.3)
optim_tilde_nelder_mead<-optim(startv_nelder_mead, f=g_tilde, 
                               control = list(fnscale=-1))
g(optim_tilde_nelder_mead$par)
round(optim_tilde_nelder_mead$par,2)

#Read the data into R. First column fertilizer (X1), second column yeild (X)
my_data <- read.delim("cressdata_masar3.txt")
fert<-my_data$X1
yeild<-my_data$X

A<-cbind(rep(1, length(fert)), fert, fert^2, fert^3)
#Define the function
g2<-function(x)
{
  A<-cbind(rep(1, length(fert)), fert, fert^2, fert^3)
  norm_value<-norm(A%*%x-yeild, type='2')^2
  norm_value
}

#test of function
g2(c(0.1,0.1,0.1, 0.2))

t<-10000
U_2<-cbind(rep(c(0),8), rep(c(-1,1),4), rep(c(-1,-1,1,1),2),c(-1,-1,-1,-1,1,1,1,1))
d_2<-t*(rep(c(-1),8))

#log barrier function
g_tilde_2<-function(s, mu=100)
{
  b<-U_2[c(1:dim(U_2)[1]),]%*%s-d_2
  g2(s)+mu*sum(log(b))
}
#test of function
#g_tilde_2(c(0.1,0.1,0.1,0.1))

#gradient function
g2_gradient<-function(x)
{
  #grad<-2*c(sum(A%*%x-yeild), fert%*%(A%*%x-yeild), fert^2%*%(A%*%x-yeild), fert^3%*%(A%*%x-yeild))
  #Gradient can also be written as 2*t(A)%*%(A%*%x-yield), better than above in efficiency 
  grad<-2*t(A)%*%(A%*%x-yeild)
  grad
}
#test of function
g2_gradient(c(1,1,1,1))

# ConstrOptim: Nelder-Mead as inner optimisation method:
startv_2<-c(1, 2, 3, 4)
res_2 <- constrOptim(startv_2, f=g2, grad=NULL, ui=U_2, ci=d_2)
round(res_2$par, 2) 
g2(res_2$par)

# ConstrOptim: BFGS as inner optimization method:
startv_2_2<-c(1, 2, 3, 4)
res_2_2 <- constrOptim(startv_2_2, f=g2, ui=U_2, ci=d_2, 
                       method = "BFGS", grad=g2_gradient)
round(res_2_2$par, 2) 
res_2_2$value

#Nelder-mead algorithm, default, testing a different algorithm for the sake of it
startv_2_3<-c(10, 2, -3, 40)
optim_tilde_nelder_mead<-optim(startv_2_3, f=g_tilde_2)
round((optim_tilde_nelder_mead$par),2)
g2(optim_tilde_nelder_mead$par)
#I assume the errors are that the points are hitting the barrier and is okey

#least square estimation, for t->+inf
sol2<-solve(t(A)%*%A)%*%t(A)%*%yeild