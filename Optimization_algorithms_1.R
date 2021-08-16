### Optimization algorithms                          
### Masar Al-Mosawi
### Packages:none

# Define the objective function
crit <- function(x, w){
  xv <- w[x==1]
  n  <- length(xv)
  X  <- cbind(rep(1, n), xv, xv^2, xv^3)  
  M  <- t(X) %*% X                        
  n*0.2 - log(det(M))
}
w <- seq(-1, 1, by=0.05)
lw <- length(w)

#Simulated annealing algorithm
Sim_ann<-function(m, alpha, tau, stages){
  #staring value
  x_0 <- rbinom(n=lw, size=1, prob=0.5)
  for (j in 1:length(m)){
    for (k in 1:stages){
      tau<-tau*(1-alpha)^(j+k-2)
      for (i in 1:m[j]){
        #Sample candidate x_star in neighborhood of x_0
        x_star<-x_0
        rand_uni<-sample(1:lw,1)
        if (x_0[rand_uni]==1){
          x_star[rand_uni]<-0
        }else{
          x_star[rand_uni]<-1
        }
        #Compute h
        h<-exp((crit(x_0,w)+crit(x_star,w))/(tau))
        #Sample x_{t+1}
        prob<-min(h,1)
        if (prob<runif(1,0,1)){
          x_t<-x_star
        }else{
          x_t<-x_0
        }
        #Update
        x_0<-x_t
      }
    }
  }
  x_t
}
design_optim<-Sim_ann(c(60,120,220), 0.1, 0.9, 5)
#b
w[design_optim==1]
crit(design_optim,w)
#c
result<-matrix(nrow = 30, ncol = 1)
for (i in 1:30) {
  design_optim<-Sim_ann(c(60,120,220), 0.1, 0.9, 5)
  result[i,1]<-crit(design_optim,w)
}
result


