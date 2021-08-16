### Optimization algorithms                           
### Masar Al-Mosawi

#PSO unidimensional
PSO_uni <- function(w,c){
  #start
  x_0<-0
  x_1<-1
  res<-0
  for (i in 1:40)
  {
    x_2<-(1+w-c)*x_1-w*x_0
    res<-c(res,x_2)
    x_0<-x_1
    x_1<-x_2
  }
  res
}
Ex_1<-PSO_uni(0.721,1.193)
Ex_2<-PSO_uni(0.9,1.193)
Ex_3<-PSO_uni(0.721,2.2)
Ex_4<-PSO_uni(0.2,3)
Ex_5<-PSO_uni(0.9,2.2)
plot(1:41,Ex_1)
plot(1:41,Ex_2)
plot(1:41,Ex_3)
plot(1:41,Ex_4)
plot(1:41,Ex_5)

#b
PSO_uni_sample <- function(w,c){
  #start
  x_0<-0
  x_1<--2
  res<-0
  for (i in 1:40)
  {
    c_random1<-runif(1,0,c)
    c_random2<-runif(1,0,c)
    x_2<-(1+w-c_random1-c_random2)*x_1-w*x_0
    res<-c(res,x_2)
    x_0<-x_1
    x_1<-x_2
  }
  res
}
PSO_uni_var<-function(w,c){
  result<-rep(0,39)
  for (j in 1:1000){
    result<-cbind(result,PSO_uni_sample(w,c)[3:41])
  }
  result
}
var_1<-apply(PSO_uni_var(0.721,1.193),1,var)
var_2<-apply(PSO_uni_var(0.9,1.193),1,var)
var_3<-apply(PSO_uni_var(0.721,2.2),1,var)
var_4<-apply(PSO_uni_var(0.2,3),1,var)
var_5<-apply(PSO_uni_var(0.9,2.2),1,var)
plot(1:39,var_1)
plot(1:39,var_2)
plot(1:39,var_3)
plot(1:39,var_4)
plot(1:39,var_5)

#c
#Order 1 stable: 1,2,3,5. 4 is not stable and it can be shown in the plot
#Order 2 stable: 1,4. 2,3 can be shown in the plot it is not stable, it is increasing
#But 3 is harder to spot, it looks to be stable but it is not since it does not fulfill
#the condition

Algo <- function(r,s,n){
  x_1<-0
  x_2<-1
  x_3<-2
  x_star<-2.5
  res<-0
  for (i in 1:n)
  {
    #r_random<-runif(1,0,r), since we have E(x), gives E(R)=r/2, E(S)=s/2
    #s_random<-runif(1,0,s)
    x_4<-x_3+r/2*(x_star-x_2)+s/2*(x_star-x_1)
    res<-c(res,x_4)
    x_0<-x_1
    x_1<-x_2
    x_2<-x_3
    x_3<-x_4
  }
  out<-c(res[n+1]-res[n], res[n]-res[n-1], res[n-1]-res[n-2])
  out
}

plot(0,0,xlim = c(0,2), ylim = c(0,2), xlab ="r", ylab="s")
r_seq<-seq(0,2,by=0.01)
s_seq<-seq(0,2,by=0.01)
for (r in r_seq){
  for (s in s_seq){
    if (abs(Algo(r,s,100)[1])<0.005 && abs(Algo(r,s,100)[2])<0.005 && 
        abs(Algo(r,s,100)[3])<0.005){
      points(r, s, col=2, pch=4, lwd=1)
    }
  }
}