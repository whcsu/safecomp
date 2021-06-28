library(MASS)
library(survival)
library(crrp)
library(cmprsk)
library(ggplot2)
library(reshape2)
library(microbenchmark)
my_data<-function(n,p,rho,c0,r){
  #r indicates the percentage of primary events
  beta1<-rep(0,p)
  beta1[1:4]<-c(0.5,0.5,-0.5,0.5)
  beta2<- -beta1
  sig<-matrix(5,p,p)
  sig<-rho^abs(row(sig)-col(sig))
  diag(sig)<-rep(1,p);
  status<-0
  i=1
  while(status==0) 
  {
    Z_temp <- mvrnorm(10*n,rep(-1,p),sig);
    haz1<- drop(Z_temp%*%beta1);
    haz2<- drop(Z_temp%*%beta2);
    u<-runif(10*n);
    a<-(exp(log(1-u)/exp(haz1))-1)/r+1;
    b<-log(1-u/(1-r)^exp(haz2))/exp(haz2)
    ind <- (1:length(a))[a>0&a<1&b<0][1:n];#取haz中大于0的指标且扩充到n
    i=i+1
    print(i)
    if(!is.na(sum(a[ind]))&!is.na(sum(b[ind])))
    {
      status <- 1;
    }
  }
  Z <- Z_temp[ind,];  
  T1<- -log(a[ind]);
  T2<- -b[ind]
  inde<-sample(1:n,n*r);
  T<-T2;
  T[inde]<-T1[inde];
  de<-rep(2,n)
  de[inde]<-1;
  C <- runif(n,min=0,max=c0);
  X<- pmin(T,C)
  status <- as.numeric(T<=C)*de
  X <- pmin(T,C);
  Y <- data.frame(X,status);
  names(Y)<-c("time","status")
  return(list(X=Z,beta=beta,de=de,Y=Y)); 
}

SER <- function(n,p,x,y,lambda,beta){
  sf<-rep(0,p)
  y<-data.frame(y)
  x<-data.frame(x)
  data<-data.frame(y,x)
  ck = apply(x[which(y$status==1),],2,FUN=sum)
  s1=rep(0,p)
  s2=rep(0,p)
  index=rep(0,p)
  data1<-data[which(data$status==1),]
  for(i in 1:nrow(data1)){
    rb1=rbind(data[which(data$time>data1[i,1]),],data[which(data$time<data1[i,1]&data$status==2),])
    rb1=na.omit(rb1)
    rb = rb1[,-c(1:2)]
    minx=apply(rb,2,FUN =min)
    s1=s1+minx
    maxx=apply(rb,2,FUN=max)
    s2=s2+maxx
  }
  s=rbind(ck-s1,s2-ck)
  sf=apply(s,2,max)
  lams=lambda/abs(beta)
  for(i in 1:p){ 
    if(!is.na(sf[i])){if(sf[i]<lams[i]){index[i]=1;m=m+1} }  
  }
  x=x[,which(index==0)]
  return(list(index=index,z=x,m=m))
}
eff_safe<-function(x,y,k,p,m,lambda,weight){
  rej<-rep(0,k)
  scre<-rep(0,k)
  tr<-rep(0,k)
  for(i in 1:k){
    print(i)
    pshfit=crrp(y$time,y$status,x,failcode = 1, cencode = 0,penalty = "LASSO",lambda=lambda[i],penalty.factor=weight, weighted=TRUE)
    tr[i]<-sum(pshfit$beta==0)
    scre[i]<-(p-m[i])/p
    if(tr[i]!=0){rej[i]<-m[i]/tr[i] }
    else{rej[i]=1}
  }
  return(list(rej=rej,scre=scre,tr=tr))
}set.seed(123) 
n=100;p=500;c0=0.2;rho=0.6;r=0.3;k=100
d<-my_data(n,p,rho,c0,r)
data<-data.frame(d$Y,d$X)
data<-data[order(data$time),]
ratio<-sum(data$status==0)/n #censoring ratio 30+%
table(data$status)
print(ratio)
y<-data[,1:2]
x1<-as.matrix(data[,-c(1:2)])
x<-scale(x1)

