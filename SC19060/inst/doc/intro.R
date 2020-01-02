## ----eval=TRUE----------------------------------------------------------------
MakeupMissingValue<-function(data)
{
  for(i in seq(ncol(data)))
    if(any(a<-is.na(data[,i])))
    {
      data[a,i]<-median(data[,i],na.rm = T)
    }
  data
}

## ----eval=TRUE----------------------------------------------------------------
Metro.Cauchy <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (1/(pi*(1+y^2))) / (1/(pi*(1+x[i-1]^2))))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      } 
  }
  return(list(x=x, k=k))
}

## -----------------------------------------------------------------------------
ctl<-c(7.82,8.41,7.39,9.02,7.62,8.27,7.93,7.26,8.53,7.74)
trt<-c(6.54,5.73,6.38,6.36,6.93,7.05,5.53,6.75,6.87,5.76)
group<-gl(2,10,20,labels=c("Ctl","Trt"))
weight<-c(ctl,trt)
lm.D9<-lm(weight~group)
summary(lm.D9)$coef

## -----------------------------------------------------------------------------
knitr::kable(head(CO2))

## -----------------------------------------------------------------------------
x<-rnorm(10)
y<-rnorm(10)
sunflowerplot(x,y)


## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))
#plot(CO2)

## -----------------------------------------------------------------------------
n<-1e5
m<-numeric(n)
u<-runif(m)
sigma<-2
x<-sqrt(-2*sigma^2*log(1-u))
hist(x,prob=TRUE,main = expression(f(x)==x/sigma^2*exp(-x^2/(2*sigma^2))))
y<-seq(0,8,.1)
lines(y,y/sigma^2*exp(-y^2/(2*sigma^2)))

## -----------------------------------------------------------------------------
n<-1e5
m<-numeric(n)
u<-runif(m)
sigma<-0.5
x<-sqrt(-2*sigma^2*log(1-u))
hist(x,prob=TRUE,main = expression(f(x)==x/sigma^2*exp(-x^2/(2*sigma^2))))
y<-seq(0,2,.1)
lines(y,y/sigma^2*exp(-y^2/(2*sigma^2)))

## -----------------------------------------------------------------------------
n<-1000
X1<-rnorm(n,0,1)
X2<-rnorm(n,3,1)
Z1<-0.75*X1+0.25*X2
Z2<-0.5*X1+0.5*X2
Z3<-0.1*X1+0.9*X2
#par(mfrow=c(1,3))
hist(Z1);hist(Z2);hist(Z3)

## -----------------------------------------------------------------------------
d<-4
n<-6
x<-c(2,1,1,1,1,2,1,1,1,1,2,1,1,1,1,2)
C<-matrix(x,nrow = 4,ncol = 4)
L<-t(chol(C))
#get a lower triangular
sigma<-L%*%t(L)
#since L is a lower triangular, so sigma can do the Choleski factorization
c1<-sqrt(rchisq(1,6))
c2<-sqrt(rchisq(1,5))
c3<-sqrt(rchisq(1,4))
c4<-sqrt(rchisq(1,3))
n21<-rnorm(1)
n31<-rnorm(1)
n32<-rnorm(1)
n41<-rnorm(1)
n42<-rnorm(1)
n43<-rnorm(1)
h<-c(c1,n21,n31,n41,0,c2,n32,n42,0,0,c3,n43,0,0,0,c4)
T<-matrix(h,nrow = 4,ncol = 4)
#T is a lower triangular 4*4 random matrix
A<-T%*%t(T)
X<-L%*%A%*%t(L)
X

## -----------------------------------------------------------------------------
m<-1e4;t<-runif(m,min = 0,max = pi/3)
theta.hat<-mean(sin(t))*(pi/3)
print(c(theta.hat,cos(0)-cos(pi/3)))

## -----------------------------------------------------------------------------
MC.Phi <- function(R = 10000, antithetic = TRUE) {
u <- runif(R/2)
if (!antithetic) v <- runif(R/2) else
v <- 1 - u
u <- c(u, v)
cdf <- numeric(length(5))
g <- exp(-u)/(1+u^2)
cdf<-mean(g)
cdf
}
set.seed(123)
MC1 <- MC.Phi(anti = FALSE)
set.seed(123)
MC2 <- MC.Phi(anti = TRUE)
print(c(MC1,MC2))

## -----------------------------------------------------------------------------
m <- 1000
MC1 <- MC2 <- numeric(m)
for (i in 1:m) {
MC1[i] <- MC.Phi(R = 1000, anti = FALSE)
MC2[i] <- MC.Phi(R = 1000)
}
print(c(sd(MC1),sd(MC2),(var(MC1) - var(MC2))/var(MC1)))

## -----------------------------------------------------------------------------
u<-c(0,.2,.4,.6,.8,1)
quintile<--log(1-u*(1-exp(-1)))#generate  quintile
round(quintile,2)

## -----------------------------------------------------------------------------
M <- 10000 #number of replicates
T2 <- numeric(5)
estimates <- matrix(0, 10, 1)
g <- function(x) {
exp(-x - log(1+x^2)) * (x > 0) * (x < 1) }
u <- runif(M) #f3, inverse transform method
x <- - log(1 - u * (1 - exp(-1)))
fg <- 5*g(x) / (exp(-x) / (1 - exp(-1)))
for (i in 1:10) {
T2[1] <- mean(g(runif(M/4, 0, .14)))
T2[2] <- mean(g(runif(M/4, .14, .29)))
T2[3] <- mean(g(runif(M/4, .29, .48)))
T2[4] <- mean(g(runif(M/4, .48, 70)))
T2[5] <- mean(g(runif(M/4, .70, 1)))
estimates[i, 1] <- mean(T2)
}
apply(estimates, 2, mean)
apply(estimates, 2, sd)

## -----------------------------------------------------------------------------
n<-20
alpha<-.05
UCL<-replicate(1000,expr = {
  x<-rchisq(n,df=2)
  mean(x)+var(x)/sqrt(n)*qt(1-alpha,df=n-1)
})
sum(UCL>2)
mean(UCL>2)
a<-abs(mean(UCL>2)/.95-1)
#calculate the difference between t-interval results and 0.95 

## -----------------------------------------------------------------------------
n<-20
alpha<-.05
UCL<-replicate(1000,expr = {
  x<-rchisq(n,df=2)
  (n-1)*var(x)/qchisq(alpha,df=n-1)
})
sum(UCL>4)
mean(UCL>4)
b<-abs(mean(UCL>4)/.95-1)
# calculate the difference between the results in example6.4 with 0.95

## -----------------------------------------------------------------------------
print(c(a,b))

## -----------------------------------------------------------------------------
n<-1e4
q<-c(.025,.05,.95,.975)
sk<-replicate(1000,expr = {
x<-rnorm(n)
xbar<-mean(x)
v3<-mean((x-xbar)^3)
v2<-var(x)
v3/(v2^(3/2)) # the skewness 
})
q1<-quantile(sk,q)
q2<-qnorm(q,0,sqrt(6/n))
# estimate the 0.025,0.05,0.95,0.975 quantiles of the skewness under normality
print(round(rbind(q1,q2),5))

## -----------------------------------------------------------------------------
q<-c(.025,.05,.95,.975)
sd<-var<-numeric(4)
for (i in 1:4) {
  var[i]<-q[i]*(1-q[i])/(n*dnorm(q1[i],0,sqrt(6*(n-2)/((n+1)*(n+3))))^2)
  sd[i]<-sqrt(var[i])
}
print(sd)

## -----------------------------------------------------------------------------
#compute the sample skewness statistic
sk<-function(x){
  xbar<-mean(x)
  m3<-mean((x-xbar)^3)
  m2<-mean((x-xbar)^2)
  return(m3/m2^1.5)
}

## -----------------------------------------------------------------------------
alpha<-.1
n<-30 #sample sizes
m<-2500 #number of replicates
a<-c(seq(.1,20,.1)) #the parameter of beta distribution
N<-length(a)
pwr<-numeric(N)
cv<-qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
#critical value  for the skewness test
for (i in 1:N) {
  sktests<-numeric(m)
  for (j in 1:m) {
    x<-rbeta(n,a[i],a[i])
    sktests[j]<-as.integer(abs(sk(x))>=cv)
  }
  pwr[i]<-mean(sktests)
}
#plot power vs a
plot(a,pwr,type = "l",
     xlab = bquote(a),ylim = c(0,1))
abline(h=.1,lty=3)
se<-sqrt(pwr*(1-pwr)/m)
lines(a,pwr+se,lty=3)
lines(a,pwr-se,lty=3)

## -----------------------------------------------------------------------------
alpha<-.1
n<-30 #sample sizes
m<-2500 #number of replicate
v<-c(seq(1,20)) #the parameter of t distribution
N<-length(v)
pwr<-numeric(N)
cv<-qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
#critical value  for the skewness test
for (i in 1:N) {
  sktests<-numeric(m)
  for (j in 1:m) {
    x<-rt(n,v[i])
    sktests[j]<-as.integer(abs(sk(x))>=cv)
  }
  pwr[i]<-mean(sktests)
}
#plot power vs v
plot(v,pwr,type = "l",
     xlab = bquote(v),ylim = c(0,1))
abline(h=.1,lty=3)
se<-sqrt(pwr*(1-pwr)/m)
lines(v,pwr+se,lty=3)
lines(v,pwr-se,lty=3)

## -----------------------------------------------------------------------------
n<-20
alpha<-.05
mu0<-1
m<-10000 #number of replicates
p<-numeric(m) #storage for p-value

## -----------------------------------------------------------------------------
for (j in 1:m) {
  x<-rchisq(n,1) #sample population is Chisq(1)
  ttest<-t.test(x,alternative = "greater",mu=mu0)
  p[j]<-ttest$p.value
}
p.hat<-mean(p<alpha)
se.hat<-sqrt(p.hat*(1-p.hat)/m)
print(c(p.hat,se.hat))

## -----------------------------------------------------------------------------
for (j in 1:m) {
  x<-runif(n,0,2) #sample population is Uniform(0,2)
  ttest<-t.test(x,alternative = "greater",mu=mu0)
  p[j]<-ttest$p.value
}
p.hat<-mean(p<alpha)
se.hat<-sqrt(p.hat*(1-p.hat)/m)
print(c(p.hat,se.hat))

## -----------------------------------------------------------------------------
for (j in 1:m) {
  x<-rexp(n,1) #sample population is Exponential(1)
  ttest<-t.test(x,alternative = "greater",mu=mu0)
  p[j]<-ttest$p.value
}
p.hat<-mean(p<alpha)
se.hat<-sqrt(p.hat*(1-p.hat)/m)
print(c(p.hat,se.hat))

## -----------------------------------------------------------------------------
library(bootstrap)
pairs(scor) 
#display the scatter plots for each pair of test scores
cor(scor) #the sample correlation matrix
library(corrplot)
corrplot(cor(scor),addCoef.col = 'yellow')

## -----------------------------------------------------------------------------
B<-200 
#number of replicates
n<-nrow(scor)
#sample sizes
R<-numeric(B)
for (b in 1:B) {
  i<-sample(1:n,size = n,replace = TRUE)
  mec<-scor$mec[i]
  vec<-scor$vec[i]
  alg<-scor$alg[i]
  ana<-scor$ana[i]
  sta<-scor$sta[i]
  R[b]<-cor(mec,vec)
}
print(se.R<-sd(R))

## -----------------------------------------------------------------------------
R[b]<-cor(alg,ana)
print(se.R<-sd(R))

## -----------------------------------------------------------------------------
R[b]<-cor(alg,sta)
print(se.R<-sd(R))

## -----------------------------------------------------------------------------
R[b]<-cor(ana,sta)
print(se.R<-sd(R))

## -----------------------------------------------------------------------------
sk<-0;n<-1e1;m<-1e3;library(boot);set.seed(12345)
boot.sk <- function(x,i){
  Ubar<-mean(U[i])
  m3<-mean((U[i]-Ubar)^3)
  m2<-mean((U[i]-Ubar)^2)
  R<-m3/m2^1.5
  return(R)
}
ci.norm<-ci.basic<-ci.perc<-matrix(NA,m,2)
for(i in 1:m){
  U<-rnorm(n)
  de <- boot(data=U,statistic=boot.sk, R = 999)
  ci <- boot.ci(de,type=c("norm","basic","perc"))
  ci.norm[i,]<-ci$norm[2:3]
  ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5]
}
cat('norm =',mean(ci.norm[,1]<=sk & ci.norm[,2]>=sk),
'basic =',mean(ci.basic[,1]<=sk & ci.basic[,2]>=sk),
'perc =',mean(ci.perc[,1]<=sk & ci.perc[,2]>=sk))

## -----------------------------------------------------------------------------
sk<-sqrt(8/5);n<-1e1;m<-1e3;library(boot);set.seed(12345)
boot.sk <- function(x,i){
  Ubar<-mean(U[i])
  m3<-mean((U[i]-Ubar)^3)
  m2<-mean((U[i]-Ubar)^2)
  R<-m3/m2^1.5
  return(R)
}
ci.norm<-ci.basic<-ci.perc<-matrix(NA,m,2)
for(i in 1:m){
  U<-rchisq(n,5)
  de <- boot(data = U,statistic=boot.sk, R = 999)
  ci <- boot.ci(de,type=c("norm","basic","perc"))
  ci.norm[i,]<-ci$norm[2:3]
  ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5]
}
cat('norm =',mean(ci.norm[,1]<=sk & ci.norm[,2]>=sk),
'basic =',mean(ci.basic[,1]<=sk & ci.basic[,2]>=sk),
'perc =',mean(ci.perc[,1]<=sk & ci.perc[,2]>=sk))

## -----------------------------------------------------------------------------
library(bootstrap)
data<-scor
n<-nrow(data)
sigma<-cov(data) #compute the covariance matrix
lamda<-eigen(sigma)$values #compute the eigenvalues
theta<-max(lamda)/sum(lamda)

#compute the jackknife replicates, leave-one-out estimates
theta.jack<-numeric(n)
for (i in 1:n) {
    sigma.jack<-cov(data[-i])
    lamda.jack<-eigen(sigma.jack)$values
    theta.jack[i]<-max(lamda.jack)/sum(lamda.jack)
bias<-(n-1)*(mean(theta.jack)-theta)
se.jack<-sqrt((n-1)*mean((theta.jack-theta)^2))
}
print(c(bias,se.jack))

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
    y <- magnetic[-k]
    x <- chemical[-k]
    
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
    e1[k] <- magnetic[k] - yhat1
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
             J2$coef[3] * chemical[k]^2
    e2[k] <- magnetic[k] - yhat2
    
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
    yhat3 <- exp(logyhat3)
    e3[k] <- magnetic[k] - yhat3
    
    J4 <- lm(y ~ x + I(x^2) + I(x^3))
    yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] +
                J4$coef[3] * chemical[k]^2 +
                J4$coef[4] * chemical[k]^3
    e4[k] <- magnetic[k] - yhat4
}

## -----------------------------------------------------------------------------
c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2))
#the estimates of prediction error

## -----------------------------------------------------------------------------
L1<-lm(magnetic ~ chemical)
L2<-lm(magnetic ~ chemical + I(chemical^2))
L3<-lm(log(magnetic) ~ chemical)
L4<-lm(log(magnetic) ~ log(chemical))
round(c(summary(L1)$adj.r.square,summary(L2)$adj.r.square,summary(L3)$adj.r.square,summary(L4)$adj.r.square),4)
#the adjusted R-square

## -----------------------------------------------------------------------------
L2

## -----------------------------------------------------------------------------
maxout <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(as.integer(max(c(outx, outy)) > 5))
}

n1 <- n2 <- 20
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 1000

#The Count 5 test for equal variance
alphahat <- mean(replicate(m, expr={
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
x <- x - mean(x) #centered by sample mean
y <- y - mean(y)
maxout(x, y)
}))

#permutation test for equal variance
R <- 999
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
x <- x - mean(x) #centered by sample mean
y <- y - mean(y)
z <- c(x, y)
n<-length(x)
K <- 1:n
reps <- numeric(R)
t0 <- alphahat
for (i in 1:R) {
k <- sample(K, size = n, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
reps[i] <- mean(maxout(x,y))
}
p <- mean(abs(c(t0, reps)) >= abs(t0))
round(c(t0,p),3)

## -----------------------------------------------------------------------------
dCov <- function(x, y) {
x <- as.matrix(x)
y <- as.matrix(y)
n <- nrow(x)
m <- nrow(y)
if (n != m || n < 2) stop("Sample sizes must agree")
if (! (all(is.finite(c(x, y)))))
stop("Data contains missing or infinite values")
Akl <- function(x) {
d <- as.matrix(dist(x))
m <- rowMeans(d)
M <- mean(d)
a <- sweep(d, 1, m)
b <- sweep(a, 2, m)
return(b + M)
}
A <- Akl(x)
B <- Akl(y)
dCov <- sqrt(mean(A * B))
}

## -----------------------------------------------------------------------------
ndCov2 <- function(z, ix, dims) {
  p <- dims[1]
  q1 <- dims[2] + 1
  d <- p + dims[2]
  x <- z[ , 1:p] 
  y <- z[ix, q1:d]
  return(nrow(z) * dCov(x, y)^2)
}

library(MASS)
sigma<-matrix(c(7,3,3,4),2,2)
x<-mvrnorm(n=100,rep(0,2),sigma)
e<-mvrnorm(n=100,rep(0,2),sigma)
y1<-x/4+e
y2<-x/4*e

library(boot)
library(Ball)
z1<-matrix(rbind(x,y1),100,4)
z2<-matrix(rbind(x,y2),100,4)
boot.obj1 <- boot(data = z1, statistic = ndCov2, R = 999,
sim = "permutation", dims = c(2, 2))
boot.obj2 <- boot(data = z2, statistic = ndCov2, R = 999,
sim = "permutation", dims = c(2, 2))
tb1 <- c(boot.obj1$t0, boot.obj1$t)
tb2 <- c(boot.obj2$t0, boot.obj2$t)

m<-1e2
p.values <- matrix(NA,m,4)
for(i in 1:m){
p.values[i,1] <- mean(tb1 >= boot.obj1$t0)
p.values[i,2] <- mean(tb2 >= boot.obj2$t0)
p.values[i,3] <- bcov.test(x=x,y=y1,R=999,seed=i*54321)$p.value
p.values[i,4] <- bcov.test(x=x,y=y2,R=999,seed=i*54321)$p.value
}
alpha <- 0.1;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
        if (u[i] <= (0.5*exp(-abs(y)) / (0.5*exp(-abs(x[i-1])))))
          x[i] <- y else {
            x[i] <- x[i-1]
            k <- k + 1
          } 
    }
  return(list(x=x, k=k))
}
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
#number of candidate points rejected
print(c(rw1$k, rw2$k, rw3$k, rw4$k))
print(c(rw1$k/N,rw2$k/N,rw3$k/N,rw4$k/N))

## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))
index<-1:N
y1<-rw1$x[index]
y2<-rw2$x[index]
y3<-rw3$x[index]
y4<-rw4$x[index]
plot(index,y1,main="sigma=0.05",type = "l", ylab = "X")
plot(index,y2,main="sigma=0.5",type = "l", ylab = "X")
plot(index,y3,main="sigma=2",type = "l", ylab = "X")
plot(index,y4,main="sigma=16",type = "l", ylab = "X")

## -----------------------------------------------------------------------------
x<-runif(10)
A<-log(exp(abs(x)))
B<-exp(log(abs(x)))
A
B
A==B
isTRUE(all.equal(A,B))

## -----------------------------------------------------------------------------
k<-c(4:25,100,500,1000)
F<-function(a,k){
  M<-2/sqrt(pi*(k-1))*exp(lgamma(k/2)-lgamma((k-1)/2))
  N<-2/sqrt(pi*k)*exp(lgamma((k+1)/2)-lgamma(k/2))
  ck1<-sqrt(a^2*(k-1)/(k-a^2))
  ck2<-sqrt(a^2*k/(k+1-a^2))
  I1<-function(u,k){
    (1+u^2/(k-1))^(-k/2)
  }
  I2<-function(u,k){
    (1+u^2/k)^(-(k+1)/2)
  }
  L<-integrate(I1,lower=0.01,upper=ck1-0.001,rel.tol=.Machine$double.eps^0.25,k=k)$value
  R<-integrate(I2,lower=0.01,upper=ck1-0.001,rel.tol=.Machine$double.eps^0.25,k=k)$value
  AK<-M*L-N*R
}
solution<-numeric(length(k))
for (i in 1:length(k)) {
  f<-function(a) F(a,k[i])
  solution[i]<-uniroot(f,c(0,1.9))$root
}
solution

## -----------------------------------------------------------------------------
blood<-function(p,n.obs){
  n<-sum(n.obs)
  na<-n.obs[1]
  nb<-n.obs[2]
  noo<-n.obs[3]
  nab<-n.obs[4]
  
  pa<-pb<-po<-rep(0,20)
  pa[1]<-p[1]
  pb[1]<-p[2]
  po[1]<-1-p[1]-p[2]
  for (i in 2:20){
    pa.old<-pa[i-1]
    pb.old<-pb[i-1]
    po.old<-po[i-1]
    
    temp1<-pa.old^2+2*pa.old*po.old
    naa<-na*pa.old^2/temp1
    nao<-2*na*pa.old*po.old/temp1
    temp2<-pb.old^2+2*pb.old*po.old
    nbb<-nb*pb.old^2/temp2
    nbo<-2*nb*pb.old*po.old/temp2
    pa[i]<-(2*naa+nao+nab)/(2*n)
    pb[i]<-(2*nbb+nbo+nab)/(2*n)
    po[i]<-(2*noo+nao+nbo)/(2*n)
  }
  return(list(pa=pa,pb=pb,po=po))
}

n.obs<-c(28,24,41,70)
p<-c(1/3,1/3)
a<-blood(p,n.obs)
p<-a$pa
q<-a$pb
r<-a$po
print(cbind(p,q,r))

## -----------------------------------------------------------------------------
plot(c(1:20),p,xlab = "iteration",ylab='p',type="b",main="The interation plot of probability of 'p'")

## -----------------------------------------------------------------------------
plot(c(1:20),q,xlab = "iteration",ylab='q',type="b",main="The interation plot of probability of 'q'")

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
la<-lapply(formulas,lm,data = mtcars) #use lapply()
la

## -----------------------------------------------------------------------------
out<-vector("list",length(formulas))
for (i in seq_along(formulas)) {
  out[[i]]<-lm(formulas[[i]], data = mtcars)
}
out #use loop

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
lapply(bootstraps,lm,formula = mpg ~ disp) #use lapply()

## -----------------------------------------------------------------------------
out<-vector("list",length(bootstraps))
for (i in seq_along(bootstraps)) {
  out[[i]]<-lm(mpg ~ disp, data = bootstraps[[i]])
}
out #use loop

## -----------------------------------------------------------------------------
# exercise 3's lapply()
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
la<-lapply(formulas,lm,data = mtcars)
rsq <- function(mod) summary(mod)$r.squared
sapply(la,rsq)

## -----------------------------------------------------------------------------
# exercise 3's loop
out<-vector("list",length(formulas))
for (i in seq_along(formulas)) {
  out[[i]]<-lm(formulas[[i]], data = mtcars)
}
sapply(out, rsq)

## -----------------------------------------------------------------------------
# exercise 4's lapply()
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
a<-lapply(bootstraps,lm,formula = mpg ~ disp)
sapply(a, rsq)

## -----------------------------------------------------------------------------
# exercise 4's loop
out<-vector("list",length(bootstraps))
for (i in seq_along(bootstraps)) {
  out[[i]]<-lm(mpg ~ disp, data = bootstraps[[i]])
}
sapply(out, rsq)

## -----------------------------------------------------------------------------
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)
#use an anonymous function
sapply(trials,function(x) x[["p.value"]]) 

## -----------------------------------------------------------------------------
# use sapply()
sapply(trials,"[[","p.value")

## -----------------------------------------------------------------------------
#library(parallel)
#boot_df <- function(x) x[sample(nrow(x), rep = T), ]
#rsquared <- function(mod) summary(mod)$r.square
#boot_lm <- function(i) {
#dat <- boot_df(mtcars)
#rsquared(lm(mpg ~ wt + disp, data = dat))
#}
#n <- 1e4
#system.time(mcsapply(1:n, boot_lm, mc.cores = 4))

