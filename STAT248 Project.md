**Disclaimer**

Some of the `R` codes in this exercise are inspired and adapted from Prof. Aditya Guntuboyina's lectures on STATISTICS 248: ANALYSIS OF TIME SERIES of UC Berkeley Department of Statistics.



# Part 1: simulation study

### simulate some data

Define the function of simulation data

```r
len=1000;mu = -5;phi = 0.8;sigma = 0.58
h0 = rnorm(1, mu, sd = sigma)
h = rep(NA,len+1)
h[1] = h0
y0= rnorm(1,0,sd=exp(h0/2))
y = rep(NA,len+1)
y[1] = y0
for( i in 2:(len+1)){
  h[i] = rnorm(1,mean=mu+phi*(h[i-1]-mu),sigma)
  y[i] = rnorm(1,0,exp(h[i]/2))
}
sv_true = exp(h[-1]/2)
y_true = y[-1]

sv_true = as.numeric(sv_true$x)
y_true = as.numeric(y_true$x)
```

Plot data

```r
par(mar=c(5, 4, 4, 6) + 0.1)
plot(sv_true, type="l", yaxt="n",
     xlab="Time", ylab="", col="blue", main='stochastic volatility v.s. Return^2')
axis(2, at= pretty(range(sv_true),10),
     col.ticks="black", col.axis="black")
points(y_true^2, type="l", col="red")
legend("top",legend=c("stochastic volatility","Return^2"),
       text.col=c("blue","red"),pch=c(15,15),col=c("blue","red"))
```



### Particle filtering and FFBS

Define bootstrap particle filter

```R
bootstrap.filter = function(sigev, phi,mu, data, N = 10000)
{
  T = (length(data) - 1)
  samples = matrix(NA, N, (T+1))
  loglikvec = rep(-1, (T+1))
  h0tilde = rnorm(N, mu, sd = sigev) #generation
  sigeps = sqrt(exp(h0tilde))
  logwts0 = ((((data[1])^2)*(-1))/(2*sigeps^2)) - ((0.5)*(log(2*pi*sigeps^2)))
  M = max(logwts0)
  logwts0 = logwts0 - M
  wts0 = exp(logwts0)
  loglikvec[1] = log(mean(wts0)) + M
  wts0 = wts0/(sum(wts0)) #weights normalization
  resmpls = sample(1:N, prob = wts0, size = N, replace = T)
  samples[,1] = h0tilde[resmpls]
  prev.samples = samples[,1]
  for(t in 1:T)
  {
    #new.samples = prev.samples + rnorm(N, mean = 0, sd = sigev) #generation
    new.samples = mu+phi*(prev.samples-mu)+rnorm(N, mean = 0, sd = sigev)
    #wts = dnorm(data[(t+1)], mean = new.samples, sd = sigeps) #weights
    sigeps = exp(new.samples/2)
    logwts = ((((data[(t+1)])^2)*(-1))/(2*sigeps^2)) - ((0.5)*(log(2*pi*sigeps^2)))
    M = max(logwts)
    logwts = logwts - M
    wts = exp(logwts)
    loglikvec[(t+1)] = log(mean(wts)) + M
    wts = wts/(sum(wts)) #normalization
    resmpls = sample(1:N, prob = wts, size = N, replace = T) #resampling
    samples[,(t+1)] = new.samples[resmpls] #resampling
    prev.samples = samples[,(t+1)]
  }
  loglik = sum(loglikvec)
  ret.list = list("samples" = samples, "loglik" = loglik)   

}

# padding function for MLE grid search
pad_seq = function(x,len,inc){
  left_vec = seq(from=x-len*inc,to=x,by=inc)
  right_vec = seq(from=x,to=x+len*inc,by=inc)
  out = append(left_vec,right_vec[-1])
  return(out)
}
```

Numerical grid search:

```r
mu_grid = pad_seq(-5.02,2,0.01);mu_grid
phi_grid = pad_seq(0.76,2,0.002);phi_grid
sigev_grid = pad_seq(0.541,2,0.002);sigev_grid


LL_res = rep(NA,length(mu_grid)*length(phi_grid)*length(sigev_grid))
param_grid = expand.grid(mu_grid,phi_grid,sigev_grid)
param_grid$LL = NA

for (i in 1:nrow(param_grid)){
  bf.samples = bootstrap.filter(sigev=param_grid[i,3], 
                                phi=param_grid[i,2],
                                mu=param_grid[i,1], 
                                data=y,
                                N=1000)
  print(i)
  param_grid$LL[i] = bf.samples$loglik
}

param_grid[which.max(param_grid$LL),]

```

FFBS of particle filter

```r
particle_FFBS = function(sigev, phi, mu,data, N = 10000, nsmoo = 50)
{
  T_ = (length(data) - 1)
  bf = bootstrap.filter(sigev = sigev, phi = phi, mu=mu,data, N = N)
  filt.samples = bf$samples
  smoo.samples = matrix(NA, nsmoo, (T_+1))
  smoo.samples[,(T_+1)] = filt.samples[sample(1:N, size = nsmoo), (T_+1)]
  for(i in 1:nsmoo)
  {
    Xtplus1 = smoo.samples[sample(1:nsmoo)[1], (T_+1)]
    for(t in T_:1)
    {
      #wts = dnorm(filt.samples[,t], mean = mu+Xtplus1, sd = sigev)
      error_t = Xtplus1-(mu+phi*(filt.samples[,t]-mu))
      logwts = (((error_t^2)*(-1))/(2*sigev^2)) - (0.5)*(log(2*pi*sigev^2)) 
      M = max(logwts)
      logwts = logwts - M
      wts = exp(logwts)
      wts = wts/(sum(wts)) #normalization
      Xtplus1 = filt.samples[sample(1:N, prob = wts, size = 1, replace = TRUE), t]
      smoo.samples[i, t] = Xtplus1
    }
    print(i)
  }
  ret.list = list("ffbs.samples" = smoo.samples)   
  return(ret.list)
}

ffbs_result = particle_FFBS(sigev=0.543, 
                            phi=0.76,
                            mu=-5.03, 
                            data=y_true,
                            N=5000, 
                            nsmoo = 500)

ffbs.smoo = ffbs_result$ffbs.samples
sv_smoothed = exp(apply(ffbs.smoo, 2, mean)/2)

```

Plot smoothed stochastic volatility as well as particle degeneracy

```R
plot(sv_true, type="l", yaxt="n",
     xlab="Time", ylab="", col="blue", main='true SV versus particle smoothed CV')
axis(2, at= pretty(range(sv_true),10),
     col.ticks="black", col.axis="black")
points(sv_smoothed, type="l", col="red")
legend("top",legend=c("true SV","particle smoothed CV"),
       text.col=c("blue","red"),pch=c(15,15),col=c("blue","red"))

T_ = length(y_true)-1
num.uniq = rep(-1, (T_+1))
for(t in 0:T_)
{
  num.uniq[(t+1)] = length(unique(bf.samples$samples[,(t+1)]))
}
plot(num.uniq, type = "l")
num.uniq
min(num.uniq)

```

### MCMC via Winbugs Gibbs sampler

Loading R2Winbugs and coda for interfacing between R and Winbugs

```R
library(R2WinBUGS)
library(coda)
```

Hand-code the probability structure of Winbugs

```r
svar1 = function()
{ # prior distribution
  mu ~ dnorm(0,0.01)
  phistar ~ dbeta(5,2)
  sigma ~ dgamma(0.5,5)
  
  phi <- 2*phistar-1
  tau <- 1/sigma
  std <- sqrt(sigma)
  
  # initialize
  var0 <- tau/(1-pow(phi,2))
  h_mu0 ~ dnorm(mu,var0)
  h_mu[1] <- mu + phi*(h_mu0-mu)
  h[1] ~ dnorm(h_mu[1], tau)
  
  for (t in 2:N)
  { 
    h_mu[t] <- mu + phi*(h[t-1]-mu)
    h[t] ~ dnorm(h_mu[t],tau)
    sigma_y[t] <- exp(h[t]/2)
    tau_y[t] <- 1/exp(h[t])
    y[t] ~ dnorm(0,tau_y[t])
  }
}

write.model(svar1, 'svar1.bug')

```

Run Winbugs with 2 chains, 10000 iterations, 5000 burn-in and thinning of 10:

```r
N = length(y_true)
data <- list(y=y_true,N=N)
inits <- function(){list(mu=0,phistar=0.9,sigma=0.1)}
param = c("mu","phi","std","sigma_y")

svar_mod = bugs(data=data,
                 inits=inits,
                 model.file="svar1.bug",
                 parameters=param,
                 n.chains=2,
                 n.iter=10000,
                 n.burnin=5000,
                 n.thin=10,
                 bugs.directory="D:/WinBUGS14/",
                 debug=TRUE,
                 DIC = FALSE,
                 codaPkg=FALSE)

svar_summary = svar_mod$summary
round(svar_summary[1:5,],4)
```

Plot traceplots and density plots

```R
par(mfrow=c(3,1))
plot(rowMeans(svar_mod$sims.array[,,"mu"]),type='l',ylab='mu')
plot(rowMeans(svar_mod$sims.array[,,"phi"]),type='l',ylab='phi')
plot(rowMeans(svar_mod$sims.array[,,"std"]),type='l',ylab='sigma')

par(mfrow=c(1,3))
plot(density(rowMeans(svar_mod$sims.array[,,"mu"])),type='l',main='mu',ylab="")
plot(density(rowMeans(svar_mod$sims.array[,,"phi"])),type='l',main='phi',ylab="")
plot(density(rowMeans(svar_mod$sims.array[,,"std"])),type='l',main='sigma',ylab="")


sv_mcmc = svar_summary[4:nrow(svar_summary),1]
plot(sv_true[-1], type="l", yaxt="n",
     xlab="Time", ylab="", col="blue", main='true SV versus MCMC smoothed SV')
axis(2, at= pretty(range(sv_true),10),
     col.ticks="black", col.axis="black")
points(sv_mcmc, type="l", col="red")
legend("top",legend=c("true SV","MCMC smoothed SV"),
       text.col=c("blue","red"),pch=c(15,15),col=c("blue","red"))

```





# Part 2: S&P 500 experiment

## Read data

```R
vix = read.csv('VIX_History.csv')
sp500 = read.csv('spy.csv')

len = length(vix)

vix_close = vix$CLOSE[len-500:len]

library(xts)
vix_close = xts(vix$CLOSE, order.by=as.Date(vix$DATE, "%m/%d/%Y"))
sp500_close = xts(sp500$Close, order.by=as.Date(sp500$Date, "%Y-%m-%d"))

data = merge(vix_close,sp500_close,all=FALSE)
data$ret = diff(log(data$sp500_close))*100
data$realvol = data$ret^2
data = na.omit(data)

data = window(data, start = "2012-01-01", end = "2021-12-31")

```

Plot VIX index data

```R
train_data = window(data, start = "2012-01-01", end = "2019-12-31")
test_data = window(data, start = "2020-01-01", end = "2021-12-31")

plot(index(data),as.numeric(data$vix_close),type='l',
     main='VIX index (2012-01-01 to 2021-12-31)',
     xlab="time",ylab='VIX index')
```



### MCMC SVAR1 model

Run model on estimating period

```R
sp500 = as.numeric(train_data$ret)
N = length(sp500)
data <- list(y=sp500,N=N)
inits <- function(){list(mu=0,phistar=0.9,sigma=0.1)}
param = c("mu","phi","std","sigma_y")

svar_mod = bugs(data=data,
                inits=inits,
                model.file="svar1.bug",
                parameters=param,
                n.chains=2,
                n.iter=5000,
                n.burnin=2000,
                n.thin=5,
                bugs.directory="D:/WinBUGS14/",
                debug=TRUE,
                DIC = FALSE,
                codaPkg=FALSE)

svar_summary = svar_mod$summary
round(svar_summary[1:5,],4)
sp500_sv_mcmc = svar_summary[4:nrow(svar_summary),1]

vix_train = as.numeric(train_data$vix_close)

dev.new()
par(mar=c(5, 4, 4, 6) + 0.1)
plot(index(train_data)[-1],vix_train[-1], type="l", yaxt="n",
     xlab="Time", ylab="", col="blue", 
     main='VIX index versus MCMC smoothed SV, estimating period (2012-01-01 to 2019-12-31)')
mtext("VIX index",side=2,col="blue",line=2) 
axis(2, at= pretty(range(vix_train[-1]),10),
     col.ticks="blue", col.axis="blue")
par(new=T)
plot(index(train_data)[-1],sp500_sv_mcmc, type="l", yaxt="n",
     xlab="Time", ylab="", col="red")
mtext("MCMC smoothed SV",side=4,col="red",line=4) 
axis(4, at = pretty(range(sp500_sv_mcmc),10),col.ticks="red", col.axis="red")
legend("topleft",legend=c("VIX index","MCMC smoothed SV"),
       text.col=c("blue","red"),pch=c(15,15),col=c("blue","red"))
```

Run particle filter FFBS on testing period

```R
sp500_test = as.numeric(test_data$ret)

ffbs_sp500_result = particle_FFBS(sigev=0.3638, 
                            phi=0.9214,
                            mu=-0.8694, 
                            data=sp500_test,
                            N=5000, 
                            nsmoo=5000)

ffbs_sp500_smoo = ffbs_sp500_result$ffbs.samples
sp500_sv_ffbs = exp(apply(ffbs_sp500_smoo, 2, mean)/2)

vix_test = as.numeric(test_data$vix_close)
plot(index(test_data),vix_test, type="l", yaxt="n",
     xlab="Time", ylab="", col="blue", 
     main='VIX index versus FFBS smoothed SV, testing period(2020-01-01 to 2021-12-31)')
mtext("VIX index",side=2,col="blue",line=2) 
axis(2, at= pretty(range(vix_test),10),
     col.ticks="blue", col.axis="blue")
par(new=T)
plot(index(test_data),sp500_sv_ffbs, type="l", yaxt="n",
     xlab="Time", ylab="", col="red")
mtext("FFBS smoothed SV",side=4,col="red",line=4) 
axis(4, at = pretty(range(sp500_sv_ffbs),10),col.ticks="red", col.axis="red")
legend("topright",legend=c("VIX index","FFBS smoothed SV"),
       text.col=c("blue","red"),pch=c(15,15),col=c("blue","red"))
```



# S&P 500 exercise extension

### SVAR1 via MCMC:

```r

datasvar <- list(y=sp500_test,N=length(sp500_test))
inits <- function(){list(mu=0,phistar=0.9,sigma=0.1)}
svar_test = bugs(data=datasvar,
                 inits=inits,
                 model.file="svar1.bug",
                 parameters=c("mu","phi","std","sigma_y"),
                 n.chains=2,
                 n.iter=10000,
                 n.burnin=5000,
                 n.thin=10,
                 bugs.directory="D:/WinBUGS14/",
                 debug=TRUE,
                 DIC = FALSE,
                 codaPkg=FALSE)

svar_test_summary = svar_test$summary
round(svar_test_summary[1:5,],)
sp500_svar = svar_test_summary[4:nrow(svar_test_summary),1]

```

### SVMIX via MCMC

Define model structure:

```R
svmix1 = function() 
{ # prior distribution
  lambda[1] ~ dnorm(0,0.001)
  theta ~ dnorm(0,0.001)I(0, )
  lambda[2] <- lambda[1] + theta
  P[1:2] ~ ddirch(alpha[])
  
  phistar ~ dbeta(10,1.5)
  sigma ~ dgamma(0.5,5)
  std <- sqrt(sigma)
  tau <- 1/sigma
  phi <- 2*phistar-1
  

  # initialize
  T[1] ~ dcat(P[])
  mu[1] <- lambda[T[1]]
  var0 <- tau/(1-pow(phi,2))
  
  h_mu0 ~ dnorm(mu[1],var0)
  h_mu[1] <- mu[1] + phi*(h_mu0-mu[1])
  h[1] ~ dnorm(h_mu[1], tau)
  
  for (t in 2:N)
  { 
    T[t] ~ dcat(P[])
    mu[t] <- lambda[T[t]]
    h_mu[t] <- mu[t] + phi*(h[t-1]-mu[t])
    h[t] ~ dnorm(h_mu[t],tau)
    sigma_y[t] <- exp(h[t]/2)
    tau_y[t] <- 1/exp(h[t])
    
    y[t] ~ dnorm(0,tau_y[t])
  }
}
write.model('svmix1.bug')

```

Run the model:

```R
svmix_mod = bugs(data=data,
                 inits=inits,
                 model.file="svmix1.bug",
                 parameters=c("phi","theta","P","std","lambda","sigma_y"),
                 n.chains=2,
                 n.iter=10000,
                 n.burnin=5000,
                 n.thin=10,
                 bugs.directory="D:/WinBUGS14/",
                 debug=TRUE,
                 DIC = FALSE,
                 codaPkg=FALSE)

svmix_summary = svmix_mod$summary
round(svmix_summary[1:9,],3)

sp500_svmix = svmix_summary[8:nrow(svmix_summary),1]
sv_mcmc = svar_summary[4:nrow(svar_summary),1]

plot(sp500_svmix, type="l", yaxt="n",
     xlab="Time", ylab="", col="blue", main='SVAR1 versus SVMIX (2020-01-01 to 2021-12-31)')
axis(2, at= pretty(range(sp500_svmix),10),
     col.ticks="black", col.axis="black")
points(sp500_svar, type="l", col="red")
legend("top",legend=c("SVMIX","SVAR1"),
       text.col=c("blue","red"),pch=c(15,15),col=c("blue","red"))
```

Plot mixed posterior mean distributions:

```r
# plot densities
plot(density(rowMeans(svar_test$sims.array[,,"mu"])),type='l',
ylab="",xlim=c(-5,5),main="estimated mean of h_t")
points(density(rowMeans(svmix_mod$sims.array[,,"lambda[1]"])),type='l',main='phi',
       ylab="",col='red')
points(density(rowMeans(svmix_mod$sims.array[,,"lambda[2]"])),type='l',
       main='phi',ylab="",col='blue')
legend("topleft",legend=c("mu","lambda1","lambda2"),
       text.col=c("black","blue","red"),
       pch=c(15,15,15),col=c("black","blue","red"))
```



### SVMS via MCMC

Markov switching stochastic volatility model

Define model structure:

```R
svmarkov1 = function()
{ # prior distribution
  lambda[1] ~ dnorm(0,0.001)
  theta ~ dnorm(0,0.001)I(0, )
  lambda[2] <- lambda[1] + theta 
  p ~ dbeta(3,2)
  q ~ dbeta(3,2)
  #P[1:2] ~ ddirch(alpha[])
  gamma0 ~ dbern(0.5)
  P[1] <- (1-p)*(1-gamma0)+q*(gamma0)
  
  phistar ~ dbeta(10,1.5)
  sigma ~ dgamma(0.5,5)
  tau <- 1/sigma
  phi <- 2*phistar-1
  std <- sqrt(sigma)
  
  # initialize
  #T[1] ~ dcat(P[])
  gamma[1] ~ dbern(P[1])
  T[1] <- gamma[1] + 1
  #proba[1] <- (1-p+q)/2
  
  mu[1] <- lambda[T[1]]
  var0 <- tau/(1-pow(phi,2))
  
  h_mu0 ~ dnorm(mu[1],var0)
  h_mu[1] <- mu[1] + phi*(h_mu0-mu[1])
  h[1] ~ dnorm(h_mu[1], tau)
  
  for (t in 2:N)
  { 
    P[t] <- (1-p)*(1-gamma[t-1])+q*(gamma[t-1])
    gamma[t] ~ dbern(P[t])
    T[t] <- gamma[t] + 1
    #proba[t] <- (1-p)*(1-proba[t-1]) + q*(proba[t-1])
    
    mu[t] <- lambda[T[t]]
    h_mu[t] <- mu[t] + phi*(h[t-1]-mu[t])
    h[t] ~ dnorm(h_mu[t],tau)
    
    sigma_y[t] <- exp(h[t]/2)
    tau_y[t] <- 1/sigma_y[t]
    
    y[t] ~ dnorm(0,tau_y[t])
  }
}
write.model('svmarkov1.bug')
```

Run model via MCMC:

```R
data <- list(y=sp500_test,N=length(sp500_test))
inits <- function(){list(phistar=rbeta(1,10,1.5),
                         sigma=0.5)}

svms_mod = bugs(data=data,
                 inits=inits,
                 model.file="svmarkov1.bug",
                 parameters=c("phi","p","q","std","lambda","theta"),
                 n.chains=2,
                 n.iter=10000,
                 n.burnin=5000,
                 n.thin=10,
                 bugs.directory="D:/WinBUGS14/",
                 debug=TRUE,
                 DIC = FALSE,
                 codaPkg=FALSE)

svms_summary = svms_mod$summary
round(svms_summary[1:9,],3)
sp500_svms = svms_summary[7:nrow(svms_summary),1]

plot(sp500_svmix, type="l", yaxt="n",
     xlab="Time", ylab="", col="blue", main='SVAR1 versus SVMIX (2020-01-01 to 2021-12-31)')
axis(2, at= pretty(range(sp500_svmix),10),
     col.ticks="black", col.axis="black")
points(sp500_svar, type="l", col="red")
points(sqrt(sp500_svms), type="l", col="darkgreen")
legend("top",legend=c("SVMIX","SVAR1","SVMS"),
       text.col=c("blue","red","darkgreen"),pch=c(15,15,15),
       col=c("blue","red","darkgreen"))

```

Traceplots:

```R

par(mfrow=c(4,2))
plot(rowMeans(svms_mod$sims.array[,,"phi"]),type='l',
     ylab="",main="phi")
plot(rowMeans(svms_mod$sims.array[,,"std"]),type='l',
     ylab="",main="sigma")
plot(rowMeans(svms_mod$sims.array[,,"lambda[1]"]),type='l',
     ylab="",main="lambda 1")
plot(rowMeans(svms_mod$sims.array[,,"lambda[2]"]),type='l',
     ylab="",main="lambda 2")
plot(rowMeans(svms_mod$sims.array[,,"theta"]),type='l',
     ylab="",main="theta")
plot(rowMeans(svms_mod$sims.array[,,"p"]),type='l',
     ylab="",main="p")
plot(rowMeans(svms_mod$sims.array[,,"q"]),type='l',
     ylab="",main="q")
```

