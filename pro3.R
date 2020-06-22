library("lme4")

nlist = list()
nlist[[1]] = c(15,1)
nlist[[2]] = c(6,1,2)
nlist[[3]] = c(6,6)
nlist[[4]] = c(7,2,3,NA,2)
nlist[[5]] = c(16,9,3,3,1)
nlist[[6]] = c(57,38,17,2,2)
nlist[[7]] = c(119,81,45,6,1,NA,NA,1)
nlist[[8]] = c(173,118,57,16,3,NA,NA,NA,1)
nlist[[9]] = c(136,103,50,13,6,1,1)
nlist[[10]] = c(54,51,32,5,1,NA,NA,NA,NA,1)
nlist[[11]] = c(13,15,12,3,1)
nlist[[12]] = c(NA,4,3,1)
nlist[[13]] = c(NA,NA,1,NA,NA,NA,NA,1)

myobs = vector()
count = 0
for (k in 1:13) {
  for (l in 1:length(nlist[[k]])) {
    if(!is.na(nlist[[k]][l]))
      for (m in 1:nlist[[k]][l]) {
        count = count + 1
        myobs = rbind(myobs,cbind(c(rep(0,k-l+1),rep(1,l-1)),count))
      }
  }
}
myobs = data.frame(y=myobs[,1],litter=factor(myobs[,2]))

model1 = glmer(y~1+(1|litter),data = myobs,family = "binomial")
summary(model1)
initial_mu = fixef(model1)
initial_sigma = sqrt(unlist(VarCorr(model1)))
####################################################################
I = nlevels(myobs$litter)
y_i = split(myobs$y,myobs$litter)
sum_y_i = sapply(y_i, sum)
n_i = sapply(y_i, length)
m = 10000
burn_out = 500

##################################################3
y_likelihood = function(mu,alpha){
  return(sum_y_i%*%(mu+alpha)-n_i%*%log(1+exp(mu+alpha)))
}

alpha_likelihood = function(alpha,sigma){
  return(sum(dnorm(alpha,mean = 0,sigma,log = T)))
}

#compute the acceptance rate
r = function(alpha0,alpha1,mu){
  p = sum_y_i*(alpha1-alpha0) - n_i*log((1+exp(mu+alpha1))/(1+exp(mu+alpha0)))
  return(pmin(exp(p),1))
}

#sample from the target distribution: f(alpha|y,mu)
metropoli = function(mu,sigma){
  iterations = m + burn_out
  initial_alpha = rnorm(I,0,sigma)
  alpha = matrix(0,I,iterations)
  alpha[,1] = initial_alpha
  for (i in 2:iterations) {
    alpha_candidate = rnorm(I,0,sigma)
    u = runif(I)
    r = r(alpha[,i-1],alpha_candidate,mu)
    jump = which(u<r)
    length(jump)
    alpha[jump,i] = alpha_candidate[jump]
    alpha[-jump,i] = alpha[-jump,i-1]
  }
  burn_out_alpha = alpha[,-(1:burn_out)]
  return(burn_out_alpha)
}

set.seed(1)
z0 = metropoli(-2.99,0.69)
write.table(z0,"z0.txt")
plot(z0[1,500:10000],type="l") 

########################optimize
optimize_mu = function(alpha,mu){
  return(-mean(apply(alpha, 2, y_likelihood,mu=mu)))
}

optimize_sigma = function(alpha,sigma){
  return(-alpha_likelihood(alpha,sigma))
}

#optim(-2.29,optimize_mu,alpha=z0,method="BFGS")
#######
MCEM = function(mu_initial,sigma_initial,times){
  mu = rep(0,times)
  mu[1] = mu_initial
  sigma = rep(0,times)
  sigma[1] = sigma_initial
  for (i in 2:times) {
    alpha = metropoli(mu[i-1],sigma[i-1])
    mu[i] = optim(mu[i-1],optimize_mu,alpha=alpha,method="L-BFGS-B",lower=-3,upper=-1)$par
    sigma[i] = optim(sigma[i-1],optimize_sigma,alpha=alpha,method="L-BFGS-B",lower=0.3,upper=0.8)$par
    dmu = abs(mu[i]-mu[i-1])
    dsigma = abs(sigma[i]-sigma[i-1])
    cat("mu is:",mu[i],"df is:",dmu,"Iteration is:",i,"\n")
    cat("sigma is:",sigma[i],"df is:",dsigma,"\n")
  }
  return(list(mu,sigma))  
}

set.seed(1)
result = MCEM(0,0.5,300)
names(result)=c("mu","sigma")
write.table(result,"result.txt")
par(mfrow=c(1,2))
plot(c(1:300),result$mu,type="l",main = "1 to 300 iterations of u",xlab="1 to 300",ylab = "value of u",cex.main=0.8)
plot(c(30:300),result$mu[30:300],type = "l",main = "30 to 300 iterations of u",xlab="30 to 300",ylab = "value of u",cex.main=0.8)

plot(c(1:300),result$sigma,type="l",main = "1 to 300 iterations of sigma",xlab="1 to 300",ylab = "value of sigma",cex.main=0.8)
plot(c(30:300),result$sigma[30:300],type = "l",main = "30 to 300 iterations of sigma",xlab="30 to 300",ylab = "value of sigma",cex.main=0.8)

mu = result$mu[300]
sigma = result$sigma[300]
alpha = z0
s1 = sum(sum_y_i)-sum(n_i%*% ((1/(1+exp(mu+alpha)))*exp(mu+alpha)) )
s2 = -10533/sigma+sum(alpha^2/sigma^3)
s = as.matrix(c(s1,s2),ncol=1)
B1 = sum(n_i %*% (exp(mu+alpha)/(1+exp(mu+alpha))^2))
B4 = 10533/sigma^2-3/(sigma^4)*sum(alpha^2)
B = matrix(c(B1,0,0,B4),ncol=2)
I = -B + s %*% t(s)
(solve(I))^(1/2)*100


