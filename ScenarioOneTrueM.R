library(zoo)
library(parallel)
library(harvestr)
library(tmle)
library(MASS)
library(clusterGeneration)
library(ggplot2)

## I'm commenting out the bootstrapping bits for now.
#source("TMLEfunctionsMScen1.R")
source("TMLEfunctionsMestScen1.R")

#detectCores()

expit<-function(p){
	exp(p)/(1+exp(p))
}

## Define seeds
#trials<-5
#seed.temp <- gather(trials,seed=123)
#Seed <- matrix(nrow=trials,ncol=6)
#for(i in 1:trials){
# Seed[i,] <- seed.temp[[i]][2:7]
#}

nsims <- 1000

#dataset with n individuals. 
set.seed(132)
n=100000

# QUESTION: I'd like to just make everything mvn for simpliciity, but this isn't exactly aligned with the case study

sigmaW1 <- as.matrix(read.csv(file="W1.csv", header=FALSE))
sigmaW2<-genPositiveDefMat(16, covMethod= "unifcorrmat", alphad=1, rangeVar=c(1,10))

w1vars<-mvrnorm(n=n, mu=c(20, -15,-10, 2.5,-1,0), Sigma=sigmaW1)
w2vars<-mvrnorm(n=n, mu=c(0,0,2,0,10 ,5,-5,10,-2,0, 0,0,0,0,0,0), Sigma=sigmaW2$Sigma)

#w1scvars<-sweep(w1vars, 2, c(.5, 1, 1, 1, 1, 1), "*")

#z1<-apply(w1scvars, 1, sum)
z1<-apply(w1vars, 1, sum)
z2<-apply(w2vars, 1, sum)
#I checked that these new covariates have similar overlap as previous simple covariates

#previous covariates
#z1<-rnorm(n)
#z2<-rnorm(n,2,1)

# make treatment variable
#P(tx ~ 1/3)
#previous
#prob.trt<-expit(-2.5 + (log(1.1)*z2) + log(1.005)*I(z2^2))
prob.trt<-expit(-2.5 + (log(1.1)*z2) )
t<-rbinom(n, 1, prob.trt)

# make selection into larger sample
# Mean P(selection)=.1
beta0 <- -1.5
beta1 <- log(1.09)
prob.sel <- expit(beta0+beta1*z1)
summary(prob.sel)

ipsvyselw<-1/prob.sel

meany0<- -3 + z2 
y0<-rnorm(n,meany0,2)
y1<-y0 + 3 + 2*z1 + rnorm(n,0,.5)

y<-ifelse(t==1, y1, y0)

# True average effect
ate <- mean(y1)-mean(y0)
print(ate)

#unadjusted/biased estimate of ate
mean(y[t==1]) - mean(y[t==0])

# Iterate drawing samples and estimating treatment effects
effects <- bootvar <- coverage.ate <- matrix(NA, ncol=7, nrow=nsims)
colnames(effects) <- colnames(bootvar) <- colnames(coverage.ate) <- c("Naive", "trueIPTW", "trueIPSW", "trueIPTSvyW",  "trueIPW", "trueTMLE2", "trueTMLE3")

# bounded y
yb1<-(y-min(y))/(max(y)-min(y))
yb2<-ifelse(yb1<0.001, 0.001, yb1)
yb<-ifelse(yb2>0.999, 0.999, yb2)

for (i in 1:nsims) {
	print(i)
	
	ss<-rbinom(n,1, prob.sel)

	pop <- data.frame(t, y0,y1, y, yb, ipsvyselw,  z1, z2, ss, prob.trt, prob.sel)
	colnames(pop)<-c("t","y0","y1", "y", "yb", "ipsvyselwt",  "z1", "z2", "insample", "prob.trt", "prob.svysel")
	pop$id<-index(pop)

	#look at characteristics of stabilized weights
	svysamp<-pop[ss==1,]
	svysamp$ssvywt<-mean(svysamp$prob.svysel)/svysamp$prob.svysel
	summary(mean(svysamp$prob.svysel)/svysamp$prob.svysel)
	sd(mean(svysamp$prob.svysel)/svysamp$prob.svysel)
	
	svysamp$stxwt[svysamp$t==1]<-mean(svysamp$prob.trt)/svysamp$prob.trt[svysamp$t==1]
	svysamp$stxwt[svysamp$t==0]<-mean(1-svysamp$prob.trt)/(1-svysamp$prob.trt[svysamp$t==0])
	summary(mean(svysamp$prob.trt)/svysamp$prob.trt[svysamp$t==1])
	sd(mean(svysamp$prob.trt)/svysamp$prob.trt[svysamp$t==1])
	summary(mean(1-svysamp$prob.trt)/(1-svysamp$prob.trt[svysamp$t==0]))
	sd(mean(1-svysamp$prob.trt)/(1-svysamp$prob.trt[svysamp$t==0]))

	# Selection model for selection into smaller sample 
	# Mean P(selection)=.25
	beta0 <- -2.4
	beta1 <- log(1.05)
	beta2 <- log(1.05)

	svysamp$prob.subsel <- expit(beta0+ beta1*svysamp$z1 + beta2*svysamp$z2 )
	svysamp$insubsample<-rbinom(nrow(svysamp), 1, svysamp$prob.subsel)
	svysamp$ssubselwt[svysamp$t==1]<-svysamp$insubsample[svysamp$t==1]*mean(svysamp$prob.subsel[svysamp$t==1])/svysamp$prob.subsel[svysamp$t==1]
	svysamp$ssubselwt[svysamp$t==0]<-svysamp$insubsample[svysamp$t==0]*mean(svysamp$prob.subsel[svysamp$t==0])/svysamp$prob.subsel[svysamp$t==0]

	svysub <- svysamp[svysamp$insubsample==1,]

	#look at distribution of stabilized weights
	summary(svysub$ssubselwt[svysub$t==1])
	sd(svysub$ssubselwt[svysub$t==1])
	summary(svysub$ssubselwt[svysub$t==0])
	sd(svysub$ssubselwt[svysub$t==0])

	summary(svysub$ssubselwt[svysub$t==1]*svysub$stxwt[svysub$t==1])
	sd(svysub$ssubselwt[svysub$t==1]*svysub$stxwt[svysub$t==1])
	summary(svysub$ssubselwt[svysub$t==0]*svysub$stxwt[svysub$t==0])
	sd(svysub$ssubselwt[svysub$t==0]*svysub$stxwt[svysub$t==0])

	summary(svysub$ssubselwt[svysub$t==1]*svysub$stxwt[svysub$t==1]*svysub$ssvywt[svysub$t==1])
	sd(svysub$ssubselwt[svysub$t==1]*svysub$stxwt[svysub$t==1]*svysub$ssvywt[svysub$t==1])
	summary(svysub$ssubselwt[svysub$t==0]*svysub$stxwt[svysub$t==0]*svysub$ssvywt[svysub$t==0])
	sd(svysub$ssubselwt[svysub$t==0]*svysub$stxwt[svysub$t==0]*svysub$ssvywt[svysub$t==0])

	sum(svysub$t)
	sum(I(svysub$t==0))
	#about 1000 more people in the simulation than in the example. 2721 v. 1608

summary(lm(y ~ t, data=pop))
summary(lm(y ~ t, data=svysamp))
summary(lm(y ~ t, data=svysub))

# First naive estimate
temp <- lm(y ~ t, data=svysub)
effects[i,"Naive"] <- summary(temp)$coef[2,1]	
sd <- summary(temp)$coef[2,2]
bootvar[i, "Naive"] <- sd^2
coverage.ate[i,"Naive"] <- ifelse(effects[i,"Naive"]-2*sd < ate & effects[i,"Naive"] +2*sd > ate, 1, 0)

###IPW
effects[i, "trueIPTW"] <-iptw(svysamp)
#res.iptw<-unlist(mclapply(1:trials, boot.iptw, mc.cores=6))
#bootvar[i, "trueIPTW"] <-var(res.iptw)
#cov<-quantile(res.iptw, probs=c(.025, .975))
#coverage.ate[i,"trueIPTW"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

effects[i, "trueIPSW"] <-ipsw(svysamp)
#res.ipsw<-unlist(mclapply(1:trials, boot.ipsw, mc.cores=6))
#bootvar[i, "trueIPSW"] <-var(res.ipsw)
#cov<-quantile(res.ipsw, probs=c(.025, .975))
#coverage.ate[i,"trueIPSW"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

effects[i, "trueIPTSvyW"] <-iptsvyw(svysamp)
#res.iptsvyw<-unlist(mclapply(1:trials, boot.iptsvyw, mc.cores=6))
#bootvar[i, "trueIPTSvyW"] <-var(res.iptsvyw)
#cov<-quantile(res.iptsvyw, probs=c(.025, .975))
#coverage.ate[i,"trueIPTSvyW"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)
	
effects[i, "trueIPW"] <-ipw(svysamp)
#res.ipw<-unlist(mclapply(1:trials, boot.ipw, mc.cores=6))
#bootvar[i, "trueIPW"] <-var(res.ipw)
#cov<-quantile(res.ipw, probs=c(.025, .975))
#coverage.ate[i,"trueIPW"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

###TMLE

## True model	
effects[i, "trueTMLE2"] <-joffe(svysamp)
#res.joffe<-unlist(mclapply(1:trials, boot.joffe, mc.cores=6))
#bootvar[i, "trueTMLE2"] <-var(res.joffe)
#cov<-quantile(res.joffe, probs=c(.025, .975))
#coverage.ate[i,"trueTMLE2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

effects[i, "trueTMLE3"] <-tmle3(svysamp)
#res.tmle3<-unlist(mclapply(1:trials, boot.tmle3, mc.cores=6))
#bootvar[i, "trueTMLE3"] <-var(res.tmle3)
#cov<-quantile(res.tmle3, probs=c(.025, .975))
#coverage.ate[i,"trueTMLE3"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

}

write.table(effects, file="EffectsS1MMTrue.csv", row.names=FALSE, col.names=TRUE, sep=",")
#write.table(bootvar, file="VarS1MMTrue.csv", row.names=FALSE, col.names=TRUE, sep=",")
#write.table(coverage.ate, file="CoverageATES1MMTrue.csv", row.names=FALSE, col.names=TRUE, sep=",")

# look at simulation results
effects<-read.csv("EffectsS1MMTrue.csv", header=TRUE)
long<-reshape(effects, varying=c("Naive", "trueIPTW", "trueIPSW", "trueIPTSvyW", "trueIPW", "trueTMLE2", "trueTMLE3"), new.row.names=(1:7000), times=c("Naive", "trueIPTW", "trueIPSW", "trueIPTSvyW", "trueIPW", "trueTMLE2", "trueTMLE3"), v.names="effects", direction="long", timevar="method")
uci<-function(x){
	uciest<-quantile(x, probs=.975)
	return(uciest)
}
lci<-function(x){
	lciest<-quantile(x, probs=.025)
	return(lciest)
}
class(long$method)
long$method<-factor(long$method)
levels(long$method)<-c("Naive",  "IPSW", "IPSvyW", "IPTW", "IPW", "DRWLS", "TMLE")
a<-summaryBy(effects~method, long, FUN=c(mean, lci, uci))
p<-ggplot(long, aes(y=effects, x=method))
p+geom_jitter(position = position_jitter(width = .1, height=.1), size=.5) + layer(data = a, mapping = aes(x = method, y = effects.mean), geom = "point", size = 2, color = "red") + layer(data = a, mapping = aes(x = method, y = effects.mean,  ymin = effects.lci, ymax = effects.uci), geom = "errorbar", size = 0.5, color = "red", width = 0.1) + geom_hline(yintercept=ate, color="Blue")