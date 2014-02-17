require(zoo)
require(parallel)

joffe <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ z2 , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, mean(sampled.data$pscore.svy)/sampled.data$pscore.svy, mean(1-sampled.data$pscore.svy)/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-mean(sampled.data$pscore.subsamp1[sampled.data$t==1])/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-mean(sampled.data$pscore.subsamp1[sampled.data$t==0])/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ssvywt, 0)
			
	model<-glm(y ~t + z2 + t:z1, data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$ssvywt
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$ssvywt)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
		}

joffeT1 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ z2, data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + poly(z2,2) + t:z1, data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
		}

joffeT2 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ z1, data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + poly(z2,2) + t:z1, data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
		}

joffeS1 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + poly(z2,2) + t:z1, data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
		}

joffeS2 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1 , data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + poly(z2,2) + t:z1, data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
		}

joffeO1 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + z2 + t:z1, data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
		}

joffeO2 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2 , data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + t:z1, data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
		}

joffeO3 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + poly(z2,2) + z1, data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
		}

tmle3<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ z2, data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-mean(sampled.data$pscore.subsamp1[sampled.data$t==1])/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-mean(sampled.data$pscore.subsamp1[sampled.data$t==0])/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ z1 + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + W2 + A:W1, g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$ipsvyselwt)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.0001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ssvywt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ssvywt)
	return(ate.tmle.fluc.svy)
}

tmle3T1<-function(svysamp){
	sampled.data<-svysamp
	sampled.data$pscore.svy <- predict(glm(t ~ z2, data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + poly(W2,2) + A:W1, g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$ipsvyselwt)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.0001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

tmle3T2<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ z1, data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + poly(W2,2) + A:W1, g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$ipsvyselwt)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.0001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

tmle3S1<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2), data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ z1 + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + poly(W2,2) + A:W1, g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$ipsvyselwt)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.0001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

tmle3S2<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2), data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ z1, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + poly(W2,2) + A:W1, g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$ipsvyselwt)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.0001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

tmle3O1<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2), data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + W2 + A:W1, g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$ipsvyselwt)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.0001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

tmle3O2<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2), data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + A:W1, g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$ipsvyselwt)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.0001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

tmle3O3<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2), data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + poly(W2,2) + W1, g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$ipsvyselwt)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.0001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

iptw <- function(svysamp){
	sampled.data<-svysamp

	#sampled.data<-svysamp
	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ z2 , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, mean(sampled.data$pscore.svy)/sampled.data$pscore.svy, mean(1-sampled.data$pscore.svy)/(1-sampled.data$pscore.svy))
	
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$t.wt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 

ipsw <- function(svysamp){
	sampled.data<-svysamp

	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-mean(sampled.data$pscore.subsamp1[sampled.data$t==1])/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~z1 + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-mean(sampled.data$pscore.subsamp1[sampled.data$t==0])/sampled.data$pscore.subsamp1[sampled.data$t==0]
	#combine weights
	sampled.data$subsvywt<-ifelse(sampled.data$insubsample==1, sampled.data$subsamp.wt*sampled.data$ssvywt, 0)
		
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$subsvywt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 

iptsvyw <- function(svysamp){
	sampled.data<-svysamp

	#sampled.data<-svysamp
	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ z2 , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, mean(sampled.data$pscore.svy)/sampled.data$pscore.svy, mean(1-sampled.data$pscore.svy)/(1-sampled.data$pscore.svy))

	#combine weights
	sampled.data$txsvywt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$ssvywt, 0)
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$txsvywt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    }     

ipw <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ z2 , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, mean(sampled.data$pscore.svy)/sampled.data$pscore.svy, mean(1-sampled.data$pscore.svy)/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-mean(sampled.data$pscore.subsamp1[sampled.data$t==1])/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-mean(sampled.data$pscore.subsamp1[sampled.data$t==0])/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ssvywt, 0)
		
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 

ipwT1 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ z2 , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
		
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 

ipwT2 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ z1 , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
		
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 

ipwS1 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1+ z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
		
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 

ipwS2 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
		
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 
