    #####################################
##### Kernel machine regression in subgroup analysis by Cho, Zhan and Ghosh, 2021, Submitted to Biometrics
#####################################

### The main R function fits least squares kernel machine regression E(Y|X)=beta+h(x), 
### using the algorithm proposed in Liu, Lin and Ghosh, Biometrics, 2007, 1079-1088


KMR = function(y, x, xnew){
## y: response nx1
## x: covariate design matrix pxn (each column is a sample)
## xnew: new observations pxn* for prediction (if NULL then no predictrion is returned)
## Output: Estimates of beta0, h(x), sig, rho, lam and predicted values (optional)  

n = length(y)
xx=matrix(rep(1,n), nrow=1) ## Intercept
dMat = sqare.diff.mat(x)
	
########### Step 0: Assigning initial values to parameters
bet.hat = matrix(.2, ncol=1, nrow=nrow(xx))
sig.hat = 1
rho.hat = 1
lam.hat = 0.01
Kmat = kernel.Gauss(dMat, rho.hat)

conv.crit = 0.00001
criter = 999
bet.old = bet.hat
lam.old = lam.hat
sig.old = sig.hat
rho.old = rho.hat
niter = 0
while( (criter > conv.crit) & (niter < 50) ){
	
########### STEP 1 :fix  sigma, lambda and rho to estimate beta and h().
Kmat = kernel.Gauss(dMat, rho.hat)
V = sig.hat* diag(1, nrow=n, ncol=n) + lam.hat*Kmat
iV = solve(V)
bet.hat = solve(xx%*%iV%*%t(xx)) %*% xx%*%iV%*%y
h.hat = lam.hat*Kmat%*%iV%*% (y - t(xx)%*%bet.hat)	

	
########### STEP 2: fix beta and h() to estimate sigma, lambda and rho
Px = solve(xx%*%iV%*%t(xx))
H.mat = (t(xx)%*%Px%*%xx)%*%iV + 
lam.hat*Kmat%*%iV%*%( diag(1, nrow=n, ncol=n) - t(xx)%*%Px%*%xx%*%iV )

sig.hat = sum( (y - t(xx)%*%bet.hat - h.hat)^2 )/ (n - Trace(H.mat))
if((n - Trace(H.mat)) < 0)	{	
sig.hat = (sum((y - t(xx)%*%bet.hat - h.hat)^2)/ (n))	
}

lam.hat = optimize(loglik.lam, interval=c(0.0,60), 
							Y=(y - t(xx)%*%bet.hat), 
							Kmat=Kmat,
							xx = xx,
							sigmasq = sig.hat)$minimum 

							
rho.hat = optimize(loglik.rho, interval=c(0.0,40), 
							Y=(y - t(xx)%*%bet.hat), 
							dMat=dMat,
							xx = xx,
							lambda = lam.hat,
							sigmasq = sig.hat)$minimum
##
criter = sum(abs(bet.hat - bet.old)) + 
					abs(lam.hat - lam.old) +
					abs(sig.hat - sig.old) +
					abs(rho.hat - rho.old)

lik.val = loglik.lam(lambda = lam.hat, Y=(y - t(xx)%*%bet.hat), Kmat=Kmat, xx=xx, sigmasq=sig.hat)			
					
bet.old = bet.hat
lam.old = lam.hat
sig.old = sig.hat
rho.old = rho.hat
niter = niter + 1
}


## sanity check
if(niter >= 100){
	print("Did not converge")
	stop("Did not converge")
}

Kmat = kernel.Gauss(dMat, rho.hat)
V = sig.hat* diag(1, nrow=n, ncol=n) + lam.hat*Kmat
iV = solve(V)
bet.hat = solve(xx%*%iV%*%t(xx)) %*% xx%*%iV%*%y
h.hat = lam.hat*Kmat%*%iV%*% (y - t(xx)%*%bet.hat)	

if (is.null(xnew)){
return(list(bet=bet.hat, h.hat = h.hat, sig=sig.hat, lam=lam.hat, rho=rho.hat))				
} else {
Kvec=matrix(NA,nrow=ncol(xnew),ncol=n)
for(jj in 1:ncol(xnew)){ xtemp=xnew[,jj]
for(ii in 1:n){ Kvec[jj,ii]= exp(-(xtemp-x[,ii])%*%(xtemp-x[,ii])/rho.hat)
}}
pred.hat= Kvec%*%iV%*% (y - t(xx)%*%bet.hat)
return(list(bet=bet.hat, h.hat = h.hat, sig=sig.hat, lam=lam.hat, rho=rho.hat, 
pred=pred.hat))
}
}


###########################
## Auxiliary functions
###########################
loglik.lam = function(lambda, Y, Kmat, xx, sigmasq){

	n = length(Y)
	V = sigmasq* diag(1, nrow=n, ncol=n) + lambda*Kmat
	iV = solve(V)

	out = -log(det(V))/2 - log( det(xx%*%iV%*%t(xx)) )/2 - t(Y)%*%iV%*%Y/2
	return(-out/n)
}


loglik.sig = function(sigmasq, Y, Kmat, xx, lambda){

	n = length(Y)
	V = sigmasq* diag(1, nrow=n, ncol=n) + lambda*Kmat
	iV = solve(V)

	out = -log(det(V))/2 - log( det(xx%*%iV%*%t(xx)) )/2 - t(Y)%*%iV%*%Y/2
	return(-out/n)
}

loglik.rho = function(rho, Y, dMat, xx, lambda, sigmasq){

	n = length(Y)
	Kmat = kernel.Gauss(dMat, rho)
	V = sigmasq* diag(1, nrow=n, ncol=n) + lambda*Kmat
	iV = solve(V)

	out = -log(det(V))/2 - log( det(xx%*%iV%*%t(xx)) )/2 - t(Y)%*%iV%*%Y/2
	return(-out/n)
}


Trace = function(M){sum(diag(M))}

## The numerator and denominator of Gaussian kernel are separate here so that 
## it is easier to estimate the denominator as a variance component parameter 
## in the main R function named KMR

sqare.diff.mat = function(x){  # Each column is a sample
Kmat = matrix(NA, nrow=ncol(x), ncol=ncol(x))
for (ii in 1:ncol(x)){Kmat[ii,] = diag(t(x-x[,ii])%*%(x-x[,ii]))}
return(Kmat)
}

kernel.Gauss = function(Kmat, rho){
return(exp(-Kmat/rho))
}


## Examples:
set.seed(2021)
n<-200
p<-10
zz<-matrix(NA,p,n)
for(i in 1:p){ zz[i,]<-runif(n,0,1) }
u<-rnorm(n,0,1)
x1<-3*cos(zz[1,])+2*u
e<-rnorm(n,0,1)
h<-function(z){2*cos(z[1])-3*z[2]^2+2*exp(-z[3])*z[4]-1.6*sin(z[5])*cos(z[3])
+4*z[1]*z[5]}
hz<-rep(0,n)
for(i in 1:n){hz[i]<-h(zz[,i])} 
x0<-rep(1,n)
yy=x0+hz+e 

## Estimation: 
fit1<-KMR(yy,zz,xnew=NULL)
fit1
RKHSerr=fit1$h.hat + as.numeric(fit1$bet)-hz-1
fit.lm=lm(yy~t(zz))
c(sum(RKHSerr^2)/n,sum(fit.lm$res^2)/n)  

## Prediction:
xnew1 <-rnorm(p,0,1)
xnew2 <-rnorm(p,0,1)
xnew=cbind(xnew1,xnew2)
fit2<-KMR(yy,zz,xnew)
fit2