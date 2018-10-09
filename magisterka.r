#PRACA MAGISTERSKA
##################################################
###########wczytanie danych i bibliotek###########
##################################################
setwd("d:/dokumenty/Rfiles/do magisterki")
#setwd("J:/04022011/Rfiles/do magisterki")
rm(list=ls())
dane <- read.csv("bst.csv", header=TRUE, sep=";", dec=".", row.names="data")
dane
library(tseries)
library(fGarch)
library(urca)
library(FinTS)
library(vars)
library(systemfit)

#b-banki, s-spoz, t-tele


##################################################
################opis szeregów#####################
##################################################

attach(dane)
ls(2)
#z_banki z_spoz z_tele
b<-z_banki
s<-z_spoz
t<-z_tele
detach()
#statystyki (zrobiæ w excelu)

st<-function(a){
min<-min(a)
max<-max(a)
sr<-mean(a)
od<-sqrt(var(a))
skew<-skewness(a)
kurt<-kurtosis(a)
s<-c(min,max,sr,od,skew,kurt)
s
}
st(b)
st(s)
st(t)

par(mfrow=c(3,1), cex=.7)
plot(b, main="Sektor bankowy", type="l")
plot(s, main="Sektor spo¿ywczy", type="l")
plot(t, main="Sektor telekomunikacyjny", type="l")

par(mfrow=c(3,1), cex=.7)
acf(b, main="Sektor bankowy", ylim=c(-.2,.2), lwd=3)
acf(s, main="Sektor spo¿ywczy",ylim=c(-.2,.2), lwd=3)
acf(t, main="Sektor telekomunikacyjny",ylim=c(-.2,.2), lwd=3)

Box.test(t, lag=1)

####stacjonarnoœæ

summary(ur.df(t,type="none",lags=2,selectlags="AIC"))
summary(ur.kpss(t,type="mu", lags="long", use.lag=NULL))

#postaæ sk³adnika losowego
#autokorelacja - test
###############################

###############################
#b_ts<-ts(b)
#t_ts<-ts(t)
#s_ts<-ts(s)
#b_stat<-FinTS.stats(b_ts)
#t_stat<-FinTS.stats(t_ts)
#s_stat<-FinTS.stats(s_ts)
#t(b_stat)

##################################################
###########wybór reprezentacji ARMA-Garch#########
##################################################

###wybór ARMA###
ARMA<- function(data, pmax=4, qmax=4){
	k <- matrix(NA, pmax+1, qmax+1)
	l<-k
	for (p in 0:pmax){
		for (q in 0:qmax){
	a<-arima(data,c(p,0,q))
	k[p+1,q+1]<-((-2)*a$log+(p+q+1)*2)/length(data)
	l[p+1,q+1]<-((-2)*a$log+(p+q+1)*log(length(data)))/length(data)
	}
		}

	#print("Kryteria Akaike'a dla poszczególnych opóŸnieñ")
	#rownames(k) <- paste('ar',0:pmax, sep="")
	#colnames(k) <- paste('ma',0:qmax, sep="")
	#return(k)

	print("", quote=FALSE)
	print("Kryterium BIC")
	rownames(l) <- paste('ar',0:pmax, sep="")
	colnames(l) <- paste('ma',0:qmax, sep="")
	return(l)
print("wybierz najlepsze, znaczy najmniejsze")
}
ARMA(b)
ARMA(s)
ARMA(t)

sink("NUL")
arma_b<-arima(b,c(0,0,1))
arma_s<-arima(s,c(1,0,0))
arma_t<-arima(t,c(0,0,1))
sink()
arma_b
arma_s
arma_t

summary(tsdiag(arma_b))
Box.test(arma_b$r, lag=1)
Box.test(arma_s$r)
Box.test(arma_t$r)


#b_r<-arma_b$r/arma_b$sigma
#s_r<-arma_s$r/arma_s$sigma
#t_r<-arma_t$r/arma_t$sigma
b_r<-arma_b$r
s_r<-arma_s$r
t_r<-arma_t$r
par(mfrow=c(3,1), cex=.7)
plot(b_r^2, main="Sektor bankowy",ylab="kwadraty reszt", type="l")
plot(s_r^2, main="Sektor spo¿ywczy",ylab="kwadraty reszt", type="l")
plot(t_r^2, main="Sektor telekomunikacyjny",ylab="kwadraty reszt", type="l")

par(mfrow=c(3,1), cex=.7)
acf(b_r, main="Sektor bankowy", ylim=c(-.2,.2), lwd=3)
acf(s_r, main="Sektor spo¿ywczy",ylim=c(-.2,.2), lwd=3)
acf(t_r, main="Sektor telekomunikacyjny",ylim=c(-.2,.2), lwd=3)



b_arma_aic<-c(0,1)
b_arma_bic<-c(0,1)
s_arma_aic<-c(1,1)
s_arma_bic<-c(1,0)
t_arma_aic<-c(0,1)
t_arma_bic<-c(0,1)

####opis reszt#### i sprawdzenie za³o¿eñ co do sk³adnika losowego ####

#########################################################################
####LM-test#### œci¹gniêty###############################################
#########################################################################
archlmtest <- function (x, lags, demean = FALSE) 
{
  x <- as.vector(x)
  if(demean) x <- scale(x, center = TRUE, scale = FALSE)
    lags <- lags + 1
    mat <- embed(x^2, lags)
    arch.lm <- summary(lm(mat[, 1] ~ mat[, -1]))
    STATISTIC <- arch.lm$r.squared * length(resid(arch.lm))
    names(STATISTIC) <- "Chi-squared"
    PARAMETER <- lags - 1
    names(PARAMETER) <- "df"
    PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
    METHOD <- "ARCH LM-test"
    result <- list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD, data.name =
deparse(substitute(x)))
    class(result) <- "htest"
    return(result)
}

archlmtest(b_r,5)
archlmtest(s_r,5)
archlmtest(t_r,5)

###wybór reprezentacji GARCH###
###nie dzia³a dok³adnie najwy¿ej ten element zrobiê rêcznie
GARCH<-function(dane, ar=0, ma=0, pmax=4, qmax=4, dist="norm"){
	k <- matrix(NA, pmax, qmax+1)
	l<-k
	for (p in 1:pmax){
		for (q in 0:qmax){
	sink("NUL")
	model <- as.formula(paste('~arma(',ar,',',ma,')+garch(',p,',',q,')', sep="" ))
	a <- garchFit(formula = model, data = dane, include.mean=TRUE, cond.dist=dist)
	k[p,q+1]<-a@fit$ics[1]
	l[p,q+1]<-a@fit$ics[2]
	sink()
		}
	}
	#print("Kryteria Akaike'a dla poszczególnych opóŸnieñ")
	#rownames(k) <- paste('arch',1:pmax, sep="")
	#colnames(k) <- paste('garch',0:qmax, sep="")
	#return(k)
	#print("", quote=FALSE)
	print("Kryterium BIC")
	rownames(l) <- paste('arch',1:pmax, sep="")
	colnames(l) <- paste('garch',0:qmax, sep="")
	return(l)
	print("wybierz najlepsze, znaczy najmniejsze")
}
GARCH(b,ar=0,ma=1,pmax=3,qmax=3)
GARCH(s,ar=1,ma=0,pmax=3,qmax=3)
GARCH(t,ar=0,ma=1,pmax=3,qmax=3)

##################################################
###########szacowanie Garchy (znane postaci)######
##################################################

attach(dane)
sink("NUL") 
b_garch <- garchFit(formula = ~arma(0,1)+garch(1,2), data = z_banki, include.mean=TRUE, cond.dist="norm")
s_garch <- garchFit(formula = ~arma(1,0)+garch(1,1), data = z_spoz, include.mean=TRUE, cond.dist="norm")
t_garch <- garchFit(formula = ~arma(0,1)+garch(1,1), data = z_tele, include.mean=TRUE, cond.dist="norm")
sink() 
summary(b_garch)
summary(s_garch)
summary(t_garch)
detach()
b_sig <- b_garch@sigma.t
s_sig <- s_garch@sigma.t
t_sig <- t_garch@sigma.t
# reszty
b_u   <- b_garch@residuals
s_u   <- s_garch@residuals
t_u   <- t_garch@residuals
# wystandardyzowane reszty
b_e   <- b_u/b_sig
s_e   <- s_u/s_sig
t_e   <- t_u/t_sig

a<-1.95
par(mfrow=c(3,1), cex=.7)
plot(b_u, main="Sektor bankowy", type="l", ylab="zwroty z 95% przedzia³ami ufnoœci")
lines(a*b_sig, col="red")
lines(-a*b_sig,col="blue")
plot(s_u, main="Sektor spo¿ywczy", type="l", ylab="zwroty z 95% przedzia³ami ufnoœci")
lines(a*s_sig, col="red")
lines(-a*s_sig,col="blue")
plot(t_u, main="Sektor telekomunikacyjny", type="l", ylab="zwroty z 95% przedzia³ami ufnoœci")
lines(a*t_sig, col="red")
lines(-a*t_sig,col="blue")



x<-runif(3448,0,1)
hist(b_e,freq=FALSE,breaks=15, main="Rozk³ad reszt modelu dla sektora bankowego", xlab=NULL)
curve(dnorm(x), col = 2, lty = 2, lwd = 3, add = TRUE)
hist(s_e,freq=FALSE,breaks=19,main="Rozk³ad reszt modelu dla sektora spo¿ywczego", xlab=NULL)
curve(dnorm(x), col = 2, lty = 2, lwd = 3, add = TRUE)
hist(t_e,freq=FALSE,breaks=15,main="Rozk³ad reszt modelu sektora telekomunikacyjnego", xlab=NULL)
curve(dnorm(x), col = 2, lty = 2, lwd = 3, add = TRUE)

qqnorm(b_e,main="Sektor bankowy")
qqline(b_e)
qqnorm(s_e,main="Sektor spo¿ywczy")
qqline(s_e)
qqnorm(t_e,main="Sektor telekomunikacyjny")
qqline(t_e)

acf(s_e, main="Sektor spo¿ywczy - reszty modelu",ylim=c(-.1,.1), lwd=3)
acf(t_e^2, main="Sektor telekomunikacyjny - wspó³czynniki funkcji autokorelacji kwadratów reszt modelu", ylim=c(-.1,.1),lwd=3)


#############################################################
######tworzenie macierzy wystandaryzowanych reszt Dt#########
#############################################################
#macierz Dt jest macierz¹ kwadratow¹ tak¹ ¿e na przek¹tnych s¹ odpowiednie odchylenia standardowe

Dt<-rep(0, times=3*3*length(b_e))
dim(Dt)<-c(3,3,3448)
Dt[1,1,]<-b_sig
Dt[2,2,]<-s_sig
Dt[3,3,]<-t_sig

############################################################
#########CCC################################################
############################################################
## R - macierz bezwarunkowych korelacji miêdzy wystandaryzowanymi resztami

R<-rep(1,9)
dim(R)<-c(3,3)
R[1,2]<-cor(b_e,s_e)->R[2,1]
R[1,3]<-cor(b_e,t_e)->R[3,1]
R[2,3]<-cor(s_e,t_e)->R[3,2]
R

CCC <-function(Dt,R){
for (i in 1:dim(Dt)[3]){
H<-Dt
H[,,i]<-Dt[,,i]%*%R%*%Dt[,,i]
}
H
}
Ht_CCC<-CCC(Dt,R)

##########################################################
############wyg³adzanie wyk³adnicze uproszczone###########
##########################################################
ewma<- function(lambda=0.06 ,u1=b_e,u2=s_e,u3=t_e){

#a-wektor parametrów
#Qt bêdzie wymiaru 3x3x3447
#dQt=odw(sqrt(dQt))
#Rt 3x3x3447
#U jest wektorem 3x3447,
U<-rep(NA, 3*(length(u1)-1))
dim(U)<-c(3,length(u1)-1)
U[1,]<-u1[1:(length(u1)-1)]
U[2,]<-u2[1:(length(u2)-1)]
U[3,]<-u3[1:(length(u3)-1)]
Qt<-rep(NA,3*3*(length(u1)-1))
dim(Qt)<-c(3,3,(length(u1)-1))

dQt<-rep(0,3*3*(length(u1)-1))
dim(dQt)<-c(3,3,(length(u1)-1))

R<-rep(1,9)
dim(R)<-c(3,3,1)
R[1,2,1]<-cor(u1,u2)->R[2,1,1]
R[1,3,1]<-cor(u1,u3)->R[3,1,1]
R[2,3,1]<-cor(u2,u3)->R[3,2,1]

Rt<-rep(0, times=3*3*(length(u1)-1))
dim(Rt)<-c(3,3,(length(u1)-1))
Rt[1,1,]<-1
Rt[2,2,]<-1
Rt[3,3,]<-1
Rt[,,1]<-R[,,1]
Qt[,,1]<-R[,,1]
for (j in 2:dim(U)[2]){
	Qt[,,j]<- lambda*(U[,j-1]%*%t(U[,j-1])) + (1-lambda)*(Qt[,,j-1])
	diag(dQt[,,j])<-sqrt(diag(Qt[,,j]))
	Rt[,,j]<-solve(dQt[,,j])%*%Qt[,,j]%*%solve(dQt[,,j])
	}
Rt
}

Rt_06<-ewma()
Rt_01<-ewma(.01)
plot(Rt_01[2,3,],type="l")

##########################################################################
############wyg³adzanie wyk³adnicze uproszczone - optymalizacja###########
##########################################################################
########spróbowaæ fapply i pokrewne
########optimize()

ewma_opt<- function(lambda=0.06 ,u1=b_e,u2=s_e,u3=t_e){

#u3<-t_e
#a-wektor parametrów
#Qt bêdzie wymiaru 3x3x3447
#dQt=odw(sqrt(dQt))
#Rt 3x3x3447
#U jest wektorem 3x3447,
U<-rep(NA, 3*(length(u1)-1))
dim(U)<-c(3,length(u1)-1)
U[1,]<-u1[1:(length(u1)-1)]
U[2,]<-u2[1:(length(u2)-1)]
U[3,]<-u3[1:(length(u3)-1)]
Qt<-rep(NA,3*3*(length(u1)-1))
dim(Qt)<-c(3,3,(length(u1)-1))

dQt<-rep(0,3*3*(length(u1)-1))
dim(dQt)<-c(3,3,(length(u1)-1))

R<-rep(1,9)
dim(R)<-c(3,3,1)
R[1,2,1]<-cor(u1,u2)->R[2,1,1]
R[1,3,1]<-cor(u1,u3)->R[3,1,1]
R[2,3,1]<-cor(u2,u3)->R[3,2,1]

Rt<-rep(0, times=3*3*(length(u1)-1))
dim(Rt)<-c(3,3,(length(u1)-1))
Rt[1,1,]<-1
Rt[2,2,]<-1
Rt[3,3,]<-1
Rt[,,1]<-R[,,1]
Qt[,,1]<-R[,,1]

Lc<-0
for (j in 2:dim(U)[2]){
	Qt[,,j]<- (1-lambda)*(U[,j-1]%*%t(U[,j-1])) + lambda*(Qt[,,j-1])
	diag(dQt[,,j])<-sqrt(diag(Qt[,,j]))
	Rt[,,j]<-solve(dQt[,,j])%*%Qt[,,j]%*%solve(dQt[,,j])
	Lc<- Lc+(log(det(Rt[,,j]))+(t(U[,j])%*%solve(Rt[,,j])%*%U[,j]))*(-0.5)
	}
Lc
}

opima<-optimize(f=ewma_opt, interval=c(0,1))
#s<-constrOptim(.06,ewma_opt,NULL,ui=rbind(-1,1),ci=c(-.1,0.000001))
opt<-optima$minimum
opt

Rt_ewma_opt<-ewma(lambda=opt)
plot(Rt_ewma_opt[1,3,], type="l")
#Rt_01<-ewma(lambda=.01)
Rt_06<-ewma()
par(mfrow=c(3,1), cex=.7)
plot(Rt_06[1,3,], type="l", ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i telekomunikacyjnym")
plot(Rt_06[1,2,], type="l", ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i spo¿ywczym")
plot(Rt_06[2,3,], type="l", ylab="korelacje warunkowe", main="korelacje miêdzy sektorem spo¿ywczym i telekomunikacyjnym")

par(mfrow=c(3,1), cex=.7)
plot(Rt_ewma_opt[1,3,], type="l", ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i telekomunikacyjnym")
plot(Rt_ewma_opt[1,2,], type="l", ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i spo¿ywczym")
plot(Rt_ewma_opt[2,3,], type="l", ylab="korelacje warunkowe", main="korelacje miêdzy sektorem spo¿ywczym i telekomunikacyjnym")

par(mfrow=c(3,1), cex=.7)
plot(Rt_DCC[1,3,], type="l", ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i telekomunikacyjnym")
plot(Rt_DCC[1,2,], type="l", ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i spo¿ywczym")
plot(Rt_DCC[2,3,], type="l", ylab="korelacje warunkowe", main="korelacje miêdzy sektorem spo¿ywczym i telekomunikacyjnym")


##########################################################
#####korelacje z wyg³adzania wyk³adniczego################
##########################################################
###Rt-korelacje z wyg³adzania wyk³adniczego
###EWMA-funkcja wyg³adzaj¹ca 2 szeregi
#a,b-szeregi
#library(forecast) - zawiera wyg³adzanie wyk³adnicze ale nie pasi mi

EWMA<-function(a,b,lambda=0.94){
#tworzê pusta wektory
ro<-rep(0,length(a))
s1<-rep(0,length(a)-1)
s2<-rep(0,length(a)-1)
s3<-rep(0,length(a)-1)
for (i in (2:length(a))){
for (j in (1:i-1)){
#liczê sk³adowe wektory
s1[j]<- (a[j]*b[j]*(lambda^(i-1-j)))
s2[j]<- (a[j]^2*(lambda^(i-1-j)))
s3[j]<- (b[j]^2*(lambda^(i-1-j)))
}
#liczê poszczególne wyrazy
ro[i]<-sum(s1)/sqrt(sum(s2)*sum(s3))
}
ro
}

Rt_b_t_ewma<-EWMA(b_e,t_e,lambda=0.99)
Rt_s_t_ewma<-EWMA(s_e,t_e,lambda=0.99)
Rt_b_s_ewma<-EWMA(b_e,s_e,lambda=0.99)
Rt_b_t_ewma150<-Rt_b_t_ewma[150:3447]
Rt_s_t_ewma150<-Rt_s_t_ewma[150:3447]
Rt_b_s_ewma150<-Rt_b_s_ewma[150:3447]

par(mfrow=c(3,1), cex=.7)
plot(Rt_b_t_ewma150, type="l", ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i telekomunikacyjnym")
plot(Rt_b_s_ewma150, type="l", ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i spo¿ywczym")
plot(Rt_s_t_ewma150, type="l", ylab="korelacje warunkowe", main="korelacje miêdzy sektorem spo¿ywczym i telekomunikacyjnym")

##################################################
###################DCC############################
##################################################

Rt<-rep(0, times=3*3*length(b_e))
dim(Rt)<-c(3,3,3448)
Rt[1,1,]<-1
Rt[2,2,]<-1
Rt[3,3,]<-1
R[1,2,]<-Rt_b_s_ewma->R[2,1,]
R[1,3,]<-Rt_b_t_ewma->R[3,1,]
R[2,3,]<-Rt_s_t_ewma->R[3,2,]

DCC <-function(Dt,R){
for (i in 1:dim(Dt)[3]-1){
H<-Dt
H[,,i]<-Dt[,,i]%*%Rt[,,i]%*%Dt[,,i]
}
H
}
DCC<-DCC(Dt,Rt)

##############################################################
#############DCC-Garch ( w³asne)##############################
##############################################################
#Funkcja wystarczy ¿eby zwraca³a Rt wtedy mogê wstawiæ do tego co jest wy¿ej

Rt_garch_opt<- function(a=c(.1,.2),u1=b_e,u2=s_e,u3=t_e){

#a-wektor parametrów
#Qt bêdzie wymiaru 3x3x3447
#dQt=odw(sqrt(dQt))
#Rt 3x3x3447
#U jest wektorem 3x3447,
U<-rep(NA, 3*(length(u1)-1))
dim(U)<-c(3,length(u1)-1)
U[1,]<-u1[1:(length(u1)-1)]
U[2,]<-u2[1:(length(u2)-1)]
U[3,]<-u3[1:(length(u3)-1)]
Qt<-rep(NA,3*3*(length(u1)-1))
dim(Qt)<-c(3,3,(length(u1)-1))

dQt<-rep(0,3*3*(length(u1)-1))
dim(dQt)<-c(3,3,(length(u1)-1))

R<-rep(1,9)
dim(R)<-c(3,3,1)
R[1,2,1]<-cor(u1,u2)->R[2,1,1]
R[1,3,1]<-cor(u1,u3)->R[3,1,1]
R[2,3,1]<-cor(u2,u3)->R[3,2,1]

Rt<-rep(0, times=3*3*(length(u1)-1))
dim(Rt)<-c(3,3,(length(u1)-1))
Rt[1,1,]<-1
Rt[2,2,]<-1
Rt[3,3,]<-1
Lc<-0
Rt[,,1]<-R[,,1]
#ustalam Q0 jako bezwarunkowa(to samo co R bo wystandaryzowane)
Qt[,,1]<-R[,,1]
#odt¹d zaczyna siê optymalizacja
#a,b s¹ parametrami optymalizowanymi
for (j in 2:dim(U)[2]){
	Qt[,,j]<-(1-a[1]-a[2])*Qt[,,1] + a[1]*(U[,j-1]%*%t(U[,j-1])) + a[2]*(Qt[,,j-1])
	diag(dQt[,,j])<-sqrt(diag(Qt[,,j]))
	Rt[,,j]<-solve(dQt[,,j])%*%Qt[,,j]%*%solve(dQt[,,j])
	#if (!(det(Rt[,,j]>0))) Rt[,,j]<-Rt[,,j-1]#daje radê, ale kiepskie wyniki wychodz¹
	V<-eigen(Rt[,,j])$vectors
	V1<-solve(V)
	D<-V1%*%Rt[,,j]%*%V
	D1<-rep(0,9)
	dim(D1)<-dim(D)
	diag(D1)<-1/diag(D)
	R_1<-V%*%D1%*%V1
	#Lc<- Lc+(log(det(Rt[,,j]))+t(U[,j])%*%solve(Rt[,,j])%*%U[,j])*(-0.5)
	if (det(Rt[,,j])>0) Lc<- Lc+(log(det(Rt[,,j]))+t(U[,j])%*%R_1%*%U[,j])*(-0.5) else Lc<- Lc+(t(U[,j])%*%R_1%*%U[,j])*(-0.5)
#docelowo zwracam Rt ale dopiero po zoptymalizowaniu
}
Lc
}
#trzeba to jakoœ zoptymalizowaæ
#optim(c(.1,.2),fn=Rt_garch, method="L-BFGS-B", lower=0, upper=1)
#nlm(Rt_garch,c(.1,.2),print.level=1)

x11<-constrOptim(c(.1,.1),Rt_garch_opt,NULL,ui=rbind(c(-1,-1),c(1,0),c(0,1)),ci=c(-1,0,0))
x11
p<-x11$par[1]
q<-x11$par[2]
#przed tym zmieniæ ¿eby funkcja zwraca³a Rt
Rt_garch<- function(a=c(p,q),u1=b_e,u2=s_e,u3=t_e){

#a-wektor parametrów
#Qt bêdzie wymiaru 3x3x3447
#dQt=odw(sqrt(dQt))
#Rt 3x3x3447
#U jest wektorem 3x3447,
U<-rep(NA, 3*(length(u1)-1))
dim(U)<-c(3,length(u1)-1)
U[1,]<-u1[1:(length(u1)-1)]
U[2,]<-u2[1:(length(u2)-1)]
U[3,]<-u3[1:(length(u3)-1)]
Qt<-rep(NA,3*3*(length(u1)-1))
dim(Qt)<-c(3,3,(length(u1)-1))

dQt<-rep(0,3*3*(length(u1)-1))
dim(dQt)<-c(3,3,(length(u1)-1))

R<-rep(1,9)
dim(R)<-c(3,3,1)
R[1,2,1]<-cor(u1,u2)->R[2,1,1]
R[1,3,1]<-cor(u1,u3)->R[3,1,1]
R[2,3,1]<-cor(u2,u3)->R[3,2,1]

Rt<-rep(0, times=3*3*(length(u1)-1))
dim(Rt)<-c(3,3,(length(u1)-1))
Rt[1,1,]<-1
Rt[2,2,]<-1
Rt[3,3,]<-1
Lc<-0
Rt[,,1]<-R[,,1]
#ustalam Q0 jako bezwarunkowa(to samo co R bo wystandaryzowane)
Qt[,,1]<-R[,,1]
#odt¹d zaczyna siê optymalizacja
#a,b s¹ parametrami optymalizowanymi
for (j in 2:dim(U)[2]){
	Qt[,,j]<-(1-a[1]-a[2])*Qt[,,1] + a[1]*(U[,j-1]%*%t(U[,j-1])) + a[2]*(Qt[,,j-1])
	diag(dQt[,,j])<-sqrt(diag(Qt[,,j]))
	Rt[,,j]<-solve(dQt[,,j])%*%Qt[,,j]%*%solve(dQt[,,j])
}
Rt
}

Rt<-Rt_garch(a=c(p,q),b_e,s_e,t_e)

##############################################################
##############test na sta³oœæ korelacji#######################
##############################################################
#stworzyc I, vechu
#co to jest R^(-1/2)
const_corr<-function(R=R,u1=b_e,u2=s_e,u3=t_e,k=3,p=5){


U<-rep(NA, 3*(length(u1)))
dim(U)<-c(3,length(u1))
U[1,]<-u1[1:(length(u1))]
U[2,]<-u2[1:(length(u2))]
U[3,]<-u3[1:(length(u3))]

####tworzê R12

V<-eigen(R)$vectors
V1<-solve(V)
D<-V1%*%R%*%V
D12<-rep(0,dim(D)[1]*dim(D)[2])
dim(D12)<-c(dim(D)[1],dim(D)[2])
for (i in 1:dim(D)[1]){
D12[i,i]<-sqrt(D[i,i])
}
R12<-V%*%D12%*%V1
R_12<-solve(R12)

###macierz jednostkowa rzêdu k
I_k<-function(l=k){
I<-c(rep(c(1,rep(0,l)),l-1),1)
dim(I)<-c(l,l)
I
}
I_k<-I_k()

####operator vechu

vechu<-function(A){
#mo¿na to zrobiæ na for(ach), ale tu ¿eby szybciej dzia³a³o zrobiê rêcznie
#k<-dim(A)[1]
c<-rep(NA,3)
c[1]<-A[1,2]
c[2]<-A[1,3]
c[3]<-A[2,3]
dim(c)<-c(3,1)
c
}


###tworzê wektory Yt

Yt<-rep(NA, 3*(length(u1)))
A<-rep(NA, 9*(length(u1)))
dim(Yt)<-c(3,length(u1))
dim(A)<-c(3,3,length(u1))
for (i in 1:dim(Yt)[2]){
A[,,i]<-(R12%*%U[,i])%*%t(R_12%*%U[,i])-I_k
Yt[,i]<-vechu(A[,,i])
Yt[,i]
}

#####organizujê dane
#Ysur<-rep(NA,p+1)
#for (i in 0:p){
#paste("Ysur",p-i,sep="")<-c(t(Yt)[(p+2-i):dim(Yt)[2]-i,1],t(Yt)[(p-i+2):dim(Yt)[2]-i,2],t(Yt)[(p+2-i):dim(Yt)[2]-i,3])
#}
#wszystko jest bardzo sztywno zrobione, w szczególnoœci wyniar Ysur, nie dim a length(trrzebaby)
Ysur_5<-c(t(Yt)[2:(dim(Yt)[2]-5),1],t(Yt)[2:(dim(Yt)[2]-5),2],t(Yt)[2:(dim(Yt)[2]-5),3])
e<-rep(1,3*(dim(Yt)[2]-6))
Ysur_4<-c(t(Yt)[3:(dim(Yt)[2]-4),1],t(Yt)[3:(dim(Yt)[2]-4),2],t(Yt)[3:(dim(Yt)[2]-4),3])
Ysur_3<-c(t(Yt)[4:(dim(Yt)[2]-3),1],t(Yt)[4:(dim(Yt)[2]-3),2],t(Yt)[4:(dim(Yt)[2]-3),3])
Ysur_2<-c(t(Yt)[5:(dim(Yt)[2]-2),1],t(Yt)[5:(dim(Yt)[2]-2),2],t(Yt)[5:(dim(Yt)[2]-2),3])
Ysur_1<-c(t(Yt)[6:(dim(Yt)[2]-1),1],t(Yt)[6:(dim(Yt)[2]-1),2],t(Yt)[6:(dim(Yt)[2]-1),3])
Ysur<-c(t(Yt)[7:(dim(Yt)[2]),1],t(Yt)[7:(dim(Yt)[2]),2],t(Yt)[7:(dim(Yt)[2]),3])
dim(Ysur)<-c(10326,1)
Xsur<-t(rbind(e,Ysur_1,Ysur_2,Ysur_3,Ysur_4,Ysur_5))

#####szacujê model SUR , a w³aœciwie szacujê KMNK i liczê reszty
b<-solve(t(Xsur)%*%Xsur)%*%t(Xsur)%*%Ysur
res<- Ysur-(Xsur%*%b)
#####wyliczam statystykê zgodnie ze wzorem

statystyka<-(t(b)%*%(t(Xsur)%*%Xsur)%*%b)/var(res)

#statystyka
#####wartoœci krytyczne 

PVAL <- 1 - pchisq(statystyka, df = p+1)
war_kryt<-qchisq(0.05,6)
nazwy<-c("statystyka", "war. krytyczna(0,05)", "p-value")
wynik<-c(statystyka, war_kryt, PVAL)
dim(wynik)<-c(1,3)
colnames(wynik)<-nazwy
print(wynik)
}

const_corr(R)
a
#################################
####DCC- CCGARCH (niew³asne)####
#################################

library(ccgarch)
U<-rep(NA, 3*(length(b_e)-1))
dim(U)<-c(3,length(b_e)-1)
U[1,]<-u1[1:(length(b_e)-1)]
U[2,]<-u2[1:(length(s_e)-1)]
U[3,]<-u3[1:(length(t_e)-1)]
trU<-t(U)
trU
opt_param<-c(.059,.94)
dim(opt_param)<-c(2,1)
x11<-dcc.estimation2(trU,opt_param)
opt_param1<-x11$par[1]
opt_param2<-x11$par[2]
opt_param<-c(opt_param1,opt_param2)
dim(opt_param)<-c(2,1)
Rt<-dcc.est(trU,opt_param)$DCC

Rt_DCC<-Rt
dim(Rt_DCC)<-c(3,3,dim(Rt)[1])
Rt_DCC[1,1,]<-Rt_DCC[3,3,]<-Rt_DCC[2,2,]<-Rt[,1]
Rt_DCC[1,2,]<-Rt_DCC[2,1,]<-Rt[,2]
Rt_DCC[1,3,]<-Rt_DCC[3,1,]<-Rt[,3]
Rt_DCC[2,3,]<-Rt_DCC[3,2,]<-Rt[,6]
Rt_DCC
loglik.dcc2(opt_param,trU)

##############################################################
#####################Ht#######################################
##############################################################

CC<-function(Dt,R){
if (dim(R)[3]==1){
for (i in 1:dim(Dt)[3]-1) Rt[,,i]<-R[,,1]
}
else Rt<-R
DCC <-function(Dt=Dt,R=Rt){
H<-Dt[,,1:dim(Dt)[3]-1]
for (i in 1:dim(Dt)[3]-1){
H[,,i]<-Dt[,,i]%*%R[,,i]%*%Dt[,,i]
}
H
}
Ht<-DCC(Dt,Rt)
Ht
}

Ht_CCC<-CC(Dt,R)
Ht_DCC_06<-CC(Dt,Rt_06)
Ht_DCC_opt<-CC(Dt,Rt_ewma_opt)
#Ht_DCC_GARCH<-CC(Dt,Rt)
Ht_DCC<-CC(Dt,Rt_DCC)
Ht_DCC

##############################################################
###################wagi#######################################
##############################################################

###wagi metod¹ minimalnej wariancji

minwar_wagi<-function(Ht){
jeden<-rep(1,dim(Ht)[1])
dim(jeden)<-c(1,dim(Ht)[1])
w<-rep(NA,dim(Ht)[1]*dim(Ht)[3])
dim(w)<-c(dim(Ht)[1],dim(Ht)[3])
for(i in 1:dim(Ht)[3]){
Ht_1<-solve(Ht[,,i])
Ct<-jeden%*%Ht_1%*%t(jeden)
w[,i]<-(Ht_1%*%t(jeden))
w[,i]<-w[,i]/Ct
}
w
}

W_CCC<-minwar_wagi(Ht_CCC)
W_DCC_06<-minwar_wagi(Ht_DCC_06)
W_DCC_opt<-minwar_wagi(Ht_DCC_opt)
W_DCC_GARCH<-minwar_wagi(Ht_DCC_GARCH)
W_DCC<-minwar_wagi(Ht_DCC)
Wcc<-minwar_wagi(Rt_06)
######UWAGA!!! niektóre wagi wychodz¹ ujemne - mo¿liwa krótka sprzeda¿

#####prezentacja wag
#### wczeœniej zrobiæ jakby skumulowane
###szybszy dzia³a³ by poprzez przemno¿enie przez wektor logiczny zamiast pêtli

rys_wagi<-function(w, title="wagi dynamiczne", ylab="udzia³y w portfelu", xlab="czas"){
	W<-rep(NA,dim(w)[2]*(dim(w)[1]+1))
	dim(W)<-c(4,dim(w)[2])
	W[1,]<-w[1,]#czerwone - krótkie spo¿ywczy, d³ugie w banku
	W[2,]<-W[1,]+w[2,]#zielone - d³ugie w spo¿ywczy
	W[3,]<-W[2,]+w[3,]#niebieskie - d³ugie w telekom
	for (i in 1:dim(W)[2]){
	W[4,i]<-min(W[1,i],W[2,i]) ###¿ó³te - d³ugie(i krótkie w banku)
	}
plot(W[1,], type="h", col="red", lwd=1, main=title,  ylab=ylab, xlab=xlab, ylim=c(min(W),max(W)))
#plot(W[3,], type="h", col="blue", lwd=1)
lines(W[3,], col="blue",type="h", lwd=1)
lines(W[2,], col="green",type="h", lwd=1)
lines(W[1,], type="h", col="red", lwd=1)
lines(W[4,], type="h", col="yellow", lwd=1)
}

par(mfrow=c(4,1), cex=.7, mar = c(3,4,3,3))
rys_wagi(W_CCC, title="portfel oparty o model CCC")
rys_wagi(W_DCC_06, title="portfel oparty o model DCC_06")
rys_wagi(W_DCC_opt, title="portfel oparty o model DCC_opt")
rys_wagi(W_DCC, title="portfel oparty o model DCC_Garch(1,1)")


##wykres nie dorobiony, to jednak trzeba przeskalowaæ, dodaæ tytu³y, podpisy,
#legendê. podpis
#daæ ¿ó³te równe czerwone ale czerwone nie mo¿

rys_wagi<-function(){
#plot(W[1,], type="h", col="red", lwd=1)
plot(W[3,], type="h", col="blue", lwd=1)
#lines(W[3,], col="blue",type="h", lwd=1)
lines(W[2,], col="green",type="h", lwd=1)
lines(W[1,], type="h", col="red", lwd=1)
lines(W[4,], type="h", col="yellow", lwd=1)
}

#wagi 1(czerwony) - banki, 2(zielony) -spoz, 3(niebieski)-telekom
#ujemne wartoœci czerwonego - krótka sprzeda¿ banki,
#powy¿ej 1 krótka sprzeda¿ telekom
#jak jest krótka sprzeda¿ spozywczego to wtedy czerwone wy¿sze od zielone (w 2489,2490,2491)

lines(rep(0,dim(W)[2]), col="black", lwd=2)

#W<-c(1/3,1/3,1/3)
#W<-c(1/4,1/2,1/4)
W_const<-rep(c(1/3,1/3,1/3),3447)
dim(W_const)<-c(3,3447)
W
t(W)

#####################################################################
#VaR na podstawie macierzy wariancji-kowariancji (VaR-short, long)###
#####################################################################
#Ht-warunkowa macierz wariancji-kowariancji; W-wektor wag, c-sta³a zale¿na od rozk³adu-tu dla poziomu ufnoœci 0,95 rozk³adu normalnego
#Zak³adam, ¿e rozk³ad jest normalny(wiêc i symetryczny wzglêdem 0), wtedy short=-long
VaR<-function(Ht,W,c=1.645){
short<-rep(NA,dim(Ht)[3])
Ht2<-(c^2)*Ht
for (i in 1:dim(Ht2)[3]){
short[i]<-sqrt(t(W[,i])%*%Ht2[,,i]%*%W[,i])
}
short
}

VaR_CCC<-VaR(Ht_CCC,W_CCC)
VaR_DCC_06<-VaR(Ht_DCC_06,W_DCC_06)
VaR_DCC_opt<-VaR(Ht_DCC_opt,W_DCC_opt)
VaR_DCC<-VaR(Ht_DCC,W_DCC)

VaR_CCC<-VaR(Ht_CCC,W_CCC,c=2.33)
VaR_DCC_06<-VaR(Ht_DCC_06,W_DCC_06,c=2.33)
VaR_DCC_opt<-VaR(Ht_DCC_opt,W_DCC_opt,c=2.33)
VaR_DCC<-VaR(Ht_DCC,W_DCC,c=2.33)

##funkcja licz¹ca zwrot z inwestycji przy czym zdecydowanie do rozbudowy na macierze
Zwrot<-function(z1,z2,z3,W){
zwrot<-rep(NA,length(z1)-1)
for (i in 1:length(z1)-1){
zwrot[i]<-W[1,i]*z1[i]+W[2,i]*z2[i]+W[3,i]*z3[i]
}
zwrot
}
zw_CCC<-Zwrot(b,s,t,W_CCC)
zw_DCC_06<-Zwrot(b,s,t,W_DCC_06)
zw_DCC_opt<-Zwrot(b,s,t,W_DCC_opt)
zw_DCC<-Zwrot(b,s,t,W_DCC)

VaR_rys<-function(zw,var,main="var"){
plot(zw, main=main, type="l", ,ylab="zwrot/VaR", xlab="czas")
lines(var, col="red")
lines(-var,col="blue")
}
par(mfrow=c(4,1), cex=.7, mar = c(3,4,3,3))
VaR_rys(zw_CCC,VaR_CCC,"Na podstawie modelu CCC")
VaR_rys(zw_DCC_06,VaR_DCC_06,"Na podstawie modelu DCC_06")
VaR_rys(zw_DCC_opt,VaR_DCC_opt,"Na podstawie modelu DCC_opt")
VaR_rys(zw_CCC,VaR_CCC,"Na podstawie modelu CCC")

###################################################################
####################testy VaR######################################
###################################################################
#tworzê funkcjê zliczaj¹c¹ przekroczenia trzeba podaæ szereg zwrotu i var
HIT<-function(zw, VaR){
spos=0  #przekroczenia pozycji krótkich
sneg=0  #przekroczenia pozycji d³ugich

	for (i in 1:length(zw)){
		if (zw[i]<(-VaR[i])) sneg<-sneg+1 else{
			if (zw[i]>VaR[i]) spos<-spos+1
			}
	}
HIT<-rep(NA,4)
dim(HIT)<-c(2,2)
HIT[1,1]<-spos
HIT[1,2]<-sneg
HIT[2,1]<-spos/length(zw)
HIT[2,2]<-sneg/length(zw)
colnames(HIT)<-c("short","long")
rownames(HIT)<-c("iloœæ","udzia³")
return(HIT)
}
HIT_CCC<-HIT(zw_CCC,VaR_CCC)
HIT_DCC_06<-HIT(zw_DCC_06,VaR_DCC_06)
HIT_DCC_opt<-HIT(zw_DCC_opt,VaR_DCC_opt)
HIT_CCC
HIT_DCC_06
HIT_DCC_opt

######################################
##dokoñczyæ testy VaR
###########################################
#################test Kupca################
###########################################
kupiec<-function(zw, VaR, p=0.05){
hit<-HIT(zw,VaR)
Tpos<-length(zw)-hit[1,1]
Tneg<-length(zw)-hit[1,2]
statystyka_pos<-(-2)*(Tpos*log((1-p)/(1-hit[2,1]))+hit[1,1]*log(p/hit[2,1]))
statystyka_neg<-(-2)*(Tneg*log((1-p)/(1-hit[2,2]))+hit[1,2]*log(p/hit[2,2]))

#coœ jest nie tak albo ze statystykami albo najprawdopodobniej z wart krytyczn¹
kryt<-qchisq(1-p,1)
pval_pos<- 1 - pchisq(statystyka_pos, df = 1)
pval_neg<- 1 - pchisq(statystyka_neg, df = 1)

wynik<-list(hit,statystyka_pos,statystyka_neg,pval_pos,pval_neg,kryt)
return(wynik)
}

kupiec250<-function(zw, VaR, p=0.05){
hit<-HIT(zw[3197:3447],VaR[3197:3447])
Tpos<-250-hit[1,1]
Tneg<-250-hit[1,2]
statystyka_pos<-(-2)*(Tpos*log((1-p)/(1-hit[2,1]))+hit[1,1]*log(p/hit[2,1]))
statystyka_neg<-(-2)*(Tneg*log((1-p)/(1-hit[2,2]))+hit[1,2]*log(p/hit[2,2]))

#coœ jest nie tak albo ze statystykami albo najprawdopodobniej z wart krytyczn¹
kryt<-qchisq(1-p,1)
pval_pos<- 1 - pchisq(statystyka_pos, df = 1)
pval_neg<- 1 - pchisq(statystyka_neg, df = 1)

wynik<-list(hit,statystyka_pos,statystyka_neg,pval_pos,pval_neg,kryt)
return(wynik)
}

kupiec(zw_CCC,VaR_CCC, p=0.01)
kupiec(zw_DCC_06,VaR_DCC_06, p=0.01)
kupiec(zw_DCC_opt,VaR_DCC_opt, p=0.01)
kupiec(zw_DCC,VaR_DCC, p=0.01)


##komitet bazylejski opiera siê tylko na liczbie przekroczeñ z ostatniego roku(250 dni), do 4 model zielony, 5-9 model ¿ó³ty, wiêcej model czerwony(z³y)(przy 1%)
kupiec250(zw_CCC,VaR_CCC)
kupiec250(zw_DCC_06,VaR_DCC_06)
kupiec250(zw_DCC_opt,VaR_DCC_opt)
kupiec250(zw_DCC,VaR_DCC)

kupiec250(zw_CCC,VaR_CCC, p=0.01)
kupiec250(zw_DCC_06,VaR_DCC_06, p=0.01)
kupiec250(zw_DCC_opt,VaR_DCC_opt, p=0.01)
kupiec250(zw_DCC,VaR_DCC, p=0.01)


christoffersen<-function(zw, VaR, p=0.05){
#dla pozycji short
T00<-T01<-T11<-T10<-0
if (zw[1]>VaR[1]) T01<-1 else T00<-1
	for (i in 2:length(zw)){
			if ((zw[i]>VaR[i])&(zw[i-1]>VaR[i-1])) T11<-T11+1 
			if ((zw[i]>VaR[i])&(!(zw[i-1]>VaR[i-1]))) T01<-T01+1
			if ((!(zw[i]>VaR[i]))&(zw[i-1]>VaR[i-1])) T10<-T10+1
			if ((!(zw[i]>VaR[i]))&(!(zw[i-1]>VaR[i-1]))) T00<-T00+1
			}
qkr<-(T01+T11)/(T01+T11+T00+T10)
q10<-T10/(T10+T11)
q01<-T01/(T00+T01)
q11<-T11/(T10+T11)
statystyka <- (-2)*((T00+T10)*log(1-qkr)+(T11+T01)*log(qkr)-T00*log(1-q01)-T01*log(q01)-T10*log(1-q11)-T11*q11)
pval<- 1 - pchisq(statystyka, df = 1)
kryt<-qchisq(1-p,1)
short<- c(statystyka, kryt, pval)


###dla pozycji d³ugich
T00<-T01<-T11<-T10<-0
if (zw[1]>VaR[1]) T01<-1 else T00<-1
	for (i in 2:length(zw)){
			if ((-zw[i]>VaR[i])&(-zw[i-1]>VaR[i-1])) T11<-T11+1 
			if ((-zw[i]>VaR[i])&(!(-zw[i-1]>VaR[i-1]))) T01<-T01+1
			if ((!(-zw[i]>VaR[i]))&(-zw[i-1]>VaR[i-1])) T10<-T10+1
			if ((!(-zw[i]>VaR[i]))&(!(-zw[i-1]>VaR[i-1]))) T00<-T00+1
			}
qkr<-(T01+T11)/(T01+T11+T00+T10)
q10<-T10/(T10+T11)
q01<-T01/(T00+T01)
q11<-T11/(T10+T11)
statystyka <- (-2)*((T00+T10)*log(1-qkr)+(T11+T01)*log(qkr)-T00*log(1-q01)-T01*log(q01)-T10*log(1-q11)-T11*q11)
pval<- 1 - pchisq(statystyka, df = 1)
long<- c(statystyka, kryt, pval)

wynik<-rbind(short, long)
colnames(wynik)<-c("statystyka testowa", "wartoœæ krytyczna", "p-value")
rownames(wynik)<-c("short","long")
return(wynik)
}

christoffersen(zw_CCC,VaR_CCC)
christoffersen(zw_DCC_06,VaR_DCC_06)
christoffersen(zw_DCC_opt,VaR_DCC_opt)
christoffersen(zw_DCC,VaR_DCC)


####zakumulowany hit ma na celu pokazanie gdzie jest przyrost przekroczeñ i porównanie z rozk³adem normalnym##
accum_hit<-function(zw,VaR,c=2.33){
	spos<-sneg<-
	x<-rnorm(3447*10,0,1)
	for (i in 2:length(zw)){
		if (zw[i]<(-VaR[i])) sneg[i]<-(sneg[i-1]+1) else sneg[i]<-sneg[i-1]
		if (zw[i]>VaR[i]) spos[i]<-(spos[i-1]+1) else spos[i]<-spos[i-1]
		if (x[i]>c) y[i]<-y[i-1]+1 else y[i]<-y[i-1]
			}
	S<-rbind(sneg,spos,y)
	return(S)
}

S<-accum_hit(zw_DCC,VaR_DCC, c=2.32)
plot(S[1,],t="l")
lines(S[2,],col="blue")
 
SIM_norm<-function(){
SIM<-rep(0,3447)
for (k in 1:1000){
y<-rep(0,3447)
x<-rnorm(3447,0,1)
for (i in 2:3447){
		if (x[i]>2.34) y[i]<-y[i-1]+1 else y[i]<-y[i-1]
		}
SIM<-rbind(SIM,y)
}
for (j in 1:3447)	k[j]<-mean(SIM[2:1000,j])
return(k)
}
##SIM_norm symulacyjnie wyznacza œrednie skumulowane przekroczenia dla 3447
##obserwacji , dalej mo¿na rozbudowaæ o przedzia³y ufnoœci

SIM<-SIM_norm()
lines(SIM,col="red")

#lines(VaR_DCC*100)
#lines(zw*100,col="blue")
#lines(sqrt(VaR_DCC*100))



