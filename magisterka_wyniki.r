st(b)
st(s)
st(t)
time1<-c(2,251,501,749,999,1248,1497,1748,2003,2254,2506,2755,3006,3257,3508)
time2<-c(1997:2011)
time2
obs_no<-c(1:length(time))
obs_no
time
attach(dane)
par(mfrow=c(3,1), cex=.7)
plot(banki, main="Sektor bankowy", type="l",ylab='wartoœæ indeksu', xlab='czas', axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(spoz, main="Sektor spo¿ywczy", type="l",ylab='wartoœæ indeksu', xlab='czas',axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(tele, main="Sektor telekomunikacyjny", type="l",ylab='wartoœæ indeksu', xlab='czas',axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
detach(dane)

par(mfrow=c(3,1), cex=.7)
plot(b, main="Sektor bankowy", type="l",ylab='zwroty', xlab='czas', axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(s, main="Sektor spo¿ywczy", type="l",ylab='zwroty', xlab='czas',axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(t, main="Sektor telekomunikacyjny", type="l",ylab='zwroty', xlab='czas',axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)


par(mfrow=c(3,1), cex=.7)
acf(b, main="Sektor bankowy", ylim=c(-.2,.2), lwd=3)
acf(s, main="Sektor spo¿ywczy",ylim=c(-.2,.2), lwd=3)
acf(t, main="Sektor telekomunikacyjny",ylim=c(-.2,.2), lwd=3)

Box.test(t, lag=1)

summary(ur.df(t,type="none",lags=2,selectlags="AIC"))
summary(ur.kpss(t,type="mu", lags="long", use.lag=NULL))

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

b_r<-arma_b$r
s_r<-arma_s$r
t_r<-arma_t$r
par(mfrow=c(3,1), cex=.7)
plot(b_r^2, main="Sektor bankowy",ylab="kwadraty reszt", type="l", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(s_r^2, main="Sektor spo¿ywczy",ylab="kwadraty reszt", type="l", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(t_r^2, main="Sektor telekomunikacyjny",ylab="kwadraty reszt", type="l",axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)

par(mfrow=c(3,1), cex=.7)
acf(b_r, main="Sektor bankowy", ylim=c(-.2,.2), lwd=3)
acf(s_r, main="Sektor spo¿ywczy",ylim=c(-.2,.2), lwd=3)
acf(t_r, main="Sektor telekomunikacyjny",ylim=c(-.2,.2), lwd=3)

archlmtest(b_r,5)
archlmtest(s_r,5)
archlmtest(t_r,5)

GARCH(b,ar=0,ma=1,pmax=3,qmax=3)
GARCH(s,ar=1,ma=0,pmax=3,qmax=3)
GARCH(t,ar=0,ma=1,pmax=3,qmax=3)

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
plot(b_u, main="Sektor bankowy", type="l",xlab='czas', ylab="zwroty z 95% przedzia³ami ufnoœci",axes=FALSE)
lines(a*b_sig, col="red",lwd=2)
lines(-a*b_sig,col="blue",lwd=2)
axis(2)
axis(1, at=time1, labels=time2)
plot(s_u, main="Sektor spo¿ywczy", type="l",xlab='czas', ylab="zwroty z 95% przedzia³ami ufnoœci",axes=FALSE)
lines(a*s_sig, col="red", lwd=2)
lines(-a*s_sig,col="blue", lwd=2)
axis(2)
axis(1, at=time1, labels=time2)
plot(t_u, main="Sektor telekomunikacyjny", type="l",xlab='czas', ylab="zwroty z 95% przedzia³ami ufnoœci", axes=FALSE)
lines(a*t_sig, col="red",lwd=2)
lines(-a*t_sig,col="blue",lwd=2)
axis(2)
axis(1,at=time1,labels=time2)

a<-1.95
par(mfrow=c(3,1), cex=.7)
plot(b_e, main="Sektor bankowy", type="l",xlab='czas', ylab="wystandaryzowane reszty",axes=FALSE)
axis(2)
axis(1, at=time1, labels=time2)
plot(s_e, main="Sektor spo¿ywczy", type="l",xlab='czas', ylab="wystandaryzowane reszty",axes=FALSE)
axis(2)
axis(1, at=time1, labels=time2)
plot(t_e, main="Sektor telekomunikacyjny", type="l",xlab='czas', ylab="wystandaryzowane reszty", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)

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

Dt<-rep(0, times=3*3*length(b_e))
dim(Dt)<-c(3,3,3448)
Dt[1,1,]<-b_sig
Dt[2,2,]<-s_sig
Dt[3,3,]<-t_sig

R<-rep(1,9)
dim(R)<-c(3,3)
R[1,2]<-cor(b_e,s_e)->R[2,1]
R[1,3]<-cor(b_e,t_e)->R[3,1]
R[2,3]<-cor(s_e,t_e)->R[3,2]
R

const_corr(R)

optima<-optimize(f=ewma_opt, interval=c(0,1))
opt<-optima$minimum
opt

############DCC#######
U<-rep(NA, 3*(length(b_e)-1))
dim(U)<-c(3,length(b_e)-1)
u1<-b_e
u2<-s_e
u3<-t_e
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
######################

Rt_06<-ewma()
Rt_ewma_opt<-ewma(lambda=opt)

par(mfrow=c(3,1), cex=.7)
plot(Rt_06[1,3,], type="l",xlab='czas', ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i telekomunikacyjnym", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(Rt_06[1,2,], type="l",xlab='czas', ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i spo¿ywczym",axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(Rt_06[2,3,], type="l",xlab='czas', ylab="korelacje warunkowe", main="korelacje miêdzy sektorem spo¿ywczym i telekomunikacyjnym", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)

par(mfrow=c(3,1), cex=.7)
plot(Rt_ewma_opt[1,3,], type="l",xlab='czas', ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i telekomunikacyjnym", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(Rt_ewma_opt[1,2,], type="l",xlab='czas', ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i spo¿ywczym", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(Rt_ewma_opt[2,3,], type="l",xlab='czas', ylab="korelacje warunkowe", main="korelacje miêdzy sektorem spo¿ywczym i telekomunikacyjnym", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)

par(mfrow=c(3,1), cex=.7)
plot(Rt_DCC[1,3,], type="l",xlab='czas', ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i telekomunikacyjnym", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(Rt_DCC[1,2,], type="l",xlab='czas', ylab="korelacje warunkowe", main="korelacje miêdzy sektorem bankowym i spo¿ywczym", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(Rt_DCC[2,3,], type="l",xlab='czas', ylab="korelacje warunkowe", main="korelacje miêdzy sektorem spo¿ywczym i telekomunikacyjnym", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)


Ht_CCC<-CCC(Dt,R)
Ht_DCC_06<-CC(Dt,Rt_06)
Ht_DCC_opt<-CC(Dt,Rt_ewma_opt)
Ht_DCC<-CC(Dt,Rt_DCC)

Ht_CCC
W_CCC<-minwar_wagi(Ht_CCC)
W_DCC_06<-minwar_wagi(Ht_DCC_06)
W_DCC_opt<-minwar_wagi(Ht_DCC_opt)
W_DCC<-minwar_wagi(Ht_DCC)

par(mfrow=c(4,1), cex=.7, mar = c(3,4,3,3))
rys_wagi(W_CCC, title="struktura portfela otrzymana na podstawie modelu CCC")
rys_wagi(W_DCC_06, title="struktura portfela otrzymana na podstawie modelu DCC_0,06")
rys_wagi(W_DCC_opt, title="struktura portfela otrzymana na podstawie modelu DCC_opt")
rys_wagi(W_DCC, title="struktura portfela otrzymana na podstawie modelu DCC_Garch(1,1)")

VaR_CCC<-VaR(Ht_CCC,W_CCC)
VaR_DCC_06<-VaR(Ht_DCC_06,W_DCC_06)
VaR_DCC_opt<-VaR(Ht_DCC_opt,W_DCC_opt)
VaR_DCC<-VaR(Ht_DCC,W_DCC)

zw_CCC<-Zwrot(b,s,t,W_CCC)
zw_DCC_06<-Zwrot(b,s,t,W_DCC_06)
zw_DCC_opt<-Zwrot(b,s,t,W_DCC_opt)
zw_DCC<-Zwrot(b,s,t,W_DCC)

par(mfrow=c(4,1), cex=.7, mar = c(3,4,3,3))
VaR_rys(zw_CCC,VaR_CCC,"VaR otrzymana na podstawie modelu CCC")
VaR_rys(zw_DCC_06,VaR_DCC_06,"VaR otrzymana na podstawie modelu DCC_0,06")
VaR_rys(zw_DCC_opt,VaR_DCC_opt,"VaR otrzymana na podstawie modelu DCC_opt")
VaR_rys(zw_DCC,VaR_DCC,"VaR otrzymana na podstawie modelu DCC_Garch(1,1)")

wart_CCC<-wartosci(zw_CCC)
wart_DCC_06<-wartosci(zw_DCC_06)
wart_DCC_opt<-wartosci(zw_DCC_opt)
wart_DCC<-wartosci(zw_DCC)

par(mfrow=c(4,1), cex=.7)
plot(wart_CCC, type="l",xlab='czas', ylab="wartoœæ portfela", main="Wartoœæ portfela, którego strukturê otrzymano na podstawie modelu CCC", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(wart_DCC_06, type="l",xlab='czas', ylab="wartoœæ portfela", main="Wartoœæ portfela, którego strukturê otrzymano na podstawie modelu DCC_0,06", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(wart_DCC_opt, type="l",xlab='czas', ylab="wartoœæ portfela", main="Wartoœæ portfela, którego strukturê otrzymano na podstawie modelu DCC_opt", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)
plot(wart_DCC, type="l",xlab='czas', ylab="wartoœæ portfela", main="Wartoœæ portfela, którego strukturê otrzymano na podstawie modelu DCC_Garc(1,1)", axes=FALSE)
axis(2)
axis(1,at=time1,labels=time2)


kupiec(zw_CCC,VaR_CCC[1:3447], p=0.05)
kupiec(zw_DCC_06,VaR_DCC_06, p=0.05)
kupiec(zw_DCC_opt,VaR_DCC_opt, p=0.05)
kupiec(zw_DCC,VaR_DCC, p=0.05)

kupiec250(zw_CCC,VaR_CCC, p=0.05)
kupiec250(zw_DCC_06,VaR_DCC_06, p=0.05)
kupiec250(zw_DCC_opt,VaR_DCC_opt, p=0.05)
kupiec250(zw_DCC,VaR_DCC, p=0.05)

christoffersen(zw_CCC,VaR_CCC)
christoffersen(zw_DCC_06,VaR_DCC_06)
christoffersen(zw_DCC_opt,VaR_DCC_opt)
christoffersen(zw_DCC,VaR_DCC)


S_CCC<-accum_hit(zw_CCC,VaR_CCC, c=1.64)
S_DCC<-accum_hit(zw_DCC,VaR_DCC, c=1.64)

SIM<-SIM_norm()

par(mfrow=c(2,1), cex=.7, mar = c(3,4,3,3))

plot(S_CCC[1,],t="l", ylab="zakumulowane przekroczenia", xlab="czas", col="blue", main="Skumulowany szereg przekroczeñ przy wykorzystaniu modelu CCC",axes=FALSE)
lines(S_CCC[2,],col="red")
lines(SIM[2,],col="black")
lines(SIM[1,],col="green")
lines(SIM[3,],col="green")
axis(2)
axis(1,at=time1,labels=time2)

plot(S_DCC[1,],t="l", ylab="zakumulowane przekroczenia", xlab="czas", col="blue",main="Skumulowany szereg przekroczeñ przy wykorzystaniu modelu DCC_GARCH(1,1)",axes=FALSE)
lines(S_DCC[2,],col="red")
lines(SIM[2,],col="black")
lines(SIM[1,],col="green")
lines(SIM[3,],col="green")
axis(2)
axis(1,at=time1,labels=time2)


###próbujê z autoregresj¹ 
VaRARMA<-arima(VaR_DCC, c(1,0,1))
VaRAR<-arima(VaR_DCC, c(1,0,0))

plot(S_DCC[1,],t="l", ylab="zakumulowane przekroczenia", xlab="czas", col="blue",main="Szereg przekroczeñ wynikaj¹cy z modelu DCC_GARCH")
lines(S_DCC[2,],col="red")
lines(SIM[2,],col="black")
lines(SIM[1,],col="green")
lines(SIM[3,],col="green")
