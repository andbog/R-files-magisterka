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

CCC <-function(Dt,R){
H<-Dt
for (i in 1:dim(Dt)[3]){
H[,,i]<-Dt[,,i]%*%R%*%Dt[,,i]
}
H
}

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
##############################################################
##############test na sta³oœæ korelacji#######################
##############################################################
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

rys_wagi<-function(w, title="wagi dynamiczne", ylab="udzia³y w portfelu", xlab="czas"){
	W<-rep(NA,dim(w)[2]*(dim(w)[1]+1))
	dim(W)<-c(4,dim(w)[2])
	W[1,]<-w[1,]#czerwone - krótkie spo¿ywczy, d³ugie w banku
	W[2,]<-W[1,]+w[2,]#zielone - d³ugie w spo¿ywczy
	W[3,]<-W[2,]+w[3,]#niebieskie - d³ugie w telekom
	for (i in 1:dim(W)[2]){
	W[4,i]<-min(W[1,i],W[2,i]) ###¿ó³te - d³ugie(i krótkie w banku)	
}
plot(W[1,], type="h", col="red", lwd=1, main=title,  ylab=ylab, xlab=xlab, ylim=c(min(W),max(W)),axes=FALSE)
lines(W[3,], col="blue",type="h", lwd=1)
lines(W[2,], col="green",type="h", lwd=1)
lines(W[1,], type="h", col="red", lwd=1)
lines(W[4,], type="h", col="yellow", lwd=1)
axis(2)
axis(1,at=time1,labels=time2)
rm(W)
}

VaR<-function(Ht,W,c=1.645){
short<-rep(NA,dim(Ht)[3])
Ht2<-(c^2)*Ht
for (i in 1:dim(Ht2)[3]){
short[i]<-sqrt(t(W[,i])%*%Ht2[,,i]%*%W[,i])
}
short
}

Zwrot<-function(z1,z2,z3,W){
zwrot<-rep(NA,length(z1)-1)
for (i in 1:length(z1)-1){
zwrot[i]<-W[1,i]*z1[i]+W[2,i]*z2[i]+W[3,i]*z3[i]
}
zwrot
}

VaR_rys<-function(zw,var,main="var"){
plot(zw, main=main, type="l", ,ylab="zwrot/VaR", xlab="czas", axes=FALSE)
lines(var, col="red",lwd=2)
lines(-var,col="blue",lwd=2)
axis(2)
axis(1,at=time1,labels=time2)
}

HIT<-function(zw, VaR){
spos=0  #przekroczenia pozycji krótkich
sneg=0  #przekroczenia pozycji d³ugich

	for (i in 1:length(VaR)){
		if (zw[i]<(-VaR[i])) sneg<-sneg+1 else{
			if (zw[i]>VaR[i]) spos<-spos+1
			}
	}
HIT<-rep(NA,4)
dim(HIT)<-c(2,2)
HIT[1,1]<-spos
HIT[1,2]<-sneg
HIT[2,1]<-spos/length(VaR)
HIT[2,2]<-sneg/length(VaR)
colnames(HIT)<-c("short","long")
rownames(HIT)<-c("iloœæ","udzia³")
return(HIT)
}

kupiec<-function(zw, VaR, p=0.05){
hit<-HIT(zw,VaR)
Tpos<-length(VaR)-hit[1,1]
Tneg<-length(VaR)-hit[1,2]
statystyka_pos<-(-2)*(Tpos*log((1-p)/(1-hit[2,1]))+hit[1,1]*log(p/hit[2,1]))
statystyka_neg<-(-2)*(Tneg*log((1-p)/(1-hit[2,2]))+hit[1,2]*log(p/hit[2,2]))
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

christoffersen<-function(zw, VaR, p=0.05){
#dla pozycji short
T00<-T01<-T11<-T10<-0
if (zw[1]>VaR[1]) T01<-1 else T00<-1
	for (i in 2:length(VaR)){
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
	for (i in 2:length(VaR)){
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

####zakumulowany hit ma na celu pokazanie gdzie jest przyrost przekroczeñ i porównanie z rozk³adem normalnym##
accum_hit<-function(zw,VaR,c=1.64){
	spos<-sneg<-rep(0,length(zw)-1)
	#x<-rnorm(3447*10,0,1)
	#y<-rep(0,length(zw)-1)
	for (i in 2:length(zw)){
		if (zw[i]<(-VaR[i])) sneg[i]<-(sneg[i-1]+1) else sneg[i]<-sneg[i-1]
		if (zw[i]>VaR[i]) spos[i]<-(spos[i-1]+1) else spos[i]<-spos[i-1]
	#	if (x[i]>c) y[i]<-y[i-1]+1 else y[i]<-y[i-1]
			}
	S<-rbind(sneg,spos)
	return(S)
}
SIM_norm<-function(c=1.64){
SIM<-rep(0,3447)
for (k in 1:1000){
y<-rep(0,3447)
x<-rnorm(3447,0,1)
for (i in 2:3447){
		if (x[i]>c) y[i]<-y[i-1]+1 else y[i]<-y[i-1]
		}
SIM<-rbind(SIM,y)
}
a<-q05<-q95<-rep(0,3447)
for (j in 1:3447){
			a[j]<-mean(SIM[2:1000,j])
			q05[j]<-quantile(SIM[2:1000,j],.05)
			q95[j]<-quantile(SIM[2:1000,j],.95)
			l<-rbind(q05,a,q95)
			}
return(l)
}

#### funkcja poni¿ej wylicza dopasowane wartoœci, póki co zrobione tak, ¿e poprostu odejmujê sk³adnik losowy od szeregu, ale mo¿na podstawiaj¹c dane
arima.fit<-function(model,szereg){
fit<-szereg - model$r
fit
}

##odwrócenie szeregu zwrotów w wartoœci; pobieram zwroty i wartosc poczatkow¹, a zwracam wartoœæ koñcow¹.

wartosci<-function(zwroty_log,wart0=1297.56){
wartosci<-rep(0,3448)
wartosci[1]<-wart0
for (i in 2:3448){
wartosci[i]<-wartosci[i-1]*exp(zwroty_log[i-1])
}
return(wartosci)
}


