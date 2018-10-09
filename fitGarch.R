setwd("D:/dokumenty/Rfiles/do magisterki")
X<-read.csv2("zeszyt1.csv")
X_ts<-ts(X)
attach(X)
sns<-X$sns
data<-data
zap<-zap
detach(X)
ls()
class(X)
sns_ts<-ts(sns)
library(fGarch)

#funkcja do automatycznego wyboru odpowiedniego modelu unigarch
fitGarch<-function(x,ar=4,ma=4,arch=3,garch=3){
for a in 0:ar{
	for b in 0:ma{
		for c in 0:arch{
			for d in 0:garch{
			if a=0&&b=0&&c=0&&d=0{
m<-..............................................powien bialy szum
}
else{
if c=0&&d=0 ......... else{

m1<-garchFit(~arma(a,b)+garch(c,d), x)
m
}}
}}}}
}
fitGarch(sns_ts)