
library(tseries)
library(fGarch)
library(urca)
library(ccgarch)

rm(list=ls())
setwd("d:/dokumenty/Rfiles/do magisterki")
getwd()
set.seed(123)
ls()
dane <- read.csv("bst.csv", header=TRUE, sep=";", dec=".")
attach(dane)
ls(2)
#z_banki z_spoz z_tele
b<-z_banki
s<-z_spoz
t<-z_tele
time<-data
detach()
