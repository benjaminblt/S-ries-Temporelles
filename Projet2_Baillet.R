# Librairies

library(ggplot2)
library(forecast)
library(caschrono)
library(lmtest)
library(urca)
library(CADFtest)
library(dplyr)
library(foreach)
library(doSNOW)
library(parallel)
library(FinTS) 
library(TSA)
library(Hmisc)

# Avant de continuer, faire:
# - "Session"
# - "Set Working Directory"
# - "Choose Directory" (choisir emplacement du fichier "pétrole.csv")

# Charger les données depuis le fichier CSV
rm(list=ls())
tab <- read.csv("pétrole.csv", header=TRUE, dec=".", sep=",")

# Sélectionner les colonnes nécessaires
tab2 <- subset(tab, select=c(LOCATION, TIME, Value))

# Extraire les données pour la Thaïlande et éliminer les valeurs manquantes
thailande <- tab2[tab2$LOCATION=="THA",]
thailande <- na.omit(thailande)

# Convertir la colonne TIME en type numérique 
thailande$TIME <- as.numeric(as.character(thailande$TIME))

Petrole<-ts(na.omit(thailande[,"Value"]),start=1971,freq=1)

###################################################################
# Visualisation des données
###################################################################

plot(Petrole,type='l',col=4)

###################################################################
# Stationnarisation de la série
###################################################################

# Selon le projet 1, on sait que le PGD est DS
# En différenciant 1 fois la série devient stationnaire

dPetrole = diff(Petrole)

# On vérifie que la tendance a été enlevée
plot(dPetrole,type='l',col=4)

###################################################################
# Chronogramme
###################################################################

# chronogramme série Petrole
petrole<-ts(na.omit(thailande[,"Value"]),start=1971,freq=2, end= 2021)
decompo<-decompose(petrole,"multiplicative", filter = NULL)
plot(decompo)

# choronogramme série différenciée
decompo<-decompose(diff(petrole),"additive", filter = NULL)
plot(decompo)

###################################################################
# ACF 
###################################################################

acf(dPetrole)
# autocorrélation ordre 3 et 15

###################################################################
# PACF 
###################################################################

pacf(dPetrole)
# autocorrélation ordre 3

###################################################################
# EACF 
###################################################################

eacf(dPetrole)

# cela nous donne une indication des valeurs de p et q
# p = 0 , q = 3

###################################################################
# Estimation ARIMA (p,d,q)  pour dPetrole
###################################################################

# On utilise  ARIMA et pas ARMA car le PGD de la série Petrole est DS

# d = 1 car on différencie une fois ici

# estimation du modèle

reg_Petrole = Arima(Petrole, order = c(0,1,3))
coeftest(reg_Petrole) 
reg_Petrole = Arima(Petrole, order = c(0,1,3), fixed=c(0,0,NA))
coeftest(reg_Petrole)
# Les coefficients sont tous significatifs désormais

###################################################################
# On vérifie que les aléas sont des bruits blancs
###################################################################

#### t.test (pvalue > 5%) ####

t.test(reg_Petrole$residuals) 
# p-value = 0,1666 > 0,05
# On accepte H0
# Donc la moyenne vaut zero

#### ArchTest (pvalue > 5%) ####

ArchTest(reg_Petrole$residuals, lag = 1)


#### Box.test (pvalue < 5%) ####

Box.test(reg_Petrole$residuals,  lag = 1, type = "Ljung-Box")

###################################################################
# Trouver le meilleur modèle
###################################################################

# L'eacf donne une indiciation des valeurs de p et q 
# mais pas les meilleures valeurs
# donc  fait varier p et q

# On garde le BIC le plus petit pour trouver le meilleur modèle
# Pour chaque modèle on fait ces 3 tests et on doit obtenir
# ce qu'on a entre parenthèses pour que ça marche

#### Coeftest ####

reg0 = Arima(Petrole, order = c(0,1,1))
coeftest(reg0)
# pas le bon modèle 
# coeff de ma1 non significatif

reg1 = Arima(Petrole, order = c(0,1,2))
coeftest(reg1)
# tous les coeffs sont significatifs


reg2 = Arima(Petrole, order = c(0,1,3))
coeftest(reg2)
reg2 = Arima(Petrole, order = c(0,1,3), fixed=c(0,0,NA))
coeftest(reg2)
# tous les coeffs sont significatifs


reg3 = Arima(Petrole, order = c(0,1,4))
coeftest(reg3)
# pas le bon modèle 
# on baisse q de 4 à 3 désormais

reg4 = Arima(Petrole, order = c(0,1,7))
coeftest(reg4)
reg4 = Arima(Petrole, order = c(0,1,7), fixed=c(0,0,NA,0,0,0,NA))
coeftest(reg4)
# tous les coeffs sont significatifs


reg5 = Arima(Petrole, order = c(1,1,1))
coeftest(reg5)
# tous les coeffs sont significatifs


reg6 = Arima(Petrole, order = c(1,1,2))
coeftest(reg6)
# tous les coeffs sont significatifs


reg7 = Arima(Petrole, order = c(1,1,3))
coeftest(reg7)
reg7 = Arima(Petrole, order = c(1,1,3), fixed=c(0,0,0,NA))
coeftest(reg7)
# tous les coeffs sont significatifs


reg8 = Arima(Petrole, order = c(1,1,4))
coeftest(reg8)
# aucun coefficient n'est significatif
# Donc ARIMA(1,1,4) n'est pas le bon modèle


reg9 = Arima(Petrole, order = c(2,1,1))
coeftest(reg9)
reg9 = Arima(Petrole, order = c(2,1,1), fixed=c(NA,0,NA))
coeftest(reg9)
# tous les coeffs sont significatifs


reg10 = Arima(Petrole, order = c(2,1,2))
coeftest(reg10)
reg10 = Arima(Petrole, order = c(2,1,2), fixed=c(NA,0,NA,NA))
coeftest(reg10)
# tous les coeffs sont significatifs


reg11 = Arima(Petrole, order = c(2,1,3))
coeftest(reg11)
reg11 = Arima(Petrole, order = c(2,1,3), fixed=c(0,0,0,0,NA))
coeftest(reg11)
# tous les coeffs sont significatifs

reg12 = Arima(Petrole, order = c(3,1,1))
coeftest(reg12)
# On baisse q de 1 à 0 car le coeff de ma1 non significatif
reg12 = Arima(Petrole, order = c(3,1,0))
coeftest(reg12)
reg12 = Arima(Petrole, order = c(3,1,0), fixed=c(0,0,NA))
coeftest(reg12)
# tous les coeffs sont significatifs


reg13 = Arima(Petrole, order = c(1,1,0))
coeftest(reg13)
# tous les coeffs sont significatifs


reg14 = Arima(Petrole, order = c(2,1,0))
coeftest(reg14)
reg14 = Arima(Petrole, order = c(2,1,0), fixed=c(0,NA))
coeftest(reg14)
# tous les coeffs sont significatifs


reg15= Arima(Petrole, order = c(3,1,2))
coeftest(reg15)
reg15= Arima(Petrole, order = c(3,1,2), fixed=c(0,0,NA,0,NA))
coeftest(reg15)
# tous les coeffs sont significatifs


reg16 = Arima(Petrole, order = c(4,1,1))
coeftest(reg16)
# on baisse q de 1 à 0
reg16 = Arima(Petrole, order = c(4,1,0))
coeftest(reg16)
# pas le bon modèle

reg17 = Arima(Petrole, order = c(5,1,0))
coeftest(reg17)
# pas le bon modèle

# On conserve reg1, reg2, reg4, reg5, reg6, reg7, reg9, reg10, reg11, reg12
# reg13, reg14, reg15


#### t.test (pvalue > 5%) ####

t.test(reg1$residuals) 
t.test(reg2$residuals) 
t.test(reg4$residuals) 
t.test(reg5$residuals) 
t.test(reg6$residuals) 
t.test(reg7$residuals) 
t.test(reg9$residuals) 
t.test(reg10$residuals) 
t.test(reg11$residuals) 
t.test(reg11$residuals) 
t.test(reg12$residuals) 
t.test(reg13$residuals) 
t.test(reg14$residuals) 
t.test(reg15$residuals) 

# toutes les p-values sont > 5% donc on accepte H0
# donc espérance des aléas est nulle

#### ArchTest (pvalue > 5%) ####

ArchTest(reg1$residuals, lag = 1)
ArchTest(reg2$residuals, lag = 1)
ArchTest(reg4$residuals, lag = 1)
ArchTest(reg5$residuals, lag = 1)
ArchTest(reg6$residuals, lag = 1)
ArchTest(reg7$residuals, lag = 1)
ArchTest(reg9$residuals, lag = 1)
ArchTest(reg10$residuals, lag = 1)
ArchTest(reg11$residuals, lag = 1)
ArchTest(reg12$residuals, lag = 1)
ArchTest(reg13$residuals, lag = 1)
ArchTest(reg14$residuals, lag = 1)
ArchTest(reg15$residuals, lag = 1)

# toutes les p-values sont > 5% donc on accepte H0
# donc absence effets arch <=> variance constante des aléas

#### Box.test (pvalue < 5%) ####

Box.test(reg1$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg2$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg4$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg5$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg6$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg7$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg9$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg10$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg11$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg12$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg13$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg14$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg15$residuals,  lag = 1, type = "Ljung-Box")

# toutes les p-values sont > 5% donc on accepte H0
# donc absence autocorrélation des aléas


#### Conclusion ####

# Tous les modèles précédents ont leurs aléas qui sont des bruits blancs
# On calcule le BIC de ces modèles 

#### BIC ####

BIC(reg_Petrole)
BIC(reg1)
BIC(reg2)
BIC(reg4)
BIC(reg5)
BIC(reg6)
BIC(reg7)
BIC(reg9)
BIC(reg10)
BIC(reg11)
BIC(reg12)
BIC(reg13)
BIC(reg14)
BIC(reg15)

# Le meilleur modèle est celui dont le BIC est le plus faible : BIC(reg15)
# On garde donc : 
# Meilleur modèle : ARIMA(3,1,2)


###################################################################
# Prévisions
###################################################################

# On passe directement à la prévision 
# Car on a pas de saisonnalité à modéliser ici
# comme les données sont annuelles <=> s = 0

# La prévision c'est toujours sur la série originale
# h = 4 car prévisions de 2022 à 2025

# On fait avec le meilleur modèle
autoplot(forecast(reg15, h = 4)) +
  labs(title = "Prévisions pour la Thailande de 2022 à 2025",
       x = "Années",
       y = "Production de pétrole non raffiné en milliers de tonnes d'équivalent pétrole")

predictions <- forecast(reg15, h = 4)

# Affichage des valeurs numériques des prévisions
print(predictions$mean)

######################################################################
# Les différences entre les valeurs réelles passées et les prévisions
######################################################################

# Selon les valeurs réelles on trouve
# 85.55 millions de tonnes d'équivalent pétrole en 2022
# donc 8 555 milliers de tonnes
# Donc on a 8 555 en 2022

# 88.25 millions de tep en 2023
# donc on a 8 825 milliers de tonnes
# Donc on a 8 825 en 2023

# C'est relativement proche des prévisions qu'on a obtenues:
# 8786.733 en 2022
# 7873.134 en 2023


######################################################################
# Même travail mais avec la donnée de 2022 ajoutée à la base initiale
######################################################################

# source: https://www.statista.com/statistics/1304238/thailand-production-volume-of-crude-oil/

# Ajout de la valeur 8555 pour l'année 2022
Petrole <- c(Petrole, 8555)
Petrole

plot(Petrole,type='l',col=4)

dPetrole = diff(Petrole)

# On vérifie que la tendance a été enlevée
plot(dPetrole,type='l',col=4)

acf(dPetrole)
# autocorrélation ordre 3 et 15

pacf(dPetrole)
# autocorrélation ordre 3

eacf(dPetrole)

# cela nous donne une indication des valeurs de p et q
# p = 0 , q = 3
# On utilise  ARIMA et pas ARMA car le PGD de la série Petrole est DS
# d = 1 car on différencie une fois ici

# estimation du modèle

reg_Petrole = Arima(Petrole, order = c(0,1,3))
coeftest(reg_Petrole) 
reg_Petrole = Arima(Petrole, order = c(0,1,3), fixed=c(0,0,NA))
coeftest(reg_Petrole)
# Les coefficients sont tous significatifs désormais


t.test(reg_Petrole$residuals) 
# p-value = 0,1964 > 0,05
# On accepte H0
# Donc la moyenne vaut zero

ArchTest(reg_Petrole$residuals, lag = 4)
# p-value = 0.05174 > 0,05
# On accepte H0
# Donc absence d'effets arch <=> Variance constante


Box.test(reg_Petrole$residuals,  lag = 4, type = "Ljung-Box")
# p-value = 0.9146 > 5% 
# On accepte H0 
# Donc pas autocorrélation des aléas.

# Donc dans ce cas les aléas sont des bruits blancs

# Coeftest

reg0 = Arima(Petrole, order = c(0,1,1))
coeftest(reg0)
# pas le bon modèle 
# coeff de ma1 non significatif

reg1 = Arima(Petrole, order = c(0,1,2))
coeftest(reg1)
# pas le bon modèle 

reg2 = Arima(Petrole, order = c(0,1,3))
coeftest(reg2)
reg2 = Arima(Petrole, order = c(0,1,3), fixed=c(0,0,NA))
coeftest(reg2)
# tous les coeffs sont significatifs


reg3 = Arima(Petrole, order = c(0,1,4))
coeftest(reg3)
# pas le bon modèle 
# on baisse q de 4 à 3 désormais

reg4 = Arima(Petrole, order = c(0,1,7))
coeftest(reg4)
# pas le bon modèle 


reg5 = Arima(Petrole, order = c(1,1,1))
coeftest(reg5)
# tous les coeffs sont significatifs


reg6 = Arima(Petrole, order = c(1,1,2))
coeftest(reg6)
# pas le bon modèle


reg7 = Arima(Petrole, order = c(1,1,3))
coeftest(reg7)
reg7 = Arima(Petrole, order = c(1,1,3), fixed=c(0,0,0,NA))
coeftest(reg7)
# tous les coeffs sont significatifs


reg8 = Arima(Petrole, order = c(1,1,4))
coeftest(reg8)
reg8 = Arima(Petrole, order = c(1,1,4), fixed=c(NA,NA,0,NA,NA))
coeftest(reg8)
# tous les coeffs sont significatifs

reg9 = Arima(Petrole, order = c(2,1,1))
coeftest(reg9)
reg9 = Arima(Petrole, order = c(2,1,1), fixed=c(NA,0,NA))
coeftest(reg9)
# tous les coeffs sont significatifs


reg10 = Arima(Petrole, order = c(2,1,2))
coeftest(reg10)
# tous les coeffs sont significatifs


reg11 = Arima(Petrole, order = c(2,1,3))
coeftest(reg11)
reg11 = Arima(Petrole, order = c(2,1,3), fixed=c(0,0,0,0,NA))
coeftest(reg11)
# tous les coeffs sont significatifs

reg12 = Arima(Petrole, order = c(3,1,1))
coeftest(reg12)
# On baisse q de 1 à 0 car le coeff de ma1 non significatif
reg12 = Arima(Petrole, order = c(3,1,0))
coeftest(reg12)
reg12 = Arima(Petrole, order = c(3,1,0), fixed=c(0,0,NA))
coeftest(reg12)
# tous les coeffs sont significatifs


reg13 = Arima(Petrole, order = c(1,1,0))
coeftest(reg13)
# pas le bon modèle


reg14 = Arima(Petrole, order = c(2,1,0))
coeftest(reg14)
reg14 = Arima(Petrole, order = c(2,1,0), fixed=c(0,NA))
coeftest(reg14)
# pas le bon modèle


reg15= Arima(Petrole, order = c(3,1,2))
coeftest(reg15)
# tous les coeffs sont significatifs


reg16 = Arima(Petrole, order = c(4,1,1))
coeftest(reg16)
# pas le bon modèle

reg17 = Arima(Petrole, order = c(5,1,0))
coeftest(reg17)
# pas le bon modèle

# On conserve reg2, reg5 , reg7, reg8, reg9, reg10, reg11, reg12, reg15

# t test

t.test(reg2$residuals) 
t.test(reg5$residuals) 
t.test(reg7$residuals) 
t.test(reg8$residuals) 
t.test(reg9$residuals) 
t.test(reg10$residuals) 
t.test(reg11$residuals) 
t.test(reg11$residuals) 
t.test(reg12$residuals) 
t.test(reg15$residuals) 

# toutes les p-values sont > 5% donc on accepte H0
# donc espérance des aléas est nulle

# Archtest

ArchTest(reg2$residuals, lag = 1)
ArchTest(reg5$residuals, lag = 1) 
ArchTest(reg7$residuals, lag = 1)
ArchTest(reg8$residuals, lag = 1)
ArchTest(reg9$residuals, lag = 1) 
ArchTest(reg10$residuals, lag = 1)
ArchTest(reg11$residuals, lag = 1)
ArchTest(reg12$residuals, lag = 1) 
ArchTest(reg15$residuals, lag = 1)

# toutes les p-values sont > 5% donc on accepte H0
# donc absence effets arch <=> variance constante des aléas

# Box test

Box.test(reg2$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg5$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg7$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg8$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg9$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg10$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg11$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg12$residuals,  lag = 1, type = "Ljung-Box")
Box.test(reg15$residuals,  lag = 1, type = "Ljung-Box")

# toutes les p-values sont > 5% donc on accepte H0
# donc absence autocorrélation des aléas


# Conclusion

# Tous les modèles précédents ont leurs aléas qui sont des bruits blancs
# On calcule le BIC 

# BIC

BIC(reg_Petrole)
BIC(reg2)
BIC(reg5)
BIC(reg7)
BIC(reg8)
BIC(reg9)
BIC(reg10)
BIC(reg11)
BIC(reg12)
BIC(reg15)

# Meilleur modèle : ARIMA(3,1,0) 

# prévisions

autoplot(forecast(reg12, h = 3)) +
  labs(title = "Prévisions pour la Thailande de 2023 à 2025",
       x = "Années",
       y = "Production de pétrole non raffiné en milliers de tonnes d'équivalent pétrole")

predictions <- forecast(reg12, h = 3)

# Affichage des valeurs numériques des prévisions
print(predictions$mean)