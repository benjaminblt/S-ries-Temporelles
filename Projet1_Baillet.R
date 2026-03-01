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

###################################################################
# Courbe tendance
###################################################################
# Tracer la courbe de tendance
ggplot(thailande, aes(x=TIME, y=Value)) +
  geom_line() +
  geom_smooth(method="lm", col="red", se=FALSE) +  # Courbe de régression (tendance)
  labs(title="Évolution des données de pétrole en Thaïlande (1971-2021)",
       x="Année",
       y="Valeur")

###################################################################
# Décomposition multiplicatve
###################################################################

# On récupère la colonne des valeurs
petrole<-ts(na.omit(thailande[,"Value"]),start=1971,freq=2, end= 2021)
decompo<-decompose(petrole,"multiplicative", filter = NULL)
plot(decompo)

###################################################################
# ACF, PACF
###################################################################

# On regarde s'il y a de l'autocorrélation des aléas
Petrole<-ts(na.omit(thailande[,"Value"]),start=1971,freq=1)
TPetrole<- length(Petrole)

graphe <- par (mfrow=c(1,2)) 
Acf(Petrole, lag=80)
Pacf(Petrole) 
par (graphe) 

###################################################################
# TEST 1: DICKEY FULLER (DF)
###################################################################

#### Etape 1: TREND ####

summary(ur.df(Petrole,type="trend",lag=0))

#### Etape 2: DRIFT ####

summary(ur.df(Petrole,type="drift",lag=0))

#### Etape 3: NONE ####

summary(ur.df(Petrole,type="none",lag=0))

# On regarde s'il y a de l'autocorrélation
plot(ur.df(Petrole,type="trend",lag=0))

###################################################################
# TEST 2: DICKEY FULLER AUGMENTE (ADF)
###################################################################

# On réalise une ADF en prenant en compte cette autocorrélation

#### Etape 1: Trend ####

# Formule de Schwert:
pmax <- as.integer(12*(TPetrole/100)^(0.25))

summary(CADFtest(Petrole,criterion="MAIC",type="trend",max.lag.y=pmax))
# Max lag of the diff. dependent variable vaut 3 donc on garde cette valeur
# au lieu de BIC
pmax <- 3

#### Etape 2: Drift ####

summary(ur.df(Petrole, type = "drift", lags = pmax))
# En effet la pvalue de z.lag1 est < 0,05
# DONC c'est le bon test

#### Etape 3: NONE ####

summary(ur.df(Petrole, type = "none", lags = pmax))

###################################################################
#  TEST 3: Zivot et Andrews avec 1 date de rupture (ZA)
###################################################################

# On veut trouver la valeur du lag
pmax <- as.integer(12*(TPetrole/100)^(0.25))

summary(ur.za(Petrole, model="both",lag=pmax))
# abs(tvalue(gamma_10))=2,074 > 1,6
# abs(tvalue(gamma_9))=1,406 < 1,6
# Donc on enlève 1 lag

summary(ur.za(Petrole, model="both",lag=pmax))

# pvalue(delta2)= 0.154774 > 0,05
# On passe à la spécification "intercept"

summary(ur.za(Petrole, model="intercept",lag=pmax))

# pvalue (delta1) = 0.0021 < 0,05

plot(ur.za(Petrole, model="intercept",lag=pmax))

###################################################################
#  TEST 4:  Lee et Strazicich avec Boolstrap et 1 date de rupture (LS)
###################################################################

# faire "Code" puis "Source File" puis emplacement de 
# "LeeStrazicichUnitRootTestParallelization.R"

#### Bootstrap 1 date de rupture

myBreaks <- 1
# Assumed break in the series, "crash" - break in intercept; "break" - break in intercept and trend
myModel <- "crash" # car intercept dans ZA
# Number of lags to be used in fixed specification or maximum number of lags, when using the GTOS method
myLags <- 5

#Define number of cores to use. By default the maximum available number minus one core is used
cl <- makeCluster(max(1, detectCores() - 1))
registerDoSNOW(cl)

myParallel_LS <- ur.ls.bootstrap(y= Petrole , model = myModel, breaks = myBreaks, lags = myLags, method = "Fixed",pn = 0.1, critval = "bootstrap", print.results = "print")

#### Bootstrap 2 dates de rupture

# car 2 dates rupture
myBreaks <- 2
# Assumed break in the series, "crash" - break in intercept; "break" - break in intercept and trend
myModel <- "crash" # car intercept dans ZA
# Number of lags to be used in fixed specification or maximum number of lags, when using the GTOS method
myLags <- 5

#Define number of cores to use. By default the maximum available number minus one core is used
cl <- makeCluster(max(1, detectCores() - 1))
registerDoSNOW(cl)

myParallel_LS <- ur.ls.bootstrap(y= Petrole , model = myModel, breaks = myBreaks, lags = myLags, method = "Fixed",pn = 0.1, critval = "bootstrap", print.results = "print")

###################################################################
#  On différencie à l'ordre 1 et on refait les 2 premiers tests
###################################################################

dPetrole=diff(Petrole)

# TEST 1: DICKEY FULLER (DF)

#### Etape 1: TREND ####

summary(ur.df(dPetrole,type="trend",lag=0))

#### Etape 2: DRIFT ####

summary(ur.df(dPetrole,type="drift",lag=0))

#### Etape 3: NONE ####

summary(ur.df(dPetrole,type="none",lag=0))


# On regarde s'il y a de l'autocorrélation
plot(ur.df(dPetrole,type="trend",lag=0))

###################################################################
# TEST 2: DICKEY FULLER AUGMENTE (ADF)
###################################################################

# Donc on réalise une ADF en prenant en compte cette autocorrélation

#### Etape 1: Trend ####

TdPetrole<- length(dPetrole)

# Formule de Schwert:
pmax <- as.integer(12*(TdPetrole/100)^(0.25))
summary(CADFtest(dPetrole,criterion="MAIC",type="trend",max.lag.y=pmax))
# Max lag of the diff. dependent variable vaut 3 donc on garde cette valeur
# au lieu de BIC
pmax <- 3

summary(ur.df(dPetrole, type = "trend", lags = pmax))

#### Etape 2: Drift ####

summary(ur.df(dPetrole, type = "drift", lags = pmax))

#### Etape 3: NONE ####

summary(ur.df(dPetrole, type = "none", lags = pmax))