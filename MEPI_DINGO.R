################################
#~~~HARDY Erwan et FRITZ Leo~~~#
################################

################################################################################
# Script pour l'étude d'une épidemie de dengue au Sri Lanka entre 2019 et 2022 #
# Prise en compte des migrations inter-urbaines et de la pluviométrie à        #
# l'échelle des districts administratifs.                                      #
################################################################################
#setwd("~/Documents/Master/M2/MEPI/projet/MEPI_dengue/dataset_dengo")
setwd("C:/Users/Nitro/Documents/Cours/MEPI_dengue/dataset_dengo")
# Import 
#devtools::install_github("GaelBn/BRREWABC")
library(BRREWABC)

# Import datasets
data19 <- read.csv(file = "dengue_2019.csv", header = TRUE, sep = ",", dec = ".",)
data20 <- read.csv(file = "dengue_2020.csv", header = TRUE, d = ",", dec = ".")
data21 <- read.csv(file = "dengue_2021.csv", header = TRUE, sep = ",", dec = ".")
############# Déterministe ###############

modele_dengue_deter=function(y,t,param){
  # Définition des paramètres à opimiser
  beta_h = param[1]
  
  # Définition des paramètres fixés
  z = 2.5
  gamma = 1/2
  beta_v = 0.375
  mu_v = 1/6
  
  # Définition des classes de population
  # Humains
  Ih = y[1]
  Rh = y[2]
  
  Nh = Sh+Ih+Rh
  
  # Moustiques
  Iv = y[3]
  Nv= z * Nh
  
  # Définition des dérivées
  dIdt = - gamma * Ih + beta_h * z * (Iv/Nv) * (Nh - Ih -Rh)/Nh
  dRdt = gamma * Ih
  dVdt = beta_v * (Nv - Iv)  - mu_v * V
  return(c(dIdt,dRdt,dVdt))
}

############# Stochastique ###############
modele_dengue_stoch=function(Sh0,Ih0,Rh0,Sv0,Iv0,param,tmax){
  Sh=Sh0
  Ih=Ih0
  Rh=Rh0
  Nv=2*Sh0 #Nombre de moustique par habitant à définir
  Iv=0.2*Nv #Proportion de moustique porteur à définir
  Sv=Nv-Iv
  betah=param[1]
  gamma=param[2]
  betav=param[3]
  for (i in 1:tmax){
    # Définition des taux de transmission
    tx_transmi_h=
    tx_transmi_v=
    tx_recovery=
    # Définition des probas de transmissions 
    pb_inf_h=1-exp(-tx_transmi_h)
    pb_inf_v=1-exp(-tx_transmi_v)
    pb_recov=1-exp(-tx_recovery)
    # Tirage des individus transitionnant
    new_inf_v=rbinom(Sh,pb_inf_h)
    new_inf_h=rbinom(Sv,pb_inf_v)
    new_recov=rbinom(Ih,pb_recov)
    # Ajout des individus transitionnant
    Sh=Sh-new_inf_h
    Ih=Ih+new_inf_h-new_recov
    Rh=Rh+new_recov
    Sv=Sv-new_inf_v
    Iv=Iv+new_inf_v
    
  }
  
}
