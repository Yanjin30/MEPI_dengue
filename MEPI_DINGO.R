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
devtools::install_github("GaelBn/BRREWABC")
library(BRREWABC)

# Import datasets
data19 <- read.csv(file = "dengue_2019.csv", header = TRUE, sep = ",", dec = ".",)
data20 <- read.csv(file = "dengue_2020.csv", header = TRUE, sep = ",", dec = ".")
data21 <- read.csv(file = "dengue_2021.csv", header = TRUE, sep = ",", dec = ".")

modele_dengue=function(Sh0,Ih0,Rh0,Sv0,Iv0,param,tmax){
  Sh=Sh0
  Ih=Ih0
  Rh=Rh0
  Sv=2*Sh0 #Nombre de moustique par habitant à définir
  Iv=0.2*Sv #Proportion de moustique porteur à définir
  betah=param[1]
  gamma=param[2]
  betav=param[3]
  for (i in 1:tmax){
    
  }
  
}