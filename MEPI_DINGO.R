################################
#~~~HARDY Erwan et FRITZ Leo~~~#
################################

################################################################################
# Script pour l'étude d'une épidemie de dengue au Sri Lanka entre 2019 et 2022 #
# Prise en compte des migrations inter-urbaines et de la pluviométrie à        #
# l'échelle des districts administratifs.                                      #
################################################################################
setwd("~/Documents/Master/M2/MEPI/projet/MEPI_dengue/datasets_dengo")

# Import 
# devtools::install_github("GaelBn/BRREWABC")
library(BRREWABC)

# Import datasets
# Rhaaa ça marche pô
data19 <- read.csv(file = "dengue_2019.csv", header = TRUE, sep = ",", dec = ".",)
data20 <- read.csv(file = "dengue_2020.csv", header = TRUE, sep = ",", dec = ".")
data21 <- read.csv(file = "dengue_2021.csv", header = TRUE, sep = ",", dec = ".")
