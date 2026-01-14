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
library(deSolve)
library(tidyr)
library(dplyr)
library(ggplot2)

# Import datasets
data19 = read.csv(
  file = "dengue_2019.csv",
  header = TRUE,
  sep = ",",
  dec = ".",
)

data20 = read.csv(
  file = "dengue_2020.csv",
  header = TRUE,
  sep = ",",
  dec = "."
)

data21 = read.csv(
  file = "dengue_2021.csv",
  header = TRUE,
  sep = ",",
  dec = "."
)

# Create a single reported cases table
data19$Year <- 2019
data20$Year <- 2020
data21$Year <- 2021
data_all <- rbind(data19, data20, data21)

# Tidy and mutate the binded dataset
reported_cases <- data_all %>%
  pivot_longer(
    cols = Jan:Dec,
    names_to = "Month",
    values_to = "Cases") %>%
  mutate(
    Month = match(Month,
                  c("Jan","Feb","Mar","Apr","May","June",
                    "July","Aug","Sept","Oct","Nov","Dec")),
    Date = as.Date(paste(Year, Month, 30, sep = "-")))

# Check for consistency
head(reported_cases)

#########
#   #   #
#  ###  # Attention à la signification de la colonne "Total" !!!
# ##### #
#########

# Compute a dataset with Sri Lankan monthly cases
Data_monthly = reported_cases %>%
  group_by(Date, Province, District)%>%
  summarise(Cases = sum(Cases))
Data_monthly$Province = as.factor(Data_monthly$Province)

SriLankan_monthly = Data_monthly %>%
  filter(Province=="TOTAL")%>%
  select(Date,Cases)

# Import weather data
weather = read.csv(
  file = "SriLanka_Weather_Dataset_2019_2021.csv",
  header = TRUE,
  sep = ",",
  dec = ".")

################ Data plotting to see the trends ############
ggplot(data = SriLankan_monthly, aes(x=Date, y=Cases))+
  geom_area(fill='orchid3', col='black',alpha=0.4)+
  labs(title =  "Sri Lankan pooled data (no district considered)", y = "Reported cases")+
  theme_minimal()

ggplot(data = Data_monthly, aes(x=Date, y=Cases))+
  geom_area(aes(fill=Province), col='black',alpha=0.3)+
  labs(title =  "Reported case per province", y = "Reported cases")+
  theme_minimal()

################ Deterministic model ########################
z_t <- function(Z,t){
  period=seq(9,162,9)
  for (j in 1:length(period)){
    if (t<period[j]){
      print(j)
      zt=Z[j]
      break
    }
  }
  #zt = runif(n = 1, min = 0, max = 1.5)
  return(zt)
}

# Modèle déterministe à base d'ODE
modele_dengue_deter = function(t, y, param) {
  # Définition des paramètres à opimiser
  #beta_h = param[1]
  Z = param
  
  # Définition des paramètres fixés
  z = z_t(Z,t)
  gamma = 1 / 2
  beta_v = 0.375
  mu_v = 1 / 6
  mu_h = 0
  beta_h = 0.75
  
  # Définition des classes de population
  # Humains
  Ih = y[1]
  Rh = y[2]
  
  Nh = 10000
  
  # Moustiques
  Iv = y[3]
  Nv = z * Nh
  
  # Définition des dérivées
  dIdt = -(gamma + mu_h) * Ih + beta_h * z * Iv  * (Nh - Ih - Rh) / Nh
  dRdt = gamma * Ih - mu_h * Rh
  dVdt = beta_v * (Nv - Iv) * Ih / Nh - mu_v * Iv
  
  return(list(c(dIdt, dRdt, dVdt)))
}

## Fonction qui réalise une simulation à partir d'un modèle (ODE)
simulation_deter = function(y, tmax, param, delta_t) {
  # Definition du pas de temps
  temps = seq(0, tmax, delta_t)
  
  # Résolution de l'équation
  result = ode(
    y = y,
    times = temps,
    func = modele_dengue_deter,
    parms = param,
    method = "rk4"
  )
  
  # Visualisation
  par(mfrow = c(1, 2))
  plot(
    result[, 1],
    result[, 3],
    type = "l",
    col = "blue",
    ylim = c(0, 10000)
  )
  lines(result[, 1], result[, 2], type = "l", col = "red")
  plot(result[, 1], result[, 4], type = "l")
  return(result)
}

Z=seq(0.1,1.8,0.1)
test = simulation_deter(c(100, 100, 75), 156, Z, 1)

################# 1st fit ####################################
# L'objectif ici est de fit z (le ratio du nombre de moustiques par rapport au nombre d'humains)
# pour l'ensemble des districts aggregés et pour plusieurs et sur plusieurs intervalles de temps 
# (on suppose que z change toutes les 2, 4, 8, 16 semaines) et on compare avec les AICs 
# du modèle fitté : on justfie comme ça le choix de la periodicité 
# (sans Fast Fourier Transform, même si en vrai c'est pas très compliqué)







## Fonction qui prend une dynamique d'infection et de recovery afin d'avoir le nombre d'incidence par mois
summary_extract = function(vect_inf, vect_recov) {
  monthly_new_case = c()
  u = 1
  for (i in seq(30, length(vect_inf), 30)) {
    if (u == 1) {
      monthly_new_case = c(monthly_new_case, sum(diff(c(0, vect_inf[u:i])) + diff(c(
        0, vect_recov[u:i]
      ))))
    }
    if (u != 1) {
      print(u)
      monthly_new_case = c(monthly_new_case, sum(diff(vect_inf[(u - 1):i]) +
                                                   diff(vect_recov[(u - 1):i])))
    }
    u = i + 1
  }
  return(monthly_new_case)
}


## Fonction qui calcule la distance à partir d'une trajectoire et de nos statistiques résumées
## Prend en entrée les statistiques résumées et un vecteur de paramètre
distance_deter = function(param, ssobs) {
  # Définition des conditions initiales
  y0 = c(100, 0, 100)
  
  # Simulation de la trajectoire
  simu = simulation_deter(
    y = y0,
    tmax = 360,
    param = param,
    delta_t = 1
  ) # On considère grossièrement que chaque mois -> 30 jours
  infected_dyna = simu[, 2]
  recover_dyna = simu[, 3]
  
  # Création d'un résumé de nos statistiques simulées
  all_mensual_case = summary_extract(infected_dyna, recover_dyna) # On récupère toutes les nouvelles infection de chaque mois
  all_mensual_case_obs = rbinom(all_mensual_case, 0.05) # Modèle d'observation, hypothèse symptome
  
  # Comparaison de nos statistiques résumées
  dist = sum((all_mensual_case - ssobs) ** 2)
  return(c(dist))
}

############# Stochastique ###############
modele_dengue_stoch = function(Sh0, Ih0, Rh0, Sv0, Iv0, param, tmax) {
  Sh = Sh0
  Ih = Ih0
  Rh = Rh0
  Nv = 2 * Sh0 #Nombre de moustique par habitant à définir
  Iv = 0.2 * Nv #Proportion de moustique porteur à définir
  Sv = Nv - Iv
  betah = param[1]
  gamma = param[2]
  betav = param[3]
  for (i in 1:tmax) {
    # Définition des taux de transmission
    tx_transmi_h =
      tx_transmi_v =
      tx_recovery =
      # Définition des probas de transmissions
      pb_inf_h = 1 - exp(-tx_transmi_h)
    pb_inf_v = 1 - exp(-tx_transmi_v)
    pb_recov = 1 - exp(-tx_recovery)
    # Tirage des individus transitionnant
    new_inf_v = rbinom(Sh, pb_inf_h)
    new_inf_h = rbinom(Sv, pb_inf_v)
    new_recov = rbinom(Ih, pb_recov)
    # Ajout des individus transitionnant
    Sh = Sh - new_inf_h
    Ih = Ih + new_inf_h - new_recov
    Rh = Rh + new_recov
    Sv = Sv - new_inf_v
    Iv = Iv + new_inf_v
  }
}

