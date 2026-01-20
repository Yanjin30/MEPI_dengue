############Authors#############
#~~~HARDY Erwan et FRITZ Leo~~~#
################################

################################################################################
# Script pour l'étude d'une épidemie de dengue au Sri Lanka entre 2019 et 2022 #
# Prise en compte des migrations inter-urbaines et de la pluviométrie à        #
# l'échelle des districts administratifs.                                      #
################################################################################

setwd("~/Documents/Master/M2/MEPI/projet/MEPI_dengue/dataset_dengo")
#setwd("C:/Users/Nitro/Documents/Cours/MEPI_dengue/dataset_dengo")

# Import
library(BRREWABC)
library(deSolve)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)

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
  pivot_longer(cols = Jan:Dec,
               names_to = "Month",
               values_to = "Cases") %>%
  mutate(Month = match(
    Month,
    c(
      "Jan",
      "Feb",
      "Mar",
      "Apr",
      "May",
      "June",
      "July",
      "Aug",
      "Sept",
      "Oct",
      "Nov",
      "Dec"
    )
  ), Date = as.Date(paste(Year, Month, 30, sep = "-"))) %>%
  filter(Province != "TOTAL")

# Check for consistency
head(reported_cases)

# Compute a dataset with Sri Lankan monthly cases
Data_monthly = reported_cases %>%
  group_by(Date, Province, District) %>%
  summarise(Cases = sum(Cases))
Data_monthly$Province = as.factor(Data_monthly$Province)

SriLankan_monthly = Data_monthly %>%
  group_by(Date) %>%
  summarise(Cases = sum(Cases))

# Import weather data
weather = read.csv(
  file = "SriLanka_Weather_Dataset_2019_2021.csv",
  header = TRUE,
  sep = ",",
  dec = "."
)

################ Data plotting to see the trends ############
ggplot(data = SriLankan_monthly, aes(x = Date, y = Cases)) +
  geom_area(fill = 'orchid3',
            col = 'black',
            alpha = 0.4) +
  labs(title =  "Sri Lankan pooled data (no district considered)", y = "Reported cases") +
  theme_minimal()

ggplot(data = Data_monthly, aes(x = Date, y = Cases)) +
  geom_area(aes(fill = Province), col = 'black', alpha = 0.3) +
  labs(title =  "Reported case per province", y = "Reported cases") +
  theme_minimal()

################ Deterministic model ########################
z_t <- function(Z, t) {
  period = seq(9, 162, 9)
  for (j in 1:length(period)) {
    if (t < period[j]) {
      print(j)
      zt = Z[j]
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
  z = z_t(Z, t)
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
solve_determistic = function(y, tmax, param, delta_t) {
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

Z = seq(0.1, 1.8, 0.1)
test = solve_determistic(c(100, 100, 75), 156, Z, 1)

################# 1st fit ################
# L'objectif ici est de fit z (le ratio du nombre de moustiques par rapport au
# nombre d'humains) pour l'ensemble des districts aggregés et pour plusieurs et
# sur plusieurs intervalles de temps (on suppose que z change toutes les 2, 4, 8
# , 16 semaines)et on compare avec les AICs du modèle fitté : on justfie comme 
# ça le choix de la periodicité (sans Fast Fourier Transform, même si en vrai 
# c'est pas très compliqué)

########### Climate data per period ##### 
colnames(weather)[1] = "dates" 
weather = separate(weather, dates, c('year', 'month', 'day'), sep = "-",remove = FALSE)

nb_month_in_period = 2 # Period = 2 months
weather_periodized = weather %>%
  mutate(
    month = as.integer(month),
    period = ceiling(month / nb_month_in_period)) %>%
  group_by(year, period) %>%
  summarise(
    minTemp = min(temperature_2m_min, na.rm = TRUE),
    maxTemp = max(temperature_2m_max, na.rm = TRUE),
    meanTemp = mean(temperature_2m_mean, na.rm = TRUE),
    meanPrecip = mean(precipitation_sum, na.rm = TRUE),
    cumulPrecip = sum(precipitation_sum, na.rm = TRUE),
    meanEvap = mean(et0_fao_evapotranspiration, na.rm = TRUE),
    cumulEvap = sum(et0_fao_evapotranspiration, na.rm = TRUE),
    meanMaxWind = mean(windspeed_10m_max, na.rm = TRUE),
    meanGusts = mean(windgusts_10m_max, na.rm = TRUE))

########## z VS climate vars ############
z_patricules_file_path = "~/Documents/Master/M2/MEPI/all_accepted_particles2.csv"
z_estim = read.csv(z_patricules_file_path)
# On recupere les particules (distributions des parametres z_t) pour la dernière
# génération de l'ABC SMC (normalement, ces distributions convergent vers les 
# distributions des paramètres)
last_gen = z_estim %>%
  filter(gen == max(gen))%>%
  select(-c(gen,model,dist1,pWeight))
# On peut recupérer brièvement la moyenne de chacun des z
summary(last_gen)

# We need to define a function to compute the max of the distribution
MAP = function(x){
  d = density(x)
  mode_value = d$x[which.max(d$y)]
  return(mode_value)
}

last_gen_mode = apply(last_gen, 2, MAP)
z_estim = as.data.frame(cbind(last_gen_mode,
                      weather_periodized$year,
                      weather_periodized$period))
colnames(z_estim) = c("Mode","Year","Period")

z_estim = z_estim %>%
  mutate(
    Month = (as.numeric(Period) - 1) * nb_month_in_period + ceiling(nb_month_in_period / 2),
    Date = as.Date(paste(Year, Month, "1", sep = "-"))
)

par(mfrow=c(1,1))
ggplot(data = z_estim, aes(x=Date,y=Mode))+
  geom_point()+
  geom_line()

plot(z_estim$Mode, ylab = "z", xlab = "Time")
lines(z_estim$Mode, col="darkorchid")

model = lm(z_estim$Mode~weather_periodized$meanPrecip)
anova(model)

