############Authors#############
#~~~HARDY Erwan et FRITZ Leo~~~#
################################

################################################################################
# Script pour l'étude d'une épidemie de dengue au Sri Lanka entre 2019 et 2022 #                                     #
################################################################################

# Erwan 
setwd("~/Documents/Master/M2/MEPI/projet/MEPI_dengue/dataset_dengo")
# Leo
#setwd("C:/Users/Nitro/Documents/Cours/MEPI_dengue/dataset_dengo")

path = "C:/Users/Nitro/Documents/Cours/Inf"

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
  dec = ".")

################ Data plotting to see the trends ################
ggplot(data = SriLankan_monthly, aes(x = Date, y = Cases)) +
  geom_area(fill = 'orchid3',
            col = 'black',
            alpha = 0.4) +
  labs(title =  "Sri Lankan pooled data (no district considered)", y = "Reported cases") +
  theme_minimal()

ggplot(data = SriLankan_monthly, aes(x = Date, y = Cases)) +
  geom_point(col = 'blue') +
  geom_line(col = 'black') +
  labs(title =  "Sri Lankan pooled data (no district considered)", y = "Reported cases") +
  theme_minimal()

################ Modèle avec beta_h variant ####################
model_determinist = function(t, y, beta_h) {
  # Définition des paramètres fixés
  gamma = 1 / 2
  beta_v = 0.375
  mu_v = 1 / 6
  mu_h = 0

  # Définition des classes de population - humains
  Ih = y[1]
  Rh = y[2]
  
  Nh = 10000
  
  # Moustiques
  Iv = y[3]
  Nv = Nh
  
  # Définition des dérivées
  dIdt = -(gamma + mu_h) * Ih + beta_h * Iv  * (Nh - Ih - Rh) / Nh
  dRdt = gamma * Ih - mu_h * Rh
  dVdt = beta_v * (Nv - Iv) * Ih / Nh - mu_v * Iv
  
  return(list(c(dIdt, dRdt, dVdt)))}

all_solutions <- list()
for (beta_h in seq(0, 1.5, 0.25)) {
  sol <- ode(
    y = c(Ih = 100, Rh = 0, Iv = 10),
    times = seq(0, 75, 0.25),
    func = model_determinist,
    parms = beta_h,
    method = "rk4"
  )
  
  sol <- as.data.frame(sol)
  sol$beta_h <- beta_h  # on garde la valeur de beta_h
  all_solutions[[as.character(beta_h)]] <- sol
  }

df <- bind_rows(all_solutions)
df$beta_h <- as.factor(df$beta_h)
df_long <- df %>%
  pivot_longer(
    cols = c(Ih, Rh, Iv),
    names_to = "compartment",
    values_to = "value")

ggplot(
  df_long %>% filter(compartment == "Ih"),
  aes(x = time, y = value, color = beta_h, group = beta_h)) +
  geom_line() +
  labs(
    title = "Trajectoires des I pour différentes valeurs de beta_h",
    x = "Temps",
    y = "Nombre d'humains infectés",
    color = "beta_h") +
  theme_minimal()

################ Deterministic model with varying z #######################
z_t <- function(Z, t) {
  period = seq(9, 162, 9)
  for (j in 1:length(period)){
    if (t < period[j]) {
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
  
  Nh = 23300000
  
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
solve_determistic = function(y, tmax, param, delta_t, plot = FALSE) {
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
  
  if(plot){# Visualisation
  par(mfrow = c(1, 2))
  plot(
    result[, 1],
    result[, 3],
    type = "l",
    col = "blue",
    ylim = c(0, max(result[,3],result[,2]))
  )
  lines(result[, 1], result[, 2], type = "l", col = "red")
  plot(result[, 1], result[, 4], type = "l")}
  return(result)
}

Z = runif(20,0.2,1.5)
test = solve_determistic(c(100, 100, 75), 156, Z, 1, plot = TRUE)

######################################################################
######################################################################
###################### Model deterministe ############################
######################################################################
######################################################################
# L'objectif ici est de fit z (le ratio du nombre de moustiques par rapport au
# nombre d'humains) pour l'ensemble des districts aggregés
distance_determinist = function(x, ssobs) {
  ### Définition de toutes les fonctions pour ABC
  library(deSolve)
  summary_extract = function(vect_inf, vect_recov) {
    monthly_new_case = c()
    u = 1
    for (i in seq(4, length(vect_inf), 4.6)) {
      if (u == 1) {
        monthly_new_case = c(monthly_new_case, sum(diff(c(0, vect_inf[u:i])) + diff(c(
          0, vect_recov[u:i]
        ))))
      }
      if (u != 1) {
        monthly_new_case = c(monthly_new_case, sum(diff(vect_inf[(u - 1):i]) +
                                                     diff(vect_recov[(u - 1):i])))
      }
      u = i + 1
    }
    return(monthly_new_case)
  }
  simulation_determinist = function(y, tmax, param, delta_t) {
    # Definition du pas de temps
    temps = seq(0, tmax, delta_t)
    
    # Résolution de l'équation
    result = ode(
      y = y,
      times = temps,
      func = modele_dengue_determinist,
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
  
  # Define z_t as a time-dependent parameter
  z_t <- function(Z,t){
    period=seq(9,162,9)
    for (j in 1:length(period)){
      if (t<period[j]){
        zt=Z[j]
        break
      }
    }
    return(zt)
  }
  
  modele_dengue_determinist = function(t, y, param) {
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
    
    Nh = 23300000
    
    # Moustiques
    Iv = y[3]
    Nv = z * Nh
    
    # Définition des dérivées
    dIdt = -(gamma + mu_h) * Ih + beta_h * z * Iv  * (Nh - Ih - Rh) / Nh
    dRdt = gamma * Ih - mu_h * Rh
    dVdt = beta_v * (Nv - Iv) * Ih / Nh - mu_v * Iv
    
    return(list(c(dIdt, dRdt, dVdt)))
  }

  ## Fonction de distance
  # Définition des conditions initiales
  y0 = c(10000, 100000, 0)
  
  param=c()
  
  for (i in 1:18){
    param=c(param,x[[paste0("z",i)]])
  }
  # Simulation de la trajectoire
  
  simu = simulation_determinist(
    y = y0,
    tmax = 156,
    param = param,
    delta_t = 1
  ) # On considère grossièrement que chaque mois -> 30 jours
  infected_dyna = simu[, 2]
  recover_dyna = simu[, 3]
  # Création d'un résumé de nos statistiques simulées
  all_mensual_case = summary_extract(infected_dyna, recover_dyna) # On récupère toutes les nouvelles infection de chaque mois
  all_mensual_case_obs = rbinom(n = length(all_mensual_case),  size = pmax(0,round(all_mensual_case)),  prob = 0.05) # Modèle d'observation, hypothèse symptome
  # Comparaison de nos statistiques résumées
  dist = sum((all_mensual_case_obs - ssobs) ** 2)
  return(c(dist))
}
## Inférence des valeurs de z
# Définition des priors
z_priors = list("m1"=lapply(1:18, function(i) c(paste0("z", i), "unif", 0, 3)))

# Définition du modèle
model_list = list("m1" = distance_determinist)

# Définition des statistiques observées
ss_obs = SriLankan_monthly$Cases

# Réalisation de l'inférence, hyper-espace en dimension 18 donc on considère une taille de génération très grande (3000 particules/gen)
res = abcsmc(model_list = model_list, prior_dist = z_priors,
             ss_obs = ss_obs, max_number_of_gen = 60, nb_acc_prtcl_per_gen = 3000,
             new_threshold_quantile = 0.8, max_attempts = 200000, experiment_folderpath = path,
             max_concurrent_jobs = 7, verbose = TRUE,progressbar = TRUE,acceptance_rate_min = 0.001)

## Exctraction des résultats
all_accepted_particles = res$particles
all_thresholds = res$thresholds

## Plot des différentes figures de l'inférence
plot_abcsmc_res(data = all_accepted_particles, prior = z_priors, colorpal = "YlOrBr", filename = file.path(path, "abcsmc_results.png"),iter=60)
plot_thresholds(data = all_thresholds, nb_threshold = 1, colorpal = "YlOrBr", filename = file.path(path, "thresholds.png"))
plot_ess(data = all_accepted_particles, colorpal = "YlOrBr", filename = file.path(path, "ess.png"))
plot_densityridges(data = all_accepted_particles, prior = z_priors, colorpal = "YlOrBr", filename = file.path(path, "densityridges.png"))

#### Trajectoire avec les particules de la génération 60, calcule de la moyenne, medianne et intervalle de confiance à 95%
y0 = c(10000, 100000, 0)
part_to_simulate=all_accepted_particles[all_accepted_particles$gen==60,seq(3,20,1)]
simulation_result_inf=matrix(nrow=length(part_to_simulate$z1),ncol=157)
simulation_result_recov=matrix(nrow=length(part_to_simulate$z1),ncol=157)
simulation_result_stik=matrix(nrow=length(part_to_simulate$z1),ncol=157)
for (line in 1:length(part_to_simulate$z1)){
  simulation=simulation_deter(y0, 156, part_to_simulate[line,], 1)
  simulation_result_inf[line,]=simulation[,2]
  simulation_result_recov[line,]=simulation[,3]
  simulation_result_stik[line,]=simulation[,4]
}

simulation_long <- as.data.frame(simulation_result_inf) %>%
  mutate(traj = row_number()) %>%
  pivot_longer(
    cols = -traj,
    values_to = "value"
  ) %>%
  group_by(traj) %>%
  mutate(
    time = row_number() - 1   # ou row_number() si 1..157
  ) %>%
  ungroup()

simulation_long = simulation_long %>%
  group_by(time) %>%
  summarise(
    mean   = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    low90  = quantile(value, 0.05, na.rm = TRUE),
    high90 = quantile(value, 0.95, na.rm = TRUE)
  )

### GGPLOT de la dynamique du modèle sans observation (100 % des cas du modèle)
ggplot(simulation_long, aes(x = time)) +
  geom_line(
    aes(y = mean, color = "Moyenne", linetype = "Moyenne"),
    linewidth = 1
  ) +
  geom_line(
    aes(y = median, color = "Médiane", linetype = "Médiane"),
    linewidth = 1
  ) +
  geom_line(
    aes(y = low90, color = "IC 90%", linetype = "IC 90%"),
    linewidth = 0.8
  ) +
  geom_line(
    aes(y = high90, color = "IC 90%", linetype = "IC 90%"),
    linewidth = 0.8
  ) +
  scale_color_manual(
    name = "Statistique",
    values = c("Moyenne" = "red", "Médiane" = "blue", "IC 90%" = "blue")
  ) +
  scale_linetype_manual(
    name = "Statistique",
    values = c("Moyenne" = "solid", "Médiane" = "solid", "IC 90%" = "dashed")
  ) +
  
  labs(
    x = "Temps (en semaine)",
    y = "Nombre d'infecté",
    title = "Dynamique d'infection chez l'homme"
  ) +
  coord_cartesian(ylim = c(0, 1e6))+
  theme_minimal() 
### A 5 % de population 
par(mfrow=c(1,2))
simulation_long <- as.data.frame(simulation_result_inf) %>%
  mutate(traj = row_number()) %>%
  pivot_longer(
    cols = -traj,
    values_to = "value"
  ) %>%
  group_by(traj) %>%
  mutate(
    time = row_number() - 1   # ou row_number() si 1..157
  ) %>%
  ungroup()

simulation_long = simulation_long %>%
  group_by(time) %>%
  summarise(
    mean   = 0.05 * mean(value, na.rm = TRUE),
    median = 0.05 * median(value, na.rm = TRUE),
    low90  = 0.05 * quantile(value, 0.05, na.rm = TRUE),
    high90 = 0.05 * quantile(value, 0.95, na.rm = TRUE)
  )

### GGPLOT à la semaine des simulations à 5% de cas observées 
ggplot(simulation_long, aes(x = time)) +
  geom_line(
    aes(y = mean, color = "Moyenne", linetype = "Moyenne"),
    linewidth = 1
  ) +
  geom_line(
    aes(y = median, color = "Médiane", linetype = "Médiane"),
    linewidth = 1
  ) +
  geom_line(
    aes(y = low90, color = "IC 90%", linetype = "IC 90%"),
    linewidth = 0.8
  ) +
  geom_line(
    aes(y = high90, color = "IC 90%", linetype = "IC 90%"),
    linewidth = 0.8
  ) +
  scale_color_manual(
    name = "Statistique",
    values = c("Moyenne" = "red", "Médiane" = "blue", "IC 90%" = "blue")
  ) +
  scale_linetype_manual(
    name = "Statistique",
    values = c("Moyenne" = "solid", "Médiane" = "solid", "IC 90%" = "dashed")
  ) +
  
  labs(
    x = "Temps (en semaine)",
    y = "Nombre d'infecté",
    title = "Dynamique d'infection chez l'homme, 5% des cas observés"
  ) +
  coord_cartesian(ylim = c(0, 1e6 * 0.05))+
  theme_minimal() 

## Calcul des cas mensuels d'infection
monthly_simu = matrix(nrow = nrow(simulation_result_inf), ncol = 34)

results_list = list()

for (line in 1:nrow(simulation_result_inf)) {
  result = summary_extract(
  vect_inf = simulation_result_inf[line, ],
  vect_recov = simulation_result_recov[line, ]
  )
  results_list[[line]] = result
}

# Convertir en matrice
weekly_simu = do.call(rbind, results_list)

# Vérifier les dimensions
dim(weekly_simu)
head(weekly_simu)

## Calcules de statistiques à représenter
simulation_long <- as.data.frame(weekly_simu) %>%
  mutate(traj = row_number()) %>%
  pivot_longer(
    cols = -traj,
    values_to = "value"
  ) %>%
  group_by(traj) %>%
  mutate(
    time = row_number() - 1   # ou row_number() si 1..157
  ) %>%
  ungroup()

simulation_long = simulation_long %>%
  group_by(time) %>%
  summarise(
    mean   = 0.05 * mean(value, na.rm = TRUE),
    median = 0.05 * median(value, na.rm = TRUE),
    low90  = 0.05 * quantile(value, 0.05, na.rm = TRUE),
    high90 = 0.05 * quantile(value, 0.95, na.rm = TRUE)
  )
### GGPLOT des simulations mensuels à 5%  
SriLankan_monthly$time=seq(0,33,1)
ggplot(simulation_long, aes(x = time)) +
  geom_line(
    aes(y = mean, color = "Moyenne", linetype = "Moyenne"),
    linewidth = 1
  ) +
  geom_line(
    aes(y = median, color = "Médiane", linetype = "Médiane"),
    linewidth = 1
  ) +
  geom_line(
    aes(y = low90, color = "IC 90%", linetype = "IC 90%"),
    linewidth = 0.8
  ) +
  geom_line(
    aes(y = high90, color = "IC 90%", linetype = "IC 90%"),
    linewidth = 0.8
  ) +
  geom_point(
    data = SriLankan_monthly,  # votre tableau avec les points
    aes(x = time, y = Cases, color = "Données observées",linetype="Données observées"),
    size = 2
  ) +
  scale_color_manual(
    name = "Statistique",
    values = c("Moyenne" = "red", "Médiane" = "blue", "IC 90%" = "blue","Données observées"="purple")
  ) +
  scale_linetype_manual(
    name = "Statistique",
    values = c("Moyenne" = "solid", "Médiane" = "solid", "IC 90%" = "dashed","Données observées"=NA)
  ) +
  
  labs(
    x = "Temps (en mois)",
    y = "Nombre d'infecté",
    title = "Dynamique d'infection chez l'homme, 5% des cas observés"
  ) +
  coord_cartesian(ylim = c(0, 40000))+
  theme_minimal()

######################################################################
######################################################################
###################### Model stochastique ############################
######################################################################
######################################################################
modele_dengue_tauleap <- function(y0, param, tmax, dt = 1) {
  # Paramètres fixes
  Z = param
  gamma = 1 / 2
  beta_v = 0.375
  mu_v = 1 / 6
  mu_h = 0
  beta_h = 0.75
  Nh = 23300000  # Population SriLanka
  
  # Initialisation
  Ih = y0[1]
  Rh = y0[2]
  Iv = y0[3]
  
  t = 0
  res = data.frame(time = t, Ih = Ih, Rh = Rh, Iv = Iv)
  
  while (t < tmax) {
    z = z_t(Z, t)
    Nv = z * Nh
    
    # Taux d’événements par unité de temps
    rate_infect_h = beta_h * z * Iv * (Nh - Ih - Rh) / Nh
    rate_recover_h = gamma * Ih
    rate_infect_v = beta_v * (Nv - Iv) * Ih / Nh
    rate_die_v = mu_v * Iv
    
    # Tirages de Poisson pour chaque événement sur l’intervalle dt
    n_infect_h = rpois(1, rate_infect_h * dt)
    n_recover_h = rpois(1, rate_recover_h * dt)
    n_infect_v = rpois(1, rate_infect_v * dt)
    n_die_v = rpois(1, rate_die_v * dt)
    
    # Mise à jour des compartiments
    Ih = Ih + n_infect_h - n_recover_h
    Rh = Rh + n_recover_h
    Iv = Iv + n_infect_v - n_die_v
    
    # Éviter les valeurs négatives
    Ih = max(Ih, 0)
    Rh = max(Rh, 0)
    Iv = max(Iv, 0)
    
    t = t + dt
    res = rbind(res, data.frame(time = t, Ih = Ih, Rh = Rh, Iv = Iv))
  }
  return(res)
}
y0 = c(10000, 100000, 0)
Z = runif(20,0.2,1.5)
stoch = modele_dengue_tauleap(y0, param = Z, tmax = 156, dt = 1)
plot(stoch$time, stoch$Ih, type='l')

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
    Date = as.Date(paste(Year, Month, "1", sep = "-")))

par(mfrow=c(1,1))
ggplot(data = z_estim, aes(x=Date,y=Mode))+
  geom_point()+
  geom_line(col="darkorchid")

plot(z_estim$Mode, ylab = "z", xlab = "Time")
lines(z_estim$Mode, col="darkorchid")

model = lm(z_estim$Mode~weather_periodized$meanPrecip)
anova(model)

model = lm(z_estim$Mode~weather_periodized$meanTemp)
anova(model)
