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