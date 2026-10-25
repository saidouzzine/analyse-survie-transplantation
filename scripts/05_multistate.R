
# =====================================================
# 05_multistate.R
# Modele multi-etats : transplantation, rechute, deces
# Cohorte adulte : n = 1835
# =====================================================

library(survival)
library(mstate)

# 1. Chargement des donnees
adults <- readRDS("data/leucemie2.rds")

# 2. Verification de la cohorte
stopifnot(nrow(adults) == 1835)
stopifnot(all(adults$age >= 18))
stopifnot(!anyDuplicated(adults$patid))

variables <- c(
  "patid", "rel", "relstat",
  "srv", "srvstat", "age", "score"
)

stopifnot(all(complete.cases(adults[, variables])))

# Verification des chronologies
stopifnot(all(adults$rel > 0))
stopifnot(all(adults$srv > 0))

stopifnot(all(
  adults$rel[adults$relstat == 1] <
    adults$srv[adults$relstat == 1]
))

stopifnot(all(
  adults$rel[adults$relstat == 0] ==
    adults$srv[adults$relstat == 0]
))

# 3. Matrice des transitions
tmat <- trans.illdeath(
  names = c("transplant", "relapse", "death")
)

cat("\n===== MATRICE DES TRANSITIONS =====\n")
print(tmat)

# 4. Preparation des temps
adults$time_rel <- adults$rel
adults$status_rel <- adults$relstat

adults$time_death <- adults$srv
adults$status_death <- adults$srvstat

# 5. Passage au format long
long <- msprep(
  time = c(NA, "time_rel", "time_death"),
  status = c(NA, "status_rel", "status_death"),
  data = adults,
  trans = tmat,
  id = "patid",
  keep = c("age", "score")
)

# 6. Controles du format long
stopifnot(all(long$Tstop >= long$Tstart))

cat("\n===== TRANSITIONS OBSERVEES =====\n")
print(table(long$trans, long$status))

cat("\n===== DIMENSIONS DU FORMAT LONG =====\n")
print(dim(long))

# 7. Modele de Cox stratifie par transition
fit_ms <- coxph(
  Surv(Tstart, Tstop, status) ~ strata(trans),
  data = long,
  method = "breslow"
)

cat("\n===== MODELE MULTI-ETATS =====\n")
print(fit_ms)

# 8. Risques cumules de transition
ms <- msfit(
  fit_ms,
  trans = tmat
)

# 9. Probabilites d'occupation des etats
pt <- probtrans(
  ms,
  predt = 0
)


# 10. Probabilites d'occupation a 1 et 5 ans
cat("\n===== PROBABILITES A 1 ET 5 ANS =====\n")

for (annee in c(1, 5)) {
  
  temps <- annee * 365.25
  
  # Derniere estimation disponible avant ou a t
  i <- max(which(pt[[1]]$time <= temps))
  
  probabilites <- pt[[1]][
    i,
    c("pstate1", "pstate2", "pstate3")
  ]
  
  cat("\nEcheance :", annee, "an(s)\n")
  cat("Temps de la derniere estimation :",
      pt[[1]]$time[i], "jours\n")
  
  print(probabilites)
  
  cat("Somme des probabilites :",
      sum(probabilites), "\n")
  
  stopifnot(
    abs(sum(probabilites) - 1) < 1e-6
  )
}


# =====================================================
# Figure - Probabilités d'occupation des états
#      Estimateur d'Aalen-Johansen
# =====================================================

# Probabilités depuis l'état initial (état 1)
occupation <- pt[[1]]

# Conversion du temps en années
temps_annees <- occupation$time / 365.25

# Vérification : les probabilités totalisent 1
sommes <- rowSums(
  occupation[, c("pstate1", "pstate2", "pstate3")]
)

stopifnot(all(abs(sommes - 1) < 1e-6))

# Création du dossier de sortie
dir.create("figures", showWarnings = FALSE)

# Fonction de tracé, pour l'écran et le fichier PNG
tracer_occupation <- function() {
  
  matplot(
    x = temps_annees,
    y = as.matrix(occupation[
      , c("pstate1", "pstate2", "pstate3")
    ]),
    type = "s",
    lty = 1,
    lwd = 2,
    col = c("darkgreen", "blue", "red"),
    xlab = "Années depuis la transplantation",
    ylab = "Probabilité d'occupation",
    main = "Probabilités d'occupation des trois états",
    ylim = c(0, 1)
  )
  
  
  legend(
    "topright",
    inset = c(0.02, 0.02),
    legend = c(
      "État 1 : vivant sans rechute",
      "État 2 : vivant après rechute",
      "État 3 : décédé"
    ),
    col = c("darkgreen", "blue", "red"),
    lty = 1,
    lwd = 2,
    bty = "n",
    cex = 0.80
  )
  
}

# Affichage dans RStudio
tracer_occupation()

# Enregistrement pour le rapport scientifique
png(
  filename = "figures/occupation_etats_AJ.png",
  width = 2400,
  height = 1500,
  res = 300
)

tracer_occupation()
dev.off()

cat("\nFigure enregistrée dans le dossier figures.\n")

