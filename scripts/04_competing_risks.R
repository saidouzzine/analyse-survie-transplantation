
# =====================================================
# 04_competing_risks.R
# Analyse des risques concurrents - Cohorte adulte
# =====================================================

library(survival)
library(cmprsk)

# 1. Chargement de la cohorte adulte
adults <- readRDS("data/leucemie2.rds")

# 2. Contrôles de cohérence
stopifnot(nrow(adults) == 1835)
stopifnot(all(adults$age >= 18))
stopifnot(all(complete.cases(
  adults[, c("time", "status", "type")]
)))
stopifnot(all(adults$time >= 0))
stopifnot(all(adults$status %in% 0:2))

cat("\n===== REPARTITION DES EVENEMENTS =====\n")
print(table(adults$status, adults$type))

# 3. Survie sans premier événement
km_event_free <- survfit(
  Surv(time / 365.25, status != 0) ~ 1,
  data = adults
)

cat("\n===== SURVIE SANS EVENEMENT =====\n")

resultats_km <- summary(
  km_event_free,
  times = c(1, 5)
)

print(resultats_km)

# 4. Incidences cumulées des événements concurrents
ci <- cuminc(
  ftime = adults$time / 365.25,
  fstatus = adults$status,
  cencode = 0
)

cat("\n===== INCIDENCES CUMULEES =====\n")

resultats_ci <- timepoints(
  ci,
  times = c(1, 5)
)

print(resultats_ci)

# 5. Vérification de l'identité fondamentale
# S(t) + F_rechute(t) + F_deces(t) = 1

survie <- resultats_km$surv
incidences <- resultats_ci$est

somme <- survie + colSums(incidences)

cat("\n===== VERIFICATION DE COHERENCE =====\n")
print(somme)

stopifnot(all(abs(somme - 1) < 0.01))

# 6. Graphique des incidences cumulées
dir.create(
  "figures",
  showWarnings = FALSE,
  recursive = TRUE
)

png(
  "figures/incidences_cumulees.png",
  width = 1200,
  height = 800,
  res = 140
)

plot(
  ci,
  xlab = "Années depuis la transplantation",
  ylab = "Incidence cumulée",
  main = "Rechute et décès sans rechute - Adultes",
  lwd = 2
)

dev.off()

cat("\nAnalyse des risques concurrents terminée.\n")
