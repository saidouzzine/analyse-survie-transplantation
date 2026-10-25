
# =====================================================
# 02_kaplan_meier.R
# Modèle de Kaplan - survie globale
# =====================================================

library(survival)

# 1. Chargement des données adultes
fichier <- "data/leucemie2.rds"

if (!file.exists(fichier)) {
  stop("Fichier introuvable : ", fichier)
}

adults <- readRDS(fichier)

# Vérifications
stopifnot(nrow(adults) == 1835)
stopifnot(all(adults$age >= 18))
stopifnot(all(!is.na(adults[, c(
  "srv", "srvstat", "score"
)])))

# 2. Survie globale : Kaplan-Meier
km_global <- survfit(
  Surv(srv / 365.25, srvstat) ~ 1,
  data = adults
)

# Estimations à 1, 5 et 8 ans
print(summary(
  km_global,
  times = c(1, 5, 8)
))

# 3. Création du dossier figures
dir.create(
  "figures",
  showWarnings = FALSE,
  recursive = TRUE
)

# 4. Graphique de survie globale
png(
  "figures/survie_globale.png",
  width = 1200,
  height = 800,
  res = 140
)

plot(
  km_global,
  xlab = "Années depuis la transplantation",
  ylab = "Probabilité de survie globale",
  main = "Survie globale - Kaplan-Meier",
  conf.int = TRUE,
  mark.time = TRUE,
  lwd = 2
)

dev.off()

# 5. Survie selon le score de Gratwohl
km_score <- survfit(
  Surv(srv / 365.25, srvstat) ~ score,
  data = adults
)

# Test du log-rank
test_logrank <- survdiff(
  Surv(srv / 365.25, srvstat) ~ score,
  data = adults
)

print(test_logrank)

p_value <- pchisq(
  test_logrank$chisq,
  df = length(test_logrank$n) - 1,
  lower.tail = FALSE
)

cat("p-value du log-rank :", p_value, "\n")

# 6. Graphique selon le score
png(
  "figures/survie_selon_score.png",
  width = 1200,
  height = 800,
  res = 140
)

plot(
  km_score,
  col = c("red", "darkgreen", "blue"),
  lwd = 2,
  xlab = "Années depuis la transplantation",
  ylab = "Probabilité de survie",
  main = "Survie selon le score de Gratwohl"
)

legend(
  "topright",
  legend = names(km_score$strata),
  col = c("red", "darkgreen", "blue"),
  lwd = 2,
  bty = "n"
)

dev.off()

cat("\nAnalyse Kaplan-Meier terminée.\n")
