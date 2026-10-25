
# =====================================================
# 03_cox.R
# Modèle de Cox - Survie globale
# =====================================================

library(survival)

# 1. Chargement de la cohorte adulte
adults <- readRDS("data/leucemie2.rds")

# Contrôles
stopifnot(nrow(adults) == 1835)
stopifnot(all(adults$age >= 18))

variables <- c("srv", "srvstat", "score", "age")

stopifnot(
  all(complete.cases(adults[, variables]))
)

# 2. Définition du groupe de référence
adults$score <- factor(
  adults$score,
  levels = c("Low risk", "Medium risk", "High risk")
)

# 3. Ajustement du modèle de Cox
mod_cox <- coxph(
  Surv(srv / 365.25, srvstat) ~ score + age,
  data = adults,
  x = TRUE
)

# 4. Résultats du modèle
cat("\n===== MODELE DE COX =====\n")
print(summary(mod_cox))

# 5. Hazard ratios et intervalles de confiance
cat("\n===== HAZARD RATIOS =====\n")

resultats_hr <- data.frame(
  Variable = names(coef(mod_cox)),
  HR = exp(coef(mod_cox)),
  exp(confint(mod_cox)),
  row.names = NULL,
  check.names = FALSE
)

print(resultats_hr)

# 6. Hypothèse des risques proportionnels
cat("\n===== TEST DE SCHOENFELD =====\n")

test_ph <- cox.zph(mod_cox)

print(test_ph)

# 7. Graphiques diagnostiques
dir.create(
  "figures",
  showWarnings = FALSE,
  recursive = TRUE
)

png(
  "figures/diagnostic_cox.png",
  width = 1200,
  height = 1400,
  res = 140
)

par(mfrow = c(3, 1))
plot(test_ph)

dev.off()

cat("\nAnalyse de Cox terminée.\n")
