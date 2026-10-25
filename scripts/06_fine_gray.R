
# =====================================================
# 06_fine_gray.R
# Modele de Fine-Gray : risque de rechute
# Cohorte adulte : n = 1835
# =====================================================

library(cmprsk)

# 1. Chargement de la cohorte adulte
adults <- readRDS("data/leucemie2.rds")

# 2. Controles de coherence
stopifnot(nrow(adults) == 1835)
stopifnot(all(adults$age >= 18))
stopifnot(all(adults$status %in% 0:2))
stopifnot(all(adults$time > 0))
stopifnot(all(complete.cases(
  adults[, c("time", "status", "score")]
)))

cat("\n===== REPARTITION DES EVENEMENTS =====\n")
print(table(adults$status, adults$type))

# 3. Codage du groupe de risque
adults$score <- factor(
  adults$score,
  levels = c("Low risk", "Medium risk", "High risk")
)

# Matrice des covariables
x <- model.matrix(~ score, data = adults)[, -1, drop = FALSE]

cat("\n===== VARIABLES EXPLICATIVES =====\n")
print(colnames(x))

# 4. Modele de Fine-Gray
# Evenement d'interet : rechute (1)
# Evenement concurrent : deces sans rechute (2)
# Censure : 0

fit_fg <- crr(
  ftime = adults$time / 365.25,
  fstatus = adults$status,
  cov1 = x,
  failcode = 1,
  cencode = 0
)

cat("\n===== MODELE DE FINE-GRAY =====\n")
print(summary(fit_fg))

# 5. Incidences cumulees par groupe de risque
ci_score <- cuminc(
  ftime = adults$time / 365.25,
  fstatus = adults$status,
  group = adults$score,
  cencode = 0
)

cat("\n===== COMPARAISON DES GROUPES =====\n")
print(ci_score$Tests)

cat("\n===== INCIDENCES CUMULEES A 1 ET 5 ANS =====\n")
print(timepoints(ci_score, times = c(1, 5)))


# =====================================================
# 6. Graphique des incidences cumulees par groupe
# =====================================================

dir.create(
  "figures",
  showWarnings = FALSE,
  recursive = TRUE
)

# Couleurs par groupe de risque
couleurs <- c(
  "Low risk" = "#009E73",     # Vert
  "Medium risk" = "#0072B2",  # Bleu
  "High risk" = "#D55E00"     # Orange
)

# Les courbes sont ordonnees selon names(ci_score)
noms_courbes <- names(ci_score)
noms_courbes <- noms_courbes[
  noms_courbes != "Tests"
]

# Identification des groupes et des evenements
groupes <- sub(" [12]$", "", noms_courbes)
evenements <- sub("^.* ([12])$", "\\1", noms_courbes)

couleurs_courbes <- unname(couleurs[groupes])

# 1 = rechute : ligne continue
# 2 = deces sans rechute : ligne pointillee
types_lignes <- ifelse(evenements == "1", 1, 2)

png(
  "figures/fine_gray_rechute.png",
  width = 1500,
  height = 1050,
  res = 150
)

# Disposition : graphique en haut, legendes en bas
layout(
  matrix(c(1, 2), nrow = 2),
  heights = c(4, 1.2)
)

# -------- Graphique principal --------

par(mar = c(4.5, 5, 3.5, 2))

plot(
  ci_score,
  col = couleurs_courbes,
  lty = types_lignes,
  lwd = 2.5,
  xlab = "Annees depuis la transplantation",
  ylab = "Incidence cumulee",
  main = paste(
    "Incidences cumulees de rechute",
    "et de deces sans rechute"
  ),
  ylim = c(0, 0.65),
  curvlab = rep("", length(noms_courbes))
)

# -------- Zone reservee aux legendes --------

par(mar = c(0, 0, 0, 0))

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

legend(
  x = 0.04,
  y = 0.95,
  legend = names(couleurs),
  col = unname(couleurs),
  lty = 1,
  lwd = 2.5,
  bty = "n",
  cex = 1.05,
  title = "Groupe de risque"
)

legend(
  x = 0.55,
  y = 0.95,
  legend = c(
    "Rechute",
    "Deces sans rechute"
  ),
  col = "grey30",
  lty = c(1, 2),
  lwd = 2.5,
  bty = "n",
  cex = 1.05,
  title = "Evenement"
)

dev.off()

cat("\nGraphique des risques concurrents en couleurs genere.\n")


# =====================================================
# 7. Diagnostic de proportionnalite - Fine-Gray
# =====================================================

# Residu de score pour chaque covariable
residus_fg <- fit_fg$res

# Temps des evenements d'interet
temps_fg <- fit_fg$uftime

stopifnot(
  nrow(residus_fg) == length(temps_fg),
  ncol(residus_fg) == ncol(x)
)

dir.create("figures", showWarnings = FALSE)

png(
  "figures/diagnostic_fine_gray.png",
  width = 1200,
  height = 650,
  res = 140
)

par(mfrow = c(1, 2))

for (j in seq_len(ncol(residus_fg))) {
  
  plot(
    temps_fg,
    residus_fg[, j],
    pch = 16,
    cex = 0.5,
    col = "grey50",
    xlab = "Temps depuis la transplantation (annees)",
    ylab = "Residu de score",
    main = colnames(x)[j]
  )
  
  abline(h = 0, lty = 2)
  
  lines(
    lowess(temps_fg, residus_fg[, j]),
    col = "blue",
    lwd = 2
  )
}

dev.off()

cat("\nDiagnostic graphique Fine-Gray termine.\n")


# ==========================================
# 8. Test des effets dependant du temps 
# ==========================================

# Deux interactions temporelles distinctes :
# Medium risk x log(1+t)
# High risk x log(1+t)

fit_fg_temps <- crr(
  ftime = adults$time / 365.25,
  fstatus = adults$status,
  cov1 = x,
  cov2 = x,
  tf = function(t) {
    cbind(
      log1p(t),
      log1p(t)
    )
  },
  failcode = 1,
  cencode = 0
)

cat("\n===== MODELE AVEC EFFETS TEMPORELS =====\n")
print(summary(fit_fg_temps))

# Identification des coefficients temporels
# Ordre : cov1 (2 coefficients), puis cov2 (2 coefficients)

indices_temps <- (
  ncol(x) + 1
):(
  2 * ncol(x)
)

beta_temps <- fit_fg_temps$coef[indices_temps]

variance_temps <- fit_fg_temps$var[
  indices_temps,
  indices_temps,
  drop = FALSE
]

# Test global de Wald
stat_wald <- as.numeric(
  t(beta_temps) %*%
    solve(variance_temps, beta_temps)
)

ddl <- length(beta_temps)

p_value <- pchisq(
  stat_wald,
  df = ddl,
  lower.tail = FALSE
)

cat("\n===== TEST GLOBAL DES EFFETS TEMPORELS =====\n")
cat("Statistique de Wald :", stat_wald, "\n")
cat("Degres de liberte :", ddl, "\n")
cat("Valeur p :", p_value, "\n")

