###############################################################################
# EMA-ANALYSE – INTERVENTIONSEFFEKTE
# Charité Masterarbeit – Cannabis, Einsamkeit & Craving
#
# Voraussetzung: EMA_cleaning_v5_final.R wurde ausgeführt und
#                EMA_cleaned_data.RData wurde gespeichert
#
# MODELLSTRATEGIE (präregistriert, auf Empfehlung Betreuer):
#   Schritt 1: Nullmodell → ICC berechnen
#   Schritt 2: Random Intercept Modell (Hauptmodell, präregistriert)
#   Schritt 3: Random Slopes Modell (explorativ, nicht präregistriert)
#              → Modellvergleich via LRT und AIC
#              → Random Intercept behalten wenn kein besserer Fit
#
# HAUPTANALYSEN (präregistriert):
#   Modell 1 (H1/H3): craving_post ~ craving_pre + intervention_type +
#                      protocol_version + (1 | Participant)
#   Modell 2 (H2/H4): loneliness_direct_post ~ loneliness_direct_pre +
#                      intervention_type + protocol_version + (1 | Participant)
#   Referenzkategorie: Control → negativer β = Hypothese bestätigt
#
# EFFEKTMASSE (auf Empfehlung beider Korrektoren):
#   Primär:   Standardisierte β + 95% CI + Bayes Faktoren (BF₁₀)
#   Sekundär: Marginal R² + Conditional R²
#   p-Werte werden zusätzlich berichtet (two-tailed, konservativ)
#   Cohen's d nur deskriptiv für within-condition Vergleiche
#
# ROBUSTNESS CHECKS:
#   1. Sensitivitätsanalyse Wortanzahl (>= 4 / >= 6 / >= 8 Wörter)
#   2. Phasenvergleich: old_8plus_prompts vs. new_4_prompts getrennt
#   3. Intention-to-treat (ITT): alle Paare ohne pair_valid Filter
#      → adressiert trial-level compliance bias (Kritik Korrektor 2)
###############################################################################

library(tidyverse)
library(lme4)
library(lmerTest)
library(effectsize)
library(performance)
library(BayesFactor)   # NEU: für Bayes Faktoren
library(ggplot2)

# =============================================================================
# 0. DATEN LADEN
# =============================================================================

load("/Users/juliabrandl/Desktop/EMA_cleaned_data.RData")

cat("=== GELADENE OBJEKTE ===\n")
cat("  emi_pairs        (>= 4 Wörter):", nrow(emi_pairs),       "Zeilen\n")
cat("  emi_pairs_sens6  (>= 6 Wörter):", nrow(emi_pairs_sens6), "Zeilen\n")
cat("  emi_pairs_sens8  (>= 8 Wörter):", nrow(emi_pairs_sens8), "Zeilen\n")

# =============================================================================
# 1. ANALYSEDATENSATZ VORBEREITEN
#
# Ausschlusskriterien (präregistriert):
#   - Nur Interventionsphase (Studientage 8–28)
#   - Keine Abbrecher (compliance < 70% ODER < 11 Tage)
#   - Keine fehlenden Pre- oder Post-Werte
#   - Nicht beide Pre-Scores gleichzeitig = 1
#   - pair_valid == TRUE (Reappraisal ausreichend / Sport bestätigt / Spiel abgeschlossen)
# =============================================================================

prepare_data <- function(pairs_df) {
  pairs_df %>%
    filter(
      phase == "Intervention",
      aborted == FALSE | is.na(aborted),
      !is.na(craving_pre),
      !is.na(craving_post),
      !is.na(loneliness_direct_pre),
      !is.na(loneliness_direct_post),
      !(craving_pre == 1 & loneliness_direct_pre == 1),
      pair_valid == TRUE
    ) %>%
    mutate(
      intervention_type = factor(intervention_type,
                                  levels = c("Control", "Reappraisal", "Movement")),
      protocol_version  = factor(protocol_version)
    )
}

emi_main  <- prepare_data(emi_pairs)
emi_sens6 <- prepare_data(emi_pairs_sens6)
emi_sens8 <- prepare_data(emi_pairs_sens8)

cat("\n=== STICHPROBENGRÖSSEN NACH FILTER ===\n")
cat("Haupt (>= 4):", nrow(emi_main),  "Paare /", n_distinct(emi_main$Participant),  "Probanden\n")
cat("Sens  (>= 6):", nrow(emi_sens6), "Paare /", n_distinct(emi_sens6$Participant), "Probanden\n")
cat("Sens  (>= 8):", nrow(emi_sens8), "Paare /", n_distinct(emi_sens8$Participant), "Probanden\n")

# =============================================================================
# 1b. Z-STANDARDISIERUNG FÜR STANDARDISIERTE β-KOEFFIZIENTEN
#
# WARUM: Unstandardisierte β sind in der Originalskala (1–7) und schwer
# zu vergleichen. Standardisierte β zeigen den Effekt in Standardabweichungen
# und sind als Effektgröße besser interpretierbar.
# Auf Empfehlung beider Korrektoren: standardisierte β als primäres Effektmaß.
# =============================================================================

emi_main <- emi_main %>%
  mutate(
    craving_post_z    = as.numeric(scale(craving_post)),
    craving_pre_z     = as.numeric(scale(craving_pre)),
    loneliness_post_z = as.numeric(scale(loneliness_direct_post)),
    loneliness_pre_z  = as.numeric(scale(loneliness_direct_pre))
  )

# =============================================================================
# 2. DESKRIPTIVE STATISTIK
# =============================================================================

cat("\n\n=== DESKRIPTIVE STATISTIK ===\n")

# Stichprobe
cat("--- Stichprobe ---\n")
cat("Probanden gesamt:", n_distinct(emi_pairs$Participant), "\n")
cat("Abbrecher (< 70% Prompts ODER < 11 Tage):", sum(compliance$aborted, na.rm = TRUE), "\n")
cat("In Hauptanalyse:", n_distinct(emi_main$Participant), "\n\n")

# Compliance
cat("--- Compliance ---\n")
compliance %>%
  summarise(
    M   = round(mean(compliance_rate, na.rm = TRUE), 3),
    SD  = round(sd(compliance_rate,   na.rm = TRUE), 3),
    Min = round(min(compliance_rate,  na.rm = TRUE), 3),
    Max = round(max(compliance_rate,  na.rm = TRUE), 3)
  ) %>%
  print()

# Pre/Post Mittelwerte pro Interventionstyp
cat("\n--- Pre/Post Mittelwerte pro Bedingung ---\n")
cat("Negatives Delta = Verbesserung\n\n")
emi_main %>%
  group_by(intervention_type) %>%
  summarise(
    n               = n(),
    craving_pre_M   = round(mean(craving_pre,             na.rm = TRUE), 2),
    craving_pre_SD  = round(sd(craving_pre,               na.rm = TRUE), 2),
    craving_post_M  = round(mean(craving_post,            na.rm = TRUE), 2),
    craving_post_SD = round(sd(craving_post,              na.rm = TRUE), 2),
    craving_delta_M = round(mean(delta_craving,           na.rm = TRUE), 2),
    lone_pre_M      = round(mean(loneliness_direct_pre,   na.rm = TRUE), 2),
    lone_pre_SD     = round(sd(loneliness_direct_pre,     na.rm = TRUE), 2),
    lone_post_M     = round(mean(loneliness_direct_post,  na.rm = TRUE), 2),
    lone_post_SD    = round(sd(loneliness_direct_post,    na.rm = TRUE), 2),
    lone_delta_M    = round(mean(delta_loneliness_direct, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  print()

# Protokollversionsvergleich
cat("\n--- Vergleich nach Protokollversion ---\n")
emi_main %>%
  group_by(protocol_version, intervention_type) %>%
  summarise(
    n               = n(),
    craving_delta_M = round(mean(delta_craving,           na.rm = TRUE), 2),
    lone_delta_M    = round(mean(delta_loneliness_direct, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  print()

# Within-Condition Cohen's d
# Zeigt ob jede Bedingung isoliert einen Effekt hat.
# Wichtig für Control: Hat das Memory-Spiel selbst Craving/Loneliness reduziert?
# Cohen's d hier nur deskriptiv – nicht als primäres Effektmaß der Regressionen
cat("\n--- Within-Condition Cohen's d (deskriptiv) ---\n")
cat("Nur deskriptiv – primäres Effektmaß sind standardisierte β\n\n")

for (itype in c("Reappraisal", "Movement", "Control")) {
  sub <- emi_main %>% filter(intervention_type == itype)
  cat(sprintf("\n=== %s (n = %d) ===\n", itype, nrow(sub)))

  d_craving <- cohens_d(sub$craving_post, sub$craving_pre, paired = TRUE)
  cat(sprintf("  Craving:   d = %.3f [%.3f, %.3f]\n",
              d_craving$Cohens_d, d_craving$CI_low, d_craving$CI_high))

  d_lone <- cohens_d(sub$loneliness_direct_post, sub$loneliness_direct_pre,
                     paired = TRUE)
  cat(sprintf("  Loneliness: d = %.3f [%.3f, %.3f]\n",
              d_lone$Cohens_d, d_lone$CI_low, d_lone$CI_high))
}

# =============================================================================
# 3. SCHRITTWEISE MODELLSTRATEGIE – CRAVING (H1 & H3)
#
# Schritt 1: Nullmodell → ICC berechnen
#   Wie viel der Varianz im Craving-Post liegt zwischen Personen?
#   Hoher ICC = Random Intercept ist wichtig
#
# Schritt 2: Random Intercept Modell (präregistriertes Hauptmodell)
#   ANCOVA-Ansatz: craving_pre als Kovariate kontrolliert Ausgangswert
#   → methodisch besser als einfacher Pre-Post-Vergleich
#
# Schritt 3: Random Slopes Modell (explorativ)
#   Testet ob Interventionseffekte sich zwischen Personen unterscheiden
#   → nur berichten wenn signifikant besser als Random Intercept
# =============================================================================

cat("\n\n### HAUPTANALYSE CRAVING (H1: Reappraisal, H3: Movement) ###\n")
cat("Referenz: Control | Negatives β = besser als Control\n\n")

# ── Schritt 1: Nullmodell Craving ──────────────────────────────────────────
cat("── Schritt 1: Nullmodell ──\n")
m1_null <- lmer(craving_post ~ 1 + (1 | Participant),
                data = emi_main, REML = FALSE)

vc_c        <- as.data.frame(VarCorr(m1_null))
icc_craving <- vc_c$vcov[vc_c$grp == "Participant"] / sum(vc_c$vcov)
cat(sprintf("ICC Craving = %.3f → %.1f%% der Varianz durch Personenunterschiede\n",
            icc_craving, icc_craving * 100))

# ── Schritt 2: Random Intercept Modell (präregistriert) ───────────────────
cat("\n── Schritt 2: Random Intercept Modell (präregistriert) ──\n")
m1_ri <- lmer(craving_post ~ craving_pre + intervention_type + protocol_version +
                (1 | Participant),
              data = emi_main, REML = FALSE)
print(summary(m1_ri))

cat("\n95% Konfidenzintervalle (Wald):\n")
print(round(confint(m1_ri, method = "Wald"), 3))

# Marginal R² und Conditional R²
cat("\nModell-Fit (R²):\n")
r2_m1 <- performance::r2(m1_ri)
cat(sprintf("  Marginal R²    = %.3f (Fixed Effects)\n",          r2_m1$R2_marginal))
cat(sprintf("  Conditional R² = %.3f (Fixed + Random Effects)\n", r2_m1$R2_conditional))

# ── Schritt 3: Random Slopes Modell (explorativ) ──────────────────────────
cat("\n── Schritt 3: Random Slopes Modell (explorativ) ──\n")
cat("Testet: Unterscheiden sich Interventionseffekte zwischen Personen?\n\n")
m1_rs <- lmer(craving_post ~ craving_pre + intervention_type + protocol_version +
                (1 + intervention_type | Participant),
              data = emi_main, REML = FALSE)

cat("Modellvergleich (LRT und AIC):\n")
print(anova(m1_ri, m1_rs))
cat(sprintf("AIC Random Intercept: %.1f\n", AIC(m1_ri)))
cat(sprintf("AIC Random Slopes:    %.1f\n", AIC(m1_rs)))

# =============================================================================
# 4. SCHRITTWEISE MODELLSTRATEGIE – LONELINESS (H2 & H4)
# =============================================================================

cat("\n\n### HAUPTANALYSE LONELINESS (H2: Reappraisal, H4: Movement) ###\n")
cat("Referenz: Control | Negatives β = besser als Control\n\n")

# ── Schritt 1: Nullmodell Loneliness ──────────────────────────────────────
cat("── Schritt 1: Nullmodell ──\n")
m2_null <- lmer(loneliness_direct_post ~ 1 + (1 | Participant),
                data = emi_main, REML = FALSE)

vc_l     <- as.data.frame(VarCorr(m2_null))
icc_lone <- vc_l$vcov[vc_l$grp == "Participant"] / sum(vc_l$vcov)
cat(sprintf("ICC Loneliness = %.3f → %.1f%% der Varianz durch Personenunterschiede\n",
            icc_lone, icc_lone * 100))

# ── Schritt 2: Random Intercept Modell (präregistriert) ───────────────────
cat("\n── Schritt 2: Random Intercept Modell (präregistriert) ──\n")
m2_ri <- lmer(loneliness_direct_post ~ loneliness_direct_pre + intervention_type +
                protocol_version + (1 | Participant),
              data = emi_main, REML = FALSE)
print(summary(m2_ri))

cat("\n95% Konfidenzintervalle (Wald):\n")
print(round(confint(m2_ri, method = "Wald"), 3))

# Marginal R² und Conditional R²
cat("\nModell-Fit (R²):\n")
r2_m2 <- performance::r2(m2_ri)
cat(sprintf("  Marginal R²    = %.3f (Fixed Effects)\n",          r2_m2$R2_marginal))
cat(sprintf("  Conditional R² = %.3f (Fixed + Random Effects)\n", r2_m2$R2_conditional))

# ── Schritt 3: Random Slopes Modell (explorativ) ──────────────────────────
cat("\n── Schritt 3: Random Slopes Modell (explorativ) ──\n")
cat("Testet: Unterscheiden sich Interventionseffekte zwischen Personen?\n\n")
m2_rs <- lmer(loneliness_direct_post ~ loneliness_direct_pre + intervention_type +
                protocol_version + (1 + intervention_type | Participant),
              data = emi_main, REML = FALSE)

cat("Modellvergleich (LRT und AIC):\n")
print(anova(m2_ri, m2_rs))
cat(sprintf("AIC Random Intercept: %.1f\n", AIC(m2_ri)))
cat(sprintf("AIC Random Slopes:    %.1f\n", AIC(m2_rs)))

# =============================================================================
# 5. ERGEBNISTABELLE HAUPTANALYSE (unstandardisiert)
# =============================================================================

cat("\n\n=== ERGEBNISTABELLE HAUPTANALYSE (unstandardisiert) ===\n")
cat("Berichtet wird das präregistrierte Random Intercept Modell.\n\n")

# Hilfsfunktion für Robustness Checks
koef_sens <- function(modell, outcome, cutoff) {
  k <- coef(summary(modell))
  tibble(
    cutoff    = cutoff,
    outcome   = outcome,
    vergleich = c("Reappraisal vs. Control", "Movement vs. Control"),
    beta      = round(k[c("intervention_typeReappraisal",
                           "intervention_typeMovement"), "Estimate"], 3),
    p         = round(k[c("intervention_typeReappraisal",
                           "intervention_typeMovement"), "Pr(>|t|)"], 3)
  )
}

k_m1 <- coef(summary(m1_ri))
k_m2 <- coef(summary(m2_ri))

ergebnisse_haupt <- bind_rows(
  tibble(
    outcome   = "Craving",
    vergleich = c("Reappraisal vs. Control", "Movement vs. Control"),
    beta      = round(k_m1[c("intervention_typeReappraisal",
                               "intervention_typeMovement"), "Estimate"], 3),
    SE        = round(k_m1[c("intervention_typeReappraisal",
                               "intervention_typeMovement"), "Std. Error"], 3),
    t         = round(k_m1[c("intervention_typeReappraisal",
                               "intervention_typeMovement"), "t value"], 3),
    p         = round(k_m1[c("intervention_typeReappraisal",
                               "intervention_typeMovement"), "Pr(>|t|)"], 3)
  ),
  tibble(
    outcome   = "Loneliness",
    vergleich = c("Reappraisal vs. Control", "Movement vs. Control"),
    beta      = round(k_m2[c("intervention_typeReappraisal",
                               "intervention_typeMovement"), "Estimate"], 3),
    SE        = round(k_m2[c("intervention_typeReappraisal",
                               "intervention_typeMovement"), "Std. Error"], 3),
    t         = round(k_m2[c("intervention_typeReappraisal",
                               "intervention_typeMovement"), "t value"], 3),
    p         = round(k_m2[c("intervention_typeReappraisal",
                               "intervention_typeMovement"), "Pr(>|t|)"], 3)
  )
) %>%
  mutate(sig = case_when(
    p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "ns"
  ))

print(ergebnisse_haupt)

# =============================================================================
# 5a. FDR-KORREKTUR FÜR MULTIPLE TESTING
#
# WARUM: Vier Hypothesen werden gleichzeitig getestet (H1-H4) was die
# Wahrscheinlichkeit eines falsch-positiven Ergebnisses erhöht.
# Auf Empfehlung von Korrektor 2 (Prof. Markett) werden FDR-korrigierte
# p-Werte zusätzlich als Sensitivity Check berichtet.
#
# FDR = False Discovery Rate: kontrolliert den Anteil falscher Positive
# unter allen signifikanten Ergebnissen – weniger streng als Bonferroni.
#
# WICHTIG: FDR-Werte sind nur ergänzend – das primäre Entscheidungskriterium
# bleiben die unkorrgierten p-Werte + standardisierte β + CI + BF₁₀.
# =============================================================================

cat("\n\n=== FDR-KORRIGIERTE P-WERTE (Sensitivity Check) ===\n")
cat("Primäres Entscheidungskriterium bleiben unkorrgiierte p-Werte + β + CI + BF\n")
cat("FDR nur ergänzend berichtet (Empfehlung Korrektor 2)\n\n")

# p-Werte der vier Hypothesen in der richtigen Reihenfolge
p_werte <- c(
  k_m1["intervention_typeReappraisal", "Pr(>|t|)"],  # H1: Reappraisal → Craving
  k_m2["intervention_typeReappraisal", "Pr(>|t|)"],  # H2: Reappraisal → Loneliness
  k_m1["intervention_typeMovement",    "Pr(>|t|)"],  # H3: Movement → Craving
  k_m2["intervention_typeMovement",    "Pr(>|t|)"]   # H4: Movement → Loneliness
)

p_fdr <- p.adjust(p_werte, method = "fdr")

cat(sprintf("H1 Reappraisal → Craving:    p = %.3f | p_fdr = %.3f\n",
            p_werte[1], p_fdr[1]))
cat(sprintf("H2 Reappraisal → Loneliness: p = %.3f | p_fdr = %.3f\n",
            p_werte[2], p_fdr[2]))
cat(sprintf("H3 Movement    → Craving:    p = %.3f | p_fdr = %.3f\n",
            p_werte[3], p_fdr[3]))
cat(sprintf("H4 Movement    → Loneliness: p = %.3f | p_fdr = %.3f\n",
            p_werte[4], p_fdr[4]))
#
# WARUM: Korrektoren empfehlen standardisierte β als primäres Effektmaß.
# Die Modelle werden mit z-standardisierten Variablen neu geschätzt.
# β in Standardabweichungen → besser vergleichbar und interpretierbar.
# =============================================================================

cat("\n\n=== STANDARDISIERTE MODELLE ===\n")
cat("β in Standardabweichungen – primäres Effektmaß\n\n")

m1_ri_std <- lmer(craving_post_z ~ craving_pre_z + intervention_type +
                    protocol_version + (1 | Participant),
                  data = emi_main, REML = FALSE)

m2_ri_std <- lmer(loneliness_post_z ~ loneliness_pre_z + intervention_type +
                    protocol_version + (1 | Participant),
                  data = emi_main, REML = FALSE)

k_m1_std <- coef(summary(m1_ri_std))
k_m2_std <- coef(summary(m2_ri_std))

ergebnisse_std <- bind_rows(
  tibble(
    outcome   = "Craving",
    vergleich = c("Reappraisal vs. Control", "Movement vs. Control"),
    beta_std  = round(k_m1_std[c("intervention_typeReappraisal",
                                   "intervention_typeMovement"), "Estimate"], 3),
    SE_std    = round(k_m1_std[c("intervention_typeReappraisal",
                                   "intervention_typeMovement"), "Std. Error"], 3),
    p_std     = round(k_m1_std[c("intervention_typeReappraisal",
                                   "intervention_typeMovement"), "Pr(>|t|)"], 3)
  ),
  tibble(
    outcome   = "Loneliness",
    vergleich = c("Reappraisal vs. Control", "Movement vs. Control"),
    beta_std  = round(k_m2_std[c("intervention_typeReappraisal",
                                   "intervention_typeMovement"), "Estimate"], 3),
    SE_std    = round(k_m2_std[c("intervention_typeReappraisal",
                                   "intervention_typeMovement"), "Std. Error"], 3),
    p_std     = round(k_m2_std[c("intervention_typeReappraisal",
                                   "intervention_typeMovement"), "Pr(>|t|)"], 3)
  )
) %>%
  mutate(sig = case_when(
    p_std < .001 ~ "***", p_std < .01 ~ "**", p_std < .05 ~ "*", TRUE ~ "ns"
  ))

cat("95% KI standardisiertes Modell Craving (Wald):\n")
ci_std_m1 <- confint(m1_ri_std, method = "Wald")
print(round(ci_std_m1, 3))

cat("\n95% KI standardisiertes Modell Loneliness (Wald):\n")
ci_std_m2 <- confint(m2_ri_std, method = "Wald")
print(round(ci_std_m2, 3))

cat("\nStandardisierte Koeffizienten:\n")
print(ergebnisse_std)

# =============================================================================
# 5c. BAYES FAKTOREN
#
# WARUM: Auf Empfehlung beider Korrektoren. BF₁₀ quantifiziert die Evidenz
# für H1 (es gibt einen Interventionseffekt) vs. H0 (kein Effekt).
#
# Besonders nützlich für nicht-signifikante Ergebnisse:
#   p = .056 sagt nur "nicht signifikant"
#   BF sagt ob das Evidenz für H0 ist oder nur insufficient power
#
# Vorgehen: Vollmodell (mit intervention_type) / Nullmodell (ohne)
#
# Interpretation:
#   BF > 3     = moderate Evidenz für H1 (Interventionseffekt)
#   BF > 10    = starke Evidenz für H1
#   BF < 1/3   = moderate Evidenz für H0 (kein Effekt)
#   1/3 bis 3  = inconclusive (Daten nicht aussagekräftig genug)
# =============================================================================

cat("\n\n=== BAYES FAKTOREN ===\n")
cat("BF > 3    = moderate Evidenz für Interventionseffekt\n")
cat("BF < 0.33 = moderate Evidenz für keinen Effekt\n")
cat("Dazwischen = inconclusive\n\n")

df_bf <- as.data.frame(emi_main) %>%
  mutate(Participant = factor(Participant))

# Nullmodelle (ohne intervention_type)
cat("Berechne Nullmodelle...\n")
bf_null_c <- lmBF(craving_post_z ~ craving_pre_z + protocol_version,
                   data = df_bf, whichRandom = "Participant")
bf_null_l <- lmBF(loneliness_post_z ~ loneliness_pre_z + protocol_version,
                   data = df_bf, whichRandom = "Participant")

# Vollmodelle (mit intervention_type)
cat("Berechne Vollmodelle...\n")
bf_full_c <- lmBF(craving_post_z ~ craving_pre_z + protocol_version +
                    intervention_type,
                   data = df_bf, whichRandom = "Participant")
bf_full_l <- lmBF(loneliness_post_z ~ loneliness_pre_z + protocol_version +
                    intervention_type,
                   data = df_bf, whichRandom = "Participant")

# Bayes Faktoren für Interventionseffekt
cat("\n--- Bayes Faktor Craving (H1/H3) ---\n")
bf_craving <- bf_full_c / bf_null_c
print(bf_craving)

cat("\n--- Bayes Faktor Loneliness (H2/H4) ---\n")
bf_loneliness <- bf_full_l / bf_null_l
print(bf_loneliness)

# =============================================================================
# 6. ROBUSTNESS CHECK 1 – WORTANZAHL-SENSITIVITÄT
# =============================================================================

cat("\n\n### ROBUSTNESS CHECK 1: WORTANZAHL-SENSITIVITÄT ###\n\n")

m1_s6 <- lmer(craving_post ~ craving_pre + intervention_type + protocol_version +
                (1 | Participant), data = emi_sens6, REML = FALSE)
m2_s6 <- lmer(loneliness_direct_post ~ loneliness_direct_pre + intervention_type +
                protocol_version + (1 | Participant), data = emi_sens6, REML = FALSE)

m1_s8 <- lmer(craving_post ~ craving_pre + intervention_type + protocol_version +
                (1 | Participant), data = emi_sens8, REML = FALSE)
m2_s8 <- lmer(loneliness_direct_post ~ loneliness_direct_pre + intervention_type +
                protocol_version + (1 | Participant), data = emi_sens8, REML = FALSE)

cat("=== VERGLEICHSTABELLE WORTANZAHL-SENSITIVITÄT ===\n")
cat("Stabile β und Signifikanz über alle Cut-offs → Haupt-Cut-off ist robust\n\n")

bind_rows(
  koef_sens(m1_ri, "Craving",    ">= 4 (Haupt)"),
  koef_sens(m1_s6, "Craving",    ">= 6"),
  koef_sens(m1_s8, "Craving",    ">= 8"),
  koef_sens(m2_ri, "Loneliness", ">= 4 (Haupt)"),
  koef_sens(m2_s6, "Loneliness", ">= 6"),
  koef_sens(m2_s8, "Loneliness", ">= 8")
) %>%
  mutate(sig = case_when(
    p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "ns"
  )) %>%
  arrange(outcome, vergleich, cutoff) %>%
  print()

# =============================================================================
# 7. ROBUSTNESS CHECK 2 – PHASENVERGLEICH
# =============================================================================

cat("\n\n### ROBUSTNESS CHECK 2: PHASENVERGLEICH ###\n\n")

df_phase1 <- emi_main %>% filter(protocol_version == "old_8plus_prompts")
cat(sprintf("Phase 1 (old_8plus_prompts): %d Paare / %d Probanden\n",
            nrow(df_phase1), n_distinct(df_phase1$Participant)))

m1_phase1 <- lmer(craving_post ~ craving_pre + intervention_type +
                    (1 | Participant), data = df_phase1, REML = FALSE)
m2_phase1 <- lmer(loneliness_direct_post ~ loneliness_direct_pre + intervention_type +
                    (1 | Participant), data = df_phase1, REML = FALSE)

df_phase2 <- emi_main %>% filter(protocol_version == "new_4_prompts")
cat(sprintf("Phase 2 (new_4_prompts):     %d Paare / %d Probanden\n\n",
            nrow(df_phase2), n_distinct(df_phase2$Participant)))

m1_phase2 <- lmer(craving_post ~ craving_pre + intervention_type +
                    (1 | Participant), data = df_phase2, REML = FALSE)
m2_phase2 <- lmer(loneliness_direct_post ~ loneliness_direct_pre + intervention_type +
                    (1 | Participant), data = df_phase2, REML = FALSE)

cat("=== ÜBERSICHTSTABELLE PHASENVERGLEICH ===\n")
cat("Ähnliche β in beiden Phasen → Ergebnis ist robust\n\n")

bind_rows(
  koef_sens(m1_phase1, "Craving",    "old_8plus_prompts"),
  koef_sens(m1_phase2, "Craving",    "new_4_prompts"),
  koef_sens(m2_phase1, "Loneliness", "old_8plus_prompts"),
  koef_sens(m2_phase2, "Loneliness", "new_4_prompts")
) %>%
  rename(protokoll = cutoff) %>%
  mutate(sig = case_when(
    p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "ns"
  )) %>%
  arrange(outcome, vergleich, protokoll) %>%
  print()

# =============================================================================
# 7b. ROBUSTNESS CHECK 3 – INTENTION-TO-TREAT (ITT)
#
# WARUM: Korrektor 2 hat darauf hingewiesen dass trial-level Ausschlüsse
# systematischen Bias erzeugen könnten:
#   - Jemand macht keinen Sport WEIL Craving besonders hoch ist
#   - Jemand füllt Reappraisal oberflächlich aus wenn sehr distressed
#   - Diese Momente verschwinden aus der Per-Protocol Analyse
#   → Interventionen könnten besser aussehen als sie wirklich sind
#
# ITT = alle Paare einschließen unabhängig ob Intervention korrekt durchgeführt
# Kein pair_valid Filter → zeigt ob Ausschlüsse systematische Muster haben
#
# Interpretation:
#   ITT ≈ Per-Protocol → kein systematischer Bias durch Ausschlüsse ✓
#   ITT ≠ Per-Protocol → Ausschlüsse sind non-random → mit Vorsicht interpretieren
#
# HINWEIS: Paare mit fehlenden Pre/Post Werten (NA) werden auch hier ausgeschlossen
# weil die Werte schlicht nicht vorhanden sind. Das ist nicht adressierbar.
# =============================================================================

cat("\n\n### ROBUSTNESS CHECK 3: INTENTION-TO-TREAT (ITT) ###\n")
cat("Adressiert trial-level compliance bias (Kritik Korrektor 2)\n")
cat("Kein pair_valid Filter → alle Paare mit vorhandenen Werten\n\n")

emi_itt <- emi_pairs %>%
  filter(
    phase == "Intervention",
    aborted == FALSE | is.na(aborted),
    !is.na(craving_pre),
    !is.na(craving_post),
    !is.na(loneliness_direct_pre),
    !is.na(loneliness_direct_post),
    !(craving_pre == 1 & loneliness_direct_pre == 1)
    # pair_valid wird NICHT gefiltert → ITT
  ) %>%
  mutate(
    intervention_type = factor(intervention_type,
                                levels = c("Control", "Reappraisal", "Movement")),
    protocol_version  = factor(protocol_version)
  )

cat(sprintf("ITT Stichprobe: %d Paare / %d Probanden\n",
            nrow(emi_itt), n_distinct(emi_itt$Participant)))
cat(sprintf("Per-Protocol:   %d Paare / %d Probanden\n\n",
            nrow(emi_main), n_distinct(emi_main$Participant)))

m1_itt <- lmer(craving_post ~ craving_pre + intervention_type + protocol_version +
                 (1 | Participant), data = emi_itt, REML = FALSE)
m2_itt <- lmer(loneliness_direct_post ~ loneliness_direct_pre + intervention_type +
                 protocol_version + (1 | Participant), data = emi_itt, REML = FALSE)

cat("=== VERGLEICH ITT vs. PER-PROTOCOL ===\n")
cat("Ähnliche β → kein systematischer Bias durch Ausschlüsse\n\n")

bind_rows(
  koef_sens(m1_itt, "Craving",    "ITT"),
  koef_sens(m1_ri,  "Craving",    "Per-Protocol"),
  koef_sens(m2_itt, "Loneliness", "ITT"),
  koef_sens(m2_ri,  "Loneliness", "Per-Protocol")
) %>%
  rename(analyse = cutoff) %>%
  mutate(sig = case_when(
    p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "ns"
  )) %>%
  arrange(outcome, vergleich, analyse) %>%
  print()

# =============================================================================
# 8. VISUALISIERUNGEN
# =============================================================================

farben <- c(
  "Control"     = "#999999",
  "Reappraisal" = "#2C7BB6",
  "Movement"    = "#D7191C"
)

# -----------------------------------------------------------------------------
# FIGURE 1: PARTICIPANT FLOWCHART
#
# WARUM: Standard in klinischen Studien – zeigt transparent wie viele
# Probanden auf welcher Stufe ausgeschlossen wurden.
# Erstellt als einfacher Text-Plot mit ggplot2.
# -----------------------------------------------------------------------------

cat("\n\n### FIGURE 1: PARTICIPANT FLOWCHART ###\n")

flowchart_data <- tibble(
  label = c(
    paste0("Enrolled\nN = ", n_distinct(emi_pairs$Participant)),
    paste0("Dropouts excluded\nn = ", sum(compliance$aborted, na.rm = TRUE),
           "\n(< 70% compliance or < 11 days)"),
    paste0("Met compliance criteria\nn = ", sum(!compliance$aborted, na.rm = TRUE)),
    paste0("Additionally excluded\nn = 1\n(no valid pre-post pairs)"),
    paste0("Final analyzed sample\nN = ", n_distinct(emi_main$Participant))
  ),
  x = c(0.5, 0.85, 0.5, 0.85, 0.5),
  y = c(1.0,  0.75, 0.5,  0.25, 0.0),
  box = c(TRUE, FALSE, TRUE, FALSE, TRUE)
)

ggplot(flowchart_data) +
  geom_rect(data = filter(flowchart_data, box),
            aes(xmin = x - 0.22, xmax = x + 0.22,
                ymin = y - 0.08, ymax = y + 0.08),
            fill = "white", color = "black", linewidth = 0.8) +
  geom_text(aes(x = x, y = y, label = label),
            size = 3.5, hjust = 0.5, vjust = 0.5) +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.92, yend = 0.58,
           arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.42, yend = 0.08,
           arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("segment", x = 0.5, xend = 0.72, y = 0.75, yend = 0.75) +
  annotate("segment", x = 0.5, xend = 0.72, y = 0.25, yend = 0.25) +
  xlim(0.2, 1.1) + ylim(-0.15, 1.15) +
  labs(title = "Figure 1: Participant Flowchart") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))

# -----------------------------------------------------------------------------
# FIGURE 3: FOREST PLOT – PRIMÄRE ERGEBNISSE
#
# WARUM: Zeigt alle vier Kontraste gleichzeitig mit β und 95% CI.
# Vertikale Linie bei 0 = Nullhypothese.
# Punkte links = Verbesserung gegenüber Control.
#
# Die Daten kommen direkt aus den bereits berechneten standardisierten
# Modellen (m1_ri_std und m2_ri_std) – keine zusätzliche Berechnung nötig.
# -----------------------------------------------------------------------------

cat("\n\n### FIGURE 3: FOREST PLOT ###\n")

# Daten direkt aus den standardisierten KI-Tabellen zusammenstellen
# ci_std_m1 und ci_std_m2 wurden bereits in Sektion 5b berechnet
forest_data <- tibble(
  outcome  = c("Craving",    "Craving",
               "Loneliness", "Loneliness"),
  contrast = c("Reappraisal vs. Control", "Movement vs. Control",
               "Reappraisal vs. Control", "Movement vs. Control"),
  beta     = c(
    k_m1_std["intervention_typeReappraisal", "Estimate"],
    k_m1_std["intervention_typeMovement",    "Estimate"],
    k_m2_std["intervention_typeReappraisal", "Estimate"],
    k_m2_std["intervention_typeMovement",    "Estimate"]
  ),
  ci_low   = c(
    ci_std_m1["intervention_typeReappraisal", "2.5 %"],
    ci_std_m1["intervention_typeMovement",    "2.5 %"],
    ci_std_m2["intervention_typeReappraisal", "2.5 %"],
    ci_std_m2["intervention_typeMovement",    "2.5 %"]
  ),
  ci_high  = c(
    ci_std_m1["intervention_typeReappraisal", "97.5 %"],
    ci_std_m1["intervention_typeMovement",    "97.5 %"],
    ci_std_m2["intervention_typeReappraisal", "97.5 %"],
    ci_std_m2["intervention_typeMovement",    "97.5 %"]
  ),
  p        = c(
    k_m1_std["intervention_typeReappraisal", "Pr(>|t|)"],
    k_m1_std["intervention_typeMovement",    "Pr(>|t|)"],
    k_m2_std["intervention_typeReappraisal", "Pr(>|t|)"],
    k_m2_std["intervention_typeMovement",    "Pr(>|t|)"]
  )
) %>%
  mutate(
    label     = paste0(contrast, "\n(", outcome, ")"),
    sig_label = ifelse(p < .05, "p < .05", "ns")
  )

ggplot(forest_data, aes(x = beta, y = label, color = outcome)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.25, linewidth = 0.9) +
  geom_point(size = 4) +
  geom_text(aes(x = ci_high + 0.02, label = sig_label),
            hjust = 0, size = 3.5, color = "black") +
  scale_color_manual(values = c("Craving" = "#2C7BB6", "Loneliness" = "#D7191C")) +
  labs(
    title    = "Figure 3: Forest Plot – Primary Results",
    subtitle = "Standardized β with 95% CI | Dashed line = null hypothesis",
    x        = "Standardized β (negative = better than Control)",
    y        = NULL,
    color    = "Outcome"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    axis.text.y     = element_text(size = 11),
    legend.position = "right"
  )

# -----------------------------------------------------------------------------
# FIGURE 4: COHEN'S D PLOT – WITHIN-CONDITION EFFECT SIZES
#
# WARUM: Zeigt ob jede Bedingung isoliert Pre-Post Veränderungen hatte.
# Wichtig für Control: Hat das Memory-Spiel selbst Craving/Loneliness reduziert?
# Nur deskriptiv – kein Hypothesentest.
#
# Die Cohen's d Werte wurden bereits in Sektion 2 berechnet.
# Hier werden sie direkt als Tabelle eingetragen – einfach und direkt.
# -----------------------------------------------------------------------------

cat("\n\n### FIGURE 4: COHEN'S D PLOT ###\n")

# Cohen's d Werte aus Sektion 2 direkt eingetragen
# d_reapp_c, d_reapp_l, d_move_c, d_move_l, d_ctrl_c, d_ctrl_l
# wurden bereits berechnet – hier werden sie nochmal explizit abgerufen

sub_reapp <- emi_main %>% filter(intervention_type == "Reappraisal")
sub_move  <- emi_main %>% filter(intervention_type == "Movement")
sub_ctrl  <- emi_main %>% filter(intervention_type == "Control")

d_rc <- cohens_d(sub_reapp$craving_post,            sub_reapp$craving_pre,            paired = TRUE)
d_rl <- cohens_d(sub_reapp$loneliness_direct_post,  sub_reapp$loneliness_direct_pre,  paired = TRUE)
d_mc <- cohens_d(sub_move$craving_post,             sub_move$craving_pre,             paired = TRUE)
d_ml <- cohens_d(sub_move$loneliness_direct_post,   sub_move$loneliness_direct_pre,   paired = TRUE)
d_cc <- cohens_d(sub_ctrl$craving_post,             sub_ctrl$craving_pre,             paired = TRUE)
d_cl <- cohens_d(sub_ctrl$loneliness_direct_post,   sub_ctrl$loneliness_direct_pre,   paired = TRUE)

cohens_d_data <- tibble(
  intervention = c("Reappraisal", "Reappraisal",
                   "Movement",    "Movement",
                   "Control",     "Control"),
  outcome      = c("Craving", "Loneliness",
                   "Craving", "Loneliness",
                   "Craving", "Loneliness"),
  d            = c(d_rc$Cohens_d, d_rl$Cohens_d,
                   d_mc$Cohens_d, d_ml$Cohens_d,
                   d_cc$Cohens_d, d_cl$Cohens_d),
  ci_low       = c(d_rc$CI_low,   d_rl$CI_low,
                   d_mc$CI_low,   d_ml$CI_low,
                   d_cc$CI_low,   d_cl$CI_low),
  ci_high      = c(d_rc$CI_high,  d_rl$CI_high,
                   d_mc$CI_high,  d_ml$CI_high,
                   d_cc$CI_high,  d_cl$CI_high)
) %>%
  mutate(intervention = factor(intervention,
                                levels = c("Control", "Reappraisal", "Movement")))

ggplot(cohens_d_data,
       aes(x = intervention, y = d, color = intervention, shape = outcome)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.2, linewidth = 0.9,
                position = position_dodge(width = 0.4)) +
  geom_point(size = 4, position = position_dodge(width = 0.4)) +
  scale_color_manual(values = farben) +
  scale_shape_manual(values = c("Craving" = 16, "Loneliness" = 17)) +
  labs(
    title    = "Figure 4: Within-Condition Cohen's d (Descriptive)",
    subtitle = "Pre-to-post change within each intervention type | 95% CI",
    x        = NULL,
    y        = "Cohen's d",
    color    = "Intervention",
    shape    = "Outcome"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# -----------------------------------------------------------------------------
# FIGURE 2: PRE-POST PLOTS
# -----------------------------------------------------------------------------

cat("\n\n### FIGURE 2: PRE-POST PLOTS ###\n")

plot_pre_post <- function(data, outcome_post, outcome_pre, titel) {
  df_long <- data %>%
    select(Participant,
           Pre  = all_of(outcome_pre),
           Post = all_of(outcome_post),
           intervention_type) %>%
    filter(!is.na(Pre), !is.na(Post)) %>%
    pivot_longer(c(Pre, Post), names_to = "Zeitpunkt", values_to = "Wert") %>%
    mutate(Zeitpunkt = factor(Zeitpunkt, levels = c("Pre", "Post")))

  mittelwerte <- df_long %>%
    group_by(intervention_type, Zeitpunkt) %>%
    summarise(M = mean(Wert, na.rm = TRUE), .groups = "drop")

  ggplot(df_long,
         aes(x = Zeitpunkt, y = Wert, color = intervention_type,
             group = interaction(Participant, intervention_type))) +
    geom_line(alpha = 0.2, linewidth = 0.4) +
    geom_point(alpha = 0.2, size = 1.2) +
    geom_line(data = mittelwerte,
              aes(x = Zeitpunkt, y = M, group = intervention_type,
                  color = intervention_type),
              linewidth = 2, inherit.aes = FALSE) +
    geom_point(data = mittelwerte,
               aes(x = Zeitpunkt, y = M, color = intervention_type),
               size = 3.5, inherit.aes = FALSE) +
    scale_color_manual(values = farben) +
    facet_wrap(~ intervention_type) +
    labs(title    = titel,
         subtitle = "Dünne Linien = individuelle Paare | Dicke Linien = Mittelwert",
         x = "Zeitpunkt", y = "Wert (1–7)") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none")
}

print(plot_pre_post(emi_main, "craving_post", "craving_pre",
                    "Craving: Pre vs. Post"))
print(plot_pre_post(emi_main, "loneliness_direct_post", "loneliness_direct_pre",
                    "Loneliness: Pre vs. Post"))

# =============================================================================
# 9. DESKRIPTIV: WITHIN-PERSON ZUSAMMENHANG CRAVING & LONELINESS
#
# Fragestellung: Treten Momente hoher Einsamkeit mit hohem Craving zusammen auf?
# Methode: Person-Mean-Centering → Abweichung vom eigenen Mittelwert
# Basis: EMI_Tagesabfrage UND EMI_Morgenabfrage in der Interventionsphase
# KEIN Hypothesentest – rein deskriptiv, präregistriert
# =============================================================================

cat("\n\n### DESKRIPTIV: WITHIN-PERSON CRAVING & LONELINESS ###\n")
cat("(Kein Hypothesentest – rein deskriptiv, präregistriert)\n\n")

within_data <- data_clean %>%
  filter(Form %in% c("EMI_Tagesabfrage", "EMI_Morgenabfrage")) %>%
  filter(phase == "Intervention") %>%
  filter(aborted == FALSE | is.na(aborted)) %>%
  filter(Participant %in% unique(emi_main$Participant)) %>%
  mutate(
    craving_pre_wp    = if_else(Form == "EMI_Tagesabfrage",
                                as.numeric(Craving_EMI_C),
                                as.numeric(Craving_EMI_morgen_C)),
    loneliness_pre_wp = if_else(Form == "EMI_Tagesabfrage",
                                as.numeric(directLoneliness_EMI_C),
                                as.numeric(directLoneliness_EMI_morgen_C))
  ) %>%
  filter(!is.na(craving_pre_wp), !is.na(loneliness_pre_wp)) %>%
  group_by(Participant) %>%
  mutate(
    craving_within    = craving_pre_wp    - mean(craving_pre_wp,    na.rm = TRUE),
    loneliness_within = loneliness_pre_wp - mean(loneliness_pre_wp, na.rm = TRUE)
  ) %>%
  ungroup()

cat("N Beobachtungen:", nrow(within_data), "\n")
cat("N Probanden:",     n_distinct(within_data$Participant), "\n\n")

mod_wp <- lmer(craving_within ~ loneliness_within + (1 | Participant),
               data = within_data, REML = FALSE)
print(summary(mod_wp))

k_wp <- coef(summary(mod_wp))
cat(sprintf("\nβ = %.3f, SE = %.3f, t = %.3f, p = %.3f\n",
            k_wp["loneliness_within", "Estimate"],
            k_wp["loneliness_within", "Std. Error"],
            k_wp["loneliness_within", "t value"],
            k_wp["loneliness_within", "Pr(>|t|)"]))

ggplot(within_data, aes(x = loneliness_within, y = craving_within)) +
  geom_point(alpha = 0.25, color = "#2C7BB6", size = 1.5) +
  geom_smooth(method = "lm", color = "#D7191C", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title    = "Within-Person: Einsamkeit & Craving",
    subtitle = "Abweichung vom Personenmittelwert | Rein deskriptiv",
    x        = "Einsamkeit (within-person zentriert)",
    y        = "Craving (within-person zentriert)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# =============================================================================
# 10. SESSION INFO
# =============================================================================

cat("\n\n=== SESSION INFO ===\n")
sessionInfo()
