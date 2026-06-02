
# Charité Masterarbeit – Cannabis, Einsamkeit & Craving
#
#   Rohdaten von movisens (Excel-Dateien pro Proband) einlesen,
#   bereinigen und Pre-Post-Paare für die Interventionsanalyse erstellen.
#
# ERGEBNIS:
#   EMA_cleaned_data.RData mit drei fertigen Datensätzen:
#   - emi_pairs        → Hauptanalyse (Reappraisal >= 4 Wörter pro Komponente)
#   - emi_pairs_sens6  → Sensitivitätsanalyse (>= 6 Wörter)
#   - emi_pairs_sens8  → Sensitivitätsanalyse (>= 8 Wörter)
#
# PRE-FORMEN (Abfrage direkt VOR der Intervention):
#   EMI_Tagesabfrage  → Craving_EMI_C          / directLoneliness_EMI_C
#   EMI_Morgenabfrage → Craving_EMI_morgen_C   / directLoneliness_EMI_morgen_C
#
# POST-FORMEN (Abfrage direkt NACH der Intervention):
#   EMI_PostCraving    → Craving_PostEMI_C
#   EMI_PostLoneliness → directLoneliness_PostEMI_C
#
# ABLAUF PRO INTERVENTION:
#   Control:     Pre → Control_Task (+1) → Post
#   Reappraisal: Pre → Direct Reappraisal → Start → Items → Ende → Feedback (+5) → Post
#   Bewegung:    Pre → PA_Startitem → [variabel] → PA_Kontrollfrage (+1) → Post
#
# ABBRUCHKRITERIUM (präregistriert):
#   Abbrecher = compliance_rate < 70% ODER tage_dabei < 11
###############################################################################

library(tidyverse)   # filter, mutate, group_by, summarise etc.
library(lubridate)   # Datum/Uhrzeit-Funktionen (ymd_hms, as_date)
library(readxl)      # Excel-Dateien einlesen

# =============================================================================

data_path <- "/Users/juliabrandl/Desktop/MA Datensätze alle bisher"
save_path <- "/Users/juliabrandl/Desktop/EMA_cleaned_data.RData"

# =============================================================================
# 1. ALLE EXCEL-DATEIEN EINLESEN UND ZUSAMMENFÜGEN
#
# WAS PASSIERT:
#   list.files() findet alle .xlsx-Dateien im Ordner.
#   Die for-Schleife liest jede Datei ein und speichert sie in einer Liste.
#   bind_rows() klebt alle Tabellen untereinander zusammen.
#
# col_types = "text": alles als Text einlesen damit R beim Import keine
# Zahlen oder Daten falsch interpretiert (z.B. Datum als Zahl).
# Wir korrigieren die Datentypen danach selbst.
# =============================================================================

alle_dateien <- list.files(
  path       = data_path,
  pattern    = "\\.xlsx$",   # nur Dateien die auf .xlsx enden
  full.names = TRUE           # vollständigen Pfad zurückgeben, nicht nur Dateiname
)

cat("Gefundene Dateien:", length(alle_dateien), "\n")

alle_tabellen <- list()
for (datei in alle_dateien) {
  cat("Lese:", basename(datei), "\n")
  alle_tabellen[[datei]] <- read_excel(datei, col_types = "text")
}
data_raw <- bind_rows(alle_tabellen)

cat("\nGesamtzeilen:", nrow(data_raw), "\n")
cat("Probanden:",    n_distinct(data_raw$Participant), "\n")

# =============================================================================
# 2. BASISBEREINIGUNG
#
# WARUM: Die Rohdaten haben Leerzeichen in Spaltennamen und alle Werte
# sind noch Text. Wir müssen Datentypen korrigieren und den Studientag
# berechnen bevor wir weiterarbeiten können.
#
# WAS PASSIERT:
#   1. Leerzeichen in Spaltennamen → Unterstriche
#   2. Participant als Text, ".0" entfernen (Excel speichert "17" als "17.0")
#   3. Zeitstempel einlesen: Excel speichert Daten als Zahlen (Tage seit
#      30.12.1899) → * 86400 rechnet Tage in Sekunden um → as.POSIXct()
#      macht daraus ein echtes Datum
#   4. date_only: nur Datum ohne Uhrzeit – brauchen wir für Studientag
#   5. study_start: frühestes Datum pro Proband = Studientag 1
#   6. study_day: Tage seit Studienstart (+1 damit Tag 1 nicht Tag 0 ist)
#   7. arrange(): nach Abschlusszeit sortieren – SEHR WICHTIG für Pre-Post-Logik!
# =============================================================================

data_raw <- data_raw %>%
  rename_with(~ str_replace_all(.x, " ", "_"))

data_clean <- data_raw %>%
  mutate(Participant      = str_remove(as.character(Participant), "\\.0$")) %>%
  mutate(form_finish      = as.POSIXct(as.numeric(Form_finish_date) * 86400,
                                       origin = "1899-12-30", tz = "UTC")) %>%
  mutate(trigger_datetime = as.POSIXct(as.numeric(Trigger_date) * 86400,
                                       origin = "1899-12-30", tz = "UTC")) %>%
  mutate(date_only        = as_date(trigger_datetime)) %>%
  group_by(Participant) %>%
  mutate(
    study_start = min(date_only, na.rm = TRUE),
    study_day   = as.numeric(date_only - study_start) + 1
  ) %>%
  ungroup() %>%
  arrange(Participant, form_finish)

# =============================================================================
# 3. STUDIENPHASEN ZUWEISEN
#
# Phasenzuweisung um später nur die Interventionsphase
# zu analysieren und Compliance korrekt zu berechnen.
#
# WAS PASSIERT:
#   case_when() ist wie ein if/else für mehrere Bedingungen gleichzeitig.
#   Jede Zeile bekommt eine Phase zugewiesen basierend auf dem Studientag.
#   "Other" ist der Auffangfall für alles was nicht in die drei Phasen passt.
#
#   Tag  1-7:  Pre-Phase          (Woche 1, nur EMA, keine Interventionen)
#   Tag  8-28: Interventionsphase (Wochen 2-4, hier passieren die JITAIs)
#   Tag 29-35: Post-Phase         (Woche 5, wieder nur EMA)
# =============================================================================

data_clean <- data_clean %>%
  mutate(phase = case_when(
    study_day >= 1  & study_day <= 7  ~ "Pre",
    study_day >= 8  & study_day <= 28 ~ "Intervention",
    study_day >= 29 & study_day <= 35 ~ "Post",
    TRUE                              ~ "Other"
  ))

cat("\nZeilen pro Studienphase:\n")
print(table(data_clean$phase))

# =============================================================================
# 4. COMPLIANCE UND ABBRECHER
#
# Probanden mit zu wenig Teilnahme liefern nicht genug Daten für
# die Analyse und werden ausgeschlossen.
#
# WAS PASSIERT:
# zwei Kriterien:
#   1. compliance_rate = Anteil beantworteter Random-Time-Prompts
#      Beantwortet = Missing-Spalte ist leer (NA oder "")
#   2. tage_dabei = wie viele Tage war der Proband in der Interventionsphase
#
#   Abbrecher (aborted = TRUE) wenn:
#     compliance_rate < 70% (weniger als 70% der Prompts beantwortet)
#     ODER tage_dabei < 11  (weniger als 50% der 21 Interventionstage dabei)
#
#   Probanden die gar keine Interventionsdaten haben bekommen
#   compliance_rate = 0 und tage_dabei = 0 → automatisch Abbrecher
# =============================================================================

# Wie viele Tage war jeder Proband in der Interventionsphase?
tage_intervention <- data_clean %>%
  filter(phase == "Intervention") %>%
  group_by(Participant) %>%
  summarise(
    tage_dabei = n_distinct(date_only),
    .groups = "drop"
  )

# Compliance = Anteil beantworteter Random-Time-Prompts
# Nur Interventionsphase, nur zufällige Prompts (nicht fixe Abendabfragen)
compliance_raw <- data_clean %>%
  filter(phase == "Intervention") %>%
  filter(str_detect(Trigger, "Random Time")) %>%
  group_by(Participant) %>%
  summarise(
    n_prompts_total = n(),
    n_responded     = sum(is.na(Missing) | Missing == "", na.rm = TRUE),
    compliance_rate = n_responded / n_prompts_total,
    .groups = "drop"
  )

# Alle Probanden zusammenführen und Abbrecher markieren
compliance <- data_clean %>%
  distinct(Participant) %>%
  left_join(compliance_raw,    by = "Participant") %>%
  left_join(tage_intervention, by = "Participant") %>%
  mutate(
    n_prompts_total = replace_na(n_prompts_total, 0),
    n_responded     = replace_na(n_responded,     0),
    compliance_rate = replace_na(compliance_rate, 0),
    tage_dabei      = replace_na(tage_dabei,      0),
    aborted         = compliance_rate < 0.7 | tage_dabei < 11
  )

cat("\n--- Compliance-Übersicht ---\n")
print(compliance)
cat("\nAbbrecher gesamt:", sum(compliance$aborted), "\n")
cat("  davon < 70% Prompts:", sum(compliance$compliance_rate < 0.7), "\n")
cat("  davon < 11 Tage:",     sum(compliance$tage_dabei < 11), "\n")

data_clean <- data_clean %>%
  left_join(compliance, by = "Participant")

# =============================================================================
# 5. PROTOKOLLVERSION ERKENNEN
#
#   old_8plus_prompts: Phase 1, bis zu 8 Prompts/Tag
#                      Control = Memory-Spiel + FruitTapGame
#                      Probanden: 6, 7, 9, 10, 11, 12, 13, 14, 16, 17,
#                                 18, 19, 20, 21, 24
#   new_4_prompts:     Phase 2, bis zu 4 Prompts/Tag
#                      Control = nur Memory-Spiel
#                      Alle anderen Probanden
#
# WAS PASSIERT:
#   Protokollversion direkt aus bekannten Probandennummern zuweisen
# =============================================================================

phase1_probanden <- c("6", "7", "9", "10", "11", "12",
                      "13", "14", "16", "17", "18", "19",
                      "20", "21", "24")

data_clean <- data_clean %>%
  mutate(
    protocol_version = if_else(
      Participant %in% phase1_probanden,
      "old_8plus_prompts",
      "new_4_prompts"
    )
  )

cat("\nProtokollversionen:\n")
print(table(data_clean$protocol_version))

# =============================================================================
# 6. REAPPRAISAL-QUALITÄT PRÜFEN (WORTANZAHL)
#
# GÜLTIG wenn ALLE DREI Komponenten den Schwellenwert erfüllen:
#   Situation:  Direct_Loneliness_C ODER Direct_Craving_C
#   Gedanke:    Direct_Loneliness_Trigger_C ODER Direct_Craving_Trigger_C
#   Strategie:  Maximum aus den vier Strategie-Spalten
#
# Drei Schwellenwerte für Sensitivitätsanalyse:
#   reapp_valid_4 = Haupt-Cut-off (>= 4 Wörter pro Komponente, präregistriert)
#   reapp_valid_6 = Sensitivität  (>= 6 Wörter, strenger)
#   reapp_valid_8 = Sensitivität  (>= 8 Wörter, noch strenger)
# =============================================================================

# Hilfsfunktion: Wörter zählen
# Gibt 0 zurück wenn Text leer, NA oder "NA" als Text ist
zaehle_woerter <- function(text) {
  if (is.na(text) || str_squish(text) == "" || text == "NA") return(0L)
  str_count(str_squish(as.character(text)), boundary("word"))
}

# Wortanzahlen pro Zeile berechnen
# sapply() wendet zaehle_woerter() auf jede Zeile einer Spalte an
data_clean <- data_clean %>%
  mutate(
    wc_situation = pmax(
      sapply(Direct_Loneliness_C,       zaehle_woerter),
      sapply(Direct_Craving_C,          zaehle_woerter)
    ),
    wc_gedanke = pmax(
      sapply(Direct_Loneliness_Trigger_C, zaehle_woerter),
      sapply(Direct_Craving_Trigger_C,    zaehle_woerter)
    ),
    # Proband wählt eine von vier Strategien – nur eine ist befüllt
    wc_strategie_max = pmax(
      sapply(Eigene_Erfahrungen_Antwort_C,           zaehle_woerter),
      sapply(Welche_Nachteile_Antwort_C,             zaehle_woerter),
      sapply(Welche_anderen_Anhaltspunkte_Antwort_C, zaehle_woerter),
      sapply(Welchen_Rat_Antwort_C,                  zaehle_woerter)
    )
  )

# Pro Reappraisal-Sequenz aggregieren
# Eine Sequenz = alle Zeilen mit demselben Trigger-Zeitstempel
# Maximum nehmen weil die Wörter über mehrere Formulare verteilt sind
reapp_wordcounts <- data_clean %>%
  filter(Form %in% c("Direct Loneliness Reappraisal",
                     "Direct Craving Reappraisal",
                     "Reappraisal Items")) %>%
  group_by(Participant, trigger_datetime) %>%
  summarise(
    seq_wc_situation = max(wc_situation,     na.rm = TRUE),
    seq_wc_gedanke   = max(wc_gedanke,       na.rm = TRUE),
    seq_wc_strategie = max(wc_strategie_max, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Inf-Werte durch 0 ersetzen (entsteht wenn alle Werte NA waren)
  mutate(
    seq_wc_situation = if_else(is.infinite(seq_wc_situation), 0L, as.integer(seq_wc_situation)),
    seq_wc_gedanke   = if_else(is.infinite(seq_wc_gedanke),   0L, as.integer(seq_wc_gedanke)),
    seq_wc_strategie = if_else(is.infinite(seq_wc_strategie), 0L, as.integer(seq_wc_strategie))
  ) %>%
  # TRUE = alle drei Komponenten erfüllen den Schwellenwert
  mutate(
    reapp_valid_4 = seq_wc_situation >= 4 & seq_wc_gedanke >= 4 & seq_wc_strategie >= 4,
    reapp_valid_6 = seq_wc_situation >= 6 & seq_wc_gedanke >= 6 & seq_wc_strategie >= 6,
    reapp_valid_8 = seq_wc_situation >= 8 & seq_wc_gedanke >= 8 & seq_wc_strategie >= 8
  )

cat("\n=== REAPPRAISAL WORTANZAHL ===\n")
cat("Gültig >= 4 Wörter:", sum(reapp_wordcounts$reapp_valid_4, na.rm = TRUE),
    "von", nrow(reapp_wordcounts), "Sequenzen\n")
cat("Gültig >= 6 Wörter:", sum(reapp_wordcounts$reapp_valid_6, na.rm = TRUE),
    "von", nrow(reapp_wordcounts), "Sequenzen\n")
cat("Gültig >= 8 Wörter:", sum(reapp_wordcounts$reapp_valid_8, na.rm = TRUE),
    "von", nrow(reapp_wordcounts), "Sequenzen\n")

data_clean <- data_clean %>%
  left_join(reapp_wordcounts, by = c("Participant", "trigger_datetime"))

# =============================================================================
# 7. PRE-POST-PAARE ERSTELLEN
#
# WARUM: Für die Hauptanalyse brauchen wir pro Intervention genau einen
# Pre-Wert (direkt davor) und einen Post-Wert (direkt danach).
#
# WAS PASSIERT:
#   Die Funktion build_emi_pairs() geht für jeden Probanden alle Zeilen durch.
#   Wenn sie eine Interventionszeile findet, sucht sie:
#     PRE:  die Zeile direkt davor (i-1) muss EMI_Tagesabfrage oder
#           EMI_Morgenabfrage sein. Wenn nicht → Paar ausgeschlossen
#     POST: je nach Interventionstyp an fixer oder variabler Position danach
#
# POST-INDEX je nach Interventionstyp:
#   Control:     +1 Zeile nach Control_Task → EMI_PostCraving
#   Reappraisal: +5 Zeilen → EMI_PostCraving
#                (Start, Items, Ende, Feedback = 4 Formulare dazwischen)
#   Bewegung:    variabel – vorwärts suchen bis PA_Kontrollfrage, dann +1
#
# GÜLTIGKEITSPRÜFUNG:
#   Reappraisal: alle 3 Wortanzahl-Komponenten >= Schwellenwert
#   Bewegung:    Kontrollfrage_PA_C == 1 (hat Sport gemacht)
#   Control:     WorkingMemoryGame_C nicht abgebrochen
#
# TECHNISCHE FORMULARE werden vorher entfernt damit die Positions-Logik
# (i-1, i+5 etc.) nicht verfälscht wird.
# =============================================================================

pre_formen <- c("EMI_Tagesabfrage", "EMI_Morgenabfrage")

build_emi_pairs <- function(data, cutoff_col) {

  technische_forms <- c(
    "GeoCoding", "Time Setting", "SST (FruitTapGame)", "Abendabfrage",
    "Reflexion_Abendabfrage", "ERQ_Momentary", "Cannabis_Tagesabfrage",
    "Cannabis_Morgenabfrage", "Zwischenziel_verpasst", "Aktivitätsziel_verpasst",
    "PA_Startbenachrichtigung", "WarnungSensor50%", "WarnungSensor30%",
    "ERQ_Abend", "Missing"
  )

  ergebnis_liste <- list()

  n_missing_pre   <- 0
  n_missing_post  <- 0
  n_invalid_reapp <- 0
  n_invalid_pa    <- 0
  n_invalid_ctrl  <- 0

  for (vp in unique(data$Participant)) {

    d <- data %>%
      filter(Participant == vp) %>%
      filter(!Form %in% technische_forms) %>%
      filter(is.na(Missing) | Missing %in% c("", "Incomplete")) %>%
      arrange(form_finish)

    for (i in seq_len(nrow(d))) {

      if (!d$Form[i] %in% c("Control_Task",
                             "Direct Craving Reappraisal",
                             "Direct Loneliness Reappraisal",
                             "PA_Startitem")) next

      intervention_typ <- case_when(
        d$Form[i] == "Control_Task"                       ~ "Control",
        d$Form[i] %in% c("Direct Craving Reappraisal",
                          "Direct Loneliness Reappraisal") ~ "Reappraisal",
        d$Form[i] == "PA_Startitem"                       ~ "Movement"
      )

      # PRE-WERT HOLEN
      # Die Zeile direkt VOR der Intervention muss eine gültige Pre-Form sein.
      # EMI_Tagesabfrage und EMI_Morgenabfrage haben verschiedene Spaltennamen!
      if (i == 1 || !d$Form[i-1] %in% pre_formen) {
        n_missing_pre  <- n_missing_pre + 1
        pre_craving    <- NA_real_
        pre_loneliness <- NA_real_
      } else {
        if (d$Form[i-1] == "EMI_Tagesabfrage") {
          pre_craving    <- as.numeric(d$Craving_EMI_C[i-1])
          pre_loneliness <- as.numeric(d$directLoneliness_EMI_C[i-1])
        } else {
          pre_craving    <- as.numeric(d$Craving_EMI_morgen_C[i-1])
          pre_loneliness <- as.numeric(d$directLoneliness_EMI_morgen_C[i-1])
        }
      }

      # POST-INDEX BESTIMMEN
      # Control: immer +1, Reappraisal: immer +5, Bewegung: variabel
      if (intervention_typ == "Control") {
        index_post_craving    <- 1L
        index_post_loneliness <- 2L
      } else if (intervention_typ == "Reappraisal") {
        index_post_craving    <- 5L
        index_post_loneliness <- 6L
      } else {
        # Bewegung: vorwärts suchen bis PA_Kontrollfrage (max. 10 Zeilen)
        idx_kf <- NA_integer_
        for (k in seq(i + 1, min(i + 10, nrow(d)))) {
          if (!is.na(d$Form[k]) && d$Form[k] == "PA_Kontrollfrage") {
            idx_kf <- k - i
            break
          }
        }
        if (is.na(idx_kf)) {
          index_post_craving    <- NA_integer_
          index_post_loneliness <- NA_integer_
        } else {
          index_post_craving    <- as.integer(idx_kf + 1L)
          index_post_loneliness <- as.integer(idx_kf + 2L)
        }
      }

      # POST-WERT HOLEN
      # Drei Bedingungen: Index nicht NA, innerhalb Datensatz, richtige Form
      if (is.na(index_post_craving) ||
          i + index_post_craving > nrow(d) ||
          d$Form[i + index_post_craving] != "EMI_PostCraving") {
        n_missing_post  <- n_missing_post + 1
        post_craving    <- NA_real_
        post_loneliness <- NA_real_
      } else {
        post_craving    <- as.numeric(d$Craving_PostEMI_C[i + index_post_craving])
        post_loneliness <- as.numeric(d$directLoneliness_PostEMI_C[i + index_post_loneliness])
      }

      # GÜLTIGKEITSPRÜFUNG
      # Reappraisal: Wortanzahl, Bewegung: Sport gemacht, Control: Spiel gespielt
      if (intervention_typ == "Reappraisal") {
        valid_col  <- d[[cutoff_col]][i]
        ist_gueltig <- !is.na(valid_col) && valid_col == TRUE
        if (!ist_gueltig) n_invalid_reapp <- n_invalid_reapp + 1

      } else if (intervention_typ == "Movement") {
        if (is.na(index_post_craving)) {
          ist_gueltig <- FALSE
          n_invalid_pa <- n_invalid_pa + 1
        } else {
          idx_kf      <- index_post_craving - 1L
          pa_val      <- as.numeric(d$Kontrollfrage_PA_C[i + idx_kf])
          ist_gueltig <- !is.na(pa_val) && pa_val == 1
          if (!ist_gueltig) n_invalid_pa <- n_invalid_pa + 1
        }

      } else {
        wm_val      <- d$WorkingMemoryGame_C[i]
        ist_gueltig <- !is.na(wm_val) && wm_val != '{"canceled": true}'
        if (!ist_gueltig) n_invalid_ctrl <- n_invalid_ctrl + 1
      }

      ergebnis_liste[[length(ergebnis_liste) + 1]] <- list(
        Participant            = vp,
        trigger_datetime       = as.character(d$trigger_datetime[i]),
        study_day              = d$study_day[i],
        phase                  = d$phase[i],
        intervention_type      = intervention_typ,
        protocol_version       = d$protocol_version[i],
        pre_form               = if (i > 1 && d$Form[i-1] %in% pre_formen)
                                   d$Form[i-1] else NA_character_,
        craving_pre            = pre_craving,
        loneliness_direct_pre  = pre_loneliness,
        craving_post           = post_craving,
        loneliness_direct_post = post_loneliness,
        pair_valid             = ist_gueltig,
        wc_situation           = d$seq_wc_situation[i],
        wc_gedanke             = d$seq_wc_gedanke[i],
        wc_strategie           = d$seq_wc_strategie[i]
      )
    }
  }

  cat("\n--- Diagnose (Cut-off:", cutoff_col, ") ---\n")
  cat("Fehlende Pre-Formen:             ", n_missing_pre,   "\n")
  cat("Fehlende Post-Formulare:         ", n_missing_post,  "\n")
  cat("Ungültige Reappraisal-Antworten: ", n_invalid_reapp, "\n")
  cat("Bewegung nicht gemacht:          ", n_invalid_pa,    "\n")
  cat("Control abgebrochen:             ", n_invalid_ctrl,  "\n")

  pairs <- bind_rows(ergebnis_liste) %>%
    mutate(
      delta_craving           = craving_post           - craving_pre,
      delta_loneliness_direct = loneliness_direct_post - loneliness_direct_pre
    ) %>%
    left_join(
      compliance %>% select(Participant, compliance_rate, aborted),
      by = "Participant"
    )

  return(pairs)
}

# =============================================================================
# 8. DREI DATENSÄTZE ERSTELLEN
#
# Wir brauchen drei Datensätze mit unterschiedlichen Reappraisal-
# Cut-offs für die Sensitivitätsanalyse. So können wir zeigen dass die
# Ergebnisse nicht davon abhängen welchen Cut-off wir gewählt haben.
# =============================================================================

cat("\n\n=== ERSTELLE PRE-POST-PAARE ===\n")

cat("\nHaupt-Cut-off >= 4 Wörter...\n")
emi_pairs       <- build_emi_pairs(data_clean, "reapp_valid_4")

cat("\nSensitivitäts-Cut-off >= 6 Wörter...\n")
emi_pairs_sens6 <- build_emi_pairs(data_clean, "reapp_valid_6")

cat("\nSensitivitäts-Cut-off >= 8 Wörter...\n")
emi_pairs_sens8 <- build_emi_pairs(data_clean, "reapp_valid_8")

# =============================================================================
# 9. ÜBERSICHT DER ERGEBNISSE
#
# Plausibilitätsprüfung – wie viele Paare haben wir pro Bedingung?
# Wie viele Probanden sind Abbrecher? Woher kommen die Pre-Messungen?
# =============================================================================

cat("\n\n=== ÜBERSICHT PRE-POST-PAARE ===\n")
cat("Probanden in Hauptanalyse (nicht abgebrochen):", n_distinct(emi_pairs$Participant[!emi_pairs$aborted]), "\n")

# Gültige Paare pro Interventionstyp und Cut-off
# nrow() würde hier irreführend gleiche Zahlen zeigen weil alle drei Datensätze
# gleich viele Zeilen haben – pair_valid entscheidet welche Paare wirklich eingehen
cat("\nGültige Paare pro Interventionstyp und Wortanzahl-Cut-off:\n")
cat("(Nur pair_valid == TRUE, nur nicht-Abbrecher)\n\n")

for (cutoff_label in c(">= 4 (Haupt)", ">= 6", ">= 8")) {

  df <- switch(cutoff_label,
    ">= 4 (Haupt)" = emi_pairs,
    ">= 6"         = emi_pairs_sens6,
    ">= 8"         = emi_pairs_sens8
  )

  zahlen <- df %>%
    filter(pair_valid == TRUE, aborted == FALSE) %>%
    count(intervention_type)

  cat(sprintf("Cut-off %s:\n", cutoff_label))
  print(zahlen)
  cat("\n")
}

cat("\nHerkunft der Pre-Messungen:\n")
emi_pairs %>%
  count(pre_form) %>%
  print()

cat("\nAbbrecher:\n")
emi_pairs %>%
  distinct(Participant, aborted) %>%
  count(aborted) %>%
  print()

cat("\nVorschau emi_pairs (erste 5 Zeilen):\n")
emi_pairs %>%
  select(Participant, study_day, intervention_type, pre_form,
         craving_pre, craving_post, delta_craving,
         loneliness_direct_pre, loneliness_direct_post, pair_valid) %>%
  head(5) %>%
  print()

# =============================================================================
# 10. DATEN SPEICHERN
#
# Analyse-Skript lädt diese eine Datei mit load() und hat
# sofort alle Datensätze parat – kein erneutes Einlesen der Rohdaten nötig.
#
# Gespeicherte Objekte:
#   emi_pairs        – Hauptanalyse, Cut-off >= 4 Wörter
#   emi_pairs_sens6  – Sensitivität, Cut-off >= 6 Wörter
#   emi_pairs_sens8  – Sensitivität, Cut-off >= 8 Wörter
#   data_clean       – vollständiger gecleanter Rohdatensatz
#   compliance       – Compliance-Übersicht pro Proband
# =============================================================================

save(
  emi_pairs,
  emi_pairs_sens6,
  emi_pairs_sens8,
  data_clean,
  compliance,
  file = save_path
)

cat("\n✓ Gespeichert unter:", save_path, "\n")
cat("  emi_pairs        – Haupt-Cut-off >= 4 Wörter\n")
cat("  emi_pairs_sens6  – Cut-off >= 6 Wörter\n")
cat("  emi_pairs_sens8  – Cut-off >= 8 Wörter\n")
cat("  data_clean       – vollständiger Datensatz\n")
cat("  compliance       – Compliance pro Proband\n")

sessionInfo()

###############################################################################
# FERTIG
###############################################################################
