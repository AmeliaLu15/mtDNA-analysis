# =============================================================================
# mtDNA Deletion Quality Control Pipeline
# =============================================================================
# Description: Quality control pipeline for mtDNA deletion data from eKLIPse.
#              Filters by BLAST support, patient metadata, missing fragments,
#              contamination, heteroplasmy level, and groups by frequency/size.
#
# Input files:
#   - eKLIPse output CSV (deletion calls)
#   - Patients metadata CSV
#
# Output files:
#   - df_all.csv              : Raw data after formatting
#   - BLAST_QC.csv            : After BLAST filter
#   - Patients.csv            : Cleaned patient metadata
#   - case_QC.csv             : After patient matching
#   - fragments_QC.csv        : After missing fragment filter
#   - contamination_QC.csv    : After contamination filter
#   - Heteroplasmy_QC.csv     : After heteroplasmy level filter
#   - final_QC.csv            : Final dataset with grouping variables
#
# =============================================================================


# --- 1. LOAD REQUIRED PACKAGES -----------------------------------------------

# Install any missing packages automatically
required_packages <- c("dplyr", "readr")
installed <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed) install.packages(pkg)
}

library(dplyr)
library(readr)


# --- 2. SET FILE PATHS -------------------------------------------------------
# No changes needed if you cloned this repository from GitHub.
# Files are read from the data/ folder and results saved to output/ folder.
# Folder structure expected:
#
#   mtDNA-analysis/
#   ├── data/
#   │   ├── eKLIPse_deletions.csv
#   │   └── Patients.csv
#   └── scripts/
#       └── mtDNA_deletion_QC.R   ← this script

eklipse_file  <- "data/eKLIPse_deletions.csv"   # eKLIPse output
patients_file <- "data/Patients.csv"             # Patient metadata


# --- 3. QC THRESHOLDS ---------------------------------------

blast_min        <- 3      # Minimum BLAST support for both 5' and 3' breakpoints
contamination_max <- 0.1   # Maximum contamination level allowed
heteroplasmy_min  <- 0.1   # Minimum heteroplasmy frequency (%)


# --- 4. LOAD DATA ------------------------------------------------------------

message("Loading data...")

df <- read_csv(eklipse_file, show_col_types = FALSE)

Patients <- read_csv(patients_file, show_col_types = FALSE)

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# --- 5. FORMAT RAW DATA ------------------------------------------------------

message("Formatting raw data...")

# Fix frequency column: remove thousand separators and convert to numeric
df$Freq <- as.numeric(gsub(",", ".", gsub("\\.", "", df$Freq)))

# Rename breakpoint columns for clarity
df <- df %>%
  rename(
    Start = `X5..breakpoint`,
    End   = `X3..breakpoint`
  )

write_csv(df, file.path(output_dir, "df_all.csv"))
message("  Saved: df_all.csv  (n = ", nrow(df), ")")


# --- 6. QC STEP 1: BLAST FILTER ----------------------------------------------

message("QC Step 1: BLAST filter (>= ", blast_min, " each end)...")

df_BLAST <- df %>%
  filter(X3..Blast >= blast_min & X5..Blast >= blast_min)

write_csv(df_BLAST, file.path(output_dir, "BLAST_QC.csv"))
message("  Saved: BLAST_QC.csv  (n = ", nrow(df_BLAST), ")")


# --- 7. CLEAN PATIENT METADATA -----------------------------------------------

message("Cleaning patient metadata...")

Patients <- Patients %>%
  rename(ID = Patients)

write_csv(Patients, file.path(output_dir, "Patients.csv"))
message("  Saved: Patients.csv  (n = ", nrow(Patients), ")")


# --- 8. QC STEP 2: MATCH TO PATIENT METADATA ---------------------------------

message("QC Step 2: Matching to patient metadata...")

df_BLAST <- left_join(df_BLAST, Patients, by = "ID") %>%
  filter(!is.na(ID))

write_csv(df_BLAST, file.path(output_dir, "case_QC.csv"))
message("  Saved: case_QC.csv  (n = ", nrow(df_BLAST), ")")


# --- 9. QC STEP 3: MISSING FRAGMENTS FILTER ----------------------------------

message("QC Step 3: Missing fragments filter (Fragments == 'Y')...")

df_MF <- df_BLAST %>%
  filter(Fragments == "Y")

write_csv(df_MF, file.path(output_dir, "fragments_QC.csv"))
message("  Saved: fragments_QC.csv  (n = ", nrow(df_MF), ")")


# --- 10. QC STEP 4: CONTAMINATION FILTER -------------------------------------

message("QC Step 4: Contamination filter (< ", contamination_max, ")...")

df_contamination <- df_MF %>%
  filter(Contamination < contamination_max)

write_csv(df_contamination, file.path(output_dir, "contamination_QC.csv"))
message("  Saved: contamination_QC.csv  (n = ", nrow(df_contamination), ")")


# --- 11. QC STEP 5: HETEROPLASMY LEVEL FILTER --------------------------------

message("QC Step 5: Heteroplasmy level filter (>= ", heteroplasmy_min, "%)...")

df_HL <- df_contamination %>%
  filter(Freq >= heteroplasmy_min)

write_csv(df_HL, file.path(output_dir, "Heteroplasmy_QC.csv"))
message("  Saved: Heteroplasmy_QC.csv  (n = ", nrow(df_HL), ")")


# --- 12. ADD GROUPING VARIABLES ----------------------------------------------

message("Adding grouping variables (Condition, Freq sector, Size density)...")

df_final <- df_HL %>%
  # Set Condition as ordered factor
  mutate(Condition = factor(Condition, levels = c("Control", "PD"))) %>%
  # Re-join patient info to include samples with 0 observations
  left_join(Patients, by = "ID") %>%
  # Heteroplasmy frequency groups
  mutate(Freq.sector = cut(
    Freq,
    breaks = c(0.1, 1, 10, 100),
    labels = c("0.1-1%", "1-10%", ">10%"),
    include.lowest = TRUE
  )) %>%
  # Deletion size groups
  mutate(Size.density = cut(
    Size,
    breaks = c(5, 140, 1400, 2400, 3400, 4400, 5400, 6400),
    labels = c("5-140bp", "140-1400bp", "1400-2400bp",
               "2400-3400bp", "3400-4400bp", "4400-5400bp", "5400-6400bp"),
    include.lowest = TRUE
  ))

write_csv(df_final, file.path(output_dir, "final_QC.csv"))
message("  Saved: final_QC.csv  (n = ", nrow(df_final), ")")
