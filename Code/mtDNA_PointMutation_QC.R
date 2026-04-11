# =============================================================================
# mtDNA Point Mutation Quality Control Pipeline
# =============================================================================
# Description: Quality control pipeline for mtDNA point mutation data from
#              mtDNA-Server. Merges per-batch variant and contamination files,
#              filters by heteroplasmy level, strand bias, variant type,
#              missing fragments, contamination, and patient metadata.
#              Identifies somatic (non-inherited) mutations per patient.
#
# Input files:
#   data/variants/ folder (7 batches, 68 samples total):
#   - data/variants/1-10.txt
#   - data/variants/11-20.txt
#   - data/variants/21-30.txt
#   - data/variants/31-40.txt
#   - data/variants/41-50.txt
#   - data/variants/51-60.txt
#   - data/variants/61-68.txt
#
#   data/contamination/ folder (7 batches):
#   - data/contamination/1-10.txt
#   - data/contamination/11-20.txt
#   - data/contamination/21-30.txt
#   - data/contamination/31-40.txt
#   - data/contamination/41-50.txt
#   - data/contamination/51-60.txt
#   - data/contamination/61-68.txt
#
#   Patient metadata:
#   - data/Patients.csv
#
# Expected column names in variant files (from mtDNA-Server):
#   ID, Filter, Pos, Ref, Variant, VariantLevel, Coverage,
#   MajorBase, MajorLevel, MinorBase, MinorLevel,
#   CoverageFWD, CoverageREV, Type
#   (ID contains values like "10.merged.sorted.bam")
#
# Output files (saved to output/):
#   - variant_combined.csv       : All variants merged across batches
#   - contamination_combined.csv : All contamination data merged
#   - variant_QC.csv             : After heteroplasmy, strand bias, type filter
#   - case_QC.csv                : After patient matching
#   - fragments_QC.csv           : After missing fragment filter
#   - contamination_QC.csv       : After contamination filter
#   - somatic_QC.csv             : After somatic mutation filter
#   - final_QC.csv               : Final dataset with frequency groups
#
# =============================================================================


# --- 1. LOAD REQUIRED PACKAGES -----------------------------------------------

required_packages <- c("dplyr", "readr")
installed <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed) install.packages(pkg)
}

library(dplyr)
library(readr)


# --- 2. SET FILE PATHS -------------------------------------------------------
# No changes needed if you cloned this repository from GitHub.
# Folder structure expected:
#
#   mtDNA-analysis/
#   ├── data/
#   │   ├── variants/              ← folder of variant batch files
#   │   │   ├── 1-10.txt
#   │   │   ├── 11-20.txt
#   │   │   ├── 21-30.txt
#   │   │   ├── 31-40.txt
#   │   │   ├── 41-50.txt
#   │   │   ├── 51-60.txt
#   │   │   └── 61-68.txt
#   │   ├── contamination/         ← folder of contamination batch files
#   │   │   ├── 1-10.txt
#   │   │   ├── 11-20.txt
#   │   │   ├── 21-30.txt
#   │   │   ├── 31-40.txt
#   │   │   ├── 41-50.txt
#   │   │   ├── 51-60.txt
#   │   │   └── 61-68.txt
#   │   └── Patients.csv
#   └── Code/
#       └── mtDNA_pointmutation_QC.R

variants_dir      <- "data/variants/"
contamination_dir <- "data/contamination/"
patients_file     <- "data/Patients.csv"


# --- 3. QC THRESHOLDS ---------------------------------------

heteroplasmy_min   <- 0.001   # Minimum variant level (0.1%)
contamination_max  <- 0.005   # Maximum contamination level allowed
variant_type       <- 2       # Variant type to keep (2 = heteroplasmic)
strand_bias_filter <- "PASS"  # Only keep variants that pass strand bias filter


# --- 4. LOAD AND MERGE VARIANT BATCH FILES -----------------------------------

message("Loading variant batch files from: ", variants_dir)

variant_files <- list.files(
  path       = variants_dir,
  pattern    = "\\.txt$",
  full.names = TRUE
)

if (length(variant_files) == 0) stop("No .txt files found in: ", variants_dir)
message("  Found ", length(variant_files), " file(s): ",
        paste(basename(variant_files), collapse = ", "))

# Read and merge all batch files.
# CI_MitoTool and DuplSeq_rCRS_pos can differ in type across batches,
# so standardise to character before merging, then convert back.
variant <- do.call("rbind", lapply(variant_files, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
  if ("CI_MitoTool" %in% names(df))      df$CI_MitoTool      <- as.character(df$CI_MitoTool)
  if ("DuplSeq_rCRS_pos" %in% names(df)) df$DuplSeq_rCRS_pos <- as.character(df$DuplSeq_rCRS_pos)
  df
}))

if ("CI_MitoTool" %in% names(variant)) variant$CI_MitoTool <- as.numeric(variant$CI_MitoTool)

# Extract numeric sample number from ID (e.g. "10.merged.sorted.bam" -> "10")
variant$Sample <- gsub("\\D", "", variant$ID)

write_csv(variant, file.path(output_dir, "variant_combined.csv"))
message("  Saved: variant_combined.csv  (n = ", nrow(variant), " variants)")


# --- 5. LOAD AND MERGE CONTAMINATION BATCH FILES -----------------------------

message("Loading contamination batch files from: ", contamination_dir)

contamination_files <- list.files(
  path       = contamination_dir,
  pattern    = "\\.txt$",
  full.names = TRUE
)

if (length(contamination_files) == 0) stop("No .txt files found in: ", contamination_dir)
message("  Found ", length(contamination_files), " file(s): ",
        paste(basename(contamination_files), collapse = ", "))

contamination <- do.call("rbind", lapply(contamination_files, function(f) {
  read.table(f, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
}))

# Extract numeric sample number
contamination$Sample <- gsub("\\D", "", contamination$Sample)

write_csv(contamination, file.path(output_dir, "contamination_combined.csv"))
message("  Saved: contamination_combined.csv  (n = ", nrow(contamination), " samples)")
message("  Mean coverage: ", round(mean(contamination$Sample.Coverage), 1),
        "  SD: ", round(sd(contamination$Sample.Coverage), 1))


# --- 6. MERGE VARIANTS WITH CONTAMINATION ------------------------------------

message("Merging variant and contamination data...")

df <- left_join(variant, contamination, by = "Sample")


# --- 7. QC STEP 1: HETEROPLASMY LEVEL, STRAND BIAS, VARIANT TYPE -------------

message("QC Step 1: Heteroplasmy (> ", heteroplasmy_min,
        "), strand bias (", strand_bias_filter,
        "), variant type (", variant_type, ")...")

df_QC <- df %>%
  filter(VariantLevel > heteroplasmy_min) %>%
  filter(Filter == strand_bias_filter) %>%
  filter(Type == variant_type)

write_csv(df_QC, file.path(output_dir, "variant_QC.csv"))
message("  Saved: variant_QC.csv  (n = ", nrow(df_QC), ")")


# --- 8. LOAD PATIENT METADATA ------------------------------------------------

message("Loading patient metadata...")

Patients <- read_csv(patients_file, show_col_types = FALSE) %>%
  rename(Sample = Title)

Patients$Sample <- as.character(Patients$Sample)


# --- 9. QC STEP 2: MATCH TO PATIENT METADATA ---------------------------------

message("QC Step 2: Matching to patient metadata...")

df_QC$Sample <- as.character(df_QC$Sample)

df_patient <- left_join(df_QC, Patients, by = "Sample") %>%
  filter(!is.na(ID))

write_csv(df_patient, file.path(output_dir, "case_QC.csv"))
message("  Saved: case_QC.csv  (n = ", nrow(df_patient), ")")


# --- 10. QC STEP 3: MISSING FRAGMENTS FILTER ---------------------------------

message("QC Step 3: Missing fragments filter (Fragments == 'Y')...")

df_MF <- df_patient %>%
  filter(Fragments == "Y")

write_csv(df_MF, file.path(output_dir, "fragments_QC.csv"))
message("  Saved: fragments_QC.csv  (n = ", nrow(df_MF), ")")


# --- 11. QC STEP 4: CONTAMINATION FILTER -------------------------------------

message("QC Step 4: Contamination filter (< ", contamination_max, ")...")

df_contamination <- df_MF %>%
  filter(Contamination < contamination_max)

write_csv(df_contamination, file.path(output_dir, "contamination_QC.csv"))
message("  Saved: contamination_QC.csv  (n = ", nrow(df_contamination), ")")


# --- 12. QC STEP 5: SOMATIC MUTATION FILTER ----------------------------------
# A mutation is considered somatic if it appears in only one cell
# per patient. Mutations appearing in multiple samples
# of the same patient are likely inherited germline variants and are excluded.

message("QC Step 5: Somatic mutation filter...")
message("  Logic: keep mutations appearing in exactly 1 sample per patient")

df_somatic <- df_contamination %>%
  group_by(ID, Pos) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  mutate(Inherited = FALSE)

write_csv(df_somatic, file.path(output_dir, "somatic_QC.csv"))
message("  Saved: somatic_QC.csv  (n = ", nrow(df_somatic), ")")


# --- 13. ADD FREQUENCY GROUPS ------------------------------------------------

message("Adding frequency groups...")

df_final <- df_somatic %>%
  mutate(Freq = cut(
    VariantLevel,
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = c("0.1-1%", "1-10%", ">10%"),
    include.lowest = TRUE
  ))

write_csv(df_final, file.path(output_dir, "final_QC.csv"))
message("  Saved: final_QC.csv  (n = ", nrow(df_final), ")")
