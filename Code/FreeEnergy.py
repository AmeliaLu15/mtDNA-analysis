"""
mtDNA Deletion Breakpoint Free Energy Calculation
==================================================
Description:
    Calculates the minimum free energy (MFE) of secondary structures formed
    by sequences flanking mtDNA deletion breakpoints, using the seqfold
    package. For each deletion, sequences are extracted symmetrically around
    both breakpoints and concatenated before folding.

    Ambiguous bases (N) in sample sequences are replaced using the rCRS
    reference sequence before folding.

Input files:
    - data/Free_energy_masterfile.csv  : Deletion metadata with columns:
                                         Title (cell/sample ID), Start, End
    - data/rCRS.fasta                  : Reference sequence (rCRS, NC_012920.1)
    - data/<CellID>.fasta              : Per-sample mtDNA FASTA files
                                         (one file per cell/sample)

Output files:
    - FreeEnergyResult.csv      : Free energy results with columns:
                                         CellID, Start, End, FreeEnergy, Windows

Parameters:
    - WINDOW_SIZES (list): Window sizes to test.
                           Default: 5-49 bp (step 1) + 50-100 bp (step 5)
    - TEMP (float)       : Folding temperature in Celsius (default: 37.0)

Usage:
    python FreeEnergy.py

Requirements:
    - Python 3.12.1
    - pandas
    - seqfold

Install requirements:
    pip install pandas seqfold


import os
import pandas
from seqfold import dg


# --- 1. SET FILE PATHS -------------------------------------------------------
# No changes needed if you cloned this repository from GitHub.
# Folder structure expected:
#
#   mtDNA-analysis/
#   ├── data/
#   │   ├── Free_energy_masterfile.csv   ← deletion metadata
#   │   ├── rCRS.fasta                   ← reference sequence
#   │   └── <CellID>.fasta               ← one FASTA file per sample
#   └── scripts/
#       └── free_energy.py

data_dir   = "data/"

masterfile = os.path.join(data_dir, "Free_energy_masterfile.csv")
rcrs_file  = os.path.join(data_dir, "rCRS.fasta")


# --- 2. PARAMETERS -----------------------------------------------------------

# Window sizes to test (bp)
# Fine resolution (5-49, step 1) + coarser resolution (50-100, step 5)
WINDOW_SIZES = list(range(5, 50, 1)) + list(range(50, 101, 5))

TEMP = 37.0   # Folding temperature (degrees Celsius)


# --- 3. LOAD REFERENCE SEQUENCE (rCRS) ---------------------------------------

print("Loading rCRS reference sequence...")

with open(rcrs_file) as f:
    lines = f.read().strip().split("\n")

# Skip the FASTA header line (starts with ">")
refseq = "".join(lines[1:])
print(f"  rCRS length: {len(refseq)} bp")


# --- 4. LOAD DELETION METADATA -----------------------------------------------

print(f"Loading deletion metadata from: {masterfile}")

meta = pandas.read_csv(masterfile)
print(f"  {len(meta)} deletions found")


# --- 5. DEFINE FREE ENERGY CALCULATION FUNCTION ------------------------------

def aligned_free_energy(fasta_path, break1, break2, windowsize):
    """
    Calculate minimum free energy (MFE) of the sequence formed by
    concatenating windows around both deletion breakpoints.

    Sequences are extracted symmetrically around each breakpoint:
        seq1 = mtDNA[ (break1 - windowsize/2) : (break1 + windowsize/2) ]
        seq2 = mtDNA[ (break2 - windowsize/2) : (break2 + windowsize/2) ]

    The two sequences are concatenated and their MFE calculated at 37°C
    using seqfold.dg().

    Ambiguous bases (N) are replaced with the corresponding rCRS base
    before folding. If rCRS itself has an N at that position, C is used
    as a standard convention.

    Parameters
    ----------
    fasta_path : str
        Path to the per-sample mtDNA FASTA file.
    break1 : int
        5' breakpoint position (0-based).
    break2 : int
        3' breakpoint position (0-based).
    windowsize : int
        Window size in bp. Half of this is taken on each side of each breakpoint.

    Returns
    -------
    float
        Minimum free energy (kcal/mol) of the concatenated sequence.
    """
    with open(fasta_path) as f:
        lines = f.read().strip().split("\n")

    # Skip FASTA header line
    mtDNA = "".join(lines[1:])

    # Replace ambiguous bases (N) with rCRS reference base
    while "N" in mtDNA:
        index = mtDNA.index("N")
        replace = "C" if refseq[index] == "N" else refseq[index]
        mtDNA = mtDNA[:index] + replace + mtDNA[index + 1:]

    half = int(windowsize / 2)

    # Extract sequence window around each breakpoint
    seq1 = mtDNA[(break1 - half) : (break1 + half)]
    seq2 = mtDNA[(break2 - half) : (break2 + half)]

    # Calculate minimum free energy of concatenated sequence at 37°C
    return dg(seq1 + seq2, temp=TEMP)


# --- 6. RUN FREE ENERGY CALCULATION ------------------------------------------

print(f"\nCalculating free energy across {len(WINDOW_SIZES)} window sizes...")
print(f"  Window range: {min(WINDOW_SIZES)}-{max(WINDOW_SIZES)} bp")
print(f"  Temperature:  {TEMP}°C")
print(f"  Total iterations: {len(WINDOW_SIZES)} windows × {len(meta)} deletions "
      f"= {len(WINDOW_SIZES) * len(meta)}")

os.makedirs(output_dir, exist_ok=True)

output = []
errors = []

for windowsize in WINDOW_SIZES:
    for _, row in meta.iterrows():
        cell_id = str(row["Title"])
        start   = int(row["Start"])
        end     = int(row["End"])

        fasta_path = os.path.join(data_dir, f"{cell_id}.fasta")

        # Skip if FASTA file is missing
        if not os.path.exists(fasta_path):
            if cell_id not in errors:
                print(f"  WARNING: FASTA file not found for {cell_id} — skipping")
                errors.append(cell_id)
            continue

        free_energy = aligned_free_energy(fasta_path, start, end, windowsize)
        output.append([cell_id, start, end, free_energy, windowsize])

    print(f"  Window {windowsize} bp done")


# --- 7. SAVE OUTPUT ----------------------------------------------------------

output_file = os.path.join(output_dir, "FreeEnergyResult.csv")
df = pandas.DataFrame(output, columns=["CellID", "Start", "End", "FreeEnergy", "Windows"])
df.to_csv(output_file, index=False)

