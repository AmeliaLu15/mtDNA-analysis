# mtDNA-analysis
Scripts used in: Mitochondrial DNA damage elicits pro-survival PINK1 transcriptional upregulation in degenerating pedunculopontine-cholinergic neurons during Parkinsons's disease

## Analysis Overview
- mtDNA deletion detection: eKLIPse
- mtDNA point mutation detection: mtDNA-Server  
- Quality control: R scripts (this repository)
- Secondary structure prediction: UNAFold, SeqFold
- SeqFold get sequence scripts : Python scripts

## Files
- `Code/mtDNA_deletion_QC.R` — Quality control for mtDNA deletion data
- `Code/mtDNA_PointMutation_QC.R` — Quality control for point mutation data
- `Code/SeqFold_GetSequence.py`

## Software Requirements
- R version 4.2.1
- R package: dplyr,readr
- Python version 3.12.1
- Python package: pandas
