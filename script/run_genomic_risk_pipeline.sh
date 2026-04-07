#!/bin/bash

echo "======================================"
echo " GENOMIC RISK ANALYSIS PIPELINE START "
echo "======================================"

# -------------------------------
# PATHS (FIXED + CONSISTENT)
# -------------------------------
BASE="/Users/neha/iResearch_Genomic_Risk/Pipeline_bash"

DATA_DIR="$BASE/Data_Files"
OUTPUT_DIR="$BASE/Results"
# GENOME_DIR="$DATA_DIR"
GENOME_DIR="$DATA_DIR/genomes" 

SCRIPT_DIR="$BASE/script"

TAXID_JSON="$DATA_DIR/taxid_to_name.json"
METADATA_CSV="$DATA_DIR/E.coli_combined_diverse_commensals.csv"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$GENOME_DIR"

## -------------------------------
## STEP 0: DOWNLOAD REFERENCE GENOMES
## -------------------------------
echo "Step 0: Downloading Reference Genomes..."

python3 "$SCRIPT_DIR/download_complete_genomes.py" \
    --input "$DATA_DIR/selected_accessions.txt" \
    --out "$DATA_DIR" \
    --log "$DATA_DIR/downloaded_accessions.txt"

## -------------------------------
## STEP 1: METADATA PROCESSING
## -------------------------------
echo "Step 1: Metadata Processing..."

python3 "$SCRIPT_DIR/genomic_metadata_processing.py" \
    --output "$OUTPUT_DIR" \
    --taxid "$TAXID_JSON" \
    --metadata "$METADATA_CSV" \
    --genomes "$GENOME_DIR"

## -------------------------------
## STEP 2: ABRICATE ANALYSIS (AMR + VIRULENCE + SUMMARY)
## -------------------------------
echo "Step 2: Running Abricate Analysis..."

ABRICATE_OUT="$OUTPUT_DIR/abricate_analysis"
mkdir -p "$ABRICATE_OUT"

## ---- AMR (ResFinder) ----
echo "Running AMR (ResFinder)..."
abricate --db resfinder "$GENOME_DIR"/*.fna \
    > "$ABRICATE_OUT/amr_results.tsv"

## ---- Virulence (VFDB) ----
echo "Running Virulence (VFDB)..."
abricate --db vfdb "$GENOME_DIR"/*.fna \
    > "$ABRICATE_OUT/virulence_results.tsv"

## ---- Summary ----
echo "Generating summaries..."

abricate --summary "$ABRICATE_OUT/amr_results.tsv" \
    > "$ABRICATE_OUT/amr_results.csv"

abricate --summary "$ABRICATE_OUT/virulence_results.tsv" \
    > "$ABRICATE_OUT/virulence_results.csv"

echo "Abricate Analysis Completed ✅"

## -------------------------------
## STEP 3: THREAT CLASSIFICATION
##-------------------------------
echo "Step 3: Running Threat Classification..."

THREAT_CLASSIFICATION_OUT="$OUTPUT_DIR/threat_classification_analysis"
mkdir -p "$THREAT_CLASSIFICATION_OUT"

python3 "$SCRIPT_DIR/threat_classification.py" \
    --vir "$ABRICATE_OUT/virulence_results.csv" \
    --amr "$ABRICATE_OUT/amr_results.csv" \
    --out "$THREAT_CLASSIFICATION_OUT/final_threat_classification.csv" \
    --logdir "$THREAT_CLASSIFICATION_OUT"

## -------------------------------
## STEP 4: CLASSIFICATION STATISTICS
## -------------------------------
echo "Step 4: Running Classification Statistics..."

CLASSIFICATION_STATISTICS_OUT="$OUTPUT_DIR/classification_statistics"
mkdir -p "$CLASSIFICATION_STATISTICS_OUT"

python3 "$SCRIPT_DIR/classification_statistics.py" \
    --classification "$THREAT_CLASSIFICATION_OUT/final_threat_classification.csv" \
    --metadata "$METADATA_CSV" \
    --outdir "$CLASSIFICATION_STATISTICS_OUT"

