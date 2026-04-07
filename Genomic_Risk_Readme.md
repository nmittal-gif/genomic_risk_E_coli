### Genomic Risk Analysis Pipeline

This pipeline automates genome download, metadata processing, AMR detection, virulence analysis, and threat classification. It produces a final genomic threat registry, along with statistical analysis and visualizations.
Execution Methods

### You can run the pipeline in two ways:

1. **Local Execution (Bash Script)**

 - Run the pipeline directly on your system.
    - `bash run_genomic_risk_pipeline.sh`
 - Executes the full pipeline automatically
 - Requires all dependencies installed locally 

2. **Containerized Execution (Docker)**

 - Ensures reproducibility and avoids dependency issues.
    - Build Docker Image:
        - `docker build -t genomic-risk .`
    - Run Container:
        - `docker run --rm -v $(pwd):/app genomic-risk`

## Benefits of Docker:
 - No dependency conflicts
 - Consistent results across systems
 - Isolated and clean environment

### Pipeline Workflow

1. **Genome Download**
 - Downloads genomes using accession lists
 - Retrieves .fna and .gff files
 - Supports batch processing and retry mechanisms
 - Downloaded genome files are stored in:
    - `/app/Data_Files/`

2. **Metadata Processing**
 - Maps TaxID to organism names
 - Processes genome metadata
 - Generates summary tables and distribution plots

3. **AMR Detection (ResFinder)**
 - Uses Abricate with ResFinder database
 - Identifies antimicrobial resistance (AMR) genes
 - Produces structured AMR reports

4. **Virulence Detection (VFDB)**
 - Uses Abricate with VFDB database
 - Detects virulence-associated genes
 - Generates virulence profiling reports

5. **Threat Classification**
 - Assigns threat categories (0–3)
 - Integrates AMR and virulence profiles
 - Classifies genomes into risk groups

6. **Statistical Analysis**
 - Computes AMR burden and virulence burden
 - Generates category-wise statistics
 - Identifies high-risk genomes
 - Final Genomic Threat Registry
 - Consolidates all results into a single registry
 - Includes classification, gene counts, and metadata

7. **Visualization Generation**
 - Boxplots (AMR vs virulence burden)
 - Scatter plots for high-risk genomes
 - Stacked bar charts for category distribution

## Input Files

**Place all required files inside**:
 `/app/Data_Files/`

**Required Inputs**:
`selected_accessions.txt` → Genome accession list
`taxid_to_name.json` → TaxID to organism mapping
`E.coli_combined_diverse_commensals.csv`→ Metadata file

### Output Files

## All outputs are saved in:
`/app/Results/`
# Generated Outputs:
- Metadata Analysis
    - `taxid_genome_summary.csv` -A table mapping genome accessions to their specific organism names and TaxIDs.
    - `taxid_distribution.png` - A chart showing the frequency of different species found in your dataset.
    - `genome_metrics.csv` - Basic assembly statistics including genome length and GC content for every file.
    - `gc_content_distribution.png` - A histogram visualizing the variation in GC percentages across your samples.
    - `genome_length_distribution.png` - A plot showing the distribution of total base pairs per genome.

- AMR & Virulence Results
    - `amr_results.tsv, amr_results.csv` - Raw and summarized reports of antimicrobial resistance genes detected by Abricate (ResFinder)
    - `virulence_results.tsv, virulence_results.csv` - Raw and summarized reports of virulence-associated factors detected by Abricate (VFDB).

- Threat Classification
    - `final_threat_classification.csv` - The primary results file where each genome is assigned a Category 0–3 based on its gene profile.
    - `genome_threat_distribution_piechart.png` - A proportional breakdown of the dataset showing the percentage of genomes in each risk level.
    - `genome_threat_distribution_barchart.png` - A bar chart providing a direct count comparison of the four threat categories.

- Final Registry & Statistics
    - `Final_Genomic_Threat_Registry.csv` - Combining all metadata, gene counts, and classification levels into one file.
    - `Threat_Category_Basic_Statistics.csv` - Statistical summary comparing gene "burden" across categories.
    - `High_Risk_Category2_3_Genome.csv` - A filtered list containing only genomes with high-threat markers (Category 2 and 3).
    - `Top_High_Risk_Genomes_Final` - The most critical threats identified for immediate priority review.

- Visualizations
    - `Threat_Category_Burden_Boxplots.png` - Statistical plots showing the density of AMR and Virulence genes within each threat category.
    - `Category2_3_Threat_Scatterplot.png` - A visualization highlighting the specific relationship between gene burden and high-risk classifications.
    - `Accession_Country_StackedBar.png`- A chart showing how different threat levels are distributed across different countries of origin.

## Requirements
**System Requirements**
 - `Linux / Ubuntu`
  - `Python 3.8+`

**Required Tools**
 - `abricate`
 - `datasets (NCBI CLI)`

**Main Libraries:**
 - `pandas`
 - `numpy`
 - `matplotlib`
 - `seaborn`
 - `scipy`

## 👩‍💻Author

Neha Mittal