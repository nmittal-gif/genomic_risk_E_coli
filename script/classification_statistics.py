#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import re
import logging
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.stats import kruskal
from matplotlib.colors import ListedColormap
import argparse


# In[2]:


# ==========================================
# 1. SURGICAL FUNCTION IMPORT
# ==========================================
def setup_logging(logs_dir, log_file_name):
    os.makedirs(logs_dir, exist_ok=True)

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        filename=os.path.join(logs_dir, log_file_name),
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
# ==========================================
# 3. FINALIZE GENOMIC THREAT REGISTRY
# ==========================================
def finalize_threat_registry(input_csv, metadata_csv, output_dir_path):
    """
    Finalizes the Genomic Threat Registry by integrating classification
    results with genome metadata and computing virulence and AMR summaries.

    Pipeline Steps:
    1. Creates output directory and initializes logging.
    2. Loads classification results.
    3. Extracts virulence gene information and calculates virulence burden.
    4. Extracts AMR gene information and calculates AMR burden.
    5. Identifies last-resort resistance genes and MDR status.
    6. Maps metadata using genome accession numbers.
    7. Exports the finalized genomic threat registry.

    Outputs Generated:
        - classification_statistics/Final_Genomic_Threat_Registry.csv
        - classification_statistics/threat_registry.log
    """
    logging.info("===== THREAT REGISTRY FINALIZATION STARTED =====")
    # ------------------------------------------
    # Ensure the output directory exists
    # ------------------------------------------
    if not os.path.exists(output_dir_path):
        print(f"Directory {output_dir_path} not found. Creating it now...")
        logging.info(f"Output directory not found. Creating: {output_dir_path}")
        os.makedirs(output_dir_path, exist_ok=True)
    # ------------------------------------------
    # STEP 1: Load classification dataset
    # ------------------------------------------
    print("Step 1: Loading Week 2 Classification Data...")
    logging.info("Loading classification data")
    df = pd.read_csv(input_csv)
    logging.info(f"Loaded {len(df)} genomes from classification file")
    cols = df.columns.tolist()
    
    # Identify where virulence and AMR gene columns start and end
    idx_x = cols.index('NUM_FOUND_x') # Virulence marker
    idx_y = cols.index('NUM_FOUND_y') # AMR marker

    vir_gene_cols = cols[idx_x + 1 : idx_y]# Extract virulence gene column names
    amr_gene_cols = cols[idx_y + 1 : ] # Extract AMR gene column names
    logging.info(f"Detected {len(vir_gene_cols)} virulence genes")
    logging.info(f"Detected {len(amr_gene_cols)} AMR genes")
    
    # ------------------------------------------
    # Define gene groups used in analysis
    # ------------------------------------------ 
    hr_stx_targets = ['stx2a', 'stx2d']  # High-risk shiga toxin genes
    last_resort_targets = ['blaCTX-M', 'mcr', 'blaNDM', 'blaKPC', 'blaOXA-48'] # Critical last-resort AMR genes
    # Helper function to determine whether a gene is present
    def is_present(val):
        return str(val).strip() != '.'

    # ------------------------------------------
    # STEP 2: Virulence analysis
    # ------------------------------------------
    print("Step 2: Extracting Virulence Weapon Names & High-Risk Counts...")
    df['Virulence_Burden'] = pd.to_numeric(df['NUM_FOUND_x'], errors='coerce').fillna(0).astype(int) # Convert virulence count column to integer
    # List all virulence gene names
    df['All_Virulence_Genes'] = df.apply(lambda r: ", ".join([g for g in vir_gene_cols if is_present(r[g])]), axis=1)  # Collect all virulence gene names detected in each genome
    
    # Get High-Risk Stx Names and Count
    def get_hr_stx_data(row):
        # We convert 'g' to lowercase so it matches 'stx2a' even if the column is 'Stx2A'
        found = [g for g in vir_gene_cols if is_present(row[g]) and any(t in g.lower() for t in hr_stx_targets)]
        # found = [g for g in vir_gene_cols if is_present(row[g]) and any(t in g for t in hr_stx_targets)] # Find virulence genes matching high-risk targets
        names = ", ".join(found) if found else "None" # Convert list to readable string
        count = len(found)
        return pd.Series([names, count])
        
    df[['High_Risk_Stx_Genes', 'High_Risk_Stx_Count']] = df.apply(get_hr_stx_data, axis=1)  # Create two new columns: gene names and counts

    # ------------------------------------------
    # STEP 3: AMR analysis
    # ------------------------------------------
    print("Step 3: Calculating AMR Counts...")
    df['AMR_Burden'] = pd.to_numeric(df['NUM_FOUND_y'], errors='coerce').fillna(0).astype(int)
    # Identify critical last-resort resistance genes
    def get_last_resort(row):
        # Use g.lower() to match 'blactx-m' regardless of original column casing
        found = [g for g in amr_gene_cols if is_present(row[g]) and any(t.lower() in g.lower() for t in last_resort_targets)]
        # found = [g for g in amr_gene_cols if is_present(row[g]) and any(t in g for t in last_resort_targets)]
        return ", ".join(found) if found else "None"
    # Define a helper function to get standard genes (excluding last-resort)
    def get_standard_amr(row):
        # Filter for genes that are present AND NOT in the last_resort_targets list
        found = [g for g in amr_gene_cols if is_present(row[g]) and not any(t.lower() in g.lower() for t in last_resort_targets)]
        return ", ".join(found) if found else "None"
        
    df['Last_Resort_AMR_Genes'] = df.apply(get_last_resort, axis=1) # Store detected last-resort gene name
    df['Last_Resort_AMR_Count'] = df['Last_Resort_AMR_Genes'].apply(lambda x: 0 if x == "None" else len(x.split(',')))  # Count last-resort genes
    df['Standard_AMR_Genes'] = df.apply(get_standard_amr, axis=1)
    df['Standard_AMR_Count'] = df['AMR_Burden'] - df['Last_Resort_AMR_Count']  # Calculate standard AMR genes
    df['Is_MDR'] = df['AMR_Burden'] >= 3 # Identify multidrug resistant genomes

     # ------------------------------------------
    # STEP 3: AMR analysis
    # ------------------------------------------
    print("Step 4: Mapping Metadata...")
    
    # Extract genome accession number from filename
    def extract_accession(filepath):
        filename = os.path.basename(filepath)
        # Finds GCA_000... or GCF_000... including the .1, .2 version
        match = re.search(r'(GC[AF]_\d+\.\d+)', filename)
        if match:
            return match.group(1)
        # return filename.split('.')[0]
        return filename.replace('.fna', '').replace('.fasta', '').replace('.fa', '')

    df['accession'] = df['#FILE'].apply(extract_accession)
    if os.path.exists(metadata_csv):
        meta_df = pd.read_csv(metadata_csv)
        
        # Ensure accession column is string and stripped of spaces
        df['accession'] = df['accession'].astype(str).str.strip()
        meta_df['accession'] = meta_df['accession'].astype(str).str.strip()
        # 🔍 DEBUG: find mismatches BEFORE merge
        missing = df[~df['accession'].isin(meta_df['accession'])]
        
        print("\nDEBUG: Mismatched accessions with hidden characters:")
        for acc in missing['accession']:
            print(repr(acc))
        keep_meta = ['accession', 'tax_id', 'strain', 'host_clean', 'source_type', 'country_clean', 'year_clean', 'diversity_group']
        existing_meta = [c for c in keep_meta if c in meta_df.columns]
        
        # Merge on accession
        df = pd.merge(df, meta_df[existing_meta], on='accession', how='left')
        
        # Check mapping success
        # mapped_count = df['strain'].notna().sum()
        mapped_count = df['accession'].isin(meta_df['accession']).sum()
        print(f"Metadata Mapping Complete: {mapped_count} out of {len(df)} genomes matched.")
        logging.info(f"Metadata mapping completed: {mapped_count}/{len(df)} genomes matched")
    else:
        print(f"⚠️ Warning: Metadata file not found at {metadata_csv}")
        logging.warning(f"Metadata file not found: {metadata_csv}")

    # ------------------------------------------
    # STEP 5: Export final registry
    # ------------------------------------------
    final_cols = [
        'accession', 'tax_id', 'strain', 'host_clean', 'source_type', 'country_clean', 'diversity_group', 'year_clean', 'Threat_Category',
        'AMR_Burden', 'Is_MDR', 'Last_Resort_AMR_Count', 'Standard_AMR_Count', 'Last_Resort_AMR_Genes', 'Standard_AMR_Genes',
        'Virulence_Burden', 'High_Risk_Stx_Count', 'High_Risk_Stx_Genes', 'All_Virulence_Genes'
    ]
    
    export_cols = [c for c in final_cols if c in df.columns]
    output_file_path = os.path.join(output_dir_path, 'Final_Genomic_Threat_Registry.csv')
    
    df[export_cols].to_csv(output_file_path, index=False)
    logging.info(f"Final registry exported: {output_file_path}")
    logging.info("===== THREAT REGISTRY FINALIZATION COMPLETED =====")
    print(f"✅ Success! Registry saved to {output_file_path}")


# In[3]:


def run_statistical_analysis(registry_csv_path, output_dir_path):
    """
    Perform statistical analysis and visualization of genomic threat categories.

    This function performs the following:
    1. Loads the final genomic threat registry dataset.
    2. Computes descriptive statistics (mean, std, quartiles) for:
        - Virulence burden
        - AMR burden
       grouped by threat category.
    3. Performs Kruskal-Wallis statistical tests to assess whether:
        - Virulence burden differs across categories
        - AMR burden differs across categories
    4. Generates boxplot visualizations with overlaid data points.
    5. Extracts and saves high-risk genomes (Category 2 and Category 3).

    Parameters:
    ----------
    registry_csv_path : str
        Path to the final genomic threat registry CSV file.

    output_dir_path : str
        Directory where output files (statistics, plots, filtered data) will be saved.

    Outputs:
    -------
    - Threat_Category_Basic_Statistics.csv : Descriptive statistics table
    - Threat_Category_Burden_Boxplots.png : Boxplot visualization
    - High_Risk_Category2_3_Genome.csv : Filtered Category 2 & 3 dataset
    """
    # ==========================================
    # INITIAL SETUP
    # ==========================================
    print("Loading Final Threat Registry...")
    df = pd.read_csv(registry_csv_path)
    
    # Ensure Threat Category is treated as a string for grouping and plotting
    df['Threat_Category'] = df['Threat_Category'].astype(int).astype(str)
    
    # ==========================================
    # 1. DESCRIPTIVE STATISTICS
    # ==========================================
    print("\n--- BASIC STATISTICS BY THREAT CATEGORY ---")
    # Compute summary statistics for virulence and AMR burden
    stats_df = df.groupby('Threat_Category')[['Virulence_Burden', 'AMR_Burden']].describe()
    
    stats_out_path = os.path.join(output_dir_path, 'Threat_Category_Basic_Statistics.csv')
    stats_df.to_csv(stats_out_path)
    print(f"Statistics saved to: {stats_out_path}")
    # ==========================================
    # 2. KRUSKAL-WALLIS STATISTICAL TEST
    # ==========================================
    print("\n--- KRUSKAL-WALLIS STATISTICAL VALIDATION ---")
    categories = sorted(df['Threat_Category'].unique()) # Get sorted list of categories (e.g., ['1', '2', '3'])
    
    vir_groups = [df[df['Threat_Category'] == cat]['Virulence_Burden'].dropna() for cat in categories] # Prepare grouped data for virulence burden
    h_stat_vir, p_val_vir = kruskal(*vir_groups) # Perform Kruskal-Wallis test (non-parametric comparison)
    print(f"Virulence Burden -> H-statistic: {h_stat_vir:.2f}, p-value: {p_val_vir:.2e}")
    # Prepare grouped data for AMR burden
    amr_groups = [df[df['Threat_Category'] == cat]['AMR_Burden'].dropna() for cat in categories]
    h_stat_amr, p_val_amr = kruskal(*amr_groups)
    print(f"AMR Burden       -> H-statistic: {h_stat_amr:.2f}, p-value: {p_val_amr:.2e}")

    # ==========================================
    # 3. BOXPLOT VISUALIZATION
    # ==========================================
    print("\n--- GENERATING BOXPLOTS ---")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # Boxplot 1: Virulence Burden
    # FIXED: Added hue='Threat_Category' and legend=False
    # ADDED: linewidth=2 for bolder, darker boundaries
    sns.boxplot(
        x='Threat_Category', y='Virulence_Burden', data=df, 
        ax=ax1, palette='Oranges', order=categories,
        hue='Threat_Category', legend=False, 
        linewidth=2.5
    )
    sns.stripplot(x='Threat_Category', y='Virulence_Burden', data=df, ax=ax1, color='black', alpha=0.3, size=4, order=categories)
    ax1.set_title(f'Virulence Burden by Threat Category\n(Kruskal-Wallis p={p_val_vir:.2e})', fontweight='bold', fontsize=14)
    ax1.set_xlabel('Threat Category', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Total Virulence Factors', fontsize=12, fontweight='bold')
    ax1.grid(axis='y', linestyle='--', alpha=0.5)

    # Boxplot 2: AMR Burden
    # FIXED: Added hue='Threat_Category' and legend=False
    # ADDED: linewidth=2 for bolder, darker boundaries
    sns.boxplot(
        x='Threat_Category', y='AMR_Burden', data=df, 
        ax=ax2, palette='Blues', order=categories,
        hue='Threat_Category', legend=False,
        linewidth=2.5
    )
    sns.stripplot(x='Threat_Category', y='AMR_Burden', data=df, ax=ax2, color='black', alpha=0.3, size=4, order=categories)
    ax2.set_title(f'AMR Burden by Threat Category\n(Kruskal-Wallis p={p_val_amr:.2e})', fontweight='bold', fontsize=14)
    ax2.set_xlabel('Threat Category', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Total Resistance Genes', fontsize=12, fontweight='bold')
    ax2.grid(axis='y', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plot_out_path = os.path.join(output_dir_path, 'Threat_Category_Burden_Boxplots.png')
    plt.savefig(plot_out_path, dpi=300)
    print(f"Plots saved to: {plot_out_path}")

    # ==========================================
    # 4. EXTRACT HIGH-RISK GENOMES (CATEGORY 2 & 3)
    # ==========================================
   
    print("\n--- EXTRACTING HIGH-RISK RESERVOIRS (CAT 2 & 3) ---")
    cat2_3_df = df[df['Threat_Category'].isin(['2', '3'])].copy()

    
    cat2_3_out_path = os.path.join(output_dir_path, 'High_Risk_Category2_3_Genome.csv')
    cat2_3_df.to_csv(cat2_3_out_path, index=False)
    
    print(f"✅ Successfully extracted {len(cat2_3_df)} Category 2 and 3 genomes.")
    print(f"✅ Category 2 and 3 data saved to: {cat2_3_out_path}")



# In[4]:


def analyze_high_risk_reservoir(cat2_3_csv_path, output_dir_path):
    """
    Analyze and visualize high-risk genomic reservoirs (Category 2 and Category 3).

    This function performs the following:
    1. Loads a pre-classified dataset containing Category 2 and Category 3 genomes.
    2. Identifies:
       - All Category 3 genomes (highest-risk group)
       - Elite Category 2 genomes based on:
            a) High last-resort AMR burden
            b) High virulence (especially stx2a/stx2d presence)
    3. Combines these into a final high-risk genome list and exports it as a CSV file.
    4. Generates a bubble scatterplot:
       - X-axis: Virulence burden
       - Y-axis: AMR burden
       - Bubble size: Number of last-resort AMR genes
       - Highlights genomes with high-risk Stx genes (stx2a/stx2d)

    Parameters:
    ----------
    cat2_3_csv_path : str
        Path to the input CSV file containing Category 2 and Category 3 genome data.

    output_dir_path : str
        Directory where output files (CSV and plots) will be saved.

    Outputs:
    -------
    - Top_High_Risk_Genomes_Final.csv : Filtered high-risk genome dataset
    - Category2_3_Threat_Scatterplot.png : Bubble scatterplot visualization
    """
    # ==========================================
    # INITIAL SETUP
    # ==========================================
    if not os.path.exists(output_dir_path):
        os.makedirs(output_dir_path, exist_ok=True)

    print("Loading Category 2 & 3 Data...")
    df = pd.read_csv(cat2_3_csv_path) 
    df['has_stx'] = df['High_Risk_Stx_Count'] > 0 # Create flag for high-risk Stx genes (stx2a / stx2d)
    df['Threat_Category'] = df['Threat_Category'].astype(int) # Ensure Threat_Category is numeric for filtering

    ## ==========================================
    # 1. HIGH-RISK GENOME SELECTION (FIXED LOGIC)
    # ==========================================
    print("Extracting Category 3 and the absolute worst Category 2 genomes...")
    # Separate Category 3 (highest-risk) and Category 2
    cat3_df = df[df['Threat_Category'] == 3].copy()
    cat2_df = df[df['Threat_Category'] == 2].copy()

    
    # Select top Category 2 genomes based on AMR burden (strongest resistance reservoirs)
    cat2_top_amr = cat2_df.sort_values(
        by=['Last_Resort_AMR_Count', 'AMR_Burden'], ascending=[False, False]
    ).head(10)

    # Select top Category 2 genomes based on virulence (especially stx-positive strains)
    cat2_top_virulence = cat2_df.sort_values(
        by=['has_stx', 'Virulence_Burden'], ascending=[False, False]
    ).head(10)
    # Combine AMR-heavy and virulence-heavy selections
    # Remove duplicates (if a genome appears in both lists)
    cat2_elite = pd.concat([cat2_top_amr, cat2_top_virulence]).drop_duplicates()

    # Combine all Category 3 genomes with Elite Category 2 genomes
    high_risk_list = pd.concat([cat3_df, cat2_elite])
    
    summary_cols = [
        'accession', 'tax_id', 'strain', 'Threat_Category',
        'host_clean', 'source_type', 'country_clean', 'year_clean', 'diversity_group',
        'Virulence_Burden', 'AMR_Burden', 'Last_Resort_AMR_Count', 'Last_Resort_AMR_Genes', 
        'has_stx', 'High_Risk_Stx_Genes'
    ]
    # Keep only available columns (prevents errors if some columns are missing)
    existing_cols = [col for col in summary_cols if col in high_risk_list.columns]
    # Save final high-risk genome dataset
    csv_out_path = os.path.join(output_dir_path, 'Top_High_Risk_Genomes_Final.csv')
    high_risk_list[existing_cols].to_csv(csv_out_path, index=False)
    print(f"Final List saved to: {csv_out_path}")

    ## ==========================================
    # 2. BUBBLE SCATTER PLOT VISUALIZATION (LIGHTER & FIXED LEGENDS)
    # ==========================================
    print("Generating Clean Theme Bubble Scatterplot...")
    cat2 = df[df['Threat_Category'] == 2]
    cat3 = df[df['Threat_Category'] == 3]
    
    plt.figure(figsize=(12, 8))
    sns.set_theme(style="whitegrid", rc={"axes.edgecolor": ".15", "axes.linewidth": 1.25})
    
    # 🔴 Category 2 (Coral Red)
    sizes_cat2 = cat2['Last_Resort_AMR_Count'] * 60 + 50
    plt.scatter(
        cat2['Virulence_Burden'], cat2['AMR_Burden'],
        s=sizes_cat2, color='#E63946', alpha=0.6, 
        edgecolor='#900C3F', linewidth=1.2, label='_nolegend_'
    )
    
    # 🔵 Category 3 (Lighter, Professional Blue)
    sizes_cat3 = cat3['Last_Resort_AMR_Count'] * 60 + 50
    plt.scatter(
        cat3['Virulence_Burden'], cat3['AMR_Burden'],
        s=sizes_cat3, color='#4A90E2', alpha=0.85, # Lightened Blue
        edgecolor='#104E8B', linewidth=1.2, label='_nolegend_'
    )
    
    # ❌ Highlight genomes with high-risk Stx genes (Vivid Purple 'X')
    stx_df = df[df['has_stx'] == True]
    plt.scatter(
        stx_df['Virulence_Burden'], stx_df['AMR_Burden'],
        s=80, color='#8E44AD', marker='X', # Vivid Purple
        linewidth=1.5, label='_nolegend_'
    )
    
    plt.xlabel('Virulence Burden (Total Genes)', fontsize=14, fontweight='500', color='#333333')
    plt.ylabel('AMR Burden (Total Genes)', fontsize=14, fontweight='500', color='#333333')

    # ==========================================
    # LEGENDS (RIGHT-HAND SIDE, COMPACT & ALIGNED)
    # ==========================================
    # 1. Threat Category Legend (Placed in the Top-Right corner)
    legend_cat = plt.legend(
        handles=[
            Line2D([0], [0], marker='o', color='w', label='Category 2 (High Risk)', markerfacecolor='#E63946', markeredgecolor='#900C3F', markersize=11),
            Line2D([0], [0], marker='o', color='w', label='Category 3 (Super-Hybrid)', markerfacecolor='#4A90E2', markeredgecolor='#104E8B', markersize=11),
            Line2D([0], [0], marker='X', color='w', label='Has stx2a/2d Toxin', markerfacecolor='#8E44AD', markersize=10)
        ],
        title="Threat Classification", 
        loc='upper right', bbox_to_anchor=(0.99, 0.99), # Top Right
        frameon=True, shadow=True, edgecolor='gray'
    )
    plt.setp(legend_cat.get_title(), fontweight='bold')
    plt.gca().add_artist(legend_cat)
    
    # 2. Bubble Size Legend (Expanded list, made compact, moved UP)
    bubble_sizes = [1, 2, 3, 4, 5] # <-- Now shows 1 through 5
    legend_bubbles = plt.legend(
        handles=[
            plt.scatter([], [], s=(size * 60 + 50), color='#A8DADC', alpha=0.7, edgecolor='#457B9D', linewidth=1.2)
            for size in bubble_sizes
        ],
        labels=[f'{size} Genes' for size in bubble_sizes],
        title="Last Resort AMR Count", 
        loc='upper right', 
        bbox_to_anchor=(0.99, 0.82), # <-- Moved UP (was 0.82)
        labelspacing=0.9,            # <-- Reduced spacing to make it smaller
        borderpad=0.8,               # <-- Reduced padding inside the box
        fontsize=10,                 # <-- Slightly smaller text
        frameon=True, shadow=True, edgecolor='gray'
    )
    plt.setp(legend_bubbles.get_title(), fontweight='bold', fontsize=11)

    plt.title('Genomic Threat Landscape: Category 2 vs Category 3\nBubble Size = Last Resort AMR Genes', 
              fontsize=16, fontweight='bold', color='#111111', pad=15)
    
    # Save FIRST, then Show
    plot_out_path = os.path.join(output_dir_path, 'Category2_3_Threat_Scatterplot.png')
    plt.savefig(plot_out_path, dpi=300, bbox_inches='tight', transparent=False)
    plt.show()
    
    print(f"✅ Final Scatterplot successfully saved to: {plot_out_path}")



# In[5]:


def generate_accession_country_stacked_bar(top15_csv_path, output_dir):
    """
    Generate a stacked horizontal bar chart visualizing last-resort AMR gene
    distribution across high-risk genomes (Category 2 & Category 3).

    This function:
    - Loads a curated dataset of high-risk genomes
    - Creates informative labels including category, strain, accession, taxonomy ID, and diversity group
    - Detects high-risk toxin genes (stx2a / stx2d) and annotates labels
    - Extracts and cleans last-resort AMR gene data
    - Constructs a gene-by-genome matrix (cross-tabulation)
    - Generates a stacked bar chart where:
        • Each bar = one genome
        • Each segment = one AMR gene
    - Highlights:
        • Category 3 strains (bold navy)
        • Toxin-positive strains (bold italic purple)
    - Saves a high-resolution plot

    Parameters:
    ----------
    top15_csv_path : str
        Path to CSV containing selected high-risk genomes

    output_dir : str
        Directory where output plot will be saved

    Output:
    -------
    Accession_Country_StackedBar_Final.png
    """

    # ==========================================
    # INITIAL SETUP
    # ==========================================
    print("Loading High-Risk Genomes (Category 2 & 3)...")
    df = pd.read_csv(top15_csv_path)
    
    # Ensure Threat Category is an integer for consistent formatting
    df['Threat_Category'] = df['Threat_Category'].astype(int)
    
    # ==========================================
    # 1. DATA PREPARATION & LABEL CREATION
    # ==========================================
    gene_data = []      # Stores gene-level records (long format)
    ordered_labels = [] # Maintains original ranking order for plotting
    
    for idx, row in df.iterrows():
        # --- Extract Base Information ---
        # country = str(row['country_clean']).strip()
        diversity = str(row['diversity_group']).strip()  # Using diversity group instead of country
        strain = str(row['strain']).strip()
        accession = str(row['accession']).strip()
        cat = row['Threat_Category']
        
        # --- Extract and clean Taxonomy ID ---
        tax_id = str(row.get('tax_id', 'Unknown')).replace('.0', '').strip()
        
        # ==========================================
        # DYNAMIC TOXIN DETECTION (stx2a / stx2d)
        # ==========================================
        # Convert gene list to lowercase for consistent matching
        vir_genes = str(row.get('High_Risk_Stx_Genes', '')).lower()
        
        found_toxins = []
        if 'stx2a' in vir_genes:
            found_toxins.append('stx2a')
        if 'stx2d' in vir_genes:
            found_toxins.append('stx2d')
            
        # Create toxin marker for visualization
        if found_toxins:
            toxin_marker = f" [+{'/'.join(found_toxins)}]"
        elif row.get('has_stx', False):
            toxin_marker = " [+stx2]"  # fallback marker
        else:
            toxin_marker = ""
        
        # ==========================================
        # FINAL LABEL (USED ON Y-AXIS)
        # ==========================================
        label = f"Cat {cat} | {strain} (Acc: {accession} | TaxID: {tax_id}) | {diversity}{toxin_marker}"
        ordered_labels.append(label)
        
        # ==========================================
        # EXTRACT LAST-RESORT AMR GENES
        # ==========================================
        genes_raw_str = str(row['Last_Resort_AMR_Genes'])
        
        # Skip genomes with no AMR genes
        if pd.isna(genes_raw_str) or genes_raw_str.strip() == "" or genes_raw_str.strip() == "nan":
            continue
            
        # Split comma-separated gene list
        genes_raw = genes_raw_str.split(',')
        genes = [g.strip() for g in genes_raw if g.strip()]
        
        # Clean gene names (remove suffix like "_1")
        for g in genes:
            clean_g = g.rsplit('_', 1)[0] if '_' in g else g
            gene_data.append({'Accession_Label': label, 'Gene': clean_g})

    # Safety check: stop if no AMR genes found
    if not gene_data:
        print("Warning: No Last Resort Genes found in this dataset to plot.")
        return

    # ==========================================
    # 2. CREATE GENE MATRIX (CROSSTAB)
    # ==========================================
    df_genes = pd.DataFrame(gene_data)

    # Create matrix: rows = genomes, columns = genes
    pivot_df = pd.crosstab(df_genes['Accession_Label'], df_genes['Gene'])

    # Align with original ranking order and fill missing values
    pivot_df = pivot_df.reindex(ordered_labels).fillna(0)

    # Reverse order so highest-risk genomes appear at the top
    pivot_df = pivot_df.iloc[::-1]

    # ==========================================
    # 3. PLOT GENERATION
    # ==========================================
    print("Generating Stacked Bar Chart...")
    plt.figure(figsize=(16, 10))  # Wider figure for long labels
    sns.set_theme(style="whitegrid", rc={"axes.edgecolor": ".15", "axes.linewidth": 1.25})
    
    # ==========================================
    # COLOR PALETTE (HIGH CONTRAST)
    # ==========================================
    # Distinct colors ensure each AMR gene is visually separable
    distinct_colors = [
        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
        '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
        '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
        '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'
    ]
    
    # Match number of colors to number of genes
    num_genes = len(pivot_df.columns)
    cmap = ListedColormap(distinct_colors[:num_genes])

    # Create stacked horizontal bar chart
    ax = pivot_df.plot(
        kind='barh', 
        stacked=True, 
        colormap=cmap, 
        ax=plt.gca(), 
        edgecolor='black', 
        linewidth=1.0
    )

    # ==========================================
    # 4. AESTHETIC FORMATTING
    # ==========================================
    plt.title(
        'The Untreatable Arsenal: Last-Resort AMR & Toxin Profiles', 
        fontsize=18, fontweight='bold', pad=15, color='#111111'
    )
    
    plt.xlabel('Total Critical "Last-Resort" AMR Genes', fontsize=14, fontweight='bold')
    plt.ylabel('Threat Category | Strain (Acc | TaxID) | Diversity_Group', fontsize=14, fontweight='bold')

    # Legend configuration
    plt.legend(
        title='Last Resort AMR Genes', 
        bbox_to_anchor=(1.02, 1), 
        loc='upper left', 
        fontsize=12, 
        title_fontsize=13,
        frameon=True, shadow=True, edgecolor='gray'
    )

    # Ensure integer ticks on X-axis
    max_genes = int(pivot_df.sum(axis=1).max())
    plt.xticks(range(0, max_genes + 2), fontsize=12)
    
    # ==========================================
    # 5. DYNAMIC LABEL HIGHLIGHTING
    # ==========================================
    plt.draw()  # Ensure labels are rendered before modifying
    
    for text in ax.get_yticklabels():
        label_text = text.get_text()
        
        # Highlight Category 3 genomes
        if "Cat 3" in label_text:
            text.set_color('#1D3557')
            text.set_fontweight('bold')
        
        # Highlight toxin-positive genomes
        elif "[+stx" in label_text:
            text.set_color('#8E44AD')
            text.set_fontweight('bold')
            text.set_fontstyle('italic')

    plt.yticks(fontsize=11)

    # ==========================================
    # 6. SAVE OUTPUT
    # ==========================================
    plt.tight_layout()
    plot_out_path = os.path.join(output_dir, 'Accession_Country_StackedBar.png')

    plt.savefig(plot_out_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"✅ Stacked Bar chart successfully saved to: {plot_out_path}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--classification", required=True, help="Path to final threat classification CSV file")
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV file containing genome information")
    parser.add_argument("--outdir", required=True, help="Output directory to save classification statistics, plots, and registry files")

    args = parser.parse_args()

    setup_logging(args.outdir, "classification_statistics.log")

    # Step 1: Registry
    finalize_threat_registry(
        args.classification,
        args.metadata,
        args.outdir
    )

    registry_path = os.path.join(args.outdir, "Final_Genomic_Threat_Registry.csv")

    # Step 2: Stats
    run_statistical_analysis(registry_path, args.outdir)

    # Step 3: High-risk analysis
    cat2_3_path = os.path.join(args.outdir, "High_Risk_Category2_3_Genome.csv")
    analyze_high_risk_reservoir(cat2_3_path, args.outdir)

    # Step 4: Final plot
    top_path = os.path.join(args.outdir, "Top_High_Risk_Genomes_Final.csv")
    generate_accession_country_stacked_bar(top_path, args.outdir)