#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


# In[5]:


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
# 3. THREAT CLASSIFICATION LOGIC
# ==========================================
def assign_threat_level(row):
    """
     Determines the threat category for a single genome based on
     virulence gene and AMR gene presence.
    """
    #################### Define AMR Gene Groups ####################
    # Last-resort antibiotic resistance genes
    last_resort_amr = ['blaCTX-M', 'mcr', 'blaNDM', 'blaKPC', 'blaOXA-48']
    # Common antimicrobial resistance genes
    common_amr = ['tetA', 'tetB', 'tetM', 'sul1', 'sul2', 'sul3', 'strA', 'strB', 'dfrA1', 'dfrA12', 
                  'blaTEM-1', 'aph(3\')-Ia', 'aadA1', 'aadA2', 'floR', 'catA1']
    
     #################### Helper Function ####################
    # To check if ANY of the genes in a list are present (partial match)
    def has_genes(gene_list):
        """
        Checks if any gene from a given list is present
        in the genome row.
        """
        for col_name in row.index:
            if str(row[col_name]).strip() != '.':
                if any(target_gene in col_name for target_gene in gene_list):
                    return True
        return False
    #################### Detect Virulence Factors ####################   
    #Check for specific conditions
    has_high_risk_stx = False
    has_stx1 = False
    for col_name in row.index:
        col_lower = col_name.lower()
        if str(row[col_name]).strip() != '.':
            if 'stx2a' in col_lower or 'stx2b' in col_lower:
                has_high_risk_stx = True
            if 'stx1' in col_lower:
                has_stx1 = True
                
    #################### Detect AMR Genes ##################### 
    has_critical_amr = has_genes(last_resort_amr) #Evaluate Critical AMR presence
    
    # Evaluate Standard AMR (True if there is AT LEAST ONE common AMR gene)
    has_standard_amr = has_genes(common_amr)

    #################### MDR Detection ####################
    # Count how many standard AMR genes are present
    amr_count = 0
    for col_name in row.index:
        if str(row[col_name]).strip() != '.':
            if any(target_gene in col_name for target_gene in common_amr):
                amr_count += 1
                
    # Define MDR: True if the bacteria has 3 or more standard resistance genes
    is_MDR = (amr_count >= 3)

    #################### Apply Threat Logic ####################
   
    # Category 3: Super-Hybrid (High-risk stx + Last-resort AMR)
    if has_high_risk_stx and has_critical_amr:
        return 3
    
    # Category 2: High-Risk 
    # (Matches the table: stx2 + MDR. OR catches edge cases: stx2 alone, or critical AMR alone)
    elif (has_high_risk_stx and is_MDR) or has_high_risk_stx or has_critical_amr:
        return 2
    
    # Category 1: Typical Pathogen 
    # (Matches the table: stx1 + Common AMR (just 1 is enough). OR catches edge case of stx1 alone)
    elif (has_stx1 and has_standard_amr) or has_stx1:
        return 1
    
    # Category 0: Baseline 
    else:
        return 0

def classify_genomes(vir_csv_path, amr_csv_path, output_csv_path):
    """
    Performs genome threat classification using virulence and AMR gene results
    generated from Abricate analysis.

    Pipeline Steps:
    1. Creates a 'threat_classification_analysis' output directory.
    2. Initializes a logging file to track the analysis and capture errors.
    3. Reads virulence gene results and antimicrobial resistance (AMR) results.
    4. Merges both datasets using genome file names.
    5. Applies the threat classification logic to each genome.
    6. Assigns a threat category score (0–3) based on virulence and AMR gene patterns.
    7. Saves the final classified genome dataset.

    Threat Categories:
        3 : Super-Hybrid Threat (High-risk virulence + Last-resort AMR)
        2 : High Risk (High-risk virulence OR Critical AMR OR MDR)
        1 : Typical Pathogen (stx1 with common AMR genes)
        0 : Baseline Risk (No major virulence or resistance detected)

    Arguments:
        output_folder_path (str): Base directory where output folder will be created.
        vir_csv_path (str): Path to the virulence results CSV file.
        amr_csv_path (str): Path to the antimicrobial resistance results CSV file.

    Outputs Generated:
        - threat_classification_analysis/threat_classification.log
        - threat_classification_analysis/final_threat_classification.csv
    """  
    try:
        ############################ START ANALYSIS ###################################
        logging.info("Starting genome threat classification...")
        
        ############################ Reading Input Files ############################## 
        vir_df = pd.read_csv(vir_csv_path, sep='\t')
        amr_df = pd.read_csv(amr_csv_path, sep='\t')
        logging.info(f"Loaded virulence data ({len(vir_df)} rows) and AMR data ({len(amr_df)} rows).")
        ############################ Dataset Merging ##################################
        # Merge data using a left join to keep all 465 genomes, and fill missing AMR hits with a dot '.'
        df = pd.merge(vir_df, amr_df, on='#FILE', how='left').fillna('.')
        logging.info("Successfully merged Virulence and AMR datasets.")
        
        ############################ Apply Classification ############################
        df['Threat_Category'] = df.apply(assign_threat_level, axis=1)
        
        # LOGGING: Record the successful classification inside the function
        logging.info(f"Successfully processed and assigned threat levels to {len(df)} genomes.") 

        # Reorder columns to put Threat_Category in the 2nd position
        cols = df.columns.tolist()
        cols.insert(1, cols.pop(cols.index('Threat_Category')))
        df = df[cols]

        # Save results
        df.to_csv(output_csv_path, index=False)
        logging.info(f"Analysis complete. Results saved to {output_csv_path}")
        
        print(f"Analysis Complete! Results saved to {output_csv_path}")
        print(df[['#FILE', 'Threat_Category']].head())
        
        return df 
    ############################ ERROR HANDLING ############################
    except Exception as e:
        # LOGGING: Log the error if something fails
        logging.error(f"ERROR in Classification Script: {str(e)}")
        print(f"An error occurred. Check the log file: {e}")
        return None




# In[4]:


def visualize_threat_distribution(plot_df, output_folder_path):
    """
    Generates visualization plots showing the distribution of genome
    threat categories.

    Pipeline Steps:
    1. Logs dataset structure and gene profile information.
    2. Calculates counts for each Threat_Category.
    3. Generates a bar chart showing the number of genomes per category.
    4. Generates a pie chart showing percentage distribution of categories.
    5. Saves the plots into the threat_classification_analysis directory.

    Arguments:
        plot_df (DataFrame): DataFrame containing genome threat classification results.
        output_folder_path (str): Directory where visualization plots will be saved.

    Outputs Generated:
        - genome_threat_distribution_barplot.png
        - genome_threat_distribution_piechart.png
    """
    try:
        ############################ START VISUALIZATION ############################
        logging.info("Starting visualization of Threat Categories...")

        # --- LOGGING: Advanced Rows and Columns Identification ---
        cols = plot_df.columns.tolist()
        
        # 1. Locate the "Markers"
        idx_x = cols.index('NUM_FOUND_x')
        idx_y = cols.index('NUM_FOUND_y')
        
        # 2. Calculate segments
        # Virulence genes are between the two NUM_FOUND markers
        num_virulence = idx_y - idx_x - 1
        
        # AMR genes are everything after the second NUM_FOUND marker
        num_amr = len(cols) - idx_y - 1
        
        # 3. Total for verification
        num_total_genes = num_virulence + num_amr
        
        # 4. Detailed Logging
        logging.info(f"Data Profile: {plot_df.shape[0]} Genomes analyzed.")
        logging.info(f"Virulence Profile: {num_virulence} genes detected (Markers: NUM_FOUND_x).")
        logging.info(f"AMR Profile: {num_amr} genes detected (Markers: NUM_FOUND_y).")
        logging.info(f"Total Unique Targets: {num_total_genes} genes.")
        
        ############################ CATEGORY DISTRIBUTION ############################
        # Calculate counts for each category
        category_counts = plot_df['Threat_Category'].value_counts().sort_index()
        
        for cat, count in category_counts.items():
            logging.info(f"Category {cat} count for visualization: {count}")

        ############################ CATEGORY LABELS & COLORS ############################
        # Map labels and colors DYNAMICALLY based on what is actually in the data
        category_labels = {
            0: 'Category 0\n(Baseline/Low Risk)',
            1: 'Category 1\n(Typical Pathogen)',
            2: 'Category 2\n(High-Risk)',
            3: 'Category 3\n(Super-Hybrid)'
        } 
        
        category_colors = {
            0: '#2ca02c', # Green
            1: '#ff7f0e', # Orange
            2: '#d62728', # Red
            3: '#9467bd'  # Purple
        }

        # Only grab the labels and colors for categories that actually exist in your dataset
        labels = [category_labels[cat] for cat in category_counts.index] 
        colors = [category_colors[cat] for cat in category_counts.index]
        counts = category_counts.values 
        
        ############################ BAR CHART ############################ 
        logging.info("Generating Threat Category Bar Chart...")
        plt.figure(figsize=(10, 6)) 
        
        # UPDATED: Added hue=labels and legend=False to fix the FutureWarning
        sns.barplot(x=labels, y=counts, hue=labels, palette=colors, legend=False, edgecolor='black', linewidth=1.5) 
        
        plt.title('Number of Genomes per Threat Category', fontsize=16, fontweight='bold') 
        plt.ylabel('Number of Genomes', fontsize=14) 
        plt.xlabel('Threat Category', fontsize=14) 
        
        # Add value labels on top of bars
        for i, count in enumerate(counts): 
            plt.text(i, count + 3, str(count), ha='center', va='bottom', fontsize=12, fontweight='bold') 
            
        plt.tight_layout()
        
        #SAVE THE FIGURE (BEFORE plt.show!)
        plot_file = os.path.join(output_folder_path, "genome_threat_distribution_barplot.png")
        plt.savefig(plot_file, dpi=300)
        
        #Display the Figure
        plt.show() 
        logging.info(f"Bar Chart generated successfully and saved to {plot_file}")


        ############################ PIE CHART ############################ 
        logging.info("Generating Threat Category Pie Chart...")
        plt.figure(figsize=(8, 8)) 
        plt.pie(counts, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors, 
                wedgeprops={'edgecolor': 'black', 'linewidth': 1.5}, textprops={'fontsize': 12}) 
        plt.title('Distribution of E. coli Threat Categories', fontsize=16, fontweight='bold')
        
        plt.tight_layout()

        #SAVE THE FIGURE (BEFORE plt.show!)
        plot_file = os.path.join(output_folder_path, "genome_threat_distribution_piechart.png")
        plt.savefig(plot_file, dpi=300)

        #Display the Figure
        plt.show()
        logging.info(f"Bar Chart generated successfully and saved to {plot_file}")
    ############################ ERROR HANDLING ############################
    except Exception as e:
        logging.error(f"ERROR during visualization: {str(e)}")
        print(f"An error occurred during plotting. Check the log file: {e}")

## -------------------------------
# MAIN (ARGPARSE)
# -------------------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--vir", required=True, help="Path to virulence results file (Abricate VFDB summary CSV)")
    parser.add_argument("--amr", required=True, help="Path to AMR results file (Abricate ResFinder summary CSV)")
    parser.add_argument("--out", required=True, help="Output file path for final threat classification CSV")
    parser.add_argument("--logdir", required=True, help="Directory to store log files for threat classification analysis")

    args = parser.parse_args()

    # ✅ ADD THIS LINE (start logging)
    setup_logging(args.logdir, "threat_classification.log")

    classify_genomes(
        args.vir,
        args.amr,
        args.out
    )

    visualize_threat_distribution(
        pd.read_csv(args.out),
        os.path.dirname(args.out)
    )