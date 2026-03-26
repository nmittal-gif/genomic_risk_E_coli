#!/usr/bin/env python
# coding: utf-8

# In[6]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import json
import os
import argparse


def setup_logging(logs_dir, log_file_name):
    """
    Creates the target directory if it doesn't exist and configures the logger.
    """
    # 1. Ensure the folder exists before creating the file
    os.makedirs(logs_dir, exist_ok=True)
    
    # 2. CLEAR HANDLERS: This is required so Jupyter lets go of the previous log file
    # when you switch to a new function!
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    # 3. Configure the logger
    logging.basicConfig(
        filename=os.path.join(logs_dir, log_file_name),
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s", 
    )
    
    logging.info(f"Logging initialized. Saving to: {logs_dir}")


# In[11]:


def analyze_assembly_metadata(output_folder_path, taxid_json_path, file_path):
    
    """
    Analyzes genomic assembly metadata, maps taxonomic IDs to organism names, 
    and generates summary statistics and distribution plots.

    Pipeline Steps:
    1. Creates a 'metadata_analysis' output directory.
    2. Initializes a logging file to track the analysis and catch errors.
    3. Reads assembly metadata (CSV) and taxonomy mapping (JSON).
    4. Calculates total assemblies and counts per unique Taxonomic ID.
    5. Saves a summary table ('taxid_genome_summary.csv').
    6. Generates and saves a bar chart distribution ('taxid_distribution.png').

    Arguments:
        output_folder_path (str): The base directory where the output folder will be created.
        taxid_json_path (str): The path to the JSON file containing TaxID to Organism mappings.
        file_path (str): The path to the input metadata CSV file.

    Outputs Generated:
        - metadata_analysis/metadata_analysis.log
        - metadata_analysis/taxid_genome_summary.csv
        - metadata_analysis/taxid_distribution.png
    """
    
    # 1. Join the target folder path with the new "metadata_analysis" folder name
    metadata_dir = os.path.join(output_folder_path, "metadata_analysis")
    
    # 2. Create the folder safely 
    os.makedirs(metadata_dir, exist_ok=True)
    
    
    # 3. Call your logging setup function
    log_file_name = 'metadata_analysis.log'
    setup_logging(metadata_dir, log_file_name)
    
    # 4. Proceed with your standard processing
    print(f"Successfully created/verified metadata_analysis directory at: {metadata_dir}")
    logging.info(f"Successfully created/verified metadata_analysis directory at: {metadata_dir}")
    
    print("===== METADATA ANALYSIS STARTED =====")
    logging.info("===== METADATA ANALYSIS STARTED =====")
    
    ############################ Reading Files ###################################
    try:
        # Load metadata csv file and read it through pandas as dataframe
        metadata_df = pd.read_csv(file_path)
        print("The metadata file is: \n", metadata_df.head())
        logging.info(f"Metadata file loaded: {file_path}")
        
        print(f"The columns in the metadata file are:\n {metadata_df.columns}")
        
        logging.info(f"The columns in the metadata file are: {list(metadata_df.columns)}")

        # Load Organisms taxid json file and read it
        with open(taxid_json_path, 'r') as file:
            taxid_json = json.load(file) 
            print("The taxid json file is: \n", taxid_json)
        logging.info("Successfully loaded JSON data.")
    
     ##############################  Metadata statistics ############################
        total_assemblies = len(metadata_df) #total rows of the metadata file;
        unique_taxids = metadata_df["tax_id"].nunique() #count of total assemblies in each unique ncbi_tax_id;
        
        # Create a summary table containing basic statistics of the metadata file
        summary_table = pd.DataFrame({"Metric": ["Genomic Assemblies", "Unique TaxIDs"],
            "Value": [total_assemblies, unique_taxids]
        })
        print("The summary of metadata: \n", summary_table)
        logging.info("\nMetadata Summary Table:\n" + summary_table.to_string(index=False))
        
     ############################## Genomic assembliess count per ncbi_tax_id ###################
        taxid_counts = metadata_df["tax_id"].value_counts() #count of genomic assemblies in each ncbi id 
        taxid_table = taxid_counts.reset_index() # Convert the value_counts result (Series) into a DataFrame
        taxid_table.columns = ["TaxID", "Genomic Assemblies"]
        # Map TaxID values to organism names using the JSON dictionary
        # TaxID is converted to string because JSON keys are stored as strings
        taxid_table["Organism"] = taxid_table["TaxID"].astype(str).map(taxid_json) 
        
        print("The ncbi taxid mapped to organism name table is: \n", taxid_table)
        logging.info("\nGenome Count per TaxID:\n" + taxid_table.to_string(index=False))
        
     ############################ Save TaxID summary #########################################
        taxid_table_file = os.path.join(metadata_dir, "taxid_genome_summary.csv")  # Define the file path for saving the TaxID genome summary table
        taxid_table.to_csv(taxid_table_file, index=False) # Save the TaxID summary table as a CSV file
        
     ############################ Plot Generation ############################################
        
        logging.info("Generating TaxID distribution plot")
        # taxid_table["Label"] = taxid_table["Organism"] + " (" + taxid_table["TaxID"].astype(str) + ")"
        sns.set_theme(style="whitegrid")# Set seaborn theme for better visualization
        plt.figure(figsize=(10,6)) # Create a figure with defined size
        # Create a bar plot showing number of genome assemblies per TaxID
        # hue="Organism" helps distinguish organisms with colors
        ax = sns.barplot(data=taxid_table, x="TaxID", y="Genomic Assemblies", hue="Organism", palette="Set2")
        # ax = sns.barplot(data=taxid_table, x="TaxID", y="Genomic Assemblies", palette="Set2")
        for bar in ax.patches:
            bar.set_width(0.5)# Reduce the width of each bar for improved visualization
        # Add values on bars
        for container in ax.containers:
            ax.bar_label(container, fontsize=11, padding=3)  # Add numeric labels on top of each bar

        plt.xlabel("Taxonomic ID", fontsize=14, fontweight="bold")
        plt.ylabel("Number of Genomes", fontsize=14, fontweight="bold")
        plt.title("Distribution of Genome Assemblies by Taxonomic ID",fontsize=16,fontweight="bold")
        plt.xticks(rotation=45, fontsize=12) # Rotate x-axis labels for better readability
        plt.yticks(fontsize=12) # Set y-axis font size
        ax.grid(True, axis="y", linestyle="--", alpha=0.6) # Add grid lines for easier comparison
        plt.legend(title="Organism name (NCBI TaxID)",loc="upper right")# Add legend showing organism names
        plt.tight_layout() # Adjust layout to prevent overlapping labels
        plot_file = os.path.join(metadata_dir, "taxid_distribution.png")
        plt.savefig(plot_file, dpi=300)
        plt.close()
        logging.info(f"Plot saved: {plot_file}")
    
    ############################ Error Handling ############################ 
    # Catch and log any errors that occur during metadata analysis
    except Exception as e:
        logging.error(f"ERROR during metadata analysis: {e}")
        print("Error occurred. Check logs/metadata_analysis.log")


# In[12]:


def analyze_genomes_length_and_gc(output_folder_path, genome_directory):
    """
    Parses a directory of FASTA-formatted genomic files (.fna) to calculate 
    sequence lengths and GC content percentages in a single, highly efficient pass.
    
    This function reads each genome sequence line-by-line to minimize memory usage, 
    making it safe for directories containing hundreds of large genomic files. It 
    calculates summary statistics, saves a combined metrics table, and generates 
    visual distributions for both metrics.

    Pipeline Steps:
    1. Creates a 'genomic_analysis' directory inside the target output folder.
    2. Initializes a custom log file ('genomic_analysis.log') to track progress and errors.
    3. Scans the provided directory for all files ending in '.fna'.
    4. Iterates through each file, skipping headers ('>'), to count total bases and G/C bases.
    5. Calculates the GC percentage: (G + C) / Total Length * 100.
    6. Compiles the results into a pandas DataFrame and calculates min/max/mean/median statistics.
    7. Exports the DataFrame to a CSV file.
    8. Generates and saves two Seaborn histograms (Length Distribution & GC% Distribution).

    Arguments:
        output_folder_path (str): The base directory where the 'genomic_analysis' 
                                  results folder will be created.
        genome_directory (str): The path to the folder containing the input .fna files.

    Outputs Generated in 'genomic_analysis/':
        - genomic_analysis.log: Execution logs and error tracking.
        - genome_metrics.csv: A table containing 'Genome_File', 'Genome_Length', 
                              and 'GC_Content_Percent' for every parsed genome.
        - genome_length_distribution.png: Histogram of genome lengths with mean/median lines.
        - gc_content_distribution.png: Histogram of GC percentages with mean/median lines.
    """
    
    # 1. Join the target folder path with the new "genomic_analysis" folder name
    genomic_dir = os.path.join(output_folder_path, "genomic_analysis")
    
    # 2. Create the folder safely 
    os.makedirs(genomic_dir, exist_ok=True)
    
    # 3. Call your logging setup function
    log_file_name = 'genomic_analysis.log'
    setup_logging(genomic_dir, log_file_name)
    
    # 4. Proceed with your standard processing
    print(f"Successfully created/verified genomic_analysis directory at: {genomic_dir}")
    logging.info(f"Successfully created/verified genomic_analysis directory at: {genomic_dir}")
    
    print("===== GENOMIC ANALYSIS STARTED =====")
    logging.info("===== GENOMIC ANALYSIS STARTED =====")

    try:
        # Lists to store our data
        genome_files = []
        genome_lengths = []
        gc_contents = []

        # Scan genome directory for FASTA genome files (.fna)
        for file in os.listdir(genome_directory):
            if file.endswith(".fna"):
                file_path = os.path.join(genome_directory, file)

                # Initialize counters
                genome_length = 0
                gc_count = 0 

                # Open the FASTA genome file
                with open(file_path, "r") as f:
                    for line in f:
                        if not line.startswith(">"):
                            seq = line.strip().upper()
                            genome_length += len(seq)
                            gc_count += seq.count('G') + seq.count('C')

                # Calculate GC percentage
                gc_percent = (gc_count / genome_length) * 100 if genome_length > 0 else 0

                # Store the results
                genome_files.append(file)
                genome_lengths.append(genome_length)
                gc_contents.append(gc_percent)

        # ------------------------------------------------------------------
        # Create a pandas DataFrame
        # ------------------------------------------------------------------
        genome_df = pd.DataFrame({
            "Genome_File": genome_files,
            "Genome_Length": genome_lengths,
            "GC_Content_Percent": gc_contents
        })

        # ------------------------------------------------------------------
        # Calculate statistics (Restored fully!)
        # ------------------------------------------------------------------
        min_len = genome_df["Genome_Length"].min()
        max_len = genome_df["Genome_Length"].max()
        mean_len = genome_df["Genome_Length"].mean()
        median_len = genome_df["Genome_Length"].median()
        
        min_gc = genome_df["GC_Content_Percent"].min()
        max_gc = genome_df["GC_Content_Percent"].max()
        mean_gc = genome_df["GC_Content_Percent"].mean()
        median_gc = genome_df["GC_Content_Percent"].median()

        # Identify genomes with extreme lengths
        min_genome = genome_df.loc[genome_df["Genome_Length"].idxmin()]
        max_genome = genome_df.loc[genome_df["Genome_Length"].idxmax()]
        
        # Identify genomes with extreme GC content
        min_genome_gc = genome_df.loc[genome_df["GC_Content_Percent"].idxmin()]
        max_genome_gc = genome_df.loc[genome_df["GC_Content_Percent"].idxmax()]

        # ------------------------------------------------------------------
        # Print results for user visibility
        # ------------------------------------------------------------------
        print("\n===== Genome Length Statistics =====")
        print(f"Minimum Genome Length : {min_len}")
        print(f"Maximum Genome Length : {max_len}")
        print(f"Mean Genome Length    : {mean_len:.2f}")
        print(f"Median Genome Length  : {median_len}")

        print("\nGenome with Minimum Length:")
        print(min_genome)

        print("\nGenome with Maximum Length:")
        print(max_genome)
        
        print("\n===== GC Content Statistics =====")
        print(f"Minimum GC Content : {min_gc:.2f}%")
        print(f"Maximum GC Content : {max_gc:.2f}%")
        print(f"Mean GC Content    : {mean_gc:.2f}%")
        print(f"Median GC Content  : {median_gc:.2f}%")

        # ------------------------------------------------------------------
        # Log the results (Restored fully!)
        # ------------------------------------------------------------------
        logging.info("--- Genome Length Statistics ---")
        logging.info(f"Minimum genome length: {min_len}")
        logging.info(f"Maximum genome length: {max_len}")
        logging.info(f"Mean genome length: {mean_len:.2f}")
        logging.info(f"Median genome length: {median_len}")
        logging.info(f"Genome with minimum length: {min_genome['Genome_File']}")
        logging.info(f"Genome with maximum length: {max_genome['Genome_File']}")
        
        logging.info("--- GC Content Statistics ---")
        logging.info(f"Minimum GC Content: {min_gc:.2f}%")
        logging.info(f"Maximum GC Content: {max_gc:.2f}%")
        logging.info(f"Mean GC Content: {mean_gc:.2f}%")
        logging.info(f"Median GC Content: {median_gc:.2f}%")
        logging.info(f"Genome with minimum GC: {min_genome_gc['Genome_File']}")
        logging.info(f"Genome with maximum GC: {max_genome_gc['Genome_File']}")

        # ------------------------------------------------------------------
        # Save genome length & GC table
        # ------------------------------------------------------------------
        genome_data_file = os.path.join(genomic_dir, "genome_metrics.csv")
        genome_df.to_csv(genome_data_file, index=False)
        logging.info(f"Genome metrics table saved: {genome_data_file}")

        # ------------------------------------------------------------------
        # Plot 1: Genome Length Distribution
        # ------------------------------------------------------------------
        sns.set(style="whitegrid")
        plt.figure(figsize=(10,6))
        
        sns.histplot(genome_df["Genome_Length"], bins=30, color="skyblue", edgecolor="black")
        plt.axvline(mean_len, color="red", linestyle="--", linewidth=2, label=f"Mean: {mean_len:.0f}")
        plt.axvline(median_len, color="green", linestyle="-.", linewidth=2, label=f"Median: {median_len:.0f}")
        
        plt.xlabel("Genome Length (bp)", fontsize=14, fontweight="bold")
        plt.ylabel("Frequency", fontsize=14, fontweight="bold")
        plt.title("Genome Length Distribution", fontsize=16, fontweight="bold")
        plt.legend()
        plt.tight_layout()
        
        plot_file_len = os.path.join(genomic_dir, "genome_length_distribution.png")
        plt.savefig(plot_file_len, dpi=300)
        plt.close()

        # ------------------------------------------------------------------
        # Plot 2: GC Content Distribution
        # ------------------------------------------------------------------
        plt.figure(figsize=(10,6))
        
        sns.histplot(genome_df["GC_Content_Percent"], bins=30, color="lightgreen", edgecolor="black")
        plt.axvline(mean_gc, color="red", linestyle="--", linewidth=2, label=f"Mean: {mean_gc:.2f}%")
        plt.axvline(median_gc, color="darkgreen", linestyle="-.", linewidth=2, label=f"Median: {median_gc:.2f}%")
        
        plt.xlabel("GC Content (%)", fontsize=14, fontweight="bold")
        plt.ylabel("Frequency", fontsize=14, fontweight="bold")
        plt.title("GC Content Distribution", fontsize=16, fontweight="bold")
        plt.legend()
        plt.tight_layout()
        
        plot_file_gc = os.path.join(genomic_dir, "gc_content_distribution.png")
        plt.savefig(plot_file_gc, dpi=300)
        plt.close()
        
        logging.info("Both distribution plots saved successfully.")

    except Exception as e:
        logging.error(f"ERROR during genomic analysis: {e}")
        print("Error occurred. Check genomic_analysis/genomic_analysis.log for details.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Genomic Metadata Processing Pipeline")

    parser.add_argument("--output", required=True, help="Output folder path")
    parser.add_argument("--taxid", required=True, help="TaxID JSON file")
    parser.add_argument("--metadata", required=True, help="Metadata CSV file")
    parser.add_argument("--genomes", required=True, help="Genome directory")

    args = parser.parse_args()

    print("Running Metadata Analysis...")
    analyze_assembly_metadata(args.output, args.taxid, args.metadata)

    print("Running Genome Analysis...")
    analyze_genomes_length_and_gc(args.output, args.genomes)

    print("All processing completed!")