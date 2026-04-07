import os
import subprocess
import zipfile
import time
import argparse

# -------------------------
# ARGUMENTS 
# -------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Accession file")
parser.add_argument("--out", required=True, help="Output directory")
parser.add_argument("--log", required=True, help="Log file path")

args = parser.parse_args()

ACCESSION_FILE = args.input
BASE_DIR = args.out
LOG_FILE = args.log

BATCH_SIZE = 5
MAX_RETRIES = 3

GENOME_DIR = os.path.join(BASE_DIR, "genomes")
GFF_DIR = os.path.join(BASE_DIR, "gff")
TMP_DIR = os.path.join(BASE_DIR, "tmp")

os.makedirs(GENOME_DIR, exist_ok=True)
os.makedirs(GFF_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)

start_time = time.time()

# -------------------------
# Read accession list
# -------------------------
with open(ACCESSION_FILE) as f:
    accessions = [x.strip() for x in f if x.strip()]

downloaded = set()

if os.path.exists(LOG_FILE):
    with open(LOG_FILE) as f:
        downloaded = set(x.strip() for x in f)

# Create list of accessions that still need to be downloaded by comparing the full list with 
# the downloaded set. This allows resuming if the script is interrupted, and also prevents 
# re-downloading genomes that were already
remaining = [a for a in accessions if a not in downloaded] 

print("Total accessions:", len(accessions))
print("Already downloaded:", len(downloaded))
print("Remaining:", len(remaining))

# -------------------------
# Calculate number of batches
# -------------------------
total_batches = (len(remaining) + BATCH_SIZE - 1) // BATCH_SIZE

completed = len(downloaded)

# -------------------------
# Process batches
# -------------------------
for i in range(0, len(remaining), BATCH_SIZE):

    batch = remaining[i:i+BATCH_SIZE]
    batch_number = i // BATCH_SIZE + 1

    print(f"\nBatch {batch_number}/{total_batches}")

    zip_file = os.path.join(TMP_DIR, "download.zip")

    cmd = [
        "datasets",
        "download",
        "genome",
        "accession",
        *batch,
        "--include",
        "genome,gff3",
        "--no-progressbar",
        "--filename",
        zip_file
    ]

    success = False

    for attempt in range(MAX_RETRIES):

        try:
            subprocess.run(cmd, check=True)

            with zipfile.ZipFile(zip_file) as z:
                z.extractall(TMP_DIR)

            success = True
            break

        except Exception:

            print("Download failed. Retrying...")

            if os.path.exists(zip_file):
                os.remove(zip_file)

    if not success:
        print("Batch failed after retries")
        continue

    data_path = os.path.join(TMP_DIR, "ncbi_dataset/data")

    for acc in batch:

        acc_dir = os.path.join(data_path, acc)

        if not os.path.exists(acc_dir):
            continue

        found_fna = False

        for file in os.listdir(acc_dir):

            src = os.path.join(acc_dir, file)

            if file.endswith(".fna"):
                os.rename(src, os.path.join(GENOME_DIR, f"{acc}.fna"))
                found_fna = True

            elif file.endswith(".gff"):
                os.rename(src, os.path.join(GFF_DIR, f"{acc}.gff"))

        if found_fna:

            with open(LOG_FILE, "a") as log:
                log.write(acc + "\n")

            completed += 1
            print(f"Downloaded {completed}/{len(accessions)}")

    # cleanup
    if os.path.exists(zip_file):
        os.remove(zip_file)

    dataset_folder = os.path.join(TMP_DIR, "ncbi_dataset")

    if os.path.exists(dataset_folder):
        subprocess.run(["rm", "-rf", dataset_folder])

# -------------------------
# Finish
# -------------------------
elapsed = time.time() - start_time

minutes = int(elapsed // 60)
seconds = int(elapsed % 60)

print("\nDownload finished")
print(f"Total time: {minutes} minutes {seconds} seconds")

