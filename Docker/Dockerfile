# Use a specific version of the base image for reproducibility
FROM continuumio/miniconda3:24.1.2-0

# Set working directory
WORKDIR /app

# Create environment and install tools in a single step to reduce layer size
# Note: Pinning versions for all major tools is best practice
RUN conda create -y -n genomic_env -c bioconda -c conda-forge \
    python=3.10 \
    abricate=1.0.1 \
    blast=2.15.0 \
    emboss=6.6.0 \
    ncbi-datasets-cli \
    && conda clean -afy

# Set the shell to use the conda environment for subsequent RUN commands
SHELL ["conda", "run", "-n", "genomic_env", "/bin/bash", "-c"]

# Install Python dependencies
# We do this separately from the project copy to leverage Docker's cache
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application
COPY . .

# Ensure scripts are executable
RUN chmod +x run_genomic_risk_pipeline.sh

# Use -v to ensure the environment is correctly activated in the entrypoint
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "genomic_env"]
CMD ["bash", "run_genomic_risk_pipeline.sh"]

#Use these commands to create image and run the docker wrapper script.
# docker build -t genomic-risk .
# docker run --rm -v $(pwd):/app genomic-risk
