#!/bin/bash

# Check if the Docker container is running
if ! docker ps | grep -q "annovar_bioinformatics"; then
    echo "The 'annovar_bioinformatics' Docker container is not running."
    exit 1
fi

# Test Trimmomatic
docker exec -it annovar_bioinformatics java -jar trimmomatic-0.39.jar -version
if [ $? -eq 0 ]; then
    echo "Trimmomatic is working." 
else
    echo "Trimmomatic encountered an error."
fi

# Test FastQC
docker exec -it annovar_bioinformatics fastqc --version
if [ $? -eq 0 ]; then
    echo "FastQC is working."
else
    echo "FastQC encountered an error."
fi

# Test ANNOVAR
docker exec -it annovar_bioinformatics annotate_variation.pl -help
if [ $? -eq 1 ]; then
    echo "ANNOVAR is working."
else
    echo "ANNOVAR encountered an error."
fi
