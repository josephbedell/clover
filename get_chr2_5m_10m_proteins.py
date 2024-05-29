#!/usr/bin/env python3

import logging
import requests
from Bio import SeqIO
from io import StringIO

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define the URL for the GenBank file
genbank_url = "https://www.ncbi.nlm.nih.gov/nuccore/CP125838.1?report=genbank&log$=seqview&format=text"

# Define the output file for the multi-sequence FASTA
output_fasta_file = "chr2_region_proteins.faa"

def fetch_genbank_file(url):
    """Fetch the GenBank file from the provided URL."""
    logging.info(f"Fetching GenBank file from {url}...")
    response = requests.get(url)
    response.raise_for_status()
    return response.text

def parse_genbank(genbank_data, start, end, output_file):
    """Parse the GenBank file and extract protein sequences within the specified range."""
    logging.info(f"Parsing GenBank file and extracting proteins between {start}-{end}...")
    records = SeqIO.parse(StringIO(genbank_data), "genbank")
    proteins = []

    for record in records:
        logging.info(f"Processing record: {record.id}")
        for feature in record.features:
            if feature.type == "CDS":
                location = feature.location
                logging.info(f"Found CDS feature from {location.start.position} to {location.end.position}")
                if (start <= location.start.position <= end) or (start <= location.end.position <= end):
                    if "translation" in feature.qualifiers:
                        protein_seq = feature.qualifiers["translation"][0]
                        protein_id = feature.qualifiers.get("protein_id", ["unknown_protein_id"])[0]
                        protein_desc = feature.qualifiers.get("product", ["unknown_product"])[0]
                        proteins.append((protein_id, protein_desc, protein_seq))
                        logging.info(f"Extracted protein: {protein_id} - {protein_desc}")

    with open(output_file, "w") as fasta_file:
        for protein_id, protein_desc, protein_seq in proteins:
            fasta_file.write(f">{protein_id} {protein_desc}\n")
            fasta_file.write(f"{protein_seq}\n")

    if not proteins:
        logging.warning("No proteins found in the specified range.")

def main():
    genbank_data = fetch_genbank_file(genbank_url)
    parse_genbank(genbank_data, 5000000, 10000000, output_fasta_file)
    logging.info(f"Protein sequences saved to {output_fasta_file}")

if __name__ == "__main__":
    logging.info("Starting protein extraction from GenBank file...")
    main()
    logging.info("Script completed.")
