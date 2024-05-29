#!/usr/bin/env python3

import logging
from Bio import SeqIO

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define the GenBank file path
genbank_file_path = "chr2.gb"

# Define the output file for the multi-sequence FASTA
output_fasta_file = "chr2_region_proteins.faa"

def parse_genbank(file_path, start, end, output_file):
    """Parse the GenBank file and extract protein sequences within the specified range."""
    logging.info(f"Parsing GenBank file {file_path} and extracting proteins between {start}-{end}...")
    records = SeqIO.parse(file_path, "genbank")
    proteins = []

    for record in records:
        logging.info(f"Processing record: {record.id}")
        for feature in record.features:
            if feature.type == "CDS":
                location = feature.location
                logging.info(f"Found CDS feature from {location.start} to {location.end}")
                if (start <= location.start <= end) or (start <= location.end <= end):
                    if "translation" in feature.qualifiers:
                        protein_seq = feature.qualifiers["translation"][0]
                        protein_id = feature.qualifiers.get("protein_id", ["unknown_protein_id"])[0]
                        protein_desc = feature.qualifiers.get("product", ["unknown_product"])[0]
                        proteins.append((protein_id, protein_desc, protein_seq, location))
                        logging.info(f"Extracted protein: {protein_id} - {protein_desc}")

    with open(output_file, "w") as fasta_file:
        for protein_id, protein_desc, protein_seq, location in proteins:
            location_str = f"chr2:{location.start}..{location.end}"
            fasta_file.write(f">{protein_id} {protein_desc} {location_str}\n")
            fasta_file.write(f"{protein_seq}\n")

    if not proteins:
        logging.warning("No proteins found in the specified range.")

def main():
    parse_genbank(genbank_file_path, 5000000, 10000000, output_fasta_file)
    logging.info(f"Protein sequences saved to {output_fasta_file}")

if __name__ == "__main__":
    logging.info("Starting protein extraction from GenBank file...")
    main()
    logging.info("Script completed.")
