#!/usr/bin/env python3

import logging
import subprocess
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define your email for NCBI
Entrez.email = "joseph.bedell@gmail.com"

# Define the GenBank accession numbers for CLV1 and BAM3
clv1_accession = "NM_001198022"
bam3_accession = "NM_127982"

# Define the SRA accession numbers for Trifolium repens samples
sra_accessions = ["ERS12548248"]

def download_sra_sequence(accession):
    """Download SRA sequence using fastq-dump and convert to FASTA format."""
    logging.info(f"Downloading SRA sequence for {accession}...")
    subprocess.run(["fastq-dump", "--fasta", accession])
    fasta_file = f"{accession}.fasta"
    return fasta_file

def create_blast_db(fasta_file):
    """Create a local BLAST database from the FASTA file."""
    logging.info(f"Creating BLAST database for {fasta_file}...")
    subprocess.run(["makeblastdb", "-in", fasta_file, "-dbtype", "nucl"])

def fetch_sequence(accession, output_file):
    """Fetch the DNA sequence for a given GenBank accession number and save it to a file."""
    logging.info(f"Fetching sequence for {accession}...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    logging.info(f"Saving sequence to {output_file}...")
    with open(output_file, "w") as file:
        SeqIO.write(record, file, "fasta")
    return record.seq

def local_blast_search(query_file, db, output_file):
    """Perform a local BLAST search using the downloaded sequences."""
    logging.info(f"Performing local BLAST search for {query_file} against {db}...")
    subprocess.run([
        "blastn", "-query", query_file, "-db", db, "-out", output_file,
        "-outfmt", "5"
    ])

def parse_blast_results(output_file):
    """Parse and print BLAST search results."""
    logging.info(f"Parsing BLAST results from {output_file}...")
    with open(output_file) as result_handle:
        blast_records = NCBIXML.read(result_handle)
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.05:
                    logging.info(f"Sequence: {alignment.title}\nLength: {alignment.length}\nE-value: {hsp.expect}\n")

def main():
    # Fetch sequences for CLV1 and BAM3
    clv1_sequence = fetch_sequence(clv1_accession, "clv1.NM_001198022.fna")
    bam3_sequence = fetch_sequence(bam3_accession, "bam3.NM_127982.fna")

    for sra_accession in sra_accessions:
        # Download SRA sequence and create a BLAST database
        fasta_file = download_sra_sequence(sra_accession)
        create_blast_db(fasta_file)

        # Perform local BLAST search for CLV1
        local_blast_search("clv1.NM_001198022.fna", fasta_file, "clv1_blast_results.xml")
        parse_blast_results("clv1_blast_results.xml")

        # Perform local BLAST search for BAM3
        local_blast_search("bam3.NM_127982.fna", fasta_file, "bam3_blast_results.xml")
        parse_blast_results("bam3_blast_results.xml")

if __name__ == "__main__":
    logging.info("Starting local BLAST search script for Trifolium repens SRA accessions...")
    main()
    logging.info("Script completed.")
