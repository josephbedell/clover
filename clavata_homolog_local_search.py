#!/usr/bin/env python3

import logging
import os
import subprocess
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define your email for NCBI
Entrez.email = "joseph.bedell@gmail.com"

# Define the GenBank accession numbers for CLV1 mRNA and protein
clv1_mrna_accessions = ["XM_019238096.1", "XM_010473414.2"]
clv1_protein_accessions = ["XP_019093641.1", "XP_010471716.1"]

# Define the GenBank accession numbers for BAM3 mRNA and protein
bam3_mrna_accession = "NM_118146.5"
bam3_protein_accession = "NP_193760.1"

# Define the Clover genome URL and file name
clover_genome_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/408/175/GCA_030408175.1_UTM_Trep_v1.0/GCA_030408175.1_UTM_Trep_v1.0_genomic.fna.gz"
clover_genome_file = "GCA_030408175.1_genomic.fna.gz"
clover_genome_unzipped_file = "GCA_030408175.1_genomic.fna"

def download_clover_genome():
    """Download the clover genome and create a BLAST database if not already present."""
    if not os.path.exists(clover_genome_file) and not os.path.exists(clover_genome_unzipped_file):
        logging.info("Downloading clover genome...")
        subprocess.run(["wget", "-O", clover_genome_file, clover_genome_url], check=True)
        logging.info("Unzipping clover genome...")
        subprocess.run(["gunzip", clover_genome_file], check=True)
    else:
        logging.info("Clover genome file already exists. Skipping download and unzip steps.")
    
    if not os.path.exists(clover_genome_unzipped_file + ".nin"):  # Check for BLAST DB index files
        logging.info("Creating BLAST database from clover genome...")
        subprocess.run(["makeblastdb", "-in", clover_genome_unzipped_file, "-dbtype", "nucl"], check=True)
    else:
        logging.info("BLAST database already exists. Skipping BLAST database creation step.")

def fetch_sequence(accession, output_file):
    """Fetch the DNA sequence for a given GenBank accession number and save it to a file."""
    logging.info(f"Fetching sequence for {accession}...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    logging.info(f"Saving sequence to {output_file}...")
    with open(output_file, "w") as file:
        SeqIO.write(record, file, "fasta")
    return output_file

def fetch_protein_sequence(accession, output_file):
    """Fetch the protein sequence for a given GenBank accession number and save it to a file."""
    logging.info(f"Fetching protein sequence for {accession}...")
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    logging.info(f"Saving protein sequence to {output_file}...")
    with open(output_file, "w") as file:
        SeqIO.write(record, file, "fasta")
    return output_file

def local_blast_search(query_file, db, output_file, blast_type="blastn"):
    """Perform a local BLAST search using the downloaded sequences."""
    logging.info(f"Performing local {blast_type} search for {query_file} against {db}...")
    subprocess.run([
        blast_type, "-query", query_file, "-db", db, "-out", output_file,
        "-outfmt", "5"
    ], check=True)

def parse_blast_results(output_file, results_file, db, blast_type):
    """Parse BLAST search results and save to a tab-delimited file."""
    logging.info(f"Parsing BLAST results from {output_file}...")
    with open(output_file) as result_handle:
        blast_records = NCBIXML.parse(result_handle)  # Use parse() instead of read() for iterable
        with open(results_file, "w") as out_handle:
            out_handle.write("query\tqlength\tsubject\tslength\tevalue\tpercent_id\taln_length\tnum_identities\tdatabase\tblast_type\n")
            for record in blast_records:
                top_hits = sorted(record.alignments, key=lambda aln: min(hsp.expect for hsp in aln.hsps))[:5]
                for alignment in top_hits:
                    for hsp in alignment.hsps:
                        if hsp.expect < 0.05:
                            query = record.query
                            qlength = record.query_length
                            subject = alignment.title
                            slength = alignment.length
                            evalue = hsp.expect
                            percent_id = (hsp.identities / hsp.align_length) * 100
                            aln_length = hsp.align_length
                            num_identities = hsp.identities
                            out_handle.write(f"{query}\t{qlength}\t{subject}\t{slength}\t{evalue}\t{percent_id:.2f}\t{aln_length}\t{num_identities}\t{db}\t{blast_type}\n")

def main():
    # Download and create BLAST database for the clover genome
    download_clover_genome()

    # Fetch sequences for CLV1 and BAM3 mRNA
    clv1_mrna_files = []
    for accession in clv1_mrna_accessions:
        clv1_mrna_files.append(fetch_sequence(accession, f"{accession}.fna"))

    bam3_mrna_file = fetch_sequence(bam3_mrna_accession, f"{bam3_mrna_accession}.fna")

    # Fetch sequences for CLV1 and BAM3 protein
    clv1_protein_files = []
    for accession in clv1_protein_accessions:
        clv1_protein_files.append(fetch_protein_sequence(accession, f"{accession}.faa"))

    bam3_protein_file = fetch_protein_sequence(bam3_protein_accession, f"{bam3_protein_accession}.faa")

    # Perform local BLAST search for CLV1 mRNA sequences
    for seq_file in clv1_mrna_files:
        local_blast_search(seq_file, clover_genome_unzipped_file, f"{seq_file}_blast.xml", blast_type="blastn")
        parse_blast_results(f"{seq_file}_blast.xml", f"{seq_file}_blast_results.tsv", clover_genome_unzipped_file, "blastn")

    # Perform local BLAST search for BAM3 mRNA sequence
    local_blast_search(bam3_mrna_file, clover_genome_unzipped_file, f"{bam3_mrna_file}_blast.xml", blast_type="blastn")
    parse_blast_results(f"{bam3_mrna_file}_blast.xml", f"{bam3_mrna_file}_blast_results.tsv", clover_genome_unzipped_file, "blastn")

    # Perform local tblastn search for CLV1 protein sequences
    for protein_file in clv1_protein_files:
        local_blast_search(protein_file, clover_genome_unzipped_file, f"{protein_file}_tblastn.xml", blast_type="tblastn")
        parse_blast_results(f"{protein_file}_tblastn.xml", f"{protein_file}_tblastn_results.tsv", clover_genome_unzipped_file, "tblastn")

    # Perform local tblastn search for BAM3 protein sequence
    local_blast_search(bam3_protein_file, clover_genome_unzipped_file, f"{bam3_protein_file}_tblastn.xml", blast_type="tblastn")
    parse_blast_results(f"{bam3_protein_file}_tblastn.xml", f"{bam3_protein_file}_tblastn_results.tsv", clover_genome_unzipped_file, "tblastn")

if __name__ == "__main__":
    logging.info("Starting local BLAST and TBLASTN search script for CLV1 and BAM3 against clover genome...")
    main()
    logging.info("Script completed.")
