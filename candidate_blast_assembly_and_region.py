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

# Define the GenBank accession numbers for the proteins to be used in the searches
protein_accessions = [
    "XP_019093641.1", "XP_010471716.1",  # CLV1 proteins
    "NP_193760.1",                       # BAM3 protein
    "XP_003556845.1",                    # GmNARK protein (Soybean)
    "XP_019447609.1",                    # LjHAR1 protein (Lotus japonicus)
    "XP_013464938.1"                     # MtSUNN protein (Medicago truncatula)
]

# Define the Clover genome URL and file name
clover_genome_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/408/175/GCA_030408175.1_UTM_Trep_v1.0/GCA_030408175.1_UTM_Trep_v1.0_genomic.fna.gz"
clover_genome_file = "GCA_030408175.1_genomic.fna.gz"
clover_genome_unzipped_file = "GCA_030408175.1_genomic.fna"

# Define the collection of proteins from the four-leaf locus on chromosome 2
four_leaf_proteins_file = "four_leaf_proteins.faa"  # Update this with the actual file path

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

def local_blast_search(query_file, db, output_file, blast_type="tblastn"):
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
                hits = []
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < 0.05:
                            hits.append((record.query, record.query_length, alignment.title, alignment.length,
                                         hsp.expect, (hsp.identities / hsp.align_length) * 100, hsp.align_length, 
                                         hsp.identities, db, blast_type))
                top_hits = sorted(hits, key=lambda x: x[4])[:5]
                for hit in top_hits:
                    out_handle.write("\t".join(map(str, hit)) + "\n")

def main():
    # Download and create BLAST database for the clover genome
    download_clover_genome()

    # Fetch sequences for the proteins
    protein_files = []
    for accession in protein_accessions:
        protein_files.append(fetch_protein_sequence(accession, f"{accession}.faa"))

    # Perform local tblastn search for the proteins against the clover genome
    for protein_file in protein_files:
        local_blast_search(protein_file, clover_genome_unzipped_file, f"{protein_file}_tblastn.xml", blast_type="tblastn")
        parse_blast_results(f"{protein_file}_tblastn.xml", f"{protein_file}_tblastn_results.tsv", clover_genome_unzipped_file, "tblastn")

    # Perform local blastp search for the proteins against the four-leaf locus protein collection
    for protein_file in protein_files:
        local_blast_search(protein_file, four_leaf_proteins_file, f"{protein_file}_blastp.xml", blast_type="blastp")
        parse_blast_results(f"{protein_file}_blastp.xml", f"{protein_file}_blastp_results.tsv", four_leaf_proteins_file, "blastp")

if __name__ == "__main__":
    logging.info("Starting local BLASTP and TBLASTN search script for homologs against clover genome and four-leaf locus proteins...")
    main()
    logging.info("Script completed.")
