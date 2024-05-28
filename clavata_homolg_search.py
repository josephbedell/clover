#!/usr/bin/env python3
from Bio import SeqIO, Entrez, SearchIO
from Bio.Blast import NCBIWWW, NCBIXML

# Define the genomic interval for the search
chromosome = '2'  # Example chromosome number
start_position = 5000000  # Start of the interval (5 Mb)
end_position = 10000000   # End of the interval (10 Mb)

# Path to the white clover genome file (FASTA format)
genome_file = 'Trifolium_repens_genome.fasta'

# Define the sequences of CLAVATA1 (CLV1) and BARELY ANY MERISTEM 3 (BAM3) for the search
clv1_sequence = '...CLV1 sequence...'  # Replace with actual sequence
bam3_sequence = '...BAM3 sequence...'  # Replace with actual sequence

def load_genome(genome_file):
    """Load the genome sequence from a FASTA file."""
    genome = {}
    with open(genome_file, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            genome[record.id] = record.seq
    return genome

def extract_interval_sequence(genome, chromosome, start, end):
    """Extract the sequence of a specified interval from the genome."""
    if chromosome in genome:
        return genome[chromosome][start:end]
    else:
        raise ValueError("Chromosome not found in the genome.")

def search_homologs(sequence, search_sequence):
    """Search for homologs of a given sequence within the target sequence using BLAST."""
    result_handle = NCBIWWW.qblast('blastn', 'nt', search_sequence)
    blast_records = NCBIXML.read(result_handle)
    return blast_records

def main():
    # Load the genome
    genome = load_genome(genome_file)

    # Extract the interval sequence
    interval_sequence = extract_interval_sequence(genome, chromosome, start_position, end_position)

    # Search for CLV1 homologs
    clv1_blast_results = search_homologs(interval_sequence, clv1_sequence)
    print("CLV1 Homologs:")
    for alignment in clv1_blast_results.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.05:
                print(f"Sequence: {alignment.title}\nLength: {alignment.length}\nE-value: {hsp.expect}\n")

    # Search for BAM3 homologs
    bam3_blast_results = search_homologs(interval_sequence, bam3_sequence)
    print("BAM3 Homologs:")
    for alignment in bam3_blast_results.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.05:
                print(f"Sequence: {alignment.title}\nLength: {alignment.length}\nE-value: {hsp.expect}\n")

if __name__ == "__main__":
    main()
