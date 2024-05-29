#!/usr/bin/env python3

from Bio import Entrez, SeqIO

# Set your email
Entrez.email = "joseph.bedell@gmail.com"

# Search for proteins in the specified region
search_handle = Entrez.esearch(db="protein", term="GCA_030408175.1[Assembly] AND 02_Occ:5,000,000-10,000,000[Location]", retmax=500)
search_results = Entrez.read(search_handle)
search_handle.close()

# Fetch protein data
protein_ids = search_results['IdList']
fetch_handle = Entrez.efetch(db="protein", id=protein_ids, rettype="gb", retmode="text")
protein_records = list(SeqIO.parse(fetch_handle, "genbank"))
fetch_handle.close()

# Save to file
with open("proteins.gb", "w") as output_handle:
    SeqIO.write(protein_records, output_handle, "genbank")

print(f"Downloaded {len(protein_records)} protein records.")
