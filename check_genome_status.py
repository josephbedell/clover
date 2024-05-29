#!/usr/bin/env python3

import requests
from bs4 import BeautifulSoup

def check_genome_status():
    url = 'https://portal.darwintreeoflife.org/data/root/details/Trifolium%20repens'
    response = requests.get(url)
    soup = BeautifulSoup(response.content, 'html.parser')
    
    status_elements = soup.find_all("div", class_="card bg-lite mb-3 filter-top")
    for element in status_elements:
        title = element.find("h3", class_="card-header filter-heading").get_text(strip=True)
        details = element.find("ul", class_="list-group").get_text(strip=True)
        print(f"Title: {title}")
        print(f"Details: {details}")
        print("-----")

check_genome_status()
