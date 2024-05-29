#!/usr/bin/env python3
import requests
from bs4 import BeautifulSoup

def check_genome_status():
    url = 'https://portal.darwintreeoflife.org/data/root/details/Trifolium%20repens'
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"Failed to retrieve the webpage. Status code: {response.status_code}")
        return
    
    soup = BeautifulSoup(response.content, 'html.parser')
    
    # Print a portion of the parsed HTML to inspect the structure
    print("HTML Content:")
    print(soup.prettify()[:1000])  # Print the first 1000 characters for inspection
    
    status_elements = soup.find_all("div", class_="card bg-lite mb-3 filter-top")
    if not status_elements:
        print("No status elements found.")
    for element in status_elements:
        title = element.find("h3", class_="card-header filter-heading")
        details = element.find("ul", class_="list-group")
        
        if title:
            title = title.get_text(strip=True)
        else:
            title = "No title found"
        
        if details:
            details = details.get_text(strip=True)
        else:
            details = "No details found"
        
        print(f"Title: {title}")
        print(f"Details: {details}")
        print("-----")

# Uncomment the call to test
check_genome_status()

