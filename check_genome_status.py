#!/usr/bin/env python3

import requests
from bs4 import BeautifulSoup

def check_genome_status():
    # Assuming this URL points directly to the search results for Trifolium repens
    url = 'https://portal.darwintreeoflife.org/data/root/details/Trifolium%20repens'
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"Failed to retrieve the webpage. Status code: {response.status_code}")
        return
    
    soup = BeautifulSoup(response.content, 'html.parser')
    
    # Print a portion of the parsed HTML to inspect the structure
    print("HTML Content:")
    print(soup.prettify()[:1000])  # Print the first 1000 characters for inspection
    
    # Update this part based on the actual structure
    status_elements = soup.find_all("div", class_="card bg-lite mb-3 filter-top")
    if not status_elements:
        print("No status elements found.")
    for element in status_elements:
        title_element = element.find("h3", class_="card-header filter-heading")
        details_element = element.find("ul", class_="list-group")
        
        if title_element:
            title = title_element.get_text(strip=True)
        else:
            title = "No title found"
        
        if details_element:
            details = details_element.get_text(strip=True)
        else:
            details = "No details found"
        
        print(f"Title: {title}")
        print(f"Details: {details}")
        print("-----")

# Call the function
check_genome_status()
