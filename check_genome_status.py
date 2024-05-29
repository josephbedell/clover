#!/usr/bin/env python3

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.options import Options
from bs4 import BeautifulSoup
import time

def check_genome_status():
    url = 'https://portal.darwintreeoflife.org/data/root/details/Trifolium%20repens'
    
    # Set up the Selenium WebDriver
    chrome_options = Options()
    chrome_options.add_argument("--headless")  # Run headless Chrome
    chrome_service = Service('/home/jbedell/bin/chromedriver')  # Update this path to the location of chromedriver

    driver = webdriver.Chrome(service=chrome_service, options=chrome_options)
    driver.get(url)
    
    # Allow time for the page to load
    time.sleep(5)

    # Get the page source after JavaScript has rendered
    html = driver.page_source
    driver.quit()
    
    # Parse the HTML with BeautifulSoup
    soup = BeautifulSoup(html, 'html.parser')
    
    # Find status elements based on the updated structure
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

    # Extract additional details
    biosample_id = soup.find("div", {"class": "biosample-id"})
    organism = soup.find("div", {"class": "organism"})
    common_name = soup.find("div", {"class": "common-name"})
    
    if biosample_id:
        biosample_id_text = biosample_id.get_text(strip=True)
        print(f"Biosample ID: {biosample_id_text}")
    else:
        print("Biosample ID not found")

    if organism:
        organism_text = organism.get_text(strip=True)
        print(f"Organism: {organism_text}")
    else:
        print("Organism not found")

    if common_name:
        common_name_text = common_name.get_text(strip=True)
        print(f"Common Name: {common_name_text}")
    else:
        print("Common Name not found")

# Call the function
check_genome_status()
