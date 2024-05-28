#!/usr/bin/env python3

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.service import Service as ChromeService
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.chrome.options import Options
import sys
from bs4 import BeautifulSoup

def search_dtol_portal(driver, query):
    url = "https://portal.darwintreeoflife.org/data"
    driver.get(url)
    
    search_box = driver.find_element(By.ID, "mat-input-0")
    search_box.clear()
    search_box.send_keys(query)
    search_box.send_keys(Keys.RETURN)
    
    driver.implicitly_wait(10)  # Wait for the results to load
    
    return driver.page_source

def check_genome_status(driver, query, output_file):
    html_content = search_dtol_portal(driver, query)
    
    with open(output_file, "w") as file:
        file.write(html_content)
    
    if html_content:
        soup = BeautifulSoup(html_content, 'html.parser')
        
        # Adjust the selector based on the actual HTML structure observed
        results = soup.find_all("div", class_="card bg-lite mb-3 filter-top")  # Example adjustment

        if not results:
            print(f"No results found in the page for {query}.")
            return

        for result in results:
            title = result.find("h3", class_="card-header filter-heading")
            details = result.find("ul", class_="list-group")
            
            if title and details:
                species = title.get_text(strip=True)
                detail_items = details.find_all("li", class_="list-group-item")
                sequencing_status = ", ".join(item.get_text(strip=True) for item in detail_items)
                
                print(f"Species: {species}")
                print(f"Sequencing Status: {sequencing_status}")
                print("---------")
            else:
                print(f"Result found but some elements are missing for {query}.")
                print(result.prettify())
    else:
        print(f"No results found or there was an error with the search for {query}.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: ./clover_genome_check.py <genome_name>")
        sys.exit(1)
    
    input_genome = sys.argv[1]
    
    chrome_options = Options()
    chrome_options.add_argument("--headless")
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-dev-shm-usage")
    
    driver = webdriver.Chrome(service=ChromeService(ChromeDriverManager().install()), options=chrome_options)
    
    print("Checking Linaria vulgaris as a control:")
    check_genome_status(driver, "Linaria vulgaris", "linaria_vulgaris_results.html")
    
    print(f"\nChecking {input_genome}:")
    check_genome_status(driver, input_genome, f"{input_genome.replace(' ', '_').lower()}_results.html")
    
    driver.quit()
