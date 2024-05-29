#!/usr/bin/env python3

import requests
from bs4 import BeautifulSoup

def search_pubmed(query):
    url = f"https://pubmed.ncbi.nlm.nih.gov/?term={query}"
    response = requests.get(url)
    soup = BeautifulSoup(response.content, 'html.parser')
    
    articles = soup.find_all("div", class_="docsum-content")
    for article in articles:
        title = article.find("a", class_="docsum-title").get_text(strip=True)
        link = "https://pubmed.ncbi.nlm.nih.gov" + article.find("a", class_="docsum-title")["href"]
        
        # Extract publication year
        publication_info = article.find("span", class_="docsum-journal-citation full-journal-citation")
        publication_year = publication_info.get_text(strip=True).split(';')[0] if publication_info else "No year found"
        
        # Extract author list
        author_list = article.find("span", class_="docsum-authors full-authors").get_text(strip=True) if article.find("span", class_="docsum-authors full-authors") else "No authors found"
        
        print(f"Title: {title}")
        print(f"Link: {link}")
        print(f"Journal: {publication_year}")
        print(f"Authors: {author_list}")
        print("-----")

search_pubmed("Trifolium repens AND (genome OR four-leaf OR multi)")
