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
        print(f"Title: {title}")
        print(f"Link: {link}")
        print("-----")

search_pubmed("Trifolium repens genome")
