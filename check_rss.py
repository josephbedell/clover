#!/usr/bin/env python3

import feedparser

def check_rss_feed(url):
    feed = feedparser.parse(url)
    for entry in feed.entries:
        print(f"Title: {entry.title}")
        print(f"Link: {entry.link}")
        print(f"Published: {entry.published}")
        print("-----")

# Example feed URL (adjust as needed)
rss_url = 'https://pubmed.ncbi.nlm.nih.gov/rss/search/1RR2b2xR_jtTovcX9bGLYpZ5H3tBbkA5AuF7hByNqk-xu5qYzZ/?limit=100&utm_campaign=pubmed-2&fc=20220514032808'
check_rss_feed(rss_url)
