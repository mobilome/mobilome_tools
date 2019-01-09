#!/usr/bin/env python3

# -*- encoding: utf-8 -*-

"""
@Author  :   Naisu
 
@Contact :   3298990@qq.com
 
"""

import requests,urllib.request,time,os
from bs4 import BeautifulSoup


url = "https://mobilednajournal.biomedcentral.com/articles"
path = "C:/Users/ASUS/Desktop/Mobile_DNA/"
def download(url):
    web_data = requests.get(url)
    
    soup = BeautifulSoup(web_data.text,"lxml")
    DOI = soup.select("article a " )
    hrefs = []
    for i in DOI:
        href = i.get("href")
        if href.startswith("/track"):
            file_url = "https://mobilednajournal.biomedcentral.com" + href
            hrefs.append(file_url)
    
    for item in hrefs:
        name = item[-17:] + ".pdf"
        file_path = os.path.join(os.getcwd(),name)
        r = requests.get(item)
        with open (file_path,"wb") as code:
            code.write(r.content)
        time.sleep(5)

if __name__  == "__main__":
    download(url)
    print("Done")
