import requests,urllib.request,time
from bs4 import BeautifulSoup


urls = ["https://weheartit.com/inspirations/taylorswift?page={}".format(str(i)) for i in range(1,20)]
path = "C://Users/ASUS/Desktop/k2p/img/"
def download(url):
    web_data = requests.get(url)
    time.sleep(2)
    soup = BeautifulSoup(web_data.text,"lxml")
    imgs = soup.select("#main-container > div.grid-responsive > div.col.span-content > div > div > div > div > div > a > img")
    down_list = []
    for img in imgs:
        down_link = img.get("src")
        down_list.append(down_link)


    for item in down_list:
        urllib.request.urlretrieve(item,path+item[-10:]+".jpg")

for a in urls:
    download(a)
    print("Done")
