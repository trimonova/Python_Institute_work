import requests
from bs4 import BeautifulSoup, SoupStrainer

url='https://www.crummy.com/software/BeautifulSoup/bs4/doc/'
r = requests.get(url)
soup = BeautifulSoup(r.content, 'html.parser', parse_only=SoupStrainer('a'))
all_links = [link['href'] for link in soup if link.has_attr('href')]
selected_links = []
for i in all_links:
    if i.startswith('http'):
        answer = requests.get(i).status_code
        if answer == 404:
            selected_links.append(i)
    else:
        answer = requests.get(url + i).status_code
        if answer == 404:
            selected_links.append(i)
print(selected_links)