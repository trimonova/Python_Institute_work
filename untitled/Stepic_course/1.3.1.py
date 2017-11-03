import requests

URL = "https://stepic.org/favicon.ico"
r = requests.get(URL)
print(r.status_code)

print(r.headers['server'])

print(r.encoding)

print(r.text)






