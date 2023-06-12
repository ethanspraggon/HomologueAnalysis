#!/usr/bin/python

from operator import ge
import httplib2 as http
import json

try:
  from urlparse import urlparse
except ImportError:
  from urllib.parse import urlparse

import requests

import multiprocessing as mp



gene_dict = {}

def getHugoData(url = "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/json/hgnc_complete_set.json"):

  response = requests.get(url)

  if response.status_code == 200:
      
      data = response.json()
      # Process the data as needed

      with open("data.json", "w") as outfile:
        json.dump(data, outfile)

  else:
      print("Error downloading data. Status code:", response.status_code)

if __name__ == "__main__":
  getHugoData()



"""
headers = {
  'Accept': 'application/json',
}

uri = 'https://rest.genenames.org'
path = '/search/symbol/*'

target = urlparse(uri+path)
method = 'GET'
body = ''

h = http.Http()

response, content = h.request(
  target.geturl(),
  method,
  body,
  headers)


print(type(content))
if response['status'] == '200':
  # assume that content is a json reply
  # parse content with the json module 
  data = json.loads(content)
  print( 'Symbol:' + data['response']['docs'][0]['symbol'])
else:
  print( 'Error detected: ' + response['status'])

"""

