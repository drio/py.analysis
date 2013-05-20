#!/usr/bin/env python

import json, requests

data = {'sites': [ {'Chrm':'1', 'Start':100}, {'Chrm':'1', 'Start':101}, {'Chrm':'2', 'Start':101} ] }
data_json = json.dumps(data)
r = requests.post("http://localhost:8080/", data=data_json)

for i in json.loads(r.text):
  print i
