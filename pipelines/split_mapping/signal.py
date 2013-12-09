#!/usr/bin/env python

import sys
import requests
import json


def update_resource(url=None, name=None, step=None):
    headers = {'Content-type': 'application/json', 'Accept': 'text/plain'}

    if url and name and step:
        r = requests.get(url).json()
        _id = None
        for s in r["samples"]:  # Pick the first one !!
            if s["name"] == name:
                _id = s["id"]
        if not _id:
            data = {'name': name}
            r = requests.post(url, data=json.dumps(data), headers=headers)
            print r
            _id = r.json()["sample"]["id"]

        data = {'step': step}
        r = requests.put(url + "/" + str(_id), data=json.dumps(data), headers=headers)

if len(sys.argv) == 4:
    url, name, step = sys.argv[1:4]
    update_resource(url, name, step)
else:
    if len(sys.argv) == 3 and sys.argv[1] == "list":
        url = sys.argv[2]
        data = requests.get(url).json()
        print json.dumps(data, sort_keys=True, indent=4)
    else:
        sys.stderr.write("Usage:\n")
        sys.stderr.write("$ tool <resource_url> <sample_name> <step>\n")
        sys.stderr.write("$ tool list <resource_url>\n")
