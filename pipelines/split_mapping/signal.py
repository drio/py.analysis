#!/usr/bin/env python

import sys
import requests
import json


def update_resource(url=None, name=None, step=None, sos=None):
    headers = {'Content-type': 'application/json', 'Accept': 'text/plain'}

    if url and name and step and sos in ["start", "stop"]:
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

        data = {'start_stop': sos, 'step': step}
        r = requests.put(url + "/" + str(_id), data=json.dumps(data), headers=headers)

if len(sys.argv) == 5:
    url, name, step, sos = sys.argv[1:5]
    update_resource(url, name, step, sos)
else:
    if len(sys.argv) == 3 and sys.argv[1] == "list":
        url = sys.argv[2]
        data = requests.get(url).json()
        print json.dumps(data, sort_keys=True, indent=4)
    else:
        sys.stderr.write("Usage:\n")
        sys.stderr.write("$ tool <resource_url> <sample_name> <step> <start or stop>\n")
        sys.stderr.write("$ tool list <resource_url>\n")
        exit(1)
