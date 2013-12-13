#!/usr/bin/env python

import sys
import requests
import json
from requests.auth import HTTPBasicAuth
import os


def check_service(url):
    try:
        requests.get(url, auth=authentication()).json()
    except:
        sys.stderr.write("Problems connecting to service (%s). Bailing out.\n" % url)
        exit(1)


def update_resource(url=None, name=None, step=None):
    headers = {'Content-type': 'application/json', 'Accept': 'text/plain'}

    if url and name and step:
        r = requests.get(url, auth=authentication()).json()
        _id = None
        for s in r["samples"]:  # Pick the first one !!
            if s["name"] == name:
                _id = s["id"]
        if not _id:
            data = {'name': name}
            r = requests.post(url, data=json.dumps(data), headers=headers, auth=authentication())
            print r
            _id = r.json()["sample"]["id"]

        data = {'step': step}
        r = requests.put(url + "/" + str(_id),
            data=json.dumps(data), headers=headers,
            auth=authentication())

def set_auth():
    try:
        user = os.environ['SAPI_USER']
        pwd = os.environ['SAPI_PWD']
    except:
        sys.stderr.write('Env vars for service credentials not setup (SAPI_USER SAPI_PWD)\n')
        exit(1)

    def a():
        return HTTPBasicAuth(user, pwd)
    return a

# Main
authentication = set_auth()

if len(sys.argv) == 4:
    url, name, step = sys.argv[1:4]
    check_service(url)
    update_resource(url, name, step)
else:
    if len(sys.argv) == 3 and sys.argv[1] == "list":
        url = sys.argv[2]
        check_service(url)
        data = requests.get(url, auth=authentication()).json()
        print json.dumps(data, sort_keys=True, indent=4)
    else:
        sys.stderr.write("Usage:\n")
        sys.stderr.write("$ tool <resource_url> <sample_name> <step>\n")
        sys.stderr.write("$ tool list <resource_url>\n")
