#!flask/bin/python

from app import app
from flask import url_for
from tornado.wsgi import WSGIContainer
from tornado.httpserver import HTTPServer
from tornado.ioloop import IOLoop

env_is = 'production'

app.setup_redis()

if env_is == 'dev':
    app.run(debug=True)
else:
    http_server = HTTPServer(WSGIContainer(app))
    http_server.listen(5000)
    IOLoop.instance().start()

#url_for('static', filename='frontend.html')
