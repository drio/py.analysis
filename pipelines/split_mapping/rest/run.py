#!flask/bin/python

from app import app
from flask import url_for
from tornado.wsgi import WSGIContainer
from tornado.httpserver import HTTPServer
from tornado.ioloop import IOLoop

app.setup_redis()

# Dev
app.run(debug=True)
# Production
"""
http_server = HTTPServer(WSGIContainer(app))
http_server.listen(5000)
IOLoop.instance().start()
"""

#url_for('static', filename='frontend.html')
