from tornado.wsgi import WSGIContainer
from tornado.httpserver import HTTPServer
from tornado.ioloop import IOLoop
from variant_mapper import app

http_server = HTTPServer(WSGIContainer(app))
http_server.listen(7000)
IOLoop.instance().start()