print('Please use pymol.pymolhttpd instead of web.pymolhttpd')

from pymol.pymolhttpd import *

if __name__ == 'pymol':
    server = PymolHttpd()
    server.start()
