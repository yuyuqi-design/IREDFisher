#!/usr/bin/env python
import sys
from rq import Connection, Worker
import os

from enzyme_evolver.app.app import create_app

production_mode = os.environ.get('PRODUCTION') or False

if __name__ == '__main__':
    app = create_app()

    app.app_context().push()

    with Connection(app.redis):
        qs = sys.argv[1:] or ['tasks']

        w = Worker(qs)
        w.work()
