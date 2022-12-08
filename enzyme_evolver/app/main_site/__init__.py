from flask import Blueprint

bp = Blueprint('main_site',
               __name__,
               template_folder='templates',
               static_folder='static',
               static_url_path='/main_site/static'
               )

from enzyme_evolver.app.main_site.routes import static_pages