from flask import Blueprint

bp = Blueprint('workflow',
               __name__,
               template_folder='templates',
               static_folder='static',
               static_url_path='/main_site/static'
               )

from enzyme_evolver.app.workflow.routes import start_workflow, jobs, download_files
from enzyme_evolver.app.workflow.routes.stand_alone_tools import evol_analysis, homo_model, docking, phylo_tree, mdpreparation
from enzyme_evolver.app.workflow.routes.workflow_tools import screen_panel_des
from enzyme_evolver.app.workflow.routes.workflow_tools import SignatureFinder
