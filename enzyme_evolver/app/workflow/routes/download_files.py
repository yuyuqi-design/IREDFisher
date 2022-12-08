from enzyme_evolver.app.workflow import bp
from enzyme_evolver.app.workflow.forms import NewWorkFlowForm
from flask import render_template, redirect, url_for, request, jsonify, send_file
from enzyme_evolver.mongo.workflow_models import Job
import uuid
from enzyme_evolver.app.workflow.functions import job_functions
from enzyme_evolver.mongo.user_models import User
from flask_security import current_user
from enzyme_evolver import database_functions
import mongoengine as db
from pathlib import Path
working_dir = Path(__file__).parents[4]

@bp.route('/download_file/<folder_id>/<file_name>', methods=['GET'])
def download_file(folder_id, file_name):
    job = Job.objects(folder_id=folder_id)[0]
    user = User.objects(id=current_user.id)[0]

    if job.owner != user:
        return redirect('/')
    elif file_name not in job.files:
        return redirect('/')

    path_to_file = f"{Path(__file__).parents[4]}/enzyme_evolver/database/{folder_id}/{file_name}"

    return send_file(path_to_file)


if __name__ == "__main__":
    working_dir = Path(__file__).parents[4]
    print(working_dir)
