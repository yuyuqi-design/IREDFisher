from enzyme_evolver.app.workflow import bp
from enzyme_evolver.app.workflow.forms import NewWorkFlowForm
from flask import render_template, redirect, url_for, request, jsonify
from enzyme_evolver.mongo.workflow_models import Job
import uuid
from enzyme_evolver.app.workflow.functions import job_functions
from enzyme_evolver.mongo.user_models import User
from flask_security import current_user
from enzyme_evolver import database_functions
import mongoengine as db

@bp.route('/go_to_job/<folder_id>', methods=['GET'])
def go_to_job(folder_id):
    job = Job.objects(folder_id=folder_id)[0]

    if job.type == "Workflow":
        return redirect(url_for("workflow.workflow_launch",
                                folder_id=folder_id))

    if job.type == 'Homology Search':
        return redirect(url_for("workflow.evolutionary_analysis",
                                folder_id=folder_id))
    if job.type == 'Phylogenetic Tree':
        return redirect(url_for("workflow.phylo_tree",
                                folder_id=folder_id))

    if job.type == 'Structure Modelling':
        return redirect(url_for("workflow.homo_model",
                                folder_id=folder_id))
    if job.type == 'Molecular Docking':
        return redirect(url_for("workflow.docking",
                        folder_id=folder_id))
    if job.type == 'MD Preparation':
        #output page
        return redirect(url_for("workflow.mdpreparation",
                        folder_id=folder_id))
    if job.type == 'Panel Design':
        return redirect(url_for("workflow.screen_panel_des", folder_id=folder_id))
    if job.type == 'SignatureFinder':
        return redirect(url_for("workflow.SignatureFinder", folder_id=folder_id))

@bp.route('/Job status', methods=['GET'])
def Job_status():
    user = User.objects(id=current_user.id)[0]

    jobs_data = list(Job.objects(db.Q(owner=user)&db.Q(type__ne="Workflow")).order_by('type').exclude('owner').as_pymongo())

    for i, job in enumerate(jobs_data):
        jobs_data[i]['_id'] = str(jobs_data[i]['_id'])

    return render_template("my_jobs.html",
                           jobs_data=jobs_data,
                           title='Job status')

@bp.route('/ComputHub Job status', methods=['GET'])
def Job_status_computhub():
    user = User.objects(id=current_user.id)[0]

    jobs_data = list(Job.objects(db.Q(owner=user)&db.Q(type__ne="Workflow")).order_by('type').exclude('owner').as_pymongo())

    for i, job in enumerate(jobs_data):
        jobs_data[i]['_id'] = str(jobs_data[i]['_id'])

    return render_template("my_jobs_computhub.html",
                           jobs_data=jobs_data,
                           title='Job status')

@bp.route('/my_workflows', methods=['GET'])
def my_workflows():
    user = User.objects(id=current_user.id)[0]

    jobs_data = list(Job.objects(db.Q(owner=user) & db.Q(type="Workflow")).order_by('type').exclude('owner').as_pymongo())

    for i, job in enumerate(jobs_data):
        jobs_data[i]['_id'] = str(jobs_data[i]['_id'])

    return render_template("my_jobs.html",
                           jobs_data=jobs_data,
                           title='My workflows')

@bp.route('/_delete_other_job', methods=['GET', 'POST'])
def delete_other_job():
    folder_id = request.form['folder_id']
    job = Job.objects(folder_id=folder_id)[0]
    user = User.objects(id=current_user.id)[0]

    if job.owner == user:
        job_functions.delete_job(folder_id)

    return jsonify(result={'status': 'success'})
