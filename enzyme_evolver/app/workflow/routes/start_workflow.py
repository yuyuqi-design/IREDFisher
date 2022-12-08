from enzyme_evolver.app.workflow import bp
from enzyme_evolver.app.workflow.forms import NewWorkFlowForm
from flask import render_template, redirect, url_for
from enzyme_evolver.mongo.workflow_models import Job
import uuid
from enzyme_evolver.app.workflow.functions import job_functions

@bp.route('/new_workflow', methods=['GET', 'POST'])
def new_workflow():
    form = NewWorkFlowForm()

    if form.validate_on_submit() == True:
        folder_id = job_functions.create_new_job(form.data['name'], "Workflow", initial_status='Idle')
        return redirect(url_for('workflow.workflow_launch', folder_id=folder_id))

    return render_template('new_workflow.html', form=form)

@bp.route('/workflow_launch/<folder_id>', methods=['GET'])
def workflow_launch(folder_id):
    workflow = Job.objects(folder_id=folder_id)[0]

    return render_template('workflow.html', name=workflow.name, folder_id=folder_id)
