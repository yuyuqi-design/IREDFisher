from flask import render_template, redirect, url_for, current_app
from enzyme_evolver.app.workflow import bp
from enzyme_evolver.app.workflow.forms import HomologySearchForm
from enzyme_evolver.app.workflow.functions import job_functions
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver.workflow_main.evolutionary_analysis_main import EvolAnalysis

@bp.route('/evolutionary_analysis_form', methods=['GET', 'POST'])
def evolutionary_analysis_form():
    form = HomologySearchForm()
    if form.validate_on_submit() == True:
        print('Creating new evolutionary analysis job..')
        protein_name = form.data['protein_name']
        protein_seq = form.data['protein_seq']
        # trim = form.data['trim']
        # ancestral = form.data['ancestral']
        folder_id = job_functions.create_new_job(protein_name, 'Homology Search')
        current_app.task_queue.enqueue(task_evol_analysis, folder_id, protein_name, protein_seq)

        return redirect(url_for("workflow.evolutionary_analysis", folder_id=folder_id))

    return render_template('stand_alone_tools/evolutionary_analysis_form.html', form=form)

@bp.route('/evolutionary_analysis/<folder_id>', methods=['GET'])
def evolutionary_analysis(folder_id):
    job = Job.objects(folder_id=folder_id)[0]

    return render_template('stand_alone_tools/evolutionary_analysis.html',
                           folder_id=folder_id,
                           protein_name=job.name,
                           status=job.status,
                           notes=job.notes,
                           files=job.files)

def task_evol_analysis(folder_id, protein_name, protein_seq):
    evol_analysis = EvolAnalysis(folder_id)
    evol_analysis.save_fasta(protein_name, protein_seq)
    evol_analysis.run()
    evol_analysis.job_finished()

    return {'status': 'success'}




