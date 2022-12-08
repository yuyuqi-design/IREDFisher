from flask import render_template, redirect, url_for, current_app, flash, send_file
from enzyme_evolver.app.workflow import bp
from enzyme_evolver.app.workflow.forms import MDForm
from werkzeug.utils import secure_filename
from enzyme_evolver.app.workflow.functions import job_functions
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver.workflow_main.md_main import MDpreparation
from pathlib import Path
import subprocess as sp
import os

@bp.route('/md_form', methods=['GET', 'POST'])
def md_form():
    form = MDForm()
    if form.validate_on_submit() == True:
        print('Creating new docking job..')
        job_name = form.data['job_name']
        folder_id = job_functions.create_new_job(job_name, 'MD Preparation')
        folder_path = f"{Path.cwd()}/enzyme_evolver/database/{folder_id}"
        print('save files at ' + folder_path)
        job_detail = form.data['job_detail']
        proteins = form.proteins.data
        # print(proteins)
        ligands = form.ligands.data
        # print(ligands)
        for file in job_detail:
            # filename = secure_filename(aln.filename)
            file.save(os.path.join(folder_path, 'file.txt'))
            sp.run(f'dos2unix {folder_path}/file.txt', shell=True)
        for pro in proteins:
            filename = secure_filename(pro.filename)
            print(filename)
            pro.save(os.path.join(folder_path, filename))
            print('file saved')
        for lig in ligands:
            if lig:
                filename = secure_filename(lig.filename)
                lig.save(os.path.join(folder_path, filename))
                print('file saved')

        current_app.task_queue.enqueue(task_docking, folder_id)
        ###go to the next route function
        return redirect(url_for("workflow.mdpreparation", folder_id=folder_id))

    return render_template('stand_alone_tools/md_form.html', form=form)

@bp.route('/md/<folder_id>', methods=['GET'])
def mdpreparation(folder_id):
    job = Job.objects(folder_id=folder_id)[0]

    return render_template('stand_alone_tools/md.html',
                           folder_id=folder_id,
                           job_name=job.name,
                           status=job.status,
                           notes=job.notes,
                           files=job.files)


def task_docking(folder_id):
    md_prep = MDpreparation(folder_id)
    md_prep.run()

    return {'status': 'success'}