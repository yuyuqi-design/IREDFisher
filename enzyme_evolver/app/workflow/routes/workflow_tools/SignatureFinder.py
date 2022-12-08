from flask import render_template, redirect, url_for, current_app, request,flash, send_file
from enzyme_evolver.app.workflow import bp
from enzyme_evolver.app.workflow.forms import SignatureFinderForm
from enzyme_evolver.app.workflow.functions import job_functions
from enzyme_evolver.mongo.workflow_models import Job
from werkzeug.utils import secure_filename
import os
import subprocess as sp
from pathlib import Path
from enzyme_evolver.workflow_main.SignatureFinder import SigFinder

def save_file(l, folder_path, filename):
    with open(f"{folder_path}/{filename}", 'w') as output:
        output.writelines(l)
def N_rename(file, N_index,folder_path):
    extract_atoms = f"grep 'HETATM\|CONECT\|ATOM' {folder_path}/{file} > {folder_path}/atom_{file}; rm {folder_path}/{file} "
    sp.run(extract_atoms, shell=True)
    rename_N = f"sh {Path.cwd()}/enzyme_evolver/workflow_main/N_rename.sh {folder_path}/atom_{file} {N_index} {folder_path}/{file}"
    sp.run(rename_N,shell=True)

def panel_sequence_name(panel_sequences):
    if panel_sequences == None:
        panel_sequences_filename = False
    else:
        panel_sequences_filename = secure_filename(panel_sequences.filename)
    return panel_sequences_filename


@bp.route("/SignatureFinder_form", methods=["GET", "POST"])
def SignatureFinder_form():
    form = SignatureFinderForm()
    if form.validate_on_submit() == True:
        print('Finding signatures for binding your interested ligand')
        job_name = form.data['job_name']
        folder_id = job_functions.create_new_job(job_name, 'Signature Finder')
        folder_path = f"{Path.cwd()}/enzyme_evolver/database/{folder_id}"
        panel_sequences = form.data['panel_sequences']

        panel_sequences.save(os.path.join(folder_path, "allseq.fasta"))
        #referene complex file
        complex = form.data['complex']
        complex.save(os.path.join(folder_path, "complex.pdb"))
        Rg = form.data['Rg']
        RMSD = form.data['RMSD']
        current_app.task_queue.enqueue(task_SignatureFinder, folder_id=folder_id, Rg=Rg, RMSD=RMSD)
        return redirect(url_for("workflow.SignatureFinder", folder_id=folder_id,Rg=Rg, RMSD=RMSD))
    return render_template('workflow_tools/SignatureFinder_form.html', form=form)

@bp.route('/SignatureFinder/<folder_id>', methods=['GET'])
def SignatureFinder(folder_id):
    job = Job.objects(folder_id=folder_id)[0]

    return render_template('workflow_tools/SignatureFinder.html',
                           folder_id=folder_id,
                           job_name=job.name,
                           status=job.status,
                           notes=job.notes,
                           files=job.files)


def task_SignatureFinder(folder_id,Rg,RMSD):
    signatureFinder = SigFinder(folder_id,sequences='allseq.fasta', complex_struct='complex.pdb', Rg=Rg, rmsd=RMSD)
    signatureFinder.run_SignatureFinder()

    return {'status': 'success'}