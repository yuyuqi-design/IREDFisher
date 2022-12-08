from flask import render_template, redirect, url_for, current_app, request,flash, send_file
from enzyme_evolver.app.workflow import bp
from enzyme_evolver.app.workflow.forms import ScreenPanelDesForm
from enzyme_evolver.app.workflow.functions import job_functions
from enzyme_evolver.mongo.workflow_models import Job
from werkzeug.utils import secure_filename
import os
import subprocess as sp
from pathlib import Path
from enzyme_evolver.workflow_main.IREDFisher_main import IREDFisher

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


@bp.route("/screen_panel_des_form", methods=["GET", "POST"])
def screen_panel_des_form():
    form = ScreenPanelDesForm()
    if form.validate_on_submit() == True:
        print('Creating new screening panel design job..')
        job_name = form.data['job_name']
        folder_id = job_functions.create_new_job(job_name, 'Panel Design')
        folder_path = f"{Path.cwd()}/enzyme_evolver/database/{folder_id}"

        panel_sequences = form.data['panel_sequences']
        check_public = form.data['check_public']
        check_inhouse = form.data['check_inhouse']
        check_screened = form.data['check_screened']
        database = ''
        if check_screened:
            database = 'Screened'
        if (not check_inhouse) and check_public:
            database = 'Public'

        panel_sequences_filename = panel_sequence_name(panel_sequences)
        if panel_sequences_filename:
            panel_sequences.save(os.path.join(folder_path, "allseq.fasta"))
        #interesed ligands (multiple files)
        ligands = form.ligands.data
        N_index = form.data['N_index']

        lig_file = []
        for lig in ligands:
            filename = secure_filename(lig.filename)
            lig_file.append(filename + '\n')
            lig.save(os.path.join(folder_path, filename))
            print('file saved')
            #N_rename(filename, N_index, folder_path)
        save_file(lig_file, folder_path, 'ligands.txt')



        ####if the protein is homodimer
                # flash('please check your input files!')

        current_app.task_queue.enqueue(task_IRED_panel_des, folder_id=folder_id,database=database,N_index=N_index,panel_sequence=panel_sequences_filename)
        return redirect(url_for("workflow.screen_panel_des", folder_id=folder_id,database=database,N_index=N_index,panel_sequence=panel_sequences_filename))
    return render_template('workflow_tools/screen_panel_des_form.html', form=form)

@bp.route('/screen_panel_des/<folder_id>', methods=['GET'])
def screen_panel_des(folder_id):
    job = Job.objects(folder_id=folder_id)[0]

    return render_template('workflow_tools/screen_panel_des.html',
                           folder_id=folder_id,
                           job_name=job.name,
                           status=job.status,
                           notes=job.notes,
                           files=job.files)


def task_IRED_panel_des(folder_id,database='Screened',N_index='3',panel_sequence=False):
    panel = IREDFisher(folder_id,database, N_index, panel_sequence)
    panel.run_IREDFisher()

    return {'status': 'success'}

####example files
@bp.route('/enzyme_fisher/tutorials/one-click/panel_ligand_example', methods=['GET'])
def panel_ligand_example():
    file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/one-click/LIG.pdb'
    return send_file(file, attachment_filename='LIG.pdb', as_attachment=True)
@bp.route('/enzyme_fisher/tutorials/one-click/panel_complex_example', methods=['GET'])
def panel_complex_example():
    file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/one-click/5ocm.pdb'
    return send_file(file, attachment_filename='5ocm.pdb', as_attachment=True)
@bp.route('/enzyme_fisher/tutorials/structure_modelling/panel_dimerT_example', methods=['GET'])
def panel_dimerT_example():
    file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/one-click/dimerT.pdb'
    return send_file(file, attachment_filename='dimerT.pdb', as_attachment=True)
@bp.route('/enzyme_fisher/tutorials/one-click/panel_singleT_example', methods=['GET'])
def panel_singleT_example():
    file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/one-click/singleT.pdb'
    return send_file(file, attachment_filename='singleT.pdb', as_attachment=True)
@bp.route('/enzyme_fisher/tutorials/one-click/panel_sequence_example', methods=['GET'])
def panel_sequence_example():
    file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/one-click/panel_sequnece.fasta'
    return send_file(file, attachment_filename='panel_sequnece.fasta', as_attachment=True)
