from flask import render_template, redirect, url_for, current_app, flash, send_file
from enzyme_evolver.app.workflow import bp
from enzyme_evolver.app.workflow.forms import PhyloTreeForm
from werkzeug.utils import secure_filename
from enzyme_evolver.app.workflow.functions import job_functions
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver.workflow_main.phylogeny_main import Phylogeny
from pathlib import Path
import os
import subprocess as sp


@bp.route('/phylo_tree_form', methods=['GET', 'POST'])
def phylo_tree_form():
    form = PhyloTreeForm()
    if form.validate_on_submit() == True:
        print('validated')
        job_name = form.data['job_name']
        folder_id = job_functions.create_new_job(job_name, 'Phylogenetic Tree')
        # folder_id = job_functions.create_new_job(job_name, 'Phylogenetic Tree')
        folder_path = f"{Path.cwd()}/enzyme_evolver/database/{folder_id}"
        selection = form.data['selection']
        print(selection)
        files = form.data['file']
        # alns.save(os.path.join(folder_path, 'input_aln.fasta'))
        print(files)
        print('save files at ' + folder_path)
        if selection=='sequence':
            for file in files:
                # filename = secure_filename(aln.filename)
                file.save(os.path.join(folder_path, 'input.fasta'))
                sp.run(f'dos2unix {folder_path}/input.fasta', shell=True)
        if selection=='alignment':
            for file in files:
                # filename = secure_filename(aln.filename)
                file.save(os.path.join(folder_path, 'input_aln.fasta'))
                sp.run(f'dos2unix {folder_path}/input_aln.fasta', shell=True)

        ancestral = form.data['ancestral']
        current_app.task_queue.enqueue(task_phylo_tree, folder_id=folder_id, aln=selection, ancestral=ancestral)


        return redirect(url_for("workflow.phylo_tree", folder_id=folder_id))

    return render_template('stand_alone_tools/phylo_tree_form.html', form=form)
####################################################################################


@bp.route('/phylo_tree/<folder_id>', methods=['GET'])
def phylo_tree(folder_id):
    job = Job.objects(folder_id=folder_id)[0]

    return render_template('stand_alone_tools/phylo_tree.html',
                           folder_id=folder_id,
                           job_name=job.name,
                           status=job.status,
                           notes=job.notes,
                           files=job.files)

def task_phylo_tree(folder_id, aln, ancestral):
    tree = Phylogeny(folder_id=folder_id, alignment=aln, bootstrap_ancestral=ancestral)
    tree.run()

    return {'status': 'success'}

# @bp.route('/enzyme_fisher/tutorials/structure_modelling/template_example', methods=['GET'])
# def template_example():
#     file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/structure_modelling/template/6l1hA.pdb'
#     return send_file(file, attachment_filename='6l1hA.pdb', as_attachment=True)
# @bp.route('/enzyme_fisher/tutorials/structure_modelling/dimerT_example', methods=['GET'])
# def dimerT_example():
#     file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/structure_modelling/homodimer/dimerT.pdb'
#     return send_file(file, attachment_filename='dimerT.pdb', as_attachment=True)
# @bp.route('/enzyme_fisher/tutorials/structure_modelling/singleT_example', methods=['GET'])
# def singleT_example():
#     file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/structure_modelling/homodimer/singleT.pdb'
#     return send_file(file, attachment_filename='singleT.pdb', as_attachment=True)




