from flask import render_template, redirect, url_for, current_app, flash, send_file
from enzyme_evolver.app.workflow import bp
from enzyme_evolver.app.workflow.forms import HomoModelForm
from werkzeug.utils import secure_filename
from enzyme_evolver.app.workflow.functions import job_functions
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver.workflow_main.HomoModelling_main import HomoModelling
from pathlib import Path
import os


@bp.route('/homo_model_form', methods=['GET', 'POST'])
def homo_model_form():
    form = HomoModelForm()
    if form.validate_on_submit() == True:
        job_name = form.data['job_name']
        folder_id = job_functions.create_new_job(job_name, 'Structure Modelling')
        folder_path = f"{Path.cwd()}/enzyme_evolver/database/{folder_id}"
        print('save files at ' + folder_path)

        fasta = form.data['fasta']

        # option 1: auto multiple template(3 tempates are used) modelling
        check_auto_multiple = form.data['check_auto_multiple']

        # option 2: use user's template structure
        check_template = form.data['check_template']
        template_single = form.template_single.data

        # option3: Homodimer modelling
        check_homodimer = form.data['check_homodimer']
        dimeT = form.data['dimeT']
        singleT = form.data['singleT']
        start_residue = form.data['start_residue']
        end_residue = form.data['end_residue']
####################################################################################
        options = check_auto_multiple + check_template + check_homodimer
        #no options are selected, run auto_single mode
        if options == 0:
            auto_single = True
            print('Creating new homology modelling job..')
            current_app.task_queue.enqueue(task_homo_model, folder_id, fasta, auto_single=True, auto_multiple=False,
                                           template_single=False, homodimer=False, start_residue='', end_residue='')
            return redirect(url_for("workflow.homo_model", folder_id=folder_id))
        # one options is selected
        if options == 1:
            auto_single = False
            if check_template == True:
                if template_single:
                    template_filename = secure_filename(template_single.filename)
                    template_single.save(os.path.join(folder_path, template_filename))
                    print('Creating new homology modelling job..')
                    current_app.task_queue.enqueue(task_homo_model, folder_id, fasta, auto_single=False, auto_multiple=False,
                                               template_single=template_filename, homodimer=False, start_residue='', end_residue='')
                    return redirect(url_for("workflow.homo_model", folder_id=folder_id))
                else:
                    flash('template structure has to be uploaded')
                    return render_template('stand_alone_tools/homo_model_form.html', form=form)
                    job_functions.delete_job(folder_id)
            if check_homodimer == True:
                if (dimeT and singleT and start_residue and end_residue):
                    dimeT.save(os.path.join(folder_path, 'dimerT.pdb'))
                    singleT.save(os.path.join(folder_path, 'singleT.pdb'))
                    print('Creating new homology modelling job..')
                    current_app.task_queue.enqueue(task_homo_model, folder_id, fasta, auto_single=False, auto_multiple=False,
                                               template_single=False, homodimer=True, start_residue=start_residue, end_residue=end_residue)
                    return redirect(url_for("workflow.homo_model", folder_id=folder_id))
                else:
                    flash('please check to fill in and upload files')
                    return render_template('stand_alone_tools/homo_model_form.html', form=form)
                    job_functions.delete_job(folder_id)
            else:
                print('Creating new homology modelling job..')
                current_app.task_queue.enqueue(task_homo_model, folder_id, fasta, auto_single=False, auto_multiple=True,
                                               template_single=False, homodimer=False, start_residue='', end_residue='')
                return redirect(url_for("workflow.homo_model", folder_id=folder_id))

        # if more than 1 options are selected
        if options >=2:
            job_functions.delete_job(folder_id)
            flash('only 1 option is allowed','error')
            return render_template('stand_alone_tools/homo_model_form.html', form=form)


    return render_template('stand_alone_tools/homo_model_form.html', form=form)

@bp.route('/homo_model/<folder_id>', methods=['GET'])
def homo_model(folder_id):
    job = Job.objects(folder_id=folder_id)[0]

    return render_template('stand_alone_tools/homo_model.html',
                           folder_id=folder_id,
                           job_name=job.name,
                           status=job.status,
                           notes=job.notes,
                           files=job.files)

def task_homo_model(folder_id, fasta_string, auto_single, auto_multiple, template_single, homodimer,start_residue, end_residue):
    hm = HomoModelling(folder_id, auto_single, auto_multiple, template_single, homodimer,start_residue, end_residue)
    hm.create_sequence_files(fasta_string)
    hm.run()

    return {'status': 'success'}

@bp.route('/enzyme_fisher/tutorials/structure_modelling/template_example', methods=['GET'])
def template_example():
    file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/structure_modelling/template/6l1hA.pdb'
    return send_file(file, attachment_filename='6l1hA.pdb', as_attachment=True)
@bp.route('/enzyme_fisher/tutorials/structure_modelling/dimerT_example', methods=['GET'])
def dimerT_example():
    file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/structure_modelling/homodimer/dimerT.pdb'
    return send_file(file, attachment_filename='dimerT.pdb', as_attachment=True)
@bp.route('/enzyme_fisher/tutorials/structure_modelling/singleT_example', methods=['GET'])
def singleT_example():
    file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/structure_modelling/homodimer/singleT.pdb'
    return send_file(file, attachment_filename='singleT.pdb', as_attachment=True)




