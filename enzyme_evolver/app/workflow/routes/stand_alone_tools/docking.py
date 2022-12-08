from flask import render_template, redirect, url_for, current_app, request,send_file
from enzyme_evolver.app.workflow import bp
from enzyme_evolver.app.workflow.forms import DockingForm
from enzyme_evolver.app.workflow.functions import job_functions
from enzyme_evolver.mongo.workflow_models import Job
from werkzeug.utils import secure_filename
import os
from pathlib import Path
from enzyme_evolver.workflow_main.dockin_main import Docking

def save_file(l, folder_path, filename):
    with open(f"{folder_path}/{filename}", 'w') as output:
        output.writelines(l)
def name_ref_lig_cof(ref_ligand, cof):
    if ref_ligand == None:
        ref_filename = ''
    else:
        ref_filename = secure_filename(ref_ligand.filename)
    if cof == None:
        cof_filename = ''
    else:
        cof_filename = 'cof_ref.pdb'
    return ref_filename, cof_filename




@bp.route("/docking_form", methods=["GET", "POST"])
def docking_form():
    form = DockingForm()
    if form.validate_on_submit() == True:
        print('Creating new docking job..')
        job_name = form.data['job_name']
        folder_id = job_functions.create_new_job(job_name, 'Molecular Docking')
        folder_path = f"{Path.cwd()}/enzyme_evolver/database/{folder_id}"
        print('save files at ' + folder_path)
        proteins = form.proteins.data
        #print(proteins)
        ligands = form.ligands.data
        #options
        add_h = form.add_h.data
        gen_3D = form.gen_3D.data
        #print(ligands)
        ref_ligand = form.ref_ligand.data
        cof = form.cof.data
        Rg = form.Rg.data
        print(Rg)

        rec_file = []
        lig_file = []
        for pro in proteins:
            filename = secure_filename(pro.filename)
            pro.save(os.path.join(folder_path, filename))
            rec_file.append(filename + '\n')
            print('file saved')
        for lig in ligands:
            filename = secure_filename(lig.filename)
            lig_file.append(filename + '\n')
            lig.save(os.path.join(folder_path, filename))
            print('file saved')
        print(rec_file)
        print(lig_file)
        save_file(rec_file,folder_path,'receptors.txt')
        save_file(lig_file, folder_path, 'ligands.txt')

        ####save files
        print(ref_ligand)
        print(cof)
        ref_filename, cof_filename = name_ref_lig_cof(ref_ligand, cof)
        if ref_filename:
            ref_ligand.save(os.path.join(folder_path, ref_filename))
        if cof_filename:
            cof.save(os.path.join(folder_path, 'cof_ref.pdb'))

        current_app.task_queue.enqueue(task_docking, folder_id,ref_filename, cof_filename,add_h,gen_3D,Rg)

        return redirect(url_for("workflow.docking", folder_id=folder_id))

    return render_template('stand_alone_tools/docking_form.html', form=form)


@bp.route('/docking/<folder_id>', methods=['GET'])
def docking(folder_id):
    job = Job.objects(folder_id=folder_id)[0]

    return render_template('stand_alone_tools/docking.html',
                           folder_id=folder_id,
                           job_name=job.name,
                           status=job.status,
                           notes=job.notes,
                           files=job.files)


def task_docking(folder_id,ref_ligand,cof,add_h,gen_3D,Rg):
    docking = Docking(folder_id, ref_ligand=ref_ligand,cof=cof,add_h=add_h,gen_3D=gen_3D,Rg=Rg)
    docking.run()

    return {'status': 'success'}


#######example files
@bp.route('/enzyme_fisher/tutorials/molecular_docking/protein_example', methods=['GET'])
def protein_example():
    file=f'{Path.cwd()}/enzyme_evolver/database/tutorial/molecular_docking/reference/5ocm_pro.pdb'
    return send_file(file, attachment_filename='5ocm_pro.pdb', as_attachment=True)

@bp.route('/enzyme_fisher/tutorials/molecular_docking/ligand_example', methods=['GET'])
def ligand_example():
    file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/molecular_docking/reference/LIG.pdb'
    return send_file(file, attachment_filename='LIG.pdb', as_attachment=True)
@bp.route('/enzyme_fisher/tutorials/molecular_docking/reflig_example', methods=['GET'])
def reflig_example():
    file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/molecular_docking/reference/ref_lig.pdb'
    return send_file(file, attachment_filename='ref_lig.pdb', as_attachment=True)
@bp.route('/enzyme_fisher/tutorials/molecular_docking/refcof_example', methods=['GET'])
def refcof_example():
    file = f'{Path.cwd()}/enzyme_evolver/database/tutorial/molecular_docking/reference/ref_cof.pdb'
    return send_file(file, attachment_filename='ref_cof.pdb', as_attachment=True)