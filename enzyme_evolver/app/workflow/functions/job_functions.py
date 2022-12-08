import uuid
from enzyme_evolver import database_functions
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver.mongo.user_models import User
from flask_security import current_user

def create_new_job(name, type, initial_status='In Queue', no_user=False):
    if no_user is False:
        user = User.objects(id=current_user.id)[0]
    else:
        user = None
    new_id = str(uuid.uuid4())
    database_functions.make_new_folder(new_id)
    job = Job(folder_id=new_id, owner=user, name=name, type=type, status=initial_status)
    job.save()

    return new_id

def delete_job(folder_id):
    job = Job.objects(folder_id=folder_id)[0]
    job.delete()
    try:
        database_functions.delete_folder(folder_id)
    except:
        print('Could not delete folder')

def save_fasta(folder_id, protein_name, protein_seq):
    path = f"enzyme_evolver/database/{folder_id}/protein.fasta"

    with open(path, 'w') as file:
        file.write('>')
        file.write(protein_name)
        file.write('\n')
        file.write(protein_seq)
        file.write('\n')

    job = Job.objects(folder_id=folder_id)[0]
    job.files.append('protein.fasta')
    job.save()

    return path
