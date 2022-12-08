import os
from pathlib import Path
import shutil
from shutil import copyfile

def make_new_folder(id):
    print('Creating new folder at:')
    path_to_folder = f"{Path(__file__).parents[0]}/database/{id}"
    print(path_to_folder)
    os.mkdir(path_to_folder)

def delete_folder(folder_id):
    path = f"{Path(__file__).parents[0]}/database/{folder_id}"
    print(f'Deleting folder - {path}')
    shutil.rmtree(path)

def test_blastxml(folder_id):
    copy_to = f"{Path(__file__).parents[0]}/database/{folder_id}/blast.xml"
    copy_from = f"{Path(__file__).parents[0]}/database/blast.xml"
    copyfile(copy_from, copy_to)