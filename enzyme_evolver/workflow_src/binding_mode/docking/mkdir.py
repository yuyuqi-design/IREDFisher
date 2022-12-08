import os
def mk_dir(rec_file,lig_file,folder):
    #make the sub directory
    name_dir = str(rec_file).split('.')[0] + '_' + str(lig_file).split('.')[0]
    try:
        os.mkdir(folder+name_dir)
    #name the input configure file as rec_lig_config.txt; for example: 4agd_b49_config.txt
    except FileExistsError:
        print('directory exists.')
    return name_dir