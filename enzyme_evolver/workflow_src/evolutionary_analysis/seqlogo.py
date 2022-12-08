import subprocess as sp
from PIL import Image

def seqlogo(file='homologous_sequences90_muscle.fas', ofolder=None):
    output_file_name = file.split(".")[0] + '_sequence_logo.eps'
    command_logo = 'weblogo   -c chemistry  <'+ ofolder + file + '> ' + ofolder + output_file_name
    #print(command_logo)
    sp.run(command_logo, shell=True)
    im = Image.open(ofolder + output_file_name)
    im.save(ofolder + file.split(".")[0] + '_sequence_logo.png',lossless=True)
    
    #command_eps2png = 'convert -colorspace sRGB -density 200 ' + ofolder+ 'sequence_logo.eps ' + '-background white -flatten ' + ofolder+'sequence_logo.png'
    #sp.run(command_eps2png,shell=True)
    sp.run('rm '+ofolder+ output_file_name,shell=True)

if __name__ == '__main__':
    seqlogo(file='homologous_sequences90_muscle.fas', ofolder='/home/g02808yy/data/webserver/EnzymeEvolver/enzyme_evolver/database/test/')
