import subprocess as sp
import os


# def pdf2tif(ifile,ofile):
#     pages = convert_from_path(ifile)
#     for page in pages:
#         page.save(ofile, 'TIFF')


def sequence_process(ofolder, file='homologous_sequences90.fas', trim=False, ofseq_ali='homologous_sequences90_muscle.fas', ofhtml='homologous_sequences90_muscle.html'):
    #function: 1. sequence alignment; 2. mview 3. Gblock triming; 4. mview the trimed seqeunce alignment file
    #output file: eg: 1. cyanobacteria_muscle.fasta (sequence alignment by using Muscle); 2. cyanobacteria_mview.pdf (a pdf to view the sequence alignment); 
    # 3. cyanobacteria_muscle.fasta-gb (Trimmed seqeunce alignment file using Gblock ); 4. cyanobacteria_mview_gblock.pdf (a pdf to view Trimmed seqeunce alignment)
    #do multiple sequence alignment
    each = file.split('.')[0]
    muscle = 'muscle -in ' + ofolder + file + ' -out ' + ofolder+ofseq_ali
    sp.run(muscle, shell=True)
    #mveiw and pdf for initial sequence alginment
    #mview cyanobacteria_muscle.fasta -threshold 100  -consensus on -con_threshold 80  -con_ignore class  -con_coloring identity -html head -bold -css on -coloring consensus  -conservation on> mview.html
    mview = 'mview ' + ofolder+ofseq_ali + ' -threshold 100  -consensus on -con_threshold 80  -con_ignore class  -con_coloring identity -html head -bold -css on -coloring consensus  -conservation on  > ' + ofolder+ofhtml
    #xml2pdf = 'wkhtmltopdf  -s B0  ' + each + '_mview.html ' + each + '_mview.pdf'
    sp.run(mview,shell=True)
    #sp.run(xml2pdf,shell=True)
    #pdf2tif(each + '_mview.pdf',each + '_mview.tif')
    ##Trim the sequence alignment by using GBlocks
    if trim == True:
        #os.chdir(ofolder)
        gblock = 'Gblocks '  + ofolder+ ofseq_ali + ' -t=p -e=-gb -b-t=p  -b4=5 -b5=a -e=-gb'
        sp.run(gblock, shell=True)
    #color the conserved motif by mview and save as pdf by wkhtmltopdf
        mview_gblock = 'mview ' + ofolder+ofseq_ali + '-gb -threshold 100  -consensus on -con_threshold 80  -con_ignore class  -con_coloring identity -html head -bold -css on -coloring consensus  -conservation on  > ' + ofolder + ofhtml + '_gblock'
       # xml2pdf_gblock = 'wkhtmltopdf  -s B0  ' + each + '_mview_gblock.html ' + each + '_mview_gblock.pdf'
        sp.run(mview_gblock,shell=True)
        # os.chdir('../../../')
        #sp.run(xml2pdf_gblock,shell=True)
        #pdf2tif(each + '_mview_gblock.pdf',each + '_mview_gblock.tif')
    #sp.run('rm *pdf *htm *htm *html', shell=True)
    
