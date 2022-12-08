from enzyme_evolver.workflow_src.homology_modelling.auto_modeller import fasta2pir, template_search

def identity2target(file='allseq.fasta'):
    fasta2pir(file)
    with open(file) as seqfile:
        target_name = seqfile.readline()[1:]
    template_search(target_name,file.split('.')[0] +'.pir')

