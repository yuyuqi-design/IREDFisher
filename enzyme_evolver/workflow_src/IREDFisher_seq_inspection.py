from enzyme_evolver.workflow_src.homology_modelling import auto_modeller
from pathlib import Path
def seq_inspection(seqcode):
    masterpath = Path.cwd()
    db = f"{masterpath}/enzyme_evolver/database/pdb_95.pir"
    auto_modeller.template_search(seqcode, db)

