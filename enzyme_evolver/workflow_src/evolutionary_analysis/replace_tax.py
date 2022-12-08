import os
import sys
sys.path.append('/home/g02808yy/data/webserver/EnzymeEvolver')
from enzyme_evolver.workflow_src.evolutionary_analysis import batch_replace
import fnmatch
file = sys.argv[1]
ofolder = sys.argv[2]
[batch_replace.batch_replace(ofolder+treefile, ofolder+'replace.txt') for treefile in os.listdir(ofolder) if fnmatch.fnmatch(treefile, str(file).strip())]
