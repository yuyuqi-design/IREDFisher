import requests
import urllib.parse
import urllib.request
from bs4 import BeautifulSoup
def pdb2uni(pdb_id):
    url = 'https://www.uniprot.org/uploadlists/'

    params = {
        'from': 'PDB_ID',
        'to': 'ACC',
        'format': 'tab',
        'query': pdb_id
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    codes = response.decode('utf-8')
    uni_id=codes.strip().split('\t')[-1]
    return uni_id


def find_catalytic_site(uni_id):

    source = requests.get('https://www.uniprot.org/uniprot/'+uni_id).text
    soup = BeautifulSoup(source,'lxml')
    function = soup.find(id='function').find("span",class_="context-help tooltipped-click html tipId-5")
    site_residues = function.find("a", class_="position tooltipped")
    if site_residues:
        site_residues=site_residues.text
        print(site_residues)
    else:
        print('no annotation')


if __name__ == '__main__':
    # pdb_id = '2AR8'
    # uni_id=pdb2uni(pdb_id)
    uni_id = 'P95480'
    find_catalytic_site(uni_id)

