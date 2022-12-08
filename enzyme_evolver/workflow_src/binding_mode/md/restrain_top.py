###############################################################################################
###############################################################################################
#include the restraint file into topol or itp file
import sys
def restraint_to_itp(file):
    with open(file) as itp:
        content = itp.readlines()
        i = content.index('; Include Position restraint file\n')
        print(i)
        content.remove('#ifdef POSRES\n')
        content.remove('#include "posre.itp"\n')
        content.remove('#endif\n')

    content.insert(i+1,'#ifdef POSRES_PROTEIN\n#include "posre_PROTEIN.itp"\n#endif\n#ifdef POSRES_Heavy\n#include "posre_Heavy.itp"\n#endif\n#ifdef POSRES_MC\n#include "posre_MC.itp"\n#endif\n#ifdef POSRES_CA\n#include "posre_CA.itp"\n#endif\n')
    with open('topol_restraint.top','w') as f2:
        f2.writelines(content)

if __name__ == '__main__':
    file=sys.argv[1]
    restraint_to_itp(file)
