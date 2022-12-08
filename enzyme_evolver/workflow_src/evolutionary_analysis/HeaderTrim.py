import subprocess as sp
import pandas as pd
def HeaderTrim(seqfile):
    with open(seqfile) as seq:
        with open(f'{seqfile}_header','w') as out:
            content = seq.readlines()
            for line in content:
                if line.startswith('>'):
                    line = line[:13]
                    out.write(line +'\n')

                else:
                    out.write(line)
    sp.run(f'mv {seqfile}_header {seqfile}',shell=True)

def GetHeader(infile,outfile):
    with open(infile) as fin:
        with open(outfile,'w') as fout:
            content = fin.readlines()
            for line in content:
                if line.startswith('>'):
                    print(line)
                    line = line[1:]
                    fout.write(line)

def Make_Header(infile,outfile):
    with open(infile) as fin:
        with open(outfile,'w') as fout:
            content = fin.readlines()
            number = len(content)
            for i in range(number):
                fout.write('seq'+str(i)+'\n')


def ReplaceHeader(infile='file.txt', A='A.txt', B='B.txt', s_r_list='replace.txt'):
    #replace A with B
    join = f'paste -d ";" {A} {B} > {s_r_list}'
    sp.run(join,shell=True)
    #s_rlist_file content: apple;banana (replace apple with banana)
    f = open(infile)
    #read the file as f_content (string)
    f_content = f.read()
    #read the replace.txt as dataframe with delimiter as ; and name the first column as stxt (searching string), the second is rtxt (replace string)
    df = pd.read_csv(s_r_list, sep=';',names=['stxt','rtxt'])
    df.head()
    #get the length of stxt then add them as a new column into dataframe
    df['length'] = df['stxt'].str.strip().str.len()
    #Rank the dataframe according to the length of stxt in decreasing order. This can avoid the substring was replaced if it is in ascending order.
    df.sort_values(by='length',ascending=False,inplace=True)
    #turn the stxt and rtxt back to list
    replace_list = list(df['stxt'] + ';' + df['rtxt'])
    #Do the replacement
    for line in replace_list:
        stxt = line.split(';')[0]
        # print(stxt)
        rtxt = line.split(';')[1]
        # print(rtxt)
        f_content = f_content.replace(stxt.strip() , rtxt.strip())
    #write to a new file
    with open(infile+'replaced', 'w') as output:
        output.write(f_content)
    f.close()
    sp.run(f'mv {infile}replaced {infile}',shell=True)


if __name__ == '__main__':
    ReplaceHeader(infile='input_aln.fasta', A='uniq_header.txt', B='original_header.txt', s_r_list='replace.txt')


