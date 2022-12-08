from Bio import SeqIO, AlignIO
import subprocess as sp
import sys

def trim_TerGap(file1,file2='temp.fasta', format1='pir',format2='fasta'):

    aln = SeqIO.parse(file1,format1)
    SeqIO.convert(file1,format1,file2,format2)
    aln = AlignIO.read(file2, format2)
    for col in range(aln.get_alignment_length()):
        if not "-" in aln[:,col]:
            position = col
            break
    for col2 in reversed(range(aln.get_alignment_length())):
        if not "-" in aln[:,col2]:
            position2 = col2+1
            break
    SeqIO.write(aln[:,position:position2],file1,format1)
    sp.run(f"sed -i 's/>XX/>P1/g' {file1}", shell=True)


if __name__ == '__main__':
    file1 = sys.argv[1]
    trim_TerGap(file1)



