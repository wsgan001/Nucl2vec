from Bio import SeqIO


def readReference(ref,format):
    dna=[]
    for seq_record in SeqIO.parse(ref, format):
        #print(seq_record.id)
        #print(repr(seq_record.seq))
        #print(len(seq_record))
        return(str(seq_record.seq),seq_record.id)
        #print(str(seq_record.seq))
        #print("\n")
        #print(seq_record.format("fasta"))
    #return dna

'''
for seq_record in SeqIO.parse("6484_snippet_1.fastq", "fastq"):
    #print(seq_record.id)
    #print(repr(seq_record.seq))
    #print(len(seq_record))
    print(str(seq_record.seq))
    print("\n")
    #print(seq_record.format("fasta"))
'''

#print(len(readReference("sequence.fasta","fasta")))


