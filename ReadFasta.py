from Bio import SeqIO

def readReference(ref,format):
    dna=[]
    for seq_record in SeqIO.parse(ref, format):
        return(str(seq_record.seq),seq_record.id)
