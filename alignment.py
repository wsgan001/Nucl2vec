from sklearn.neighbors import NearestNeighbors
from gensim.models import Word2Vec
import time
import sys
import numpy as np
import pysam as  ps
import ReadFasta as fa


# Dictionary for Base pair to index and index to base pair conversions
genomeDict = {'A': 0, 'C': 1, 'T': 2, 'G': 3}
genomeRevDict = {0: 'A', 1: 'C', 2: 'T', 3: 'G'}


encodingLength = 4
matchingPieces = 4

# This function returns the value corresponding to the match or mismatch from 'pt'
def Diagonal(n1,n2,pt):
    if(n1 == n2):
        return pt['MATCH']
    else:
        return pt['MISMATCH']


# This function gets the optional elements of the aligment matrix and returns the elements for the pointers matrix.
def Pointers(di,ho,ve):

    # Based on python default maximum(return the first element).
    pointer = max(di,ho,ve)

    if(di == pointer):
        return 'D' # Diagonal
    elif(ho == pointer):
        return 'H' # Horizontal
    else:
         return 'V' # Vertical


# This function creates the aligment and pointers matrices
def NW(s1,s2,match = 1,mismatch = -1, gap = -2):
    penalty = {'MATCH': match, 'MISMATCH': mismatch, 'GAP': gap} #A dictionary for all the penalty valuse.
    m = len(s1) + 1 #The dimension of the matrix columns.
    n = len(s2) + 1 #The dimension of the matrix rows.
    al_mat = np.zeros((m,n),dtype = int) #Initializes the alignment matrix with zeros.
    p_mat = np.zeros((m,n),dtype = str) #Initializes the alignment matrix with zeros.
    #Scans all the first rows element in the matrix and fill it with "gap penalty"
    for i in range(m):
        al_mat[i][0] = penalty['GAP'] * i
        p_mat[i][0] = 'V'
    #Scans all the first columns element in the matrix and fill it with "gap penalty"
    for j in range(n):
        al_mat[0][j] = penalty['GAP'] * j
        p_mat [0][j] = 'H'
    #Fill the matrix with the correct values.

    p_mat [0][0] = 0 #Return the first element of the pointer matrix back to 0.
    for i in range(1,m):
        for j in range(1,n):
            di = al_mat[i-1][j-1] + Diagonal(s2[j-1],s1[i-1],penalty) #The value for match/mismatch -  diagonal.
            ho = al_mat[i][j-1] + penalty['GAP'] #The value for gap - horizontal.(from the left cell)
            ve = al_mat[i-1][j] + penalty['GAP'] #The value for gap - vertical.(from the upper cell)
            al_mat[i][j] = max(di,ho,ve) #Fill the matrix with the maximal value.(based on the python default maximum)
            p_mat[i][j] = Pointers(di,ho,ve)

    # Create a score for the best alignment
    indi = m-1
    indj = n-1
    while(indi!=0 and indj!=0):
        if(p_mat[indi][indj]=='D'):
            indi = indi-1
            indj = indj-1
        elif(p_mat[indi][indj]=='H'):
            indj = indj-1
        elif(p_mat[indi][indj]=='V'):
            indi = indi-1

    ie = m-1
    je = n-1
    while(p_mat[ie][je]=='H' and je>0):
        je = je-1

    # Normalise the score by the read length and return
    return al_mat[ie][je]/len(s1), indj


# Encode the string into a distributed representation of encoded kmers
def encodeKmer(inp_str, model, encodingLength):

    ans_arr = []
    i = 0
    n = len(inp_str)

    # Iterate through the string
    # taking 'encodingLength' number of base pairs at a time
    # and converting them into encoding using the input model
    # Finally append it to an array for the final distributed representation
    while (i < n - encodingLength + 1):
        enc = model[inp_str[i:i + encodingLength]]
        ans_arr.extend(enc)
        i += encodingLength

    return ans_arr


if(len(sys.argv)!=4):
    print("ERROR : Required format -> python alignment.py <reads_fileName> <reference_fileName> <model_filename>")
    print("reads_fileName -> .bam file expected")
    print("reference_fileName -> .fasta file expected")
    print("model_filename -> The same file saved during pre_processing")
    exit()

bam_file = sys.argv[1]
print("Reading file: ", bam_file)

pmerLength = matchingPieces*encodingLength

reads_found = 0
base2num = dict(zip("ACTG", (0, 1, 2, 3)))

bamfile = ps.AlignmentFile(bam_file, "r")

kmer_arr = []
kmer_ind = []

# Reading the file contains the reads to be aligned.
# Any kmer of length less than '32' will be skipped (no aligned).
for SeqRead in bamfile:
    query_sequence = SeqRead.query_sequence
    leng = len(query_sequence)
    if (leng >= 32):
        reference_start = SeqRead.reference_start
        flag = 0

        # Check that all the base pairs are one of the four -> A, C, T or G
        # That is, it is not 'N'
        for x in range(0, len(query_sequence)):
            if (query_sequence[x] is not "N") and (query_sequence[x] is not "n"):
                flag = 0
            else:
                flag = 1
                break

        if flag == 0:
            kmer_arr.append(query_sequence)
            kmer_ind.append([reference_start])

            reads_found += 1
            # Keep track of the progress
            if reads_found % 100000 == 0:
                print("Reads Found : ", reads_found)

print("Total Reads Found : ", reads_found)


# Reading reference DNA
# The final reference DNA is stored in ref_dna as a string
dna = fa.readReference(sys.argv[2], "fasta")
ref_dna = dna[0]

model = Word2Vec.load(sys.argv[3])

startTime = time.clock()

print("Creating KNN Tree")
# Encode the reference DNA into the KNN tree
X = []
for i in range(0, len(ref_dna) - pmerLength):
    X.append(encodeKmer(ref_dna[i:i + pmerLength]))

# Create the Nearestneighbors tree
nbrs = NearestNeighbors(n_neighbors=1, algorithm='auto').fit(X)

endTime = time.clock()

print("Preprocessing Done")
print("Preprocessing time : ", endTime - startTime)


# For each of the four possible pieces, extract the string of length pmerLength frm the relevant position
Y = []
for kmr, ind in zip(kmer_arr, kmer_ind):
    Y.append(encodeKmer(kmr[:pmerLength]))

Y2 = []
kmer_ind2 = []
for kmr, inde in zip(kmer_arr, kmer_ind):
    Y2.append(encodeKmer(kmr[-pmerLength:]))
    kmer_ind2.append([inde[0] + len(kmr) - pmerLength])

Y3 = []
kmer_ind3 = []
for kmr, inde in zip(kmer_arr, kmer_ind):
    Y3.append(encodeKmer(kmr[int(len(kmr)/2)-pmerLength:int(len(kmr)/2)]))
    kmer_ind3.append([inde[0] + int(len(kmr)/2) - pmerLength])

Y4 = []
kmer_ind4 = []
for kmr, inde in zip(kmer_arr, kmer_ind):
    Y4.append(encodeKmer(kmr[int(len(kmr)/2):int(len(kmr)/2)+pmerLength]))
    kmer_ind4.append([inde[0] + int(len(kmr)/2)])

startTime = time.clock()

dist, ind = nbrs.kneighbors(Y)
dist2, ind2 = nbrs.kneighbors(Y2)
dist3, ind3 = nbrs.kneighbors(Y3)
dist4, ind4 = nbrs.kneighbors(Y4)


pred = []

# Check for the minimum distance out of the four
# Now check if the minimum distance is less than the threshold otherwise pred the index to be -1
for i1, i2, i3, i4, d1, d2, d3, d4, ele in zip(ind, ind2, ind3, ind4, dist, dist2, dist3, dist4, kmer_arr):
    if (d1[0] <= d2[0] and d1[0] <= d3[0] and d1[0] <= d4[0] and d1[0]<=0.0002):
        pred.append(i1[0])
    elif (d1[0] >= d2[0] and d3[0] >= d2[0] and d2[0] <= d4[0] and d2[0]<=0.0002):
        pred.append(i2[0]-len(ele)+pmerLength)
    elif (d1[0] >= d3[0] and d2[0] >= d3[0] and d3[0] <= d4[0] and d3[0]<=0.0002):
        pred.append(i3[0]-int(len(ele)/2)+pmerLength)
    elif (d1[0] >= d4[0] and d2[0] >= d4[0] and d4[0] <= d3[0] and d4[0]<=0.0002):
        pred.append(i4[0]-int(len(ele)/2))
    else:
        pred.append((-1))


endTime = time.clock()

print("Prediction Done")
print("Prediction time : ", endTime - startTime)

# The prediction has been done in the code above. \
# The code below is just the analysis of the predicted results and is not the part of the pipeline


# Initialising the number of examples in each category to 0
# All the six types of categories are defined in the paper
category_values = [0, 0, 0, 0, 0, 0]


for ele1, ele2, kmer in zip(pred, kmer_ind, kmer_arr):
    if(ele1==-1 and ele2[0]==-1):
        category_values[0] += 1
    elif(ele1==-1):
        category_values[5] += 1
    elif(ele2[0]==-1):
        a = NW(kmer, ref_dna[ele1: ele1 + len(kmer)])
        if(a[0]>0.8):
            category_values[1] += 1
        category_values[3] += 1
    elif(ele1 - ele2[0] < 3 and ele1 - ele2[0] > -3):
        category_values[0] += 1
    else:
        a = NW(kmer, ref_dna[ele1: ele1 + len(kmer)])
        b = NW(kmer, ref_dna[ele2[0]: ele2[0] + len(kmer)])
        if(b[0] - a[0] > 0.1):
            category_values[4] += 1
        else:
            category_values[2] += 1

# Print percentages for all 6 categories
print(category_values/len(kmer))
