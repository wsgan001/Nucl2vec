from sklearn.neighbors import NearestNeighbors
from gensim.models import Word2Vec
import pickle
import time
import sys
from keras.models import model_from_json
import numpy as np
import pysam as  ps
import ReadFasta as fa
import csv

#-------------------------------------------------------
#This function returns to values for cae of match or mismatch
def Diagonal(n1,n2,pt):
    if(n1 == n2):
        return pt['MATCH']
    else:
        return pt['MISMATCH']

#------------------------------------------------------------
#This function gets the optional elements of the aligment matrix and returns the elements for the pointers matrix.
def Pointers(di,ho,ve):

    pointer = max(di,ho,ve) #based on python default maximum(return the first element).

    if(di == pointer):
        return 'D'
    elif(ho == pointer):
        return 'H'
    else:
         return 'V'

#--------------------------------------------------------
#This function creates the aligment and pointers matrices
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


    id = m-1
    jd = n-1
    while(id!=0 and jd!=0):
        if(p_mat[id][jd]=='D'):
            id = id-1
            jd = jd-1
        elif(p_mat[id][jd]=='H'):
            jd = jd-1
        elif(p_mat[id][jd]=='V'):
            id = id-1

    ie = m-1
    je = n-1
    while(p_mat[ie][je]=='H' and je>0):
        je = je-1

    return al_mat[ie][je]/len(s1), jd


def enc(inp_str, ind):
    if (ind == 1):
        return encd1(inp_str)
    elif (ind == 2):
        return encd2(inp_str)
    # elif(ind==3):
    #     return encd3(inp_str)
    elif (ind == 4):
        return encd4(inp_str)


def encd1(inp_str):
    global model
    global enc_len

    ans_arr = []
    i = 0
    n = len(inp_str)
    while (i < n - enc_len + 1):
        enc = model[inp_str[i:i + enc_len]]
        ans_arr.extend(enc)
        i += enc_len

    return ans_arr


def encd2(inp_str):
    global model

    ans_arr = []
    i = 0
    n = len(inp_str)
    while (i < n - 3):
        enc = model[inp_str[i:i + 4]]
        ans_arr.extend(enc)
        i += 4

    return ans_arr


def encd3(inp_str):
    global model

    ans_arr = []
    i = 0
    n = len(inp_str)
    while (i < n - 3):
        temp_arr = []
        for j in range(0, 4):
            enc = model[inp_str[i + j]]
            temp_arr.extend(enc)
        ans_arr.extend(temp_arr)
        i += 4

    return ans_arr


gen_rev_dict = {0: 'A', 1: 'C', 2: 'T', 3: 'G'}
rev_dict = {'A': 0, 'C': 1, 'T': 2, 'G': 3}


def get_kmer_arr(inp_str):
    global gen_rev_dict
    global rev_dict

    value = 0
    for char in inp_str:
        value = 4 * value
        value += rev_dict[char]

    value = int(value)

    ans = []
    for i in range(0, 4):
        temp = [0., 0., 0., 0.]
        temp[value % 4] += 1
        ans.extend(temp)
        value = int(value / 4)

    return ans


def encd4(inp_str):
    global model

    ans_arr = []
    i = 0
    n = len(inp_str)
    arr = []
    while (i < n - 3):
        arr.append(get_kmer_arr(inp_str[i:i + 4]))
        i += 4

    enc = model.predict(np.array(arr))

    for ele in enc:
        ans_arr.extend(ele)

    return ans_arr


# Reference DNA and kmer matrix. Default length of kmer = 32.
# Any kmer of length less than that should not be included. Any kmer of length more than that, only first 32 characters are taken.

# The pre-processing part, to create encoding of the refernce DNA and array of kmers.

# Output files -> nbrs + str(i) and kmers + str(i)
# Here 'i' represents the encoding used

bam_file = "SRR2724099_1_addRG.bam"
print("reading file: ", bam_file)

enc_len = 4

kmer_len = 4*enc_len
reads_found = 0

base2num = dict(zip("ACTG", (0, 1, 2, 3)))

bamfile = ps.AlignmentFile(bam_file, "r")

st_time = time.clock()

kmer_arr = []
kmer_ind = []

for SeqRead in bamfile:
    query_sequence = SeqRead.query_sequence
    leng = len(query_sequence)
    if (leng >= 32):
        reference_start = SeqRead.reference_start
        flag = 0

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
            if reads_found % 100000 == 0:
                print(reads_found)

print(reads_found)



# some refernce
dna = fa.readReference("sequence.fasta", "fasta")
ref_dna = dna[0]
print(len(ref_dna))
print(ref_dna[:100])
# array of all the kmers. Length should be > 32

# value of 'i'
enc_type = 1

if (enc_type == 1):
    model = Word2Vec.load('model.bin')
elif (enc_type == 2):
    model = Word2Vec.load('model2.bin')
elif (enc_type == 3):
    model = {'A': [1, 0], 'C': [-1, 0], 'T': [0, 1], 'G': [0, -1]}
elif (enc_type == 4):
    json_file = open('model3.json', 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    model = model_from_json(loaded_model_json)
    model.load_weights("model3.h5")
    model.compile(loss='mean_squared_error', optimizer='adam')

mid_time = time.clock()

X = []
for i in range(0, len(ref_dna) - kmer_len):
    if (i % 100000 == 0):
        print(i)
    X.append(enc(ref_dna[i:i + kmer_len], enc_type))

print(np.shape(X))

nbrs = NearestNeighbors(n_neighbors=1, algorithm='auto').fit(X)

Y = []
for kmr, ind in zip(kmer_arr, kmer_ind):
    Y.append(enc(kmr[:kmer_len], enc_type))

end_time = time.clock()

Y2 = []
kmer_ind2 = []
for kmr, inde in zip(kmer_arr, kmer_ind):
    Y2.append(enc(kmr[-kmer_len:], enc_type))
    kmer_ind2.append([inde[0] + len(kmr) - kmer_len])

Y3 = []
kmer_ind3 = []
for kmr, inde in zip(kmer_arr, kmer_ind):
    Y3.append(enc(kmr[int(len(kmr)/2)-kmer_len:int(len(kmr)/2)], enc_type))
    kmer_ind3.append([inde[0] + int(len(kmr)/2) - kmer_len])

Y4 = []
kmer_ind4 = []
for kmr, inde in zip(kmer_arr, kmer_ind):
    Y4.append(enc(kmr[int(len(kmr)/2):int(len(kmr)/2)+kmer_len], enc_type))
    kmer_ind4.append([inde[0] + int(len(kmr)/2)])


print("File Loading Time Taken : ", mid_time - st_time)
print("Preprocessing Time Taken : ", end_time - mid_time)
print("Total Time Taken : ", end_time - st_time)

st_time = time.clock()
dist, ind = nbrs.kneighbors(Y)
dist2, ind2 = nbrs.kneighbors(Y2)
dist3, ind3 = nbrs.kneighbors(Y3)
dist4, ind4 = nbrs.kneighbors(Y4)
end_time = time.clock()

mat_dist = 0.0
unm_dist = 0.0


mat_ttl = 0
unm_ttl = 0

corr = 0
ocorr = 0
mism = 0
inc = 0
cat3 = 0
nomatch = 0
ttl = 0


pred = []

iterr = 0

for i1, i2, i3, i4, d1, d2, d3, d4, ele in zip(ind, ind2, ind3, ind4, dist, dist2, dist3, dist4, kmer_arr):
    if (d1[0] <= d2[0] and d1[0] <= d3[0] and d1[0] <= d4[0] and d1[0]<=0.0002):
        # scr, loc = NW(ele, ref_dna[i1[0]-10:i1[0]+len(ele)+20])
        # if(scr>0.5):
        #     pred.append(i1[0]-10 + loc)
        # else:
        #     pred.append(-1)
        pred.append(i1[0])
    elif (d1[0] >= d2[0] and d3[0] >= d2[0] and d2[0] <= d4[0] and d2[0]<=0.0002):
        # scr, loc = NW(ele, ref_dna[i2[0]-len(ele)+kmer_len-10:i2[0]+kmer_len+20])
        # if(scr>0.5):
        #     pred.append(i2[0]-len(ele)+kmer_len-10 + loc)
        # else:
        #     pred.append(-1)
        pred.append(i2[0]-len(ele)+kmer_len)
    elif (d1[0] >= d3[0] and d2[0] >= d3[0] and d3[0] <= d4[0] and d3[0]<=0.0002):
        # scr, loc = NW(ele, ref_dna[i3[0]-int(len(ele)/2)-10:i3[0]-int(len(ele)/2)+len(ele)+20])
        # if(scr>0.5):
        #     pred.append(i3[0]-int(len(ele)/2)-10 + loc)
        # else:
        #     pred.append(-1)
        pred.append(i3[0]-int(len(ele)/2)+kmer_len)
    elif (d1[0] >= d4[0] and d2[0] >= d4[0] and d4[0] <= d3[0] and d4[0]<=0.0002):
        # scr, loc = NW(ele, ref_dna[i3[0]-int(len(ele)/2)-10:i3[0]-int(len(ele)/2)+len(ele)+20])
        # if(scr>0.5):
        #     pred.append(i3[0]-int(len(ele)/2)-10 + loc)
        # else:
        #     pred.append(-1)
        pred.append(i4[0]-int(len(ele)/2))
    else:
        pred.append((-1))
    # if(iterr%10000==0):
    #     print(iterr, len(kmer_arr))
    # iterr += 1



for ele1, ele2, kmer in zip(pred, kmer_ind, kmer_arr):
    # print("new")
    ttl += 1
    if(ele1==-1 and ele2[0]==-1):
        corr += 1
    elif(ele1==-1):
        nomatch += 1
        # print(NW(kmer, ref_dna[ele2[0]: ele2[0] + len(kmer)]), "a")
    elif(ele2[0]==-1):
        a = NW(kmer, ref_dna[ele1: ele1 + len(kmer)])
        if(a[0]>0.5):
            ocorr += 1
        inc += 1
    elif(ele1 - ele2[0] < 5 and ele1 - ele2[0] > -5):
        corr += 1
    else:
        a = NW(kmer, ref_dna[ele1: ele1 + len(kmer)])
        b = NW(kmer, ref_dna[ele2[0]: ele2[0] + len(kmer)])
        if(b[0] - a[0] > 0.1):
            mism += 1
        else:
            cat3 += 1
    if(ttl%10000==0):
        print(ttl)


print("Prediction time  : ", end_time - st_time)

# print(mat_dist/mat_ttl)

print("Prediction Accuracy 1  : ", corr / ttl)
print("Prediction Accuracy 1  : ", ocorr / ttl)
print("Prediction Accuracy 1  : ", mism / ttl)
print("Prediction Accuracy 1  : ", inc / ttl)
print("Prediction Accuracy 1  : ", cat3 / ttl)
print("Prediction Accuracy 1  : ", nomatch / ttl)
