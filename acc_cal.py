from sklearn.neighbors import NearestNeighbors
from gensim.models import Word2Vec
import pickle
import time
import sys
from keras.models import model_from_json
import numpy as np
import pysam as  ps
import ReadFasta as fa

def encd1(inp_str):
    global model

    ans_arr = []
    i = 0
    n = len(inp_str)
    while (i < n - 3):
        enc = model[inp_str[i:i + 4]]
        ans_arr.extend(enc)
        i += 4

    return ans_arr

def get_val(inp_str):
    a = encd1(inp_str[:16])
    b = encd1(inp_str[-16:])

    print(inp_str[:16], inp_str[-16:])
    print(a, b)

model = Word2Vec.load('model.bin')

str1 = "GGGTTGGATCTGGTGATGGATACAGAAGGGGAAATCAAGGGCGATGATCGACAATCTTGTTGTCAATATTGAACCATTTTGAGGTCACACATAG"
str2 = "GGTTTGGCTCTGGTGATGGCTACAGAAGGGCAAATCAAGGGCGGTGATCGACAATTTTGTTGTCAATATTGAACCATTTTGAGGTCACACATAT"
str3 = "AACGCAAACACAACAACCGAGGAATGCCCATGAGTATGTCATCCATACCGTCGTCCTCCCAATCCGGGAAGCTCTATGGCTGGGTCGAAAGAAT"

get_val(str1)
get_val(str2)
get_val(str3)
