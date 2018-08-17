import numpy as np

from gensim.models import Word2Vec
from sklearn.decomposition import PCA
from matplotlib import pyplot

from numpy.random import seed
from tensorflow import set_random_seed

import matplotlib.pyplot as plt

gen_dict = {'A':0, 'C':1, 'T':2, 'G':3}
gen_rev_dict = {0:'A', 1:'C', 2:'T', 3:'G'}

def get_value(gen_str):
	global gen_dict
	global enc_len

	mult = 4**(enc_len-1)
	ans = 0
	for c in gen_str:
		ans += mult*gen_dict[c]
		mult = mult/4
	return ans

def get_kmer(value):
	global gen_rev_dict
	global enc_len

	ans = ""
	for i in range(0, enc_len):
		ans = gen_rev_dict[value%4] + ans
		value = int(value/4)

	return ans

def create_sent(ind):
	global gen_dict
	global gen_rev_dict
	global enc_len

	ref_gen = list(get_kmer(ind))
	full_data = []
	for i in range(0, enc_len):
		sent = []
		sent.append("".join(ref_gen))
		for j in range(0, 4):
			temp_c = ref_gen[i]
			ref_gen[i] = gen_rev_dict[j]
			sent.append("".join(ref_gen))
			ref_gen[i] = temp_c
		full_data.append(sent)

	sent = []
	pre_ref_gen = ref_gen[1:]
	for i in range(0, 4):
		pre_ref_gen.append(gen_rev_dict[i])
		sent.append("".join(pre_ref_gen))
		pre_ref_gen = pre_ref_gen[:-1]
	sent.append("".join(ref_gen))
	suf_ref_gen = ref_gen[:-1]
	for i in range(0, 4):
		suf_ref_gen = [gen_rev_dict[i]] + suf_ref_gen
		sent.append("".join(suf_ref_gen))
		suf_ref_gen = suf_ref_gen[1:]

	full_data.append(sent)

	return full_data

def collect_data(vocab_size):
	final_data = []
	for i in range(0, vocab_size):
		final_data.extend(create_sent(i))

	return final_data

enc_len = 4

vocab_size = 4**enc_len
data = collect_data(vocab_size)
fn_data = []

for sent in data:
	fn_sent = []
	for word in sent:
		fn_sent.append(get_value(word))
	fn_data.append(fn_sent)

print(len(fn_data))

print("Data Made")

model = Word2Vec(data, size=1, window=9, hs=1, iter=20, sg=1, min_count=1)
# model = Word2Vec.load('model.bin')
# summarize the loaded model
print(model)
# summarize vocabulary
# words = ["AATG", "AATC", "CGAC", "GACT", "CGAT", "CGTT", "GACC", "AAGG"]
# val = 0
# x = []
# y = []
# for ele in words:
# 	x.append(model[ele][0])
# 	y.append(val)
# x.append(1)
# y.append(5)


# print(x)
# print(y)

# fig, ax = plt.subplots()
# ax.scatter(x, y)

# # for i, txt in enumerate(words):
# #     ax.annotate(txt, (words[i], y[i]))

# plt.savefig("abc.png")

# save model
model.save('model.bin')
# load model
new_model = Word2Vec.load('model.bin')
print(new_model)