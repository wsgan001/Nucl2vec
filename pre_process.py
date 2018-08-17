import numpy as np
import sys
from gensim.models import Word2Vec

# Dictionary for Base pair to index and index to base pair conversions
genomeDict = {'A':0, 'C':1, 'T':2, 'G':3}
genomeRevDict = {0:'A', 1:'C', 2:'T', 3:'G'}

# Length of the input k-mer for creating encodings. (The value of 'k' in simpler terms)
encodingLength = 4
# Dimension of the output encoding from the input k-mer
encodingDimension = 1


# Convert a given k-mer into its corresponding index
def convertToIndex(gen_str):
	global genomeDict
	global encodingLength

	mult = 4**(encodingLength-1)
	ans = 0
	for c in gen_str:
		ans += mult*genomeDict[c]
		mult = mult/4
	return ans

# Convert the index value to the corresponding kmer
def convertToKmer(value):
	global genomeRevDict
	global encodingLength

	ans = ""
	for i in range(0, encodingLength):
		ans = genomeRevDict[value%4] + ans
		value = int(value/4)

	return ans

# Create a sentence for a given index value.
# Places the k-mer corresponding to this value at the center
# and all the k-mers at edit distance = 1 around it to form the sentence
def createSentence(ind):
	global genomeDict
	global genomeRevDict
	global encodingLength

	# Get the k-mer corresponding to the given value
	# and break it into an character array
	referenceKmer = list(convertToKmer(ind))

	full_data = []

	# Create sentences with substitutions having edit distance = 1
	for i in range(0, encodingLength):
		sent = []
		# Sentence comprising of 5 words (referenceKmer repeated), all of them one edit distance from each other
		sent.append("".join(referenceKmer))
		for j in range(0, 4):
			temp_c = referenceKmer[i]
			referenceKmer[i] = genomeRevDict[j]
			sent.append("".join(referenceKmer))
			referenceKmer[i] = temp_c
		full_data.append(sent)

	# Create sentences with indels having edit distance = 1
	sent = []
	kmerSuffix = referenceKmer[1:]
	for i in range(0, 4):
		kmerSuffix.append(genomeRevDict[i])
		sent.append("".join(kmerSuffix))
		kmerSuffix = kmerSuffix[:-1]
	sent.append("".join(referenceKmer))
	kmerPrefix = referenceKmer[:-1]
	for i in range(0, 4):
		kmerPrefix = [genomeRevDict[i]] + kmerPrefix
		sent.append("".join(kmerPrefix))
		kmerPrefix = kmerPrefix[1:]

	# Final data containing one sentence for the indels, and 'k' sentences for the substitutions
	full_data.append(sent)

	return full_data

# Create the data for training with number of sentences = (1+encodingLength)*vocabSize
def collectData(vocabSize):
	final_data = []
	# For every index, create a sentence for the corresponding k-mer
	for i in range(0, vocabSize):
		final_data.extend(createSentence(i))

	return final_data


if(len(sys.argv)!=2):
	print("ERROR : Required format -> python pre_process.py <model_filename>")
	exit()

vocabSize = 4**encodingLength

# Create data for training of the Skip-Gram model
data = collectData(vocabSize)

# Represent every word by a unique index, as required for skip-gram training
fn_data = []

for sent in data:
	fn_sent = []
	for word in sent:
		fn_sent.append(convertToIndex(word))
	fn_data.append(fn_sent)

print("Data Made")

# Training of Gensim Word2Vec model using the dataset generated above
# Window = 9 -> Since this is the maximum size of our sentence. Will vary for different encoding legnths
# hs = 1 -> Using heirarchial softmax
# iter = 20 -> Accuracy saturates at 20 iterations, 
# may vary for a different encoding length and encoding dimension combination
# sg = 1 -> Skip-Gram is used
model = Word2Vec(data, size=encodingDimension, window=9, hs=1, iter=20, sg=1, min_count=1)

print("Model made : ")
print(model)

# Save model into a file for later use
model.save(sys.argv[1])
