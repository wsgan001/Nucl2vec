# Nucl2Vec : Local alignment of DNA sequences using Distributed Vector Representation
Prakhar Ganesh, Gaurav Gupta, Shubhi Saini, Kolin Paul

*This repository contains the code for the novel encoding proposed, Nucl2Vec. To know more, please refer to the paper here. (Link not added yet)*

--------------
Files included
> pre_process.py -> To traing encodings before the actual prediction start. Separately done beforehand. Need to decide the encoding length and encoding dimension. Learns and saves the encoding values.

> alignment.py -> Creating KNN Tree using the learned encodings and thus predicting alignments. Reads the reference DNA and creates KNN Tree using encodings for every k-mer present in the DNA. Segments are cut strategically from the reads to be aligned and encoded. Use these encodings corresponding to every read and the already made KNN Tree to predict alignment index.
