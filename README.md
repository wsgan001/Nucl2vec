# Nucl2Vec : Local alignment of DNA sequences using Distributed Vector Representation
Prakhar Ganesh, Gaurav Gupta, Shubhi Saini, Kolin Paul

*This repository contains the code for the novel encoding proposed, Nucl2Vec. To know more, please refer to the paper.*

Biorxiv -> [preprint](https://www.biorxiv.org/content/early/2018/08/29/401851)

----


## Files included
> **pre_process.py** -> To train encodings before the actual prediction. Separately done beforehand. Need to decide the encoding length and encoding dimension. Learns and saves the encoding values.

> **alignment.py** -> Creating KNN Tree using the learned encodings and thus predicting alignments. Reads the reference DNA and creates KNN Tree using encodings for every k-mer present in the DNA. Segments are cut strategically from the reads to be aligned and encoded. Use these encodings corresponding to every read and the already made KNN Tree to predict alignment index.


## Execution Instructions
> python pre_process.py <model_name>

> python alignment.py <reads_fileName.bam> <reference_fileName.fasta> <model_name>


----
