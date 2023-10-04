#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 18:07:46 2023

@author: ijulca
"""
# from nltk.tokenize import sent_tokenize, word_tokenize
# import warnings
 
# warnings.filterwarnings(action = 'ignore')
 
import bs4 as bs
import urllib.request
import re
import nltk
import ssl
import pandas as pd
import glob
import gensim
import numpy as np
from gensim.models import Word2Vec
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt
import seaborn as sns

def get_sentences(files):
    sentences = []
    keys = []
    for file in files:
        for line in open(file):
            line = line.strip()
            data = line.split("\t")
            keys.append(data[0])
            sentence = [data[0]]+data[1].split(";")
            sentences.append(sentence)
    return sentences, keys

def model_word2vec(sentences, keys):
    # train model
    model = Word2Vec(sentences, min_count=1)
    ## window_size = 2 ### Takes up to two words in the right and left
    dades = []
    for k in keys:
        values = []
        for e in keys:
            d = model.wv.similarity(k, e)
            values.append(d)
        dades.append(values)
    df = pd.DataFrame(np.array(dades),index=keys, columns=keys)
    sns.clustermap(df, cmap="Blues")
    plt.savefig(path +"cluster_traits.png", bbox_inches='tight')
    print(df)


### main
path = "/home/ijulca/projects/QTL/"
files = glob.glob(path+"*/*.trait*")

# sentences, keys = get_sentences(files)
# model_word2vec(sentences,keys)
        
# # Create Skip Gram model
# model2 = gensim.models.Word2Vec(data, min_count = 1, vector_size = 100,
#                                              window = 5, sg = 1)
 
# # Print results
# print("Cosine similarity between 'alice' " +
#           "and 'wonderland' - Skip Gram : ",
#     model2.wv.similarity('alice', 'wonderland'))
     
# print("Cosine similarity between 'alice' " +
#             "and 'machines' - Skip Gram : ",
#       model2.wv.similarity('alice', 'machines'))

# print("End")

import io
import re
import string
import tqdm

import numpy as np

import tensorflow as tf
from tensorflow.keras import layers

# Load the TensorBoard notebook extension
## load_ext tensorboard
SEED = 42
AUTOTUNE = tf.data.AUTOTUNE
print(AUTOTUNE)

sentence = "The wide road shimmered in the hot sun"
tokens = list(sentence.lower().split())
print(len(tokens))













