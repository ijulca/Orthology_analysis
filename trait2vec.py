#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 18:07:46 2023

@author: ijulca
"""
# from nltk.tokenize import sent_tokenize, word_tokenize
# import warnings
 
# warnings.filterwarnings(action = 'ignore')
 
import pandas as pd
import glob
import numpy as np
from time import time
from collections import defaultdict
import multiprocessing
import gensim
from gensim.models import Word2Vec
from gensim.models.phrases import Phrases, Phraser
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt
import seaborn as sns

def get_sentences(files):
    sentences = {}
    for file in files:
        for line in open(file):
            line = line.strip()
            data = line.split("\t")
            if data[0] not in sentences:
                sentences[data[0]] = []
            for hog in data[1].split(";"):
                sentences[data[0]].append(hog)
    keys = sorted([x+"-"+str(len(sentences[x])) for x in sentences], key=lambda x:int(x.split("-")[1]))
    n = keys[0].split("-")[1]
    keys = [x.split("-")[0] for x in keys]
    new_sentences = []
    for k in sentences:
        sen = [k]+sentences[k]
        new_sentences.append(sen)
    return new_sentences, keys, n

def clena_df(df):
    df = df.dropna().reset_index(drop=True) # removing the missing values
    print(df.shape)
    sent = [row.split(";") for row in df['HOGS']]
    phrases = Phrases(sent, min_count=30, progress_per=10000)
    bigram = Phraser(phrases)
    sentences = bigram[sent]
    word_freq = defaultdict(int)
    for sent in sentences:
        for i in sent:
            word_freq[i] += 1
    print(sorted(word_freq, key=word_freq.get, reverse=True)[:10])
    
    w2v_model = Word2Vec(min_count=1, # Ignores all words with total absolute frequency lower than this
                          window=int(4), #The maximum distance between the current and predicted word within a sentence
                          vector_size=300, # Dimensionality of the feature vectors - (50, 300)
                          sample=6e-5, # The threshold for configuring which higher-frequency words are randomly downsampled. Highly influencial. - (0, 1e-5)
                          alpha=0.03, # The initial learning rate - (0.01, 0.05)
                          min_alpha=0.0007, # Learning rate will linearly drop to min_alpha as training progresses. To set it: alpha - (min_alpha * epochs) ~ 0.00
                          negative=20, # If > 0, negative sampling will be used, the int for negative specifies how many "noise words" should be drown. If set to 0, no negative sampling is used. - (5, 20)
                          workers=2) # Use these many worker threads to train the model 
    
    t = time() 
    w2v_model.build_vocab(sentences, progress_per=10000) 
    print('Time to build vocab: {} mins'.format(round((time() - t) / 60, 2)))
    
    t = time()
    w2v_model.train(sentences, total_examples=w2v_model.corpus_count, epochs=30, report_delay=1)
    print('Time to train the model: {} mins'.format(round((time() - t) / 60, 2)))
    w2v_model.init_sims(replace=True)
    w2v_model.wv.most_similar(positive=["HD_ORYSJ"])
    return df

def model_word2vec(sentences, keys,n):
    # train model
    model = Word2Vec(sentences,
                          min_count=1, # Ignores all words with total absolute frequency lower than this
                          window=int(n), #The maximum distance between the current and predicted word within a sentence
                          vector_size=300, # Dimensionality of the feature vectors - (50, 300)
                          sample=6e-5, # The threshold for configuring which higher-frequency words are randomly downsampled. Highly influencial. - (0, 1e-5)
                          alpha=0.03, # The initial learning rate - (0.01, 0.05)
                          min_alpha=0.0007, # Learning rate will linearly drop to min_alpha as training progresses. To set it: alpha - (min_alpha * epochs) ~ 0.00
                          negative=5, # If > 0, negative sampling will be used, the int for negative specifies how many "noise words" should be drown. If set to 0, no negative sampling is used. - (5, 20)
                          workers=2) # Use these many worker threads to train the model 
   
    # model = Word2Vec(sentences, min_count=1)
    # Create Skip Gram model
    # model = Word2Vec(sentences, min_count = 1, vector_size = 100, window = 5, sg = 1)
    dades = []
    for k in keys:
        values = []
        for e in keys:
            d = model.wv.similarity(k, e)
            values.append(d)
        dades.append(values)
    df = pd.DataFrame(np.array(dades),index=keys, columns=keys)  
    ax = sns.clustermap(df, cmap="Blues", figsize=(10, 10), yticklabels=True)
    plt.savefig(path +"cluster_traits.png", bbox_inches='tight')
    print(df)


### main
path = "/home/ijulca/projects/QTL/"
files = glob.glob(path+"*/*.trait*")

sentences, keys, n = get_sentences(files)
model_word2vec(sentences,keys,int(n))
        
 
# Print results
# print("Cosine similarity between 'alice' " +
#           "and 'wonderland' - Skip Gram : ",
#     model2.wv.similarity('alice', 'wonderland'))
     
# print("Cosine similarity between 'alice' " +
#             "and 'machines' - Skip Gram : ",
#       model2.wv.similarity('alice', 'machines'))

# print("End")

# import io
# import re
# import string
# import tqdm

# import numpy as np

# import tensorflow as tf
# from tensorflow.keras import layers

# # Load the TensorBoard notebook extension
# ## load_ext tensorboard
# SEED = 42
# # AUTOTUNE = tf.data.AUTOTUNE
# # print(AUTOTUNE)

# sentence = "The wide road shimmered in the hot sun"
# tokens = list(sentence.lower().split())
# print(len(tokens))

# vocab, index = {}, 1  # start indexing from 1
# vocab['<pad>'] = 0  # add a padding token
# for token in tokens:
#   if token not in vocab:
#     vocab[token] = index
#     index += 1
# vocab_size = len(vocab)
# print(vocab)

# inverse_vocab = {index: token for token, index in vocab.items()}
# print(inverse_vocab)

# example_sequence = [vocab[word] for word in tokens]
# print(example_sequence)

# window_size = 2
# positive_skip_grams, _ = tf.keras.preprocessing.sequence.skipgrams(
#       example_sequence,
#       vocabulary_size=vocab_size,
#       window_size=window_size,
#       negative_samples=0)
# print(len(positive_skip_grams))


# for target, context in positive_skip_grams[:5]:
#   print(f"({target}, {context}): ({inverse_vocab[target]}, {inverse_vocab[context]})")


# # Get target and context words for one positive skip-gram.
# target_word, context_word = positive_skip_grams[0]

# # Set the number of negative samples per positive context.
# num_ns = 4

# context_class = tf.reshape(tf.constant(context_word, dtype="int64"), (1, 1))
# negative_sampling_candidates, _, _ = tf.random.log_uniform_candidate_sampler(
#     true_classes=context_class,  # class that should be sampled as 'positive'
#     num_true=1,  # each positive skip-gram has 1 positive context class
#     num_sampled=num_ns,  # number of negative context words to sample
#     unique=True,  # all the negative samples should be unique
#     range_max=vocab_size,  # pick index of the samples from [0, vocab_size]
#     seed=SEED,  # seed for reproducibility
#     name="negative_sampling"  # name of this operation
# )
# print(negative_sampling_candidates)
# print([inverse_vocab[index.numpy()] for index in negative_sampling_candidates])

# # Reduce a dimension so you can use concatenation (in the next step).
# squeezed_context_class = tf.squeeze(context_class, 1)

# # Concatenate a positive context word with negative sampled words.
# context = tf.concat([squeezed_context_class, negative_sampling_candidates], 0)

# # Label the first context word as `1` (positive) followed by `num_ns` `0`s (negative).
# label = tf.constant([1] + [0]*num_ns, dtype="int64")
# target = target_word


# # Generates skip-gram pairs with negative sampling for a list of sequences
# # (int-encoded sentences) based on window size, number of negative samples
# # and vocabulary size.
# def generate_training_data(sequences, window_size, num_ns, vocab_size, seed):
#   # Elements of each training example are appended to these lists.
#   targets, contexts, labels = [], [], []

#   # Build the sampling table for `vocab_size` tokens.
#   sampling_table = tf.keras.preprocessing.sequence.make_sampling_table(vocab_size)

#   # Iterate over all sequences (sentences) in the dataset.
#   for sequence in tqdm.tqdm(sequences):

#     # Generate positive skip-gram pairs for a sequence (sentence).
#     positive_skip_grams, _ = tf.keras.preprocessing.sequence.skipgrams(
#           sequence,
#           vocabulary_size=vocab_size,
#           sampling_table=sampling_table,
#           window_size=window_size,
#           negative_samples=0)

#     # Iterate over each positive skip-gram pair to produce training examples
#     # with a positive context word and negative samples.
#     for target_word, context_word in positive_skip_grams:
#       context_class = tf.expand_dims(
#           tf.constant([context_word], dtype="int64"), 1)
#       negative_sampling_candidates, _, _ = tf.random.log_uniform_candidate_sampler(
#           true_classes=context_class,
#           num_true=1,
#           num_sampled=num_ns,
#           unique=True,
#           range_max=vocab_size,
#           seed=seed,
#           name="negative_sampling")

#       # Build context and label vectors (for one target word)
#       context = tf.concat([tf.squeeze(context_class,1), negative_sampling_candidates], 0)
#       label = tf.constant([1] + [0]*num_ns, dtype="int64")

#       # Append each element from the training example to global lists.
#       targets.append(target_word)
#       contexts.append(context)
#       labels.append(label)

#   return targets, contexts, labels


