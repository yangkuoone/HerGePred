import argparse
import numpy as np
import networkx as nx
import node2vec
from gensim.models import Word2Vec
from gensim.models.keyedvectors import KeyedVectors
import os
import pickle
import time
import gzip
from datetime import datetime

def parse_args():
    '''
	Parses the node2vec arguments.
	'''
    parser = argparse.ArgumentParser(description="Run learn embeddings.")
    parser.add_argument('--input', nargs='?', default=' ',
                        help='Input walks path')
    parser.add_argument('--output', nargs='?', default='emb/karate.emb',
                        help='Embeddings path')
    parser.add_argument('--model', nargs='?', default='emb/karate.emb',
                        help='Embeddings path')
    parser.add_argument('--dimensions', type=int, default=128,
                        help='Number of dimensions. Default is 128.')
    parser.add_argument('--window-size', type=int, default=10,
                        help='Context size for optimization. Default is 10.')
    parser.add_argument('--iter', default=5, type=int,
                        help='Number of epochs in SGD')
    parser.add_argument('--workers', type=int, default=8,
                        help='Number of parallel workers. Default is 8.')
    return parser.parse_args()

class Sentences(object):
    def __init__(self, fname):
        self.fname = fname
    def __iter__(self):
        try:
            file_loc = self.fname
            for line in open(file_loc, mode='r'):
                line = line.rstrip('\n')
                words = line.split(" ")
                yield words
        except Exception:
            print("Failed reading file:")
            print(self.fname)

def learn_embeddings(walks):

    # Learn embeddings by optimizing the Skipgram objective using SGD.
    sentences = Sentences(walks)
    model = Word2Vec(sentences, size=args.dimensions, window=args.window_size, min_count=0,
                     workers=args.workers, iter=args.iter, negative=25, sg=1)
    print("defined model using w2v")
    # model.save_word2vec_format(args.output)
    model.wv.save_word2vec_format(args.output)
    model.wv.save(args.model)
    print("saved model in word2vec format")
    return

def learnEMB(args):
    # Pipeline for representational learning for all nodes in a graph.
    start_time = datetime.now()
    print('---- input parameters -----')
    print('iterations=', args.iter)
    print('window size=', args.window_size)
    print('dimensions=', args.dimensions)
    print('workers=', args.workers)
    print('--- run program ---')
    learn_embeddings(args.input)
    end_time = datetime.now()
    print('run time: ', (end_time-start_time))

def emb_main(args):
    bp = 'D:\\exp\\DisGePred\\pred_dis_gene\\n2v\\embSDGG_all\\'
    args.input = bp + 'walks.txt'
    args.output = bp + 'emb.txt'
    args.model = bp + 'emb.model'
    args.dimensions = 128
    args.iter = 5
    args.window_size = 10
    args.workers = 6
    learnEMB(args)


if __name__ == '__main__':
    args = parse_args()
    emb_main(args)
