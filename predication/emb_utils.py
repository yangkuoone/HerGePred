from gensim.models.keyedvectors import KeyedVectors
from itertools import islice

# load vectors file.
def load_vectors(vectors_file):
    global nodes_name_dic, nodes_type_dic
    vecs = KeyedVectors.load(vectors_file)
    return vecs


# load vectors file from LINE or GenSim.
def get_emb(vecs_file):
    emb_dic = {}
    with open(vecs_file, 'r') as fr:
        for line in fr:
            arr = line.strip().split(' ')
            if len(arr) == 2: continue
            id = arr[0]
            arr1 = arr
            del arr1[0]
            arr2 = []
            for a in arr1: arr2.append(float(a))
            emb_dic[id] = arr2
    return emb_dic
