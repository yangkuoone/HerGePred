from gensim.models import Word2Vec
from gensim.models.keyedvectors import KeyedVectors
import cos_sim as cs
import emb_utils as eu
import dataset_partition as dp

nodes_name_dic = {}   # key:node id, value:node name
nodes_NI_dic = {}   # key: node name, value:node id
nodes_type_dic = {}   # key:node id, value:node type
nodes_dic = {}

def load_node_type(nodes_file):
    global nodes_name_dic, nodes_type_dic, nodes_dic, nodes_NI_dic
    with open(nodes_file, 'r') as fr:
        for line in fr:
            arr = line.strip().split('\t')
            name = arr[0]
            id = arr[1]
            type = arr[2]
            nodes_name_dic[id] = name
            nodes_NI_dic[name] = id
            nodes_type_dic[id] = type
            if type not in nodes_dic:
                s = set()
                s.add(id)
                nodes_dic[type] = s
            else:
                nodes_dic[type].add(id)

def get_pre_config(pre_config):
    global nodes_name_dic, nodes_type_dic, nodes_dic
    for key, value in pre_config.items():
        pre_node_type = key
        related_node_type = value
    pre_node_set = nodes_dic[pre_node_type]
    related_node_set = set()
    for type in related_node_type:
        node_set = nodes_dic[type]
        for node in node_set:
            related_node_set.add(node)
    return pre_node_set, related_node_set
# get most similar topn nodes of pre_node_id.
def get_most_similar(vec_dic, pre_node_id, topn):
    most_similar_nodes = []
    if pre_node_id not in vec_dic: return None
    pre_node_vec = vec_dic[pre_node_id]
    sim_nodes = {}
    for node_id, node_vec in vec_dic.items():
        if node_id == pre_node_id: continue
        sim = cs.cos(pre_node_vec, node_vec)
        sim_nodes[node_id] = sim
    # sort sim_nodes
    sorted_nodes = sorted(sim_nodes.items(), key=lambda a: a[1], reverse=True)
    counter = 0
    for node, sim in sorted_nodes:
        most_similar_nodes.append([node, sim])
        counter += 1
        if counter >= topn: break
    return most_similar_nodes
def load_dis_gene(dis_gene_file):
    dis_set = set()
    gene_set = set()
    known_edge_set = set()
    dis_dic = {}    # key： dis name, value: gene set of the dis
    with open(dis_gene_file, 'r') as fr:
        for line in fr:
            dis, gene = line.strip().split('\t')
            known_edge_set.add(line.strip())
            dis_set.add(dis)
            gene_set.add(gene)
            if dis not in dis_dic:
                dis_dic[dis] = set([gene])
            else:
                dis_dic[dis].add(gene)
    return known_edge_set, dis_set, gene_set, dis_dic


# value normalization
def norm(sim_list):
    sim_norm_list = []
    if len(sim_list) == 1:
        sim_norm_list.append(1)
        return sim_norm_list
    max_val = max(sim_list)
    min_val = min(sim_list)
    for sim in sim_list:
        if max_val == min_val:
            sim_norm = 0
        else:
            sim_norm = (sim-min_val)/(max_val-min_val)
        sim_norm_list.append(sim_norm)
    return sim_norm_list

# load train and test edges from edge_train_test_file
def load_edge_train_test(edge_train_test_file):
    dis_set = set()
    gene_set = set()
    train_edge_set = set()
    test_edge_set = set()
    train_dis_dic  = {}  # key: test dis_name, value: all the test genes of the dis
    test_dis_dic = {}  # key: train dis_name, value: all the train genes of the dis
    with open(edge_train_test_file, 'r') as fr:
        for line in fr:
            dis, gene, type = line.strip().split('\t')
            dis_set.add(dis)
            gene_set.add(gene)
            if type == 'train':
                if dis not in train_dis_dic:
                    a = set()
                    a.add(gene)
                    train_dis_dic[dis] = a
                else:
                    train_dis_dic[dis].add(gene)
                train_edge_set.add(dis+'\t'+gene)
            elif type == 'test':
                test_edge_set.add(dis+'\t'+gene)
                if dis not in test_dis_dic:
                    a = set()
                    a.add(gene)
                    test_dis_dic[dis] = a
                else:
                    test_dis_dic[dis].add(gene)
    return train_dis_dic, test_dis_dic, dis_set, gene_set, train_edge_set, test_edge_set

# similarity-based prediction algorithm
# 2. vec_type = 1: vectors file is from gensim; vec_type = 0: vectors file is from line.
# 3. sim_thld: similarity threshold, to speed up the program, and avoid out of memory.
def edges_pre_TT(nodes_file, edge_train_test_file,
                 vectors_file, vec_type, pre_file, pre_times, topN):
    global nodes_name_dic, nodes_type_dic, nodes_dic, nodes_NI_dic
    fw = open(pre_file, 'w')
    fw.truncate()
    load_node_type(nodes_file)
    train_dic, test_dic = dp.load_train_test_data(edge_train_test_file)
    if vec_type == 1: vecs = eu.load_vectors(vectors_file)
    else: vec_dic = eu.get_emb(vectors_file)
    pre_counter = 0
    all_node_num = len(nodes_name_dic)
    test_dis_dic = test_dic['dis_dic']
    test_edge_set = test_dic['edge']
    train_gene_set = train_dic['gene']
    train_edge_set = train_dic['edge']
    test_dis_num = len(test_dis_dic)
    for test_dis_name in test_dis_dic:
        pre_counter += 1
        if test_dis_name not in nodes_NI_dic: continue
        pre_node_id = nodes_NI_dic[test_dis_name]
        if pre_counter % 100 == 0: print(pre_counter, '/', test_dis_num, ':', test_dis_name)
        if vec_type == 1:
            sims = vecs.most_similar(pre_node_id, topn=all_node_num)
        else:
            sims = get_most_similar(vec_dic, pre_node_id, all_node_num)
        if sims == None: continue
        counter = 0
        test_gene_num = len(test_dis_dic[test_dis_name])   # known gene num of test dis
        for temp in sims:
            node_id = temp[0]
            sim = temp[1]
            node_name = nodes_name_dic[node_id]
            pre_edge = test_dis_name + '\t' + node_name
            if node_name in train_gene_set and pre_edge not in train_edge_set:
                counter += 1
                fw.write(test_dis_name + '\t' + node_name + '\t' + str(round(sim, 4)) + '\n')
                if counter >= test_gene_num * pre_times and counter >= topN: break
                # if counter <= test_gene_num * pre_times or counter <= topN:

    fw.flush()
    fw.close()

# test: get_most_similar
def get_most_similar_test():
    bp1 = 'G:\\Work\\Rech\\DisGenePre\\LINE\\'
    vectors_file = bp1 + 'vecs1.emb'
    vec_dic = eu.get_emb(vectors_file)
    pre_node_id = '1'
    topn = 2
    sims = get_most_similar(vec_dic, pre_node_id, topn)
    print(sims)


def edges_pre_TT_main():
    # only_dg, with_sym, with_ppi, no_pg, with_pg
    cv_num = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    for i in cv_num:
        num = str(i)
        bp1 = 'G:\\Work\\rech\\DisGePred\\pred_gene\\n2v\\embDGSG\\cv10_of' + num + '\\'
        bp2 = 'G:\\Work\\rech\\DisGePred\\pred_gene\\n2v\\embDGSG\\cv10_of' + num + '\\input_file\\'
        edge_train_test_file = bp2 + 'cv10_of' + num + '.txt'
        nodes_file = bp2 + 'nodes.txt'
        vectors_file = bp2 + 'emb.model'
        pre_file = bp1 + 'out_top_g_100.txt'
        vec_type = 1  # vec_type=1: vectors file from gensim; vec_type=0: from LINE.
        pre_times = 2
        topN = 100
        edges_pre_TT(nodes_file, edge_train_test_file, vectors_file,
                     vec_type, pre_file, pre_times, topN)

def test():
    a = [1, 2, 3]
    print(max(a))
    print(min(a))

if __name__ == '__main__':
    print('sim prediction')
    edges_pre_TT_main()
