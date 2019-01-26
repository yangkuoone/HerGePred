import cos_sim as cs
import emb_utils as eu

nodes_name_dic = {}   # key:node id, value:node name
nodes_NI_dic = {}   # key: node name, value:node id
nodes_type_dic = {}   # key:node id, value:node type
nodes_dic = {}

def load_node_type(nodes_file):
    global nodes_name_dic, nodes_type_dic, nodes_dic, nodes_NI_dic
    with open(nodes_file, 'r') as fr:
        for line in fr:
            name, id, type = line.strip().split('\t')
            nodes_name_dic[id] = name
            nodes_NI_dic[name] = id
            nodes_type_dic[id] = type
            if type not in nodes_dic:
                s = set()
                s.add(id)
                nodes_dic[type] = s
            else:
                nodes_dic[type].add(id)


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
            dis_dic.setdefault(dis, set())
            dis_dic[dis].add(gene)

    return known_edge_set, dis_set, gene_set, dis_dic


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


# load train and test edges and nodes from edge_train_test_file
def load_train_test_data(edge_train_test_file):
    train_dic = {}  # key: 'dis', 'gene', 'edge', value: the set
    test_dic = {}  # key: 'dis', 'gene', 'edge', value: the set
    train_edge_set = set()
    test_edge_set = set()
    train_dis_set = set()
    train_gene_set = set()
    test_dis_set = set()
    test_gene_set = set()
    train_dis_dic = {}  # key: test dis_name, value: all the test genes of the dis
    test_dis_dic = {}  # key: train dis_name, value: all the train genes of the dis
    with open(edge_train_test_file, 'r') as fr:
        for line in fr:
            dis, gene, type = line.strip().split('\t')
            if type == 'train':
                if dis not in train_dis_dic:
                    a = set()
                    a.add(gene)
                    train_dis_dic[dis] = a
                else:
                    train_dis_dic[dis].add(gene)
                train_edge_set.add(dis+'\t'+gene)
                train_dis_set.add(dis)
                train_gene_set.add(gene)
            elif type == 'test':
                if dis not in test_dis_dic:
                    a = set()
                    a.add(gene)
                    test_dis_dic[dis] = a
                else:
                    test_dis_dic[dis].add(gene)
                test_edge_set.add(dis + '\t' + gene)
                test_dis_set.add(dis)
                test_gene_set.add(gene)
    train_dic['dis'] = train_dis_set
    train_dic['gene'] = train_gene_set
    train_dic['edge'] = train_edge_set
    train_dic['dis_dic'] = train_dis_dic
    test_dic['dis'] = test_dis_set
    test_dic['gene'] = test_gene_set
    test_dic['edge'] = test_edge_set
    test_dic['dis_dic'] = test_dis_dic
    return train_dic, test_dic


# similarity-based prediction algorithm
# 2. vec_type = 1: vectors file is from gensim; vec_type = 0: vectors file is from line.
# 3. sim_thld: similarity threshold, to speed up the program, and avoid out of memory.
def edges_pre_TT(nodes_file, edge_train_test_file,
                 vectors_file, vec_type, pre_file, pre_times, topN):
    global nodes_name_dic, nodes_type_dic, nodes_dic, nodes_NI_dic
    fw = open(pre_file, 'w')
    fw.truncate()
    load_node_type(nodes_file)

    train_dic, test_dic = load_train_test_data(edge_train_test_file)
    vecs = None
    vec_dic = None
    if vec_type == 1:
        vecs = eu.load_vectors(vectors_file)
    else:
        vec_dic = eu.get_emb(vectors_file)
    pre_counter = 0
    all_node_num = len(nodes_name_dic)
    test_dis_dic = test_dic['dis_dic']
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
        if sims is None: continue
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

    fw.flush()
    fw.close()


def edges_pre_TT_main():

    edge_train_test_file = 'data/cv10_of0.txt'
    nodes_file = 'data/nodes_example.txt'
    vectors_file = 'data/emb.model'
    pre_file = 'data/prediction_results.txt'
    vec_type = 1  # vec_type=1: vectors file from gensim; vec_type=0: from LINE.
    pre_times = 2
    top_n = 100
    edges_pre_TT(nodes_file, edge_train_test_file, vectors_file,
                 vec_type, pre_file, pre_times, top_n)


if __name__ == '__main__':
    print('sim prediction')
    edges_pre_TT_main()
