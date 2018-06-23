import numpy as np
import math
import random
import json

'''
# Citation
 Erten S, Bebek G, Ewing R M, et al. D A, D A : Degree-Aware Algorithms for Network-Based
  Disease Gene Prioritization[J]. BioData Mining, 2011, 4(1):19.
'''

class DaDa:
    def __init__(self):
        self.gene_NI = {}   # key: gene name, value: gene id
        self.gene_IN = {}  # key: gene id, value: gene name
        self. dis_sim_set = {}    # key: dis, value: a dic, key: similar dis value sim score
        self.ppi_net = np.zeros([1, 1])
        self.ppi_net_norm = np.zeros([1, 1])
        self.net_size = 0
        self.train_data = {}   # train data of dis-gene
        self.test_data = {}     # test data of dis-gene
        self.sig = 'CENT'
        self.hyb = 'SEED'
        self.dis_sim_thld = 0.4     # similarity threshold
        self.c = 0.3    # crosstalk restart probability. 0.3 is the optimal value
        self.dis_top_thld = 100   # top disease threshold
        self.test_dis_seed = {}
        self.ppi_net_deg = []
        self.ppi_net_avg_deg = 0

        # self.dis_NI = {}  # key: dis name, value: dis id
        # self.dis_IN = {}  # key: dis id, value: dis name
        # self. dis_counter = 0   # disease counter

        """
         self.train_data['gene'] = train_gene_set
        self.train_data['edge'] = train_edge_set
        self.train_data['dis_dic'] = train_dis_dic
        self.test_data['gene'] = test_gene_set
        self.test_data['edge'] = test_edge_set
        self.test_data['dis_dic'] = test_dis_dic
        """

    def run(self, ppi_file, dis_sim_file, cv_file, cv_i, param, out_file):
        # load ppi data
        self.load_ppi(ppi_file)

        # load disease similarity file
        self.load_dis_sim_file(dis_sim_file)

        # load cv file
        self.load_cv_file(cv_file, cv_i)

        # set the parameters of algorithms
        self.set_param(param)

        # get the normalized weighted network
        self.column_normalize()

        # compute all degrees and average degree of ppi network
        self.get_net_deg()

        self.pred_gene(out_file)

    # 预测给定疾病的基因
    def pred_gene(self, out_file):

        fw = open(out_file, 'w')
        fw.truncate()

        dis_counter = 0
        self.test_dis_seed = self.get_seed_gene()
        test_dis_dic = self.test_data['dis_dic']

        train_gene_set = self.train_data['gene']
        train_dis_dic = self.train_data['dis_dic']

        test_dis_num = len(test_dis_dic)
        for test_dis, test_genes in test_dis_dic.items():
            # get seed genes of test_dis
            dis_counter += 1
            print(dis_counter, '/', test_dis_num, ':', test_dis)
            # if dis_counter == 200: break
            train_genes = train_dis_dic[test_dis]
            if test_dis not in self.test_dis_seed:
                print('no seed genes')
                continue
            else:
                predicted_genes = self.pred_g(test_dis)
                if predicted_genes is None: continue
                counter = 0
                for score, gid in predicted_genes:
                    gene_name = self.gene_IN[gid]
                    score_std = (500-score)/500
                    if gene_name in train_gene_set and gene_name not in train_genes:
                        fw.write('\t'.join([test_dis, self.gene_IN[gid], str(score_std)]) + '\n')
                        counter += 1
                        if counter >= 100 and counter >= len(test_genes) * 2: break
        fw.flush()
        fw.close()
    def pred_g(self, test_dis):
        pass
        seed_genes = self.test_dis_seed[test_dis]
        if len(seed_genes) == 0:
            return None
        seed_gid = []
        for gene, score in seed_genes.items():
            if gene in self.gene_NI:
                gid = self.gene_NI.get(gene)
                seed_gid.append([gid, score])
        # get average degree of seed genes
        seed_avg_deg = self.get_avg_seed_deg(seed_gid)

        # calculate proximities using random walk with restarts
        c_vector, _ = self.seed_crosstalkers_nodisease(self.ppi_net_norm, self.c, seed_gid)
        # find the rank of the target gene
        c_vector_ss = [c_vector[i, 0] for i in range(len(c_vector))]

        m, rwr = self.my_sort(c_vector_ss, self.net_size - 1, 'descend')

        if self.sig == 'CENT':
            # significance based on centrality
            # print('statistical adjustment method: based on network centrality')
            # calculate pagerank vector
            p_vector, _ = self.seed_crosstalkers_nodisease(self.ppi_net_norm, 0, seed_gid)
            combined_vector = c_vector / p_vector   # 若 / 前后都是矩阵，则 / 代表点除。
            # sig_vector = combined_vector
            sig_vector = [combined_vector[i, 0] for i in range(len(combined_vector)) ]
            # remove non value of sig_vector
        """
        elif self.sig == 'CAND':
            # significance based on candidate degree
            print('statistical adjustment method: based on candidate degree')
            # sig_vector = np.zeros([1, candidate_size])
            sig_vector = []
            # for i in range(self.net_size):
                # curr_gene_id = candidate_gid[i]
                # s = self.candidate_significance(c_vector, curr_gene_id, weighted_degs)
                # sig_vector.append(s)
        elif self.sig == 'SEED':
            # significance based on seed size
            print('statistical adjustment method: based on seed degrees')
            # generate buckets based on degree distribution of seed set
            print('generating random seed sets...')
            buckets = self.generate_buckets(seed_gene)
            sig_vector = self.seed_significance(c, seed_gene, c_vector, buckets)
            sig_vector = [sig_vector[gid, 0] for gid in candidate_gid]
            # sig_vector = sig_vector(candidateSet)
            # remove non value of sig_vector
        else:
            print('Error. No such option!')
        """
        # UNIFORM PRIORITIZATION
        # results for RWR and statistical adjustment models are merged here
        # c_vector_subset = [c_vector[gid, 0] for gid in candidate_gid]
        if self.hyb == 'SEED':
            dada = self.hybrid_seeddegree_sort(c_vector_ss, sig_vector, 'descend', seed_avg_deg)
            # print('uniform prioritization method: based on average seed degree')
        elif self.hyb == 'OPT':
            m, DADA = self.hybrid_best_sort(c_vector_ss, sig_vector, 'descend', self.net_size)
            print('uniform prioritization method: best ranking for each gene (optimistic)')
        elif self.hyb == 'CAND':
            m, DADA = self.hybrid_candidatedegree_sort(c_vector_ss, sig_vector, 'descend')
            print('uniform prioritization method: based on average candidate degree')
        else:
            print('Error. No such option!')

        # print('rankings are computed')

        return dada

    # get seed genes of every test disease
    def get_seed_gene(self):
        print('get_seed_gene...')
        # speed up program and coding
        test_dis_dic = self.test_data['dis_dic']
        train_dis_dic = self.train_data['dis_dic']

        test_dis_seed = {}   # key: test dis, value: seed genes, a list
        for test_dis in test_dis_dic:
            if test_dis not in self.dis_sim_set:
                continue
            similar_diseases = self.dis_sim_set.get(test_dis)
            # get sorted dis list
            sorted_diseases = sorted(similar_diseases.items(), key=lambda x: x[1], reverse=True)
            # get first dis as seed gene's score
            # sim = sorted_diseases[0][1]
            for dis, sim in sorted_diseases:
                # if sim <= self.dis_sim_thld: break
                test_dis_seed.setdefault(test_dis, {})
                # get all the genes of dis
                genes = train_dis_dic.get(dis)
                if len(test_dis_seed[test_dis]) == self.dis_top_thld: break
                if genes is not None:
                    for g in genes:
                        if len(test_dis_seed[test_dis]) == self.dis_top_thld: break
                        if g not in test_dis_seed[test_dis]:
                            test_dis_seed[test_dis][g] = sim

        return test_dis_seed

    def hybrid_seeddegree_sort(self, c_vector, sig_vector, mode, seed_avg_deg):
        merged_vect = []

        c_vector_temp = self.sort2(c_vector, mode)
        sig_vector_temp = self.sort2(sig_vector, mode)
        # c_vector_temp = [a for a, _ in c_vector_sorted]
        # sig_vector_temp = [a for a, _ in sig_vector_sorted]

        for i in range(self.net_size):
            # print(i)
            c_number = c_vector[i]
            m1 = self.get_m(c_vector_temp, c_number)
            # m1, _ = self.my_sort(c_vector_subset, i, mode)
            # m2, _ = self.my_sort(sig_vector, i, mode)
            sig_number = sig_vector[i]
            m2 = self.get_m(sig_vector_temp, sig_number)
            # d = self.ppi_net_deg[i]
            if seed_avg_deg < self.ppi_net_avg_deg:
                merged_vect.append(m2)
                # merged_vect = merged_vect + [m2]
            else:
                merged_vect.append(m1)
                # merged_vect = merged_vect + [m1]
        m, I = self.my_sort(merged_vect, self.net_size-1, 'ascend')
        # I = self.sort2(merged_vect, self.net_size - 1, 'ascend')
        return I

    def get_m(self, vect, number):
        # vect    # key: number, value: pos list
        pos_list = vect[number]
        min_value = min(pos_list)
        max_value = max(pos_list)
        if min_value == max_value:
            m = min_value
            return m
        else:
            m = math.floor((min_value + max_value)/2)
            return m
        """
        first_pos = vect.index(number)  # 第一出现number的位置
        s = len(vect)
        for i in range(first_pos, s):
            if vect[i] != number:
                l = i - 1
                m = math.floor((first_pos + l) / 2)
                return m
        l = s
        m = math.floor((first_pos + l) / 2)
        return m
        """
    def sort2(self, vect, mode):

        vect_size = len(vect)
        value_pos = {}
        vect_temp = [[vect[i], i] for i in range(vect_size)]
        if mode == 'descend':
            vect_temp.sort(reverse=True)
        elif mode == 'ascend':
            vect_temp.sort(reverse=False)
        # sorted_vect = [a for a, _ in vect_temp]
        I = vect_temp
        for i in range(vect_size):
            pos = i
            value = vect_temp[i][0]
            value_pos.setdefault(value, [])
            value_pos[value].append(pos)
        return value_pos


    # index is the index of the number to be searched
    # score_gid:   [score, gid]
    def my_sort(self, vect, index, mode):
        # print('my_sort...')
        # number to be searched in the sorted vector
        number = vect[index]
        s = len(vect)
        vect_temp = [[vect[i], i] for i in range(len(vect))]
        if mode == 'descend':
            vect_temp.sort(reverse=True)
        elif mode == 'ascend':
            vect_temp.sort(reverse=False)
        sorted_vect = [a for a, _ in vect_temp]
        I = vect_temp

        first_pos = sorted_vect.index(number)  # 第一出现number的位置
        for i in range(first_pos, s):
            if sorted_vect[i] != number:
                l = i - 1
                m = math.floor((first_pos + l) / 2)
                return m, I
        l = s
        m = math.floor((first_pos + l) / 2)
        return m, I

        '''
        m = 0
        if mode == 'descend':
            for i in range(s):
                if sorted_vect[i] == number:  # number is the last value of vect
                    h = i
                    for j in range(i+1, s):
                        if sorted_vect[j] != number:
                            l = j - 1
                            m = math.floor((h+l)/2)
                            return m, I
                    l = s
                    m = math.floor((h+l)/2)
                    return m, I
        elif mode == 'ascend':
            for i in range(s):
                if sorted_vect[i] == number:
                    h = i
                    for j in range(i+1, s):
                        if sorted_vect[j] != number:
                            l = j - 1
                            m = math.floor((h+l)/2)
                            return m, I
                    l = s
                    m = math.floor((h+l)/2)
                    return m, I
        return m, I
        '''

    # 设置算法参数,若没有，则使用默认参数
    def set_param(self, param):
        if param.get('SIG') is not None:
            self.sig = param.get('SIG')
        if param.get('HYB') is not None:
            self.hyb = param.get('HYB')
        if param.get('dis_sim_thld') is not None:
            self.dis_sim_thld = param.get('dis_sim_thld')
        if param.get('dis_top_thld') is not None:
            self.dis_top_thld = param.get('dis_top_thld')

    # load disease similarities data
    def load_dis_sim_file(self, dis_sim_file):
        with open(dis_sim_file, 'r') as fr:
            for line in fr:
                d1, d2, sim = line.strip().split('\t')
                sim = float(sim)
                self.dis_sim_set.setdefault(d1, {d2: sim})
                self.dis_sim_set.setdefault(d2, {d1: sim})
                self.dis_sim_set[d1][d2] = sim
                self.dis_sim_set[d2][d1] = sim

    # load ppi data
    def load_ppi(self, ppi_file):
        print('read ppi data ...')
        edge = []
        counter = 0
        with open(ppi_file, 'r') as fr:
            for line in fr:
                p1, p2, w = line.strip().split('\t')
                if p1 not in self.gene_NI:
                    self.gene_NI[p1] = counter
                    self.gene_IN[counter] = p1
                    counter += 1
                if p2 not in self.gene_NI:
                    self.gene_NI[p2] = counter
                    self.gene_IN[counter] = p2
                    counter += 1
                p1_id = self.gene_NI[p1]
                p2_id = self.gene_NI[p2]
                edge.append([p1_id, p2_id])
        gene_num = len(self.gene_NI)
        self.net_size = gene_num
        self.ppi_net = np.zeros([gene_num, gene_num])
        for p1_id, p2_id in edge:
            self.ppi_net[p1_id, p2_id] = 1
            self.ppi_net[p2_id, p1_id] = 1
        # return ppi_net

    # simulate seed genes and candidate genes
    def simulate_seed_candidate(self):
        candidate_gene = ['2887', '6624', '85477', '5879', '60', '3624', '107', '117', '273', '196',
                          '1124', '2', '54205', '3486', '358', '5898', '84725', '2768', '9844', '5683',
                          '55698', '7975', '9734', '3205', '25928', '168667', '10268', '51608',
                          '10392', '1956']
        seed_gene = [['673', 1], ['2064', 1], ['5071', 1]]
        return candidate_gene, seed_gene

    # get seed genes
    def get_seed(seed_gene):
        seed_set = [g for g, _ in seed_gene]
        seed_scores = [s for _, s in seed_gene]
        return seed_set, seed_scores

    # column normalize
    def column_normalize(self):
        print('column_normalize...')
        # gene_num = len(self.ppi_net)
        self.ppi_net_norm = np.zeros([self.net_size, self.net_size])
        for i in range(self.net_size):
            deg = sum(self.ppi_net[:, i])
            if deg == 0:
                self.ppi_net_norm[:, i] = self.ppi_net[:, i]
            else:
                self.ppi_net_norm[:, i] = self.ppi_net[:, i]/np.linalg.norm(self.ppi_net[:, i], 1)
        return self.ppi_net_norm

    # load train and test edges and nodes from edge_train_test_file
    def load_cv_file(self, cv_file, i):
            train_edge_set = set()
            test_edge_set = set()
            train_gene_set = set()
            test_gene_set = set()
            train_dis_dic = {}  # key: test dis_name, value: all the test genes of the dis
            test_dis_dic = {}  # key: train dis_name, value: all the train genes of the dis
            with open(cv_file, 'r') as fr:
                json_dict = json.load(fr)
            data = json_dict[i]
            for dis, gene in data['d_train']:
                if dis not in train_dis_dic:
                    train_dis_dic[dis] = set([gene])
                else:
                    train_dis_dic[dis].add(gene)
                train_edge_set.add(dis + '\t' + gene)
                train_gene_set.add(gene)
            for dis, gene in data['e_test']:
                if dis not in test_dis_dic:
                    a = set()
                    a.add(gene)
                    test_dis_dic[dis] = a
                else:
                    test_dis_dic[dis].add(gene)
                test_edge_set.add(dis + '\t' + gene)
                test_gene_set.add(gene)

            self.train_data['gene'] = train_gene_set
            self.train_data['edge'] = train_edge_set
            self.train_data['dis_dic'] = train_dis_dic
            self.test_data['gene'] = test_gene_set
            self.test_data['edge'] = test_edge_set
            self.test_data['dis_dic'] = test_dis_dic

    # get average degree of seed genes
    def get_avg_seed_deg(self, seed_gid):
        deg_sum = 0
        for gid, _ in seed_gid:
            deg_sum += sum(self.ppi_net[:, gid])
        deg_avg = deg_sum/len(seed_gid)
        return deg_avg

    # compute all degrees of ppi network
    def get_net_deg(self):
        for i in range(self.net_size):
            self.ppi_net_deg.append(sum(self.ppi_net[:, i]))
        self.ppi_net_avg_deg = sum(self.ppi_net_deg)/len(self.ppi_net_deg)

    # finds crosstalk scores of all nodes wrt to the seed set
    # network is already column_normalized
    def seed_crosstalkers_nodisease(self, net, c, seed_gid):
        threshold = 1e-10
        maxit = 100
        residue = 1
        iter = 1
        prox_vector = np.zeros([self.net_size, 1])
        for gid, score in seed_gid:
            prox_vector[gid, 0] = score
        prox_vector = prox_vector/np.linalg.norm(prox_vector, 1)
        restart_vector = prox_vector
        while residue > threshold and iter < maxit:
            old_prox_vector = prox_vector
            prox_vector = (1 - c) * np.dot(net, prox_vector) + c * restart_vector
            residue = np.linalg. norm(prox_vector - old_prox_vector, 1)
            iter = iter + 1
        return prox_vector, residue

    def candidate_significance(self, c_vector, curr_gene_id, weighted_degs):
        pass
        curr_deg = weighted_degs[curr_gene_id]
        abs_degs = []
        for i in range(self.net_size):
            a = abs(weighted_degs[i] - curr_deg)
            abs_degs.append([a, i])
        abs_degs.sort()
        mean_cr = 0
        for i in range(1000):
            ind = abs_degs[i][0]
            cr = c_vector(ind)
            mean_cr = mean_cr + cr
        mean_cr = mean_cr / 1000
        std_dev = 0
        for i in range(1000):
            ind = abs_degs[i][0]
            cr = c_vector(ind)
        std_dev = std_dev + (cr - mean_cr) ** 2

        std_dev = math.sqrt(std_dev / 100)
        curr_cr = c_vector(curr_gene_id)
        sig = (curr_cr - mean_cr) / std_dev
        return sig

    def generate_buckets(self, seed_gene):
        print('generate buckets for the random seed sets')
        # 对于每个基因种子基因都有一个bucket，每个bucket的大小等于网络中的所有基因数目
        seed_size = len(seed_gene)
        buckets = []
        for i in range(seed_size):
            buckets.append([])
        # buckets{seedSize} = [];

        # weighted degrees
        seed_degs = []
        for i in range(seed_size):
            seed = seed_gene[i][0]
            seed_id = self.gene_NI[seed]

            weighted_degree = sum(self.ppi_net_norm[seed_id, :])
            seed_degs = seed_degs + [weighted_degree]

        for i in range(self.net_size):
            weighted_degree = sum(self.ppi_net_norm[i, :])
            min = 1000
            min_index = 0
            for j in range(seed_size):
                if abs(weighted_degree-seed_degs[j]) < min:
                    min = abs(weighted_degree-seed_degs[j])
                    min_index = j
            buckets[min_index] = buckets[min_index] + [i]
        return buckets

    # generates a random instance of a given seed set by applying the bucket
    # method explained in Nibbe et al.
    def generate_rand_seed_set(seed_gene, buckets):
        new_seed_set = []
        seed_size = len(seed_gene)
        for i in range(seed_size):
            bucket_size = len(buckets[i])
            r = math.ceil(bucket_size * random.random())
            if r > 0:
                new_seed_set = new_seed_set + [buckets[i][r]]
        new_seed_set = [[i, 1] for i in new_seed_set]
        return new_seed_set

    def seed_significance(self, c, seed_gene, c_vector, buckets):
        pass
        # to store all random vectors
        all_crosstalk_vectors = {}

        iter = 10
        mean_vector_crosstalk = np.zeros([self.net_size, 1])
        for z in range(iter):
            # just get random elements from the buckets
            rand_seed_set = self.generate_rand_seed_set(seed_gene, buckets)
            # crosstalk vector wrt random seed set
            c_vect, _ = self.seed_crosstalkers_nodisease(self.ppi_net_norm, c, rand_seed_set)
            all_crosstalk_vectors[z] = c_vect
            mean_vector_crosstalk = mean_vector_crosstalk + c_vect
        mean_vector_crosstalk = mean_vector_crosstalk / iter
        dev_vector_crosstalk = np.zeros(self.net_size, 1)
        for v in range(iter):
            temp = all_crosstalk_vectors[v] - mean_vector_crosstalk
            dev_vector_crosstalk = dev_vector_crosstalk + temp ** 2
        dev_vector_crosstalk = dev_vector_crosstalk / iter
        z_crosstalk = (c_vector - mean_vector_crosstalk) / math.sqrt(dev_vector_crosstalk)
        # remove non value of z_crosstalk
        return z_crosstalk

    def hybrid_best_sort(self, c_vector_subset, sig_vector, mode, candidate_size):
        merged_vect = []
        for i in range(candidate_size):
            m1, _ = self.my_sort(c_vector_subset, i, mode)
            m2, _ = self.my_sort(sig_vector, i, mode)
            merged_vect = merged_vect + min(m1, m2)
        m, I = self.my_sort(merged_vect, candidate_size, 'ascend')
        return m, I

    def hybrid_candidatedegree_sort(self, c_vector_subset, sig_vector, mode):
        merged_vect = []
        for i in range(self.net_size):
            m1, _ = self.my_sort(c_vector_subset, i, mode)
            m2, _ = self.my_sort(sig_vector, i, mode)
            d = self.ppi_net_deg[i]
            if d < self.ppi_net_avg_deg:
                merged_vect = merged_vect + m2
            else:
                merged_vect = merged_vect + m1
        m, I = self.my_sort(merged_vect, self.net_size-1, 'ascend')
        return m, I

def dada_run():
    # 设置参数r=0.4，依据是论文van Driel M A, Bruggeman J, Vriend G, et al. A text-mining analysis of the human phenome.[J].
    # European Journal of Human Genetics, 2006, 14(5):535-542.
    # 加载疾病-症状数据，做疾病-基因预测。
    bp1 = 'D:\\exp\\pred_dis_gene\\DaDa\\data\\'
    bp2 = 'D:\\exp\\pred_dis_gene\\DaDa\\'
    ppi_file = 'D:\\exp\\pred_dis_gene\\data\\SP\\blab_ppi2016.txt'
    dis_sim_file = bp1 + 'hpo&orpha_dis_sim.txt'
    cv_file = 'D:\\exp\\pred_dis_gene\\data\\dis_gene_cv\\dis_gene_cv10.json'
    # cv_list = ['0']
    cv_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
    for cv_i in cv_list:
        param = {'SIG': 'CENT', 'HYB': 'SEED', 'dis_sim_thld': 0.4, 'dis_top_thld': 50}  # 算法默认参数
        # param = {'SIG': 'SEED', 'HYB': 'CAND'}
        out_file = bp2 + 'out_top_cv' + cv_i + '.txt'
        dada = DaDa()
        dada.run(ppi_file, dis_sim_file, cv_file, cv_i, param, out_file)


if __name__ == '__main__':
    pass
    dada_run()