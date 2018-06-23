import numpy as np
import math
import json
import subprocess as sp
import create_random_networks_for_netzcore as crnn

'''
# Citation
 Guney E, Oliva B. Exploiting Protein-Protein Interaction Networks for Genome-Wide
 Disease-Gene Prioritization. PLoS ONE 7(9): e43557 (2012). 
 [link](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0043557)
'''

class Guild:
    def __init__(self):
        self.gene_NI = {}   # key: gene name, value: gene id
        self.gene_IN = {}  # key: gene id, value: gene name
        self. dis_sim_set = {}    # key: dis, value: a dic, key: similar dis value sim score
        self.ppi_net = np.zeros([1, 1])
        self.ppi_net_norm = np.zeros([1, 1])
        self.net_size = 0
        self.train_data = {}   # train data of dis-gene
        self.test_data = {}     # test data of dis-gene
        self.dis_sim_thld = 0.4     # similarity threshold
        self.dis_top_thld = 100   # top disease threshold
        self.test_dis_seed = {}
        self.ppi_net_deg = []
        self.ppi_net_avg_deg = 0
        self.n_sample = 100

        # using NetZcore
        self.command = "./guild -s z -n dg_data/NetZcore/can_gene.txt" \
                       " -e dg_data/NetZcore/subnet&rand_net/ppi_subnet.txt" \
                       " -o dg_data/NetZcore/pred_gene.txt " \
                       "-i 5 -d dg_data/NetZcore/subnet&rand_net/ppi_subnet.txt. -x 100"
        self.can_gene_file = 'C:\\cygwin64\\home\\yk\\guild\\DisGePred\\dg_data\\NetZcore\\can_gene.txt'
        self.pred_gene_file = 'C:\\cygwin64\\home\\yk\\guild\\DisGePred\\dg_data\\NetZcore\\pred_gene.txt'

        self.network_file = 'D:\\exp\\pred_dis_gene\\guild\\NetShort\\blab_ppi_2016.txt'
        self.subnet_file = 'C:\\cygwin64\\home\\yk\\guild\\DisGePred\\dg_data\\' \
                           'NetZcore\\subnet&rand_net\\ppi_subnet.txt'


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
        # self.column_normalize()

        # compute all degrees and average degree of ppi network
        # self.get_net_deg()

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
                # 预测疾病基因。
                predicted_genes = self.pred_g(test_dis)
                if predicted_genes is None: continue
                counter = 0
                for score, gene_name in predicted_genes:
                    # gene_name = self.gene_IN[gid]
                    # score_std = (500-score)/500
                    if gene_name in train_gene_set and gene_name not in train_genes:
                        fw.write('\t'.join([test_dis, gene_name, str(score)]) + '\n')
                        counter += 1
                        if counter >= 100 and counter >= len(test_genes) * 2: break
        fw.flush()
        fw.close()

    def pred_g(self, test_dis):

        data = []
        g_list = []
        # 先获得种子基因,输出到文件。
        with open(self.can_gene_file, 'w') as fw:
            fw.truncate()
            g_set = set()
            for g, score in self.test_dis_seed.get(test_dis).items():
                fw.write('\t'.join([g, str(score)]) + '\n')
                g_list.append(g)

        # 根据基因列表g_set抽取基因子网络
        count = self.extract_subnet(g_list)

        # 生成随机网络
        if count == 0:
            return None
        crnn.sample_network_preserving_topology(self.subnet_file, self.n_sample, self.subnet_file + ".")

        # 执行cpp程序，得到结果。
        sp.call(self.command)

        # 读取结果
        with open(self.pred_gene_file, 'r') as fr:
            for line in fr:
                g, score = line.split('\t')
                data.append([float(score), g])
        data.sort(reverse=True)
        return data

    # 设置算法参数,若没有，则使用默认参数
    def set_param(self, param):
        if param.get('dis_sim_thld') is not None:
            self.dis_sim_thld = param.get('dis_sim_thld')
        if param.get('dis_top_thld') is not None:
            self.dis_top_thld = param.get('dis_top_thld')

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

    def extract_subnet(self,gene_list):
        count = 0
        with open(self.subnet_file, 'w') as fw:
            fw.truncate()
            gene_num = len(gene_list)
            for i in range(gene_num-1):
                for j in range(i+1, gene_num):
                    g1 = gene_list[i]
                    g2 = gene_list[j]
                    g1_id = self.gene_NI.get(g1)
                    g2_id = self.gene_NI.get(g2)
                    if g1_id is not None and g2_id is not None:
                        if self.ppi_net[g1_id, g2_id] == 1:
                            count += 1
                            fw.write('\t'.join([g1, '1', g2]) + '\n')
        return count

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



def guild_run():

    # guild NetZcore
    # 加载疾病-症状数据，做疾病-基因预测。
    bp1 = 'D:\\exp\\pred_dis_gene\\guild\\data\\'
    bp2 = 'D:\\exp\\pred_dis_gene\\guild\\NetZcore\\result\\'
    ppi_file = 'D:\\exp\\pred_dis_gene\\data\\SP\\blab_ppi2016.txt'
    dis_sim_file = bp1 + 'hpo&orpha_dis_sim.txt'
    cv_file = 'D:\\exp\\pred_dis_gene\\data\\dis_gene_cv\\dis_gene_cv10.json'
    # cv_list = ['0']
    cv_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
    for cv_i in cv_list:
        # param = {'dis_sim_thld': 0.4, 'dis_top_thld': 100}  # 算法默认参数
        param = {'dis_sim_thld': 0.4, 'dis_top_thld': 200}
        out_file = bp2 + 'out_top_cv' + cv_i + '.txt'
        guild = Guild()
        guild.run(ppi_file, dis_sim_file, cv_file, cv_i, param, out_file)


if __name__ == '__main__':
    guild_run()