# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 19:12:51 2016

@author: yk
"""
import numpy as np
# import Dic
from itertools import islice
import pickle

# global variable
disIN = {}          # key: dis id, value: dis name
disNI = {}          # key: dis name, value: dis id
geneIN = {}         # key: gene id, value: gene name
geneNI = {}         # key: gene name, value: gene id

#####################################
#  function: read gene-sim file
#  input:
#      geneSimFile: with title line
#  output:
#      geneMat: gene-sim matrix
#  time: 2016.11.11
#####################################
def getGeneMat(geneSimFile):
    print('load gene-gene sim file...')
    global geneIN, geneNI
    geneMat = np.eye(len(geneNI))
    with open(geneSimFile, 'r') as fr:
        for line in islice(fr, 0, None):
            ge1, ge2, sim = line.strip().split('\t')
            if ge1 not in geneNI or ge2 not in geneNI: continue
            gen1Id = geneNI[ge1]
            gen2Id = geneNI[ge2]
            sim = round(float(sim), 4)
            geneMat[gen1Id, gen2Id] = sim
            geneMat[gen2Id, gen1Id] = sim
    return geneMat

#  function: read dis-sim file
#  input:
#      disSimFile: with title line
#  output:
#     disMat: dis-sim matrix
#  time: 2016.11.11
def build_dic(dis_ge_file, dis_dis_file, ge_ge_file):
    global disNI, disIN, geneIN, geneNI
    dis_counter = 0
    gene_counter = 0
    with open(dis_ge_file, 'r') as fr:
        for line in fr:
            dis, gene, type = line.strip().split('\t')
            if type == 'train':
                if dis not in disNI:
                    disNI[dis] = dis_counter
                    disIN[dis_counter] = dis
                    dis_counter += 1
                if gene not in geneNI:
                    geneNI[gene] = gene_counter
                    geneIN[gene_counter] = gene
                    gene_counter += 1
    with open(dis_dis_file, 'r') as fr:
        for line in fr:
            dis1, dis2, _ = line.strip().split('\t')
            if dis1 not in disNI:
                disNI[dis1] = dis_counter
                disIN[dis_counter] = dis1
                dis_counter += 1
            if dis2 not in disNI:
                disNI[dis2] = dis_counter
                disIN[dis_counter] = dis2
                dis_counter += 1
    with open(ge_ge_file, 'r') as fr:
        for line in fr:
            ge1, ge2, _ = line.strip().split('\t')
            if ge1 not in geneNI:
                geneNI[ge1] = gene_counter
                geneIN[gene_counter] = ge1
                gene_counter += 1
            if ge2 not in geneNI:
                geneNI[ge2] = gene_counter
                geneIN[gene_counter] = ge2
                gene_counter += 1

def build_dic1(dis_ge_file, dis_dis_file, ge_ge_file):
    global disNI, disIN, geneIN, geneNI
    dis_counter = 0
    gene_counter = 0
    with open(dis_ge_file, 'r') as fr:
        for line in fr:
            dis, gene, type = line.strip().split('\t')
            if type == 'train':
                if dis not in disNI:
                    disNI[dis] = dis_counter
                    disIN[dis_counter] = dis
                    dis_counter += 1
                if gene not in geneNI:
                    geneNI[gene] = gene_counter
                    geneIN[gene_counter] = gene
                    gene_counter += 1

def getDisMat(disSimFile):
    print('load dis-dis sim file...')
    global disIN, disNI
    disMat = np.eye(len(disNI))
    with open(disSimFile, 'r') as fr:
        for line in islice(fr, 0, None):
            dis1, dis2, sim = line.strip().split('\t')
            if dis1 not in disNI or dis2 not in disNI: continue
            dis1Id = disNI[dis1]
            dis2Id = disNI[dis2]
            sim = round(float(sim), 4)
            disMat[dis1Id, dis2Id] = sim
            disMat[dis2Id, dis1Id] = sim
    return disMat

def load_train_test(train_test_file):
    train_dic = {}; test_dic = {}
    train_edge_set = set()
    train_dis_set = set()
    train_gene_set = set()
    test_edge_set = set()
    test_dis_set = set()
    test_gene_set = set()
    train_dis_dic = {}   # key: dis, value: the gene set
    test_dis_dic = {}  # key: dis, value: the gene set
    with open(train_test_file, 'r') as fr:
        for line in fr:
            dis, gene, type = line.strip().split('\t')
            if type == 'train':
                train_dis_set.add(dis)
                train_gene_set.add(gene)
                train_edge_set.add(dis+'\t'+gene)
                if dis not in train_dis_dic:
                    train_dis_dic[dis] = set([gene])
                else:
                    train_dis_dic[dis].add(gene)
            elif type == 'test':
                test_dis_set.add(dis)
                test_gene_set.add(gene)
                test_edge_set.add(dis + '\t' + gene)
                if dis not in test_dis_dic:
                    test_dis_dic[dis] = set([gene])
                else:
                    test_dis_dic[dis].add(gene)
    train_dic['edge'] = train_edge_set
    train_dic['dis'] = train_dis_set
    train_dic['gene'] = train_gene_set
    train_dic['dis_dic'] = train_dis_dic
    test_dic['edge'] = test_edge_set
    test_dic['dis'] = test_dis_set
    test_dic['gene'] = test_gene_set
    test_dic['dis_dic'] = test_dis_dic
    return train_dic, test_dic

#  function: read dis-gene relation file
#  input:
#      disGeneFile: with title line
#  output:
#      disGeneMat: dis-sim matrix
#  time: 2016.11.11
def getDisGeneMat(dis_ge_file):
    print('load dis-gene rel file...')
    global geneIN, geneNI, disIN, disNI
    disGeneMat = np.zeros((len(disNI), len(geneNI)))
    with open(dis_ge_file, 'r') as fr:
        for line in islice(fr, 0, None):
            dis, ge, type = line.strip().split('\t')
            if type == 'train':
                disId = disNI[dis]
                geneId = geneNI[ge]
                disGeneMat[disId, geneId] = 1
    return disGeneMat

# function: row normalization
# modified by yk, at 2017.12.5
#
def row_norm(mat):
    m = mat.copy()
    for i in range(0, len(m)):
        min_val = np.min(m[i])
        max_val = np.max(m[i])
        m[i] = (mat[i]-min_val)/(max_val-min_val)
    return m

############################
# Function: build matrix
# Parameter:
#    n:select top n in every row
# Return:
#    mat
###########################
def selTop(matrix, n):
    print('selTop ...')
    rowNum = len(matrix)
    retMat = np.zeros((rowNum, rowNum))
    if rowNum <= n:
        n = rowNum - 1
    for i in range(0, rowNum):
        sortArr = np.argsort(-matrix[i])
        #sortArr = su.sort(arr)
        for j in range(1, n+1):
            pos = sortArr[j]
            value = matrix[i][pos]
            retMat[i][pos] = value
    return retMat

#  Function: cal transition matrix
def getTransMat1(disMat1, geneMat1, disGeneMat, tao):
    global geneIN, geneNI
    UMat = row_norm(disMat1)
    VMat = row_norm(geneMat1)
    RMat = row_norm(disGeneMat)
    SMat = row_norm(disGeneMat.T)
    mat1 = np.concatenate(((1-tao)*UMat, tao*RMat), axis=1)
    mat2 = np.concatenate((tao*SMat, (1-tao)*VMat), axis=1)
    mat3 = np.concatenate((mat1, mat2), axis=0)
    W = row_norm(mat3)
    return W
    
# Function: cal L1 norm
# Parameter:
#    old and new is (n*1) matrix
# Return:
#    L1 norm
def getL1(old, new):
    L1 = 0 
    for i in range(0, len(old)):
        L1 = L1 + (new[i]-old[i])*(new[i]-old[i])
    return L1

def storeVarb(alpha, beta, dis_ge_file, dis_dis_file, ge_ge_file, varb_file):
    global disNI, disIN, geneIN, geneNI
    build_dic1(dis_ge_file, dis_dis_file, ge_ge_file)
    print('dis num:', len(disNI))
    print('gene num:', len(geneIN))
    disMat = getDisMat(dis_dis_file)
    geneMat = getGeneMat(ge_ge_file)
    disGeneMat = getDisGeneMat(dis_ge_file)
    disMat1 = selTop(disMat, alpha)
    geneMat1 = selTop(geneMat, beta)
    with open(varb_file, 'wb') as fw:
        fw.truncate()
        pickle.dump(geneIN, fw, True)
        pickle.dump(geneNI, fw, True)
        pickle.dump(disIN, fw, True)
        pickle.dump(disNI, fw, True)
        pickle.dump(disGeneMat, fw, True)
        pickle.dump(disMat1, fw, True)
        pickle.dump(geneMat1, fw, True)

def getV0(disGenesTrain, disId):
    global geneIN, geneNI, disIN, disNI
    #print('getV0')
    v0 = np.zeros((1, len(geneIN)))  # (1*n)
    genes = disGenesTrain[disId]
    for gene in genes:
        v0[0, gene] = 1
    return v0
#########################################
# Function: the pgWalk algorithm,
# Parameter:
#   alpha: number of top neighboring diseases of the highest sim score
#   beta:  number of top neighboring genes of the highest sim score
#   tao is transition probility
#   pai is transition probility
#   C is changed value of iteration
# Note:
#   LOOCV,set all genes of a disease to zeros
##########################################    
def pgWalk2(pai, C, tao,
           varbFile, train_test_file,
           outPreFile, outIterFile):
    global geneIN, geneNI, disIN, disNI
    outPreBW = open(outPreFile, 'w')
    outIterW = open(outIterFile, 'w')
    outPreBW.truncate()
    outIterW.truncate()
    
    with open(varbFile, 'rb') as fr:
        geneIN = dict(pickle.load(fr))
        geneNI = dict(pickle.load(fr))
        disIN = dict(pickle.load(fr))
        disNI = dict(pickle.load(fr))
        disGeneMat = np.matrix(pickle.load(fr))
        disMat1 = np.matrix(pickle.load(fr))
        geneMat1 = np.matrix(pickle.load(fr))
    print('dis num', len(disIN))
    print('gene num', len(geneIN))
    train_dic, test_dic = load_train_test(train_test_file)
    test_dis_dic = test_dic['dis_dic']
    train_dis_dic = train_dic['dis_dic']
    train_gene_set = train_dic['gene']
    test_edge_set = test_dic['edge']

    W = getTransMat1(disMat1, geneMat1, disGeneMat, tao)
    counter = 0
    test_dis_num = len(test_dis_dic)
    sim_dic = {}  # key: pre_node_id, node_id, value : the similarity of nodes
    for test_dis, test_g_set in test_dis_dic.items():
        counter += 1
        print(counter, '/', test_dis_num, ':', test_dis)
        train_g_set = train_dis_dic[test_dis]   # train target set of test_dis
        test_dis_id = disNI[test_dis]
        u0 = np.matrix(disMat1[test_dis_id])  # (1*m)
        v0 = np.matrix(disGeneMat[test_dis_id])   # (1*n)
        Y = np.concatenate((u0, v0), axis=1)  # 1*(m+n)
        P0 = Y.T
        P = Y.T
        iterNum = 0
        L1 = 100
        while (L1 > C or iterNum < 5):
            old = P  # i-th results
            P = (1 - pai) * np.dot(W.T, P) + pai * P0
            iterNum = iterNum + 1
            new = P  # (i+1)-th results
            L1 = getL1(old, new)
            if iterNum % 5 == 0: print('iterNum; L1:', iterNum, L1)
        # iter finish, write result to file
        outIterW.write(test_dis + '\t' + str(iterNum) + '\n')
        t_counter = 0

        score_dic = {}   # key: gene_id, value: score
        for j in range(len(disIN), len(P)):
            score = round(P[j, 0], 4)
            gene_id = j - len(disIN)
            if score > 0: score_dic[gene_id] = score
        sorted_score = sorted(score_dic.items(), key=lambda a: a[1], reverse=True)
        for temp in sorted_score:
            gene_name = geneIN[temp[0]]
            if gene_name not in train_g_set and gene_name in train_gene_set:
                t_counter += 1
                outPreBW.write(test_dis + '\t' + gene_name + '\t' + str(round(temp[1], 4)) + '\n')
                if t_counter >= len(test_g_set) * 2 and t_counter >= 100: break
    outPreBW.flush()
    outIterW.flush()
    outPreBW.close()
    outIterW.close()

def init_net():
    bp1 = 'G:\\Work\\rech\\DisGePred\\pred_gene\\n2v+hnrw\\cv10_of0\\input_file\\'
    bp2 = 'G:\\Work\\rech\\DisGePred\\pred_gene\\n2v+hnrw\\cv10_of0\\input_file\\'
    dis_ge_file = bp1 + 'cv10_of0.txt'
    dis_dis_file = bp2 + 'dis_sim_on_vecs.txt'
    ge_ge_file = bp2 + 'gene_sim_on_vecs.txt'
    varb_file = bp1 + 'out_varb.txt'
    alpha = 20
    beta = 100
    storeVarb(alpha, beta, dis_ge_file, dis_dis_file, ge_ge_file, varb_file)

def pgWalkMain():
    print('pgWalkMain...')
    bp1 = 'D:\\exp\\n2v+hnrw\\cv10_of0\\input_file\\'
    bp2 = 'D:\\exp\\n2v+hnrw\\cv10_of0\\'
    varb_file = bp1 + 'out_varb.txt'
    train_test_file = bp1 + 'cv10_of0.txt'
    out_pre_file = bp2 + 'outPre.txt'
    out_iter_file = bp2 + 'outIter.txt'
    # pai = 0.7
    # C = 1e-10
    # tao = 0.5
    pai = 0.7
    C = 1e-10
    tao = 0.8
    pgWalk2(pai, C, tao,
            varb_file, train_test_file,
            out_pre_file, out_iter_file)

if __name__ == '__main__':
    # init_net()
    pgWalkMain()