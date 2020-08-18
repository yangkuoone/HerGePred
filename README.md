# HerGePred

## Basic Usage
### 1. Generate embedding vectors of nodes in heterogeneous network
+ Step (1): run n2v/e2v_walks.py
   + Input file: 
   data/edges_example.txt (i.e., all edges of a heterogeneous network)
   + Output file: 
   data/walks.txt (i.e., the result of random walks)
+ Step (2): run n2v/learn_vecs.py
   + Input file: 
   data/walks.txt (i.e., the result of random walks)
   + Output files: 
      + data/emb.txt (i.e., embedding vectors of nodes, text format)
      + data/emb.model (i.e., embedding vectors of nodes, binary format)


### 2. Predict disease genes based on embedding vectors
Run prediction/dis_gene_pred.py
   + Input files: 
      + data/cv10_of0.txt (i.e., train and test data)
      + data/nodes.txt (i.e., nodes in network)
      + data/emb.model (i.e., embedding vectors of nodes, binary format)
   + Output file:
      + data/prediction_results.txt (i.e., prediction results for disease genes)


## Citing
If you find HerGePred useful for your research, please consider citing the following paper:
```
@article{Yang2018HerGePred,
   author = {Yang, Kuo and Wang, Ruyu and Liu, Guangming and Shu, Zixin and Wang, Ning and Zhang, Runshun and Yu, Jian and Chen, Jianxin and Li, Xiaodong and Zhou, Xuezhong},
   title = {HerGePred: heterogeneous network embedding representation for disease gene prediction},
   journal = {IEEE Journal of Biomedical and Health Informatics},
   volume = {23},
   number = {4},
   pages = {1805-1815},
   year = {2018},
   type = {Journal Article}
}
```
K. Yang, R. Wang, G. Liu, Z. Shu, N. Wang, R. Zhang, J. Yu, J. Chen, X. Li, X. Zhou\*. HerGePred: Heterogeneous Network Embedding Representation for Disease Gene Prediction, IEEE Journal of Biomedical and Health Informatics, 2018, 23(4): 1805-1815.
