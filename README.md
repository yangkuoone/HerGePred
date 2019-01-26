# HerGePred
Heterogeneous disease-gene-related network embedding representation for disease gene predication

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
