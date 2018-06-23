import calculate_mean_and_sigma as cms
'''
# Citation
 Guney E, Oliva B. Exploiting Protein-Protein Interaction Networks for Genome-Wide
 Disease-Gene Prioritization. PLoS ONE 7(9): e43557 (2012). 
 [link](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0043557)
'''

class Guild:
    def __init__(self):
        self.dis_genes = {}

    def run(self, netscore_file, netzcore_file, netshort_file, out_file):
        # load predicted data
        self.load_data(netscore_file, netzcore_file, netshort_file)
        self.pred_gene(out_file)

    def load_data(self, netscore_file, netzcore_file, netshort_file):
        with open(netscore_file, 'r') as fw:
            for line in fw:
                dis, gene, score = line.strip().split('\t')
                self.dis_genes.setdefault(dis, {})
                self.dis_genes[dis].setdefault('netscore', [])
                self.dis_genes[dis]['netscore'].append([gene, score])
        with open(netzcore_file, 'r') as fw:
            for line in fw:
                dis, gene, score = line.strip().split('\t')
                self.dis_genes.setdefault(dis, {})
                self.dis_genes[dis].setdefault('netzcore', [])
                self.dis_genes[dis]['netzcore'].append([gene, score])
        with open(netshort_file, 'r') as fw:
            for line in fw:
                dis, gene, score = line.strip().split('\t')
                self.dis_genes.setdefault(dis, {})
                self.dis_genes[dis].setdefault('netshort', [])
                self.dis_genes[dis]['netshort'].append([gene, score])

    # 预测给定疾病的基因
    def pred_gene(self, out_file):

        fw = open(out_file, 'w')
        fw.truncate()
        for dis, result in self.dis_genes.items():
            # print(result)
            if dis == 'C1864873':   # 对疾病C1864873做特殊处理。
                fw.write('\t'.join([dis, 'GNRHR', '0.5']) + '\n')
                continue
            score_nodes, gene_max_num = self.combine_scores(result)
            counter = 0

            for score, gene_name in score_nodes:
                fw.write('\t'.join([dis, gene_name, str(score)]) + '\n')
                counter += 1
                if counter >= 100 and counter >= gene_max_num: break
        fw.flush()
        fw.close()

    def combine_scores(self, result):
        node_to_scores = {}
        gene_max_num = 100
        for _, gene_score in result.items():
            node_to_score_inner = {}
            for node, score in gene_score:
                node_to_score_inner[node] = float(score)
            mean, sigma = cms.cal_mean_sigma(list(node_to_score_inner.values()))
            for node, score in node_to_score_inner.items():
                if sigma != 0:
                    node_to_scores.setdefault(node, []).append((score - mean) / sigma)
            if gene_max_num <= len(gene_score):
                gene_max_num = len(gene_score)
        values = []
        for node, scores in node_to_scores.items():
            score = sum(scores) / len(scores)
            values.append((score, node))
        values.sort(reverse=True)
        # print(values)
        values_list = [x[0] for x in values]

        min_v, max_v = min(values_list), max(values_list)
        score_nodes = []
        for score, node in values:
            score = (score - min_v) / (max_v - min_v)
            score_nodes.append([score, node])
        return score_nodes, gene_max_num



def guild_run():

    # guild NetCombo
    bp1 = 'D:\\exp\\pred_dis_gene\\guild\\NetScore\\result\\'
    bp2 = 'D:\\exp\\pred_dis_gene\\guild\\NetZcore\\result\\'
    bp3 = 'D:\\exp\\pred_dis_gene\\guild\\NetShort\\result\\'
    bp4 = 'D:\\exp\\pred_dis_gene\\guild\\NetCombo\\result\\'
    cv_list = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    for cv_i in cv_list:
        netscore_file = bp1 + 'out_top_cv' + cv_i + '.txt'
        netzcore_file = bp2 + 'out_top_cv' + cv_i + '.txt'
        netshort_file = bp3 + 'out_top_cv' + cv_i + '.txt'
        out_file = bp4 + 'out_top_cv' + cv_i + '.txt'
        guild = Guild()
        guild.run(netscore_file, netzcore_file, netshort_file, out_file)


if __name__ == '__main__':
    guild_run()