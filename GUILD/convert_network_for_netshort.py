#############################################################################
#    GUILD (Genes Underlying Inheritance Linked Disorders) implements several 
#    graph based algorithms for scoring relevance of a node in the network in 
#    terms of a phenotype using known associations in the node's neighborhood 
#    for that phenotype. GUILD has been applied to the prioritization of genes 
#    for several human disorders. 2011 - Emre Guney (Unviersitat Pompeu Fabra)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#############################################################################

import sys

def main():
    # create_edge_scores_as_node_scores_file("../data/test_proteins.txt", "../data/test_interactions.txt", "../data/test_interactions_for_netshort.txt")

    '''
    # source code
    node_file = sys.argv[1]
    network_file = sys.argv[2]
    output_file = sys.argv[3]
    '''
    bp = 'C:\\cygwin64\\home\\yk\\guild\\data\\'
    node_file = bp + 'test_proteins.txt'
    network_file = bp + 'test_interactions.txt'
    output_file = bp + 'test_interactions_for_netshort_test.txt'

    create_edge_scores_as_node_scores_file(node_file, network_file, output_file)

def create_edge_scores_as_node_scores_file(node_scores_file, edge_scores_file, output_file):
    """
    Creates edge score file from node association scores, intended comparing netshort with other algorithms without using other edge reliability/relevance score
    """
    nodes, edges, dictDummy, edge_to_score = get_nodes_and_edges_from_sif_file(edge_scores_file, store_edge_type=True, delim=None)
    nodes, setDummy, node_to_score, dictDummy = get_nodes_and_edges_from_sif_file(node_scores_file, store_edge_type=False, delim=None)
    f = open(output_file, 'w')
    for e, w in edge_to_score.items():
        u, v = e
        try:
            score_u =  node_to_score[u]
        except:
            raise ValueError("Missing node score for " + u)
        try:
            score_v =  node_to_score[v]
        except:
            raise ValueError("Missing node score for " + v)
        # w = float(w)
        f.write("%s %f %s\n" % (u, w*(score_u + score_v) / 2, v))
    f.close()

def get_nodes_and_edges_from_sif_file(file_name, store_edge_type=False, delim=None, data_to_float=True):
    """Parse sif file into node and edge sets and dictionaries
    returns setNode, setEdge, dictNode, dictEdge
    store_edge_type: if True, dictEdge[(u,v)] = edge_value
    delim: delimiter between elements in sif file, if None all whitespaces between letters are considered as delim
    """
    setNode = set()
    setEdge = set()
    dictNode = {}
    dictEdge = {}
    f = open(file_name, 'r')
    for line in f:
        if delim is None:
            words = line.split()
        else:
            words = line.split(delim)
        id1 = words[0]
        setNode.add(id1)
        if len(words) == 2:
            if data_to_float:
                score = float(words[1])
            else:
                score = words[1]
            dictNode[id1] = score
        elif len(words) == 3:
            id2 = words[2]
            setNode.add(id2)
            setEdge.add((id1, id2))
            if store_edge_type:
                if data_to_float:
                    dictEdge[(id1, id2)] = float(words[1])
                else:
                    dictEdge[(id1, id2)] = words[1]
    f.close()
    if len(setEdge) == 0:
        setEdge = None
    if len(dictNode) == 0:
        dictNode = None
    if len(dictEdge) == 0:
        dictEdge = None
    return setNode, setEdge, dictNode, dictEdge


if __name__ == "__main__":
    main()


