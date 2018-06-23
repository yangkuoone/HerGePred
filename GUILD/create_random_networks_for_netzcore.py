"""
    BIANA: Biologic Interactions and Network Analysis
    Copyright (C) 2009  Javier Garcia-Garcia, Emre Guney, Baldo Oliva

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import networkx, random, sys

def main():
    # sample_network_preserving_topology("../data/test_interactions.txt", 10, "../data/test_interactions.txt.")
    # network_file = sys.argv[1]
    # n_sample = int(sys.argv[2])

    '''
    # test function
    bp = 'D:\\exp\\pred_dis_gene\\guild\\NetZcore\\test\\'
    network_file = bp + 'test_interactions.txt'
    n_sample = 100
    '''

    bp = 'D:\\exp\\pred_dis_gene\\guild\\NetZcore\\ppi_rand\\'
    network_file = bp + 'blab_ppi_2016.txt'
    n_sample = 100


    '''
    network_file = 'C:\\cygwin64\\home\\yk\\guild\\DisGePred\\dg_data\\test1\\edges.txt'
    n_sample = 100
    '''
    sample_network_preserving_topology(network_file, n_sample, network_file + ".")


def sample_network_preserving_topology(network_sif_file, n_sample, output_prefix):
    g = create_network_from_sif_file(network_file_in_sif=network_sif_file, use_edge_data=True)
    for i in range(1, n_sample+1):
        g_sampled = randomize_graph(graph=g, randomization_type="preserve_topology_and_node_degree")
        output_network_in_sif(g_sampled, output_prefix+"%s"%i)


def create_network_from_sif_file(network_file_in_sif, use_edge_data=False, delim=None, include_unconnected=True):
    setNode, setEdge, dictDummy, dictEdge = get_nodes_and_edges_from_sif_file(network_file_in_sif, store_edge_type=use_edge_data, delim=delim)
    g = networkx.Graph()
    if include_unconnected:
        g.add_nodes_from(setNode)
    if use_edge_data:
        for e, w in dictEdge.items():
            u, v = e
            g.add_edge(u, v, {"w": w})
    else:
        g.add_edges_from(setEdge)
    return g

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


def output_network_in_sif(g, output_file_name, delim=" ", include_unconnected=False):
    f = open(output_file_name, 'w')
    included_nodes = set()
    for u, v in g.edges_iter():
        f.write("%s%s%s%s%s\n" % (u, delim, g.get_edge_data(u, v)["w"], delim, v))
        included_nodes.add(u)
        included_nodes.add(v)
    if include_unconnected:
        for u in g.nodes_iter():
            if u not in included_nodes:
                f.write("%s\n" % u)
    f.close()


def randomize_graph(graph, randomization_type, allow_self_edges=True):
    """
    Creates a random network from given network as a networkx graph
    randomization_type: 
        - "random": add same number of edges randomly between nodes of original graph
        - "preserve_topology": keep edges, shuffle nodes of original graph
        - "preserve_topology_and_node_degree": keep edges, shuffle nodes of original graph with the nodes of same degree
        - "preserve_degree_distribution": remove an edge between two random nodes with degrees k, l then add to two nodes with degrees k-1 & l-1, then shuffle nodes
        - "preserve_degree_distribution_and_node_degree": remove 2 random edges between a-b and c-d where degree(a)=degree(c) and degree(b)=degree(d) then add 2 edges between a-d and b-c, then shuffle nodes with the same degree
    - "erdos_renyi": creates a graph where edges are redistributed based on erdos renyi random model.
    - "barabasi_albert": creates a graph where edges are redistributed based on barabasi albert model (preferential attachment).
    """

    debug = False

    n_node = graph.number_of_nodes()
    n_edge = graph.number_of_edges()

    if randomization_type == "erdos_renyi":
        # raise Exception("Work in progress")
        p = float(2 * n_edge) / (n_node*n_node - 2*n_node)
        # Chooses each of the possible [n(n-1)]/2 edges with probability p
        new_graph = networkx.erdos_renyi_graph(n_node, p)
        mapping = dict(zip(new_graph.nodes(), graph.nodes()))
        new_graph = networkx.relabel_nodes(new_graph, mapping)
        available_edges = graph.edges()

        # Map graph from random model to new graph
        for edge in new_graph.edges():
            if len(available_edges) > 0:
                edge_org = available_edges.pop()
                if debug:
                    print("From random:", (edge[0], edge[1]))
                new_graph.add_edge(edge[0], edge[1], graph.get_edge_data(edge_org[0], edge_org[1]))
            # If the random model added too many edges
            else:
                if debug:
                    print("Removing:", edge)
            new_graph.remove_edge(edge[0], edge[1])

        # If the random model failed to add enough edges
        nodes = new_graph.nodes()
        for edge_org in available_edges:
            source_id = random.choice(nodes)
            target_id = random.choice(nodes)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
                source_id = random.choice(nodes)
                target_id = random.choice(nodes)
            if debug:
                print("Adding:", (source_id, target_id))
            new_graph.add_edge(source_id, target_id, graph.get_edge_data(edge_org[0], edge_org[1]))
        return new_graph

    if randomization_type == "barabasi_albert":
        # raise Exception("Work in progress")
        if n_edge >= n_node:
            # A graph of n nodes is grown by attaching new nodes each with m edges that are preferentially attached to
            # existing nodes with high degree
            new_graph = networkx.barabasi_albert_graph(n_node, n_edge / n_node)
            mapping = dict(zip(new_graph.nodes(), graph.nodes()))
            new_graph = networkx.relabel_nodes(new_graph, mapping)
        else:
            new_graph = networkx.create_empty_copy(graph)

        available_edges = graph.edges()
        degree_map = new_graph.degree()
        nodes = new_graph.nodes()

        # Map graph from random model to new graph
        for edge in new_graph.edges():
            if len(available_edges) > 0:
                edge_org = available_edges.pop()
                if debug:
                    print("From random:", (edge[0], edge[1]))
                new_graph.add_edge(edge[0], edge[1], graph.get_edge_data(edge_org[0], edge_org[1]))
            # If the random model added too many edges
            else:
                nodes_to_select = [ id for id, d in degree_map.items() for j in range(d+1) ]
                source_id = random.choice(nodes())
                target_id = random.choice(nodes_to_select)
                if debug:
                    print("Removing:", (source_id, target_id))
                new_graph.remove_edge(source_id, target_id)
                degree_map[source_id] -= 1
                degree_map[target_id] -= 1

        # If the random model failed to add enough edges
        for edge_org in available_edges:
            nodes_to_select = [ id for id, d in degree_map.items() for j in range(d+1) ]
            source_id = random.choice(nodes)
            target_id = random.choice(nodes_to_select)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
                source_id = random.choice(nodes)
                target_id = random.choice(nodes_to_select)
            if debug:
                print("Adding:", (source_id, target_id))
            new_graph.add_edge(source_id, target_id, graph.get_edge_data(edge_org[0], edge_org[1]))
            degree_map[source_id] += 1
            degree_map[target_id] += 1

        return new_graph

    new_graph = networkx.create_empty_copy(graph)
    # new_graph.add_nodes_from(graph.nodes())

    if randomization_type == "random":
        nodes = new_graph.nodes()
        for edge in graph.edges():
            source_id = random.choice(nodes)
            target_id = random.choice(nodes)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
                source_id = random.choice(nodes)
                target_id = random.choice(nodes)
            new_graph.add_edge(source_id, target_id, graph.get_edge_data(edge[0], edge[1]))
        
    elif randomization_type == "preserve_topology":    # shuffle_nodes
        nodes = graph.nodes()
        random_nodes = graph.nodes()
        random.shuffle(random_nodes)
        equivalences = dict([(nodes[i], random_nodes[i]) for i in range(len(nodes))])
        new_graph.add_edges_from([(equivalences[current_edge[0]],
                                   equivalences[current_edge[1]],
                                   graph.get_edge_data(current_edge[0], current_edge[1])) for current_edge in graph.edges()])

    # 默认执行这块代码
    elif randomization_type == "preserve_topology_and_node_degree":     # shuffle_nodes_within_same_degree
        nodes_by_degree = dict((degree, []) for degree in graph.degree().values())
        graph_degree = graph.degree()
        [nodes_by_degree[graph_degree[node]].append(node) for node in graph_degree]
        equivalences = {}
        for current_degree in nodes_by_degree.keys():
            nodes = nodes_by_degree[current_degree]
            random_nodes = list(nodes)
            random.shuffle(random_nodes)
            equivalences.update(dict([(nodes[i], random_nodes[i]) for i in range(len(nodes))]))
        new_graph.add_edges_from([(equivalences[current_edge[0]],
                                   equivalences[current_edge[1]],
                                   graph.get_edge_data(current_edge[0], current_edge[1])) for current_edge in graph.edges()])
        
    elif randomization_type == "preserve_degree_distribution":
        # add edges as well
        for current_node1, current_node2 in graph.edges():
            new_graph.add_edge(current_node1, current_node2, graph.get_edge_data(current_node1, current_node2))
        max_degree = sorted(graph.degree().values())[-1]
        # nodes_by_degree = dict( (degree,{}) for degree in graph.degree() )
        nodes_by_degree = dict((degree, {}) for degree in range(max_degree+1))
        graph_degree = graph.degree()
        [nodes_by_degree[graph_degree[node]].setdefault(node) for node in graph_degree]
        # print new_graph.nodes(), new_graph.edges()
        # print nodes_by_degree
        # if n_edge < MIN_NUMBER_OF_PERTURBATION:
        #    n_perturbation = random.randint(n_edge/2, n_edge)
        # else:
        #    n_perturbation = random.randint(MIN_NUMBER_OF_PERTURBATION, n_edge)
        n_perturbation = random.randint(n_edge/2, n_edge)
        for i in range(n_perturbation):
            n_trial = 0
            while True:
                n_trial += 1
                if n_trial > MAX_NUMBER_OF_TRIAL:
                    if debug:
                        print("Warning: Max number of trials exceeded in perturbation ", i)
                    break
                source_id = random.choice(new_graph.nodes())
                source_degree = new_graph.degree(source_id)
                while source_degree < 1:  
                    source_id = random.choice(new_graph.nodes())
                    source_degree = new_graph.degree(source_id)
                target_id = random.choice(new_graph.neighbors(source_id))
                target_degree = new_graph.degree(target_id)
                del nodes_by_degree[source_degree][source_id] 
                nodes_by_degree[source_degree-1].setdefault(source_id)
                del nodes_by_degree[target_degree][target_id] 
                nodes_by_degree[target_degree-1].setdefault(target_id)
                # not very important to check for cases where new_source = source (v.v. for targets)
                new_source_id = random.choice(nodes_by_degree[source_degree-1].keys())
                new_target_id = random.choice(nodes_by_degree[target_degree-1].keys())
            if debug:
                print(source_id, target_id, " / ", new_source_id, new_target_id)
            # check if going to add an existing edge or self edge
            if new_graph.has_edge(new_source_id, new_target_id) or new_source_id == new_target_id:
                del nodes_by_degree[source_degree-1][source_id]
                nodes_by_degree[source_degree].setdefault(source_id)
                del nodes_by_degree[target_degree-1][target_id]
                nodes_by_degree[target_degree].setdefault(target_id)
                continue
            if debug:
                print("rm %d %d" % (source_id, target_id))
            edge_data = new_graph.get_edge_data(source_id, target_id)
            new_graph.delete_edge(source_id, target_id)
            if debug:
                print("add %d %d" % (new_source_id, new_target_id))
            new_graph.add_edge(new_source_id, new_target_id, edge_data)
            del nodes_by_degree[source_degree-1][new_source_id]
            nodes_by_degree[source_degree].setdefault(new_source_id)
            del nodes_by_degree[target_degree-1][new_target_id]
            nodes_by_degree[target_degree].setdefault(new_target_id)
            break
        # self.randomize_graph(new_graph, "preserve_topology")

    elif randomization_type == "preserve_degree_distribution_and_node_degree":
        # add edges as well
        for current_node1, current_node2 in graph.edges():
            new_graph.add_edge(current_node1, current_node2, graph.get_edge_data(current_node1, current_node2))
        nodes_by_degree = dict((degree, {}) for degree in graph.degree().values())
        graph_degree = graph.degree()
        [nodes_by_degree[graph_degree[node]].setdefault(node) for node in graph_degree]
        
        # if n_edge < MIN_NUMBER_OF_PERTURBATION:
        #    n_perturbation = random.randint(1, n_edge)
        # else:
        #    n_perturbation = random.randint(MIN_NUMBER_OF_PERTURBATION, n_edge)
        n_perturbation = random.randint(n_edge/2, n_edge)
        for i in range(n_perturbation):
            source_id = random.choice(new_graph.nodes())
            source_degree = new_graph.degree(source_id)
            # find a node for which another node with the same degree exists
            # available_neighbors = []
            n_trial = 0
            while True:    # (len(nodes_by_degree[source_degree]) < 2 or len(available_neighbors) < 1):
                n_trial += 1
                if n_trial > MAX_NUMBER_OF_TRIAL:
                    if debug:
                        print("Warning: Max number of trials exceeded in perturbation ", i)
                    break
                source_id = random.choice(new_graph.nodes())
                source_degree = new_graph.degree(source_id)
                if len(nodes_by_degree[source_degree]) < 2:
                    continue
                available_neighbors = []
                # find a neighbor for which another node with the same degree exists
                for neighbor_id in new_graph.neighbors_iter(source_id):
                    if source_degree == new_graph.degree(neighbor_id): 
                        if len(nodes_by_degree[new_graph.degree(neighbor_id)]) > 2:
                            available_neighbors.append(neighbor_id)
                    else:
                        if len(nodes_by_degree[new_graph.degree(neighbor_id)]) > 1:
                            available_neighbors.append(neighbor_id)
                if len(available_neighbors) < 1:
                    continue
                target_id = random.choice(available_neighbors)
                target_degree = new_graph.degree(target_id)
                # select a new source node with different id
        n_trial2 = 0
        inner_break = False
        while True:
            n_trial2 += 1
            if n_trial2 > MAX_NUMBER_OF_TRIAL:
                if debug:
                    print("Warning: Max number of trials exceeded in perturbation ", i)
                inner_break = True
                break
            new_source_id = random.choice(nodes_by_degree[source_degree].keys())
            while new_source_id == source_id:
                new_source_id = random.choice(nodes_by_degree[source_degree].keys())
                new_available_neighbors = []
                # find a neighbor as new target node for which id is different
                # from target and has an id equivalent to target
                for neighbor_id in new_graph.neighbors_iter(new_source_id):
                    if target_degree == new_graph.degree(neighbor_id):
                            new_available_neighbors.append(neighbor_id)
                    if len(new_available_neighbors) < 1:
                        continue
                    new_target_id = random.choice(new_available_neighbors)
                    if len(new_available_neighbors) > 1:
                        while new_target_id == target_id:
                            new_target_id = random.choice(new_available_neighbors)
                            # print new_available_neighbors, new_target_id
                    else:
                        new_target_id = new_available_neighbors[0]
                    break
            if inner_break:
                break
            if debug:
                print(source_id, target_id, " / ", new_source_id, new_target_id)
            if source_id == new_target_id or new_source_id == target_id:
                continue
            if new_graph.has_edge(source_id, new_target_id) or new_graph.has_edge(new_source_id, target_id):
                continue
        if debug:
            print("rm %d %d" % (source_id, target_id))
            print("rm %d %d" % (new_source_id, new_target_id))
        edge_data_1 = new_graph.get_edge_data(source_id, target_id)
        edge_data_2 = new_graph.get_edge_data(new_source_id, new_target_id)
        new_graph.delete_edge(source_id, target_id)
        new_graph.delete_edge(new_source_id, new_target_id)
        if debug:
            print("add %d %d" % (source_id, new_target_id))
            print("add %d %d" % (new_source_id, target_id))
            new_graph.add_edge(source_id, new_target_id, edge_data_1)
            new_graph.add_edge(new_source_id, target_id, edge_data_2)

    else:
        raise Exception("Unknown randomization type %s" % randomization_type)

    return new_graph


if __name__ == "__main__":
    main()
