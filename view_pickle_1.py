
import pickle
import networkx as nx
input_data = 'progress.pkl'

with open(input_data, 'rb') as fh:
    data = pickle.load(fh)

print('data',data.keys())
graph = data['thermograph']['last_solution']['best_solution']
molecules_dict = data['ligands_data']
print('molecules_dict',molecules_dict)

#mcs_data = data['mcs_dict']
#for k, v in mcs_data.items():
#    print('k',k)
#    print('v',v.smartsString)

#print(data['mcs_dict'])
#print(graph.__dict__)
print('########')
#print(graph._succ)
#edge_data =  graph.edges.data()
#node_data = graph.nodes.data()
for node_i, node_j, edge_data in graph.edges.data():
    print(node_i,node_j)
    print('edge_data',edge_data)
#exit()
#aprint('edge_data',edge_data)
edges_info = nx.get_edge_attributes(graph,'perturbed_atoms')
#print(graph.edges(['iFXR_12','FXR_74']))
print(edges_info)
print(edges_info[('FXR_74','FXR_12')])
#print('node_data',node_data)
