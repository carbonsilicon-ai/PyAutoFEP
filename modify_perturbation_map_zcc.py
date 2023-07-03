import rdkit.Chem
import rdkit.Chem.PropertyMol
from merge_topologies_zcc import find_mcs, constrained_two_mol
import pickle
import os
import savestate_util
import networkx
from generate_perturbation_map_atom_map_wangwei import fill_thermograph
import json
import argparse
def Modify_Map(data):
    try:
        add_edge = data['add_edge']
    except:
        add_edge = []
    try:
        del_edge = data['del_edge']
    except:
        del_edge = []
    # read the progress.pkl and modify the perturbation map and corresponding atom maps
    progress_data = savestate_util.SavableState('progress.pkl')

    new_graph = progress_data['thermograph']['last_solution']['best_solution']
    molecules_dict = progress_data['ligands_data']
    bias = progress_data['thermograph']['last_solution']['bias']

    perturbation_map = {(node_i, node_j): edge_data for node_i, node_j, edge_data in new_graph.edges.data()}

    # add edges
    add_edge_mcs = []
    if len(add_edge)> 0 or add_edge !=None:
        for a_edge in add_edge:
            a_edge = tuple(a_edge)
            mcs = find_mcs([molecules_dict[a_edge[0]]['molecule'], molecules_dict[a_edge[1]]['molecule']],
                           matchValences=True, ringMatchesRingOnly=True, completeRingsOnly=True,
                           savestate=None).smartsString
            print('mcs: ', a_edge, mcs)
            if mcs:
                add_edge_mcs.append(a_edge)
            else:
                print('There is not MCS between ligand {} and ligand {}'.format(molecules_dict[a_edge[0]]['molecule'],
                                                                                molecules_dict[a_edge[1]]['molecule']))

        print('add_edge_mcs: ', add_edge_mcs)

        # calculate the thermgraph of added edges
        # full_thermograph = networkx.DiGraph()

        molecules_dict_mol = {k: v['molecule'] for k, v in molecules_dict.items()}
        print('molecules_dict_mol: ', molecules_dict_mol)

        new_graph = fill_thermograph(new_graph, molecules_dict_mol, pairlist=add_edge_mcs, use_hs=False,
                                     threads=1, savestate=progress_data, custom_mcs=None,
                                     verbosity=0)

        # print('full_thermograph', new_graph)
        print('perturbation_graph_edges:', new_graph.edges.data())
        print('perturbation_graph_nodes:', new_graph.nodes.data())
        perturbation_map_add = {(node_i, node_j): edge_data for node_i, node_j, edge_data in new_graph.edges.data()}

        # merge previous perturbation map and added perturbation map
        perturbation_map = {**perturbation_map, **perturbation_map_add}
        print('merged_perturbation_map_added: ', perturbation_map)


    # delete edges
    if len(del_edge) > 0 or del_edge != None:
        for de in del_edge:
            de = tuple(de)
            #print('de',de)
            #print('del_edge:',del_edge)
            if de not in perturbation_map.keys():
                print("Warning: edge {} does not exist.".format(de))
            else:
                del perturbation_map[de]
                new_graph.remove_edge(de[0], de[1])

        continue_delete = True  # the value should be decided by outer
        # test for the presence of disconnected ligands
        disconnected_ligands = [lig for lig in molecules_dict if
                                lig not in [k for tup in perturbation_map.keys() for k in tup]]
        if disconnected_ligands:
            print("Warning: ligands {} are not connect to any other ligands.".format(disconnected_ligands))
            if continue_delete:
                for dis_lig in disconnected_ligands:
                    new_graph.remove_node(dis_lig)
        print('merged_perturbation_map_deleted: ', perturbation_map)
        # save modified perturbation map to progress.pkl
        # perturbation_graph = networkx.DiGraph()
        # [perturbation_graph.add_edge(i, j, **data) for (i, j), data in perturbation_map.items()]
        print('perturbation_graph_edges:', new_graph.edges.data())
        print('perturbation_graph_nodes:', new_graph.nodes.data())

    progress_data['perturbation_map'] = perturbation_map.copy()
    progress_data['thermograph']['last_solution']['best_solution'] = new_graph.copy()
    progress_data.save_data()

    import matplotlib

    matplotlib.use('svg')
    import matplotlib.pyplot
    import networkx.drawing
    from copy import deepcopy

    bias = progress_data['thermograph']['last_solution']['bias']
    center_molecule = bias
    outer_edges = deepcopy(new_graph)
    outer_edges.remove_node(center_molecule)
    node_position = networkx.drawing.circular_layout(outer_edges, center=[0.0, 0.0])
    node_position[center_molecule] = [0.0, 0.0]
    print('full_thermograph:', new_graph)
    print('pos:', node_position)
    networkx.drawing.draw(new_graph, with_labels=True, pos=node_position)
    matplotlib.pyplot.savefig('best_graph_modified.svg')

def output_info_json():
    import pickle
    from rdkit import Chem

    input_data = 'progress.pkl'
    file_path = 'modified_map_info.json'

    with open(input_data, 'rb') as fh:
        data = pickle.load(fh)
    graph = data['thermograph']['last_solution']['best_solution']
    for node_i, node_j, edge_data in graph.edges.data():
        edge_data['mcs'] = Chem.MolToSmiles(edge_data['mcs'])

    map_info = {node_i + ',' + node_j: edge_data for node_i, node_j, edge_data in graph.edges.data()}
    data = json.dumps(map_info, indent=2, separators=(',', ': '))
    with open(file_path, 'w') as f:
        f.write(data)

def load_json(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data
def write_json(file_path):
    path_info = {}
    path_info['add_edge'] = [('FXR_76', 'FXR_84'), ('FXR_84', 'FXR_85')]
    path_info['del_edge'] = [('FXR_12', 'FXR_84'), ('FXR_12', 'FXR_74'), ('FXR_12', 'FXR_88')]
    data = json.dumps(path_info, indent=2, separators=(',', ': '))
    with open(file_path, 'w') as f:
        f.write(data)
if __name__ == '__main__':

    #file_path = 'modify_map_input.json'

    #write_json(file_path)
    parser = argparse.ArgumentParser(description='the path of check_ligands json file')
    parser.add_argument('--file_path', '-f', help='where is your json file?',
                        default='docs/tutorial01/check_ligs_input.json')
    args = parser.parse_args()
    try:
        data = load_json(args.file_path)
    except:
        print('cannot find your json file')
        exit()
    #add_edge = [('FXR_76', 'FXR_84'), ('FXR_84', 'FXR_85')]
    #del_edge = [('FXR_12', 'FXR_84'), ('FXR_12', 'FXR_74'), ('FXR_12', 'FXR_88')]
    Modify_Map(data)
    output_info_json()
