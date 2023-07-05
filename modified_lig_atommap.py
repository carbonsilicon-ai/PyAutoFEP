
import savestate_util
import json
import argparse
import networkx as nx

#add_atom_pair = {('FXR_76','FXR_84'):[(0,3)],('FXR_84','FXR_85'):[(0,3),(4,4)]}
#del_atom_pair = {('FXR_76','FXR_84'):[(34,35),(38,39)],('FXR_84','FXR_85'):[(52,59),(56,55)]}

def Modified_AtomMap(data):
    try:
        add_atom_pair = data['add_atom_pair']
    except:
        add_atom_pair = []
    try:
        del_atom_pair = data['del_atom_pair']
    except:
        del_atom_pair = []

    progress_data = savestate_util.SavableState('progress.pkl')
    new_graph = progress_data['thermograph']['last_solution']['best_solution']
    molecules_dict = progress_data['ligands_data']

    all_perturbed_atoms = nx.get_edge_attributes(new_graph, 'perturbed_atoms')
    # delete atom_map
    map_info = nx.get_edge_attributes(new_graph, 'atom_map')
    for pair, maps in del_atom_pair.items():
        pair = tuple(pair.split(','))
        ori_map = map_info[pair]
        ori_per_atoms = all_perturbed_atoms[pair]
        maps = [tuple(mm) for mm in maps]
        add_per_num = len(maps)
        for map in maps:
            ori_map.remove(map)
        print(pair, ori_map, ori_per_atoms + add_per_num)
        new_graph.edges[pair[0], pair[1]]['atom_map'] = ori_map
        new_graph.edges[pair[0], pair[1]]['perturbed_atoms'] = ori_per_atoms + add_per_num

    # add atom_map
    map_info = nx.get_edge_attributes(new_graph, 'atom_map')
    all_perturbed_atoms = nx.get_edge_attributes(new_graph, 'perturbed_atoms')
    for pair, maps in add_atom_pair.items():
        pair = tuple(pair.split(','))
        ori_map = map_info[pair]
        maps = [tuple(mm) for mm in maps]
        del_per_atom = len(maps)
        ori_per_atoms = all_perturbed_atoms[pair]
        if ori_per_atoms - del_per_atom <= 0:
            print(pair, 'all perturbation atoms is', ori_per_atoms - del_per_atom, ' less than 1')
            continue
        ori_map.extend(maps)
        new_graph.edges[pair[0], pair[1]]['atom_map'] = ori_map
        new_graph.edges[pair[0], pair[1]]['perturbed_atoms'] = ori_per_atoms - del_per_atom

    progress_data['thermograph']['last_solution']['best_solution'] = new_graph.copy()
    progress_data.save_data()
    print('finished modified_lig_atommap')

def load_json(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data
def write_json(file_path):
    path_info = {}
    path_info['add_atom_pair'] = {('FXR_76,FXR_84'): [(0, 3)], ('FXR_84,FXR_85'): [(0, 3), (4, 4)]}
    path_info['del_atom_pair'] = {('FXR_76,FXR_84'): [(34, 35), (38, 39)], ('FXR_84,FXR_85'): [(52, 59), (56, 55)]}
    data = json.dumps(path_info, indent=2, separators=(',', ': '))
    with open(file_path, 'w') as f:
        f.write(data)

def output_info_json():
    import pickle
    from rdkit import Chem

    input_data = 'progress.pkl'
    file_path = 'modified_atommap_info.json'

    with open(input_data, 'rb') as fh:
        data = pickle.load(fh)
    graph = data['thermograph']['last_solution']['best_solution']
    for node_i, node_j, edge_data in graph.edges.data():
        edge_data['mcs'] = Chem.MolToSmiles(edge_data['mcs'])

    map_info = {node_i + ',' + node_j: edge_data for node_i, node_j, edge_data in graph.edges.data()}
    data = json.dumps(map_info, indent=2, separators=(',', ': '))
    with open(file_path, 'w') as f:
        f.write(data)

if __name__ == '__main__':

    '''
    file_path, write_json should be supplied by user or fronthead
    '''
    file_path = 'modify_atommap_input.json'
    write_json(file_path)

    parser = argparse.ArgumentParser(description='the path of modified atommap json file')
    parser.add_argument('--file_path', '-f', help='where is your json file?',
                        default='modify_atommap_input.json')
    args = parser.parse_args()
    try:
        data = load_json(args.file_path)
    except:
        print('cannot find your json file')
        exit()
    #add_edge = [('FXR_76', 'FXR_84'), ('FXR_84', 'FXR_85')]
    #del_edge = [('FXR_12', 'FXR_84'), ('FXR_12', 'FXR_74'), ('FXR_12', 'FXR_88')]
    Modified_AtomMap(data)
    output_info_json()
