
import savestate_util
import json
import argparse
import networkx as nx
import process_user_input


def Modified_AtomMap(modify_data, progress_file):
    try:
        add_atom_pair = modify_data['add_atom_pair']
    except:
        add_atom_pair = []
    try:
        del_atom_pair = modify_data['del_atom_pair']
    except:
        del_atom_pair = []

    progress_data = savestate_util.SavableState(progress_file)
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
        new_graph.edges[pair[0], pair[1]
                        ]['perturbed_atoms'] = ori_per_atoms + add_per_num

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
            print(pair, 'all perturbation atoms is',
                  ori_per_atoms - del_per_atom, ' less than 1')
            continue
        ori_map.extend(maps)
        new_graph.edges[pair[0], pair[1]]['atom_map'] = ori_map
        new_graph.edges[pair[0], pair[1]
                        ]['perturbed_atoms'] = ori_per_atoms - del_per_atom

    progress_data['thermograph']['last_solution']['best_solution'] = new_graph.copy()
    progress_data.save_data()
    print('finished modified_lig_atommap')


def output_info_json(progress_file, map_info_file):
    import pickle
    from rdkit import Chem

    with open(progress_file, 'rb') as fh:
        data = pickle.load(fh)

    graph = data['thermograph']['last_solution']['best_solution']
    for node_i, node_j, edge_data in graph.edges.data():
        edge_data['mcs'] = Chem.MolToSmiles(edge_data['mcs'])
    map_info = {node_i+','+node_j: edge_data for node_i,
                node_j, edge_data in graph.edges.data()}
    datas = json.dumps(map_info, indent=2, separators=(',', ': '))
    with open(map_info_file, 'w') as f:
        f.write(datas)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='the path of modified atommap json file')
    parser.add_argument('--file_path', '-f', help='json file defined the atom map to be modified',
                        default='modify_atom_map_input.json')
    process_user_input.add_argparse_global_args(parser)
    arguments = process_user_input.read_options(
        parser, unpack_section='modify_atom_map')
    args = parser.parse_args()

    try:
        with open(args.file_path, 'r') as f:
            modify_data = json.load(f)
    except:
        print('cannot find your modify_atom_map_input json file')
        exit()

    Modified_AtomMap(modify_data, arguments.progress_file)
    output_info_json(arguments.progress_file, arguments.map_info_file)
