from ete3 import Tree
import re
import os
import argparse
import pandas as pd


# Returns a dictonary, where the key is the node, and the values is a dictionary of
# different transition types and their respective expectations
def get_expected_number_of_events(expectation_file_path):
    dict_node_expectations = {}
    dict_types = {1:"GAIN", 2:"LOSS", 3:"DUPLICATION", 4:"DEMI-DUPLICATION", 5:"BASE-NUMBER"}
    pattern = re.compile("#ALL EVENTS EXPECTATIONS PER NODE\nNODE\tGAIN\tLOSS\tDUPLICATION\tDEMI-DUPLICATION\tBASE-NUMBER\tTOMAX(.*?)#", re.DOTALL)
    file = open(expectation_file_path, 'r')
    content = file.read()
    file.close()
    section = pattern.findall(content)[0]
    expectations_per_nodes = re.split("\n", section)
    for line in expectations_per_nodes:
        if line == '':
            continue
        node_entries = re.split("[\s]+", line)
        node = node_entries[0]
        dict_node_expectations[node] = {}
        for i in dict_types:
            dict_node_expectations[node][dict_types[i]] = float(node_entries[i])
    return dict_node_expectations


# Returns a dictionary, where #model is the key and the respective nodes are the values
def get_branches_per_model(tree_with_shifts_file, tree_lengths):
    tree = Tree(tree_with_shifts_file, format=1)
    pattern = re.compile("([\S]+)-([\d]+)")
    dict_model_nodes = {}
    node_names = [(n.name, n) for n in tree.traverse()]
    for node_name, node in node_names:
        name, model = pattern.findall(node_name)[0]
        if node.is_root():
            continue
        father = node.up
        father_name, model_father = pattern.findall(father.name)[0]
        if int(model_father) != int(model):
            assigned_model = int(model_father)
        else:
            assigned_model = int(model)
        if assigned_model in dict_model_nodes:
            dict_model_nodes[assigned_model].append(name)
        else:
            dict_model_nodes[assigned_model] = [name]
        if not (assigned_model in tree_lengths):
            tree_lengths[assigned_model] = node.dist
        else:
            tree_lengths[assigned_model] += node.dist

    return dict_model_nodes


def get_rates_per_model(expectation_file_path, tree_with_shifts_file, output_file):
    expected_number_of_events_per_node = get_expected_number_of_events(expectation_file_path)
    models_tree_lengths = {}
    nodes_per_model = get_branches_per_model(tree_with_shifts_file, models_tree_lengths)
    tree_length = 0
    for model in models_tree_lengths:
        tree_length += models_tree_lengths[model]
    dict_df = {"model":[], "transition type":[], "rate":[]}
    models = sorted(list((nodes_per_model.keys())))
    for model in models:
        rates = {"GAIN":0, "LOSS":0, "DUPLICATION":0,
               "DEMI-DUPLICATION":0, "BASE-NUMBER":0}
        model_nodes = nodes_per_model[model]
        for node in model_nodes:
            expectations = expected_number_of_events_per_node[node]
            for transition_type in expectations:
                rates[transition_type] += expectations[transition_type]
        for transition_type in rates:
            rates[transition_type] /= models_tree_lengths[model]
            dict_df["model"].append(model)
            dict_df["transition type"].append(transition_type)
            dict_df["rate"].append(rates[transition_type])
    df = pd.DataFrame(dict_df)
    df.to_csv(output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='calculate rates')
    parser.add_argument('--res_dir_path', '-r', help='results files path')
    parser.add_argument('--output_file_path', '-o', help='table output path')
    # parse arguments
    args = parser.parse_args()
    res_dir_path = args.res_dir_path
    output_file_path = args.output_file_path
    shifts_tree_path = os.path.join(res_dir_path, "treeWithShifts.tree")
    expectation_file_path = os.path.join(res_dir_path, "expectations_second_round.txt")
    get_rates_per_model(expectation_file_path, shifts_tree_path, output_file_path)

