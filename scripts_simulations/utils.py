from ete3 import Tree
import re
import os


def get_tree_length(tree):
    sum_of_branch_lengths = 0
    for node in tree.traverse():
        sum_of_branch_lengths += node.dist
    return sum_of_branch_lengths


def rescale_tree(tree, factor):
    for node in tree.traverse():
        node.dist *= factor
    return


def read_counts_data(src_counts_path):
    file = open(src_counts_path, 'r')
    content = file.read()
    file.close()
    pattern = re.compile(">([\S]+)[\s]+([\S]+)", re.MULTILINE)
    return pattern.findall(content)


def write_heterogeneous_control_file(counts_file, tree_path, results_dir, model, control_file_path, nodes_file, number_of_models):
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    parameters = set_paths(counts_file, tree_path, results_dir)
    parameters += set_fixed_parameters()
    parameters += "_heterogeneousModel = true\n"
    parameters += "_parallelization = true\n"
    parameters += "_numOfModels = "+ str(number_of_models)+"\n"
    parameters += "_maxNumOfModels = "+ str(number_of_models+1)+"\n"
    if number_of_models > 1:
        parameters += "_minCladeSize = 20\n"
    else:
        parameters += "_minCladeSize = 5\n"
    if not (nodes_file is None):
        parameters += "_nodeIdsFilePath = "+ nodes_file +"\n"

    parameters += set_model_parameters(model, number_of_models)
    parameters += "_backwardPhase = false\n"
    control_file = open(control_file_path, 'w')
    control_file.write(parameters)
    control_file.close()



def write_homogeneous_model_control_file(working_dir, control_file_path, counts_path, tree_path, results_path, model):
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    if not os.path.exists(results_path):
        os.makedirs(results_path)
    parameters = set_paths(counts_path, tree_path, results_path)
    parameters += set_fixed_parameters()
    parameters += "_heterogeneousModel = false\n"
    parameters += "_parallelization = false\n"
    parameters += "_numOfModels = 1\n"
    parameters += "_maxNumOfModels = 1\n"
    parameters += set_model_parameters(model, 1)
    control_file = open(control_file_path, 'w')
    control_file.write(parameters)
    control_file.close()

def set_model_parameters(model, num_of_models):
    """
    Adds model specific parameters
    Given a model, where uppercase represents a function,
    lowercase represents rate type- the function sets the
    relevant parameters. Functions are as following:
    L (linear), E (exponent), C (constant), I (ignore).
    Rates: g (gain), l (loss), du (dupl), de (demi), b (base number)
    An example of a model: g+L_l+L_du+C_de+C_b+C
    """
    dic_rates_functions = get_functions(model)
    dic_rates = {"_gainFunc":"_gain_", "_lossFunc":"_loss_", "_duplFunc":"_dupl_", "_demiDuplFunc":"_demiPloidyR_",
                 "_baseNumRFunc":"_baseNumR_"}
    parameters = ""
    for func in dic_rates_functions:
        parameters += func + " = "+dic_rates_functions[func] +"\n"
    cat_counter = 1
    for i in range(num_of_models):
        if dic_rates_functions["_baseNumRFunc"] != "IGNORE":
            parameters += "_baseNum_"+ str(i+1)+ " = "+ str(cat_counter)+";4\n"
            cat_counter += 1
        for func in dic_rates:
            rate = dic_rates[func]
            function = dic_rates_functions[func]
            if function == "IGNORE":
                continue
            default_parameters = get_default_parameters_for_func(function)
            parameters += rate + str(i + 1) + " = " + str(cat_counter) + ";"+default_parameters+"\n"
            cat_counter += 1
    return parameters


def get_default_parameters_for_func(function):
    parameters = ""
    if function == "CONST":
        parameters += "2"
    elif function == "LINEAR":
        parameters += "2,0.1"
    elif function == "EXP":
        parameters += "2,0.01"
    else:
        raise Exception("get_default_parameters_for_func(): Not implemented yet!!")
    return parameters


def write_general_params(parameters):
    content = ""
    for param in parameters:
        content += param + " = "+ str(parameters[param]) + "\n"
    return content


def write_rate_parameters(rate_parameters):
    content = ""
    category = 1
    for model in rate_parameters:
        for param in rate_parameters[model]:
            content += "_"+param+"_"+ str(model)+" = "+ str(category)+";"
            for i in range(len(rate_parameters[model][param])):
                value = rate_parameters[model][param][i][1]
                if i != len(rate_parameters[model][param])-1:
                    content += str(value)+","
                else:
                    content += str(value)+"\n"
            category += 1
    return content


def test_max_on_tree(counts_file, tree_file):
    """
    calculates the maximal transition on the inferred tree
    :param counts_file: counts input file
    :param tree_file: ml inferred tree
    :return: if the maximal transitions if equal or larger than the original range -
                return 0
            otherwise
                return a list of the new base number and the max base transition
    """
    counts = get_counts(counts_file)
    counts_range = max(counts)-min(counts)
    max_base_on_tree = get_max_transition(tree_file)
    if max_base_on_tree >= counts_range:
        return 0
    max_base = max(max_base_on_tree, 3)
    return max_base


def get_max_transition(tree_file):
    """
    searches for the largest transition that was made on the phylogeny itself, between internal node and tip, or between internal nodes. Used in base_num_models.py
    :param tree_file: phylogeny file
    :return: a number representing the maximal transition on the tree
    """
    t = Tree(tree_file, format=1)
    max_transition = 0
    for node in t.traverse():
        if node.name == "":
            continue
        if not node.is_leaf():
            num1 = regex_internal(node.name)
            for child in node.get_children():
                if child.is_leaf():  # if the child is a tip - parse number from tip label
                    num2 = regex_tip(child.name)
                else:  # if the child is an internal node - take number
                    num2 = regex_internal(child.name)
                tmp_score = abs(num1 - num2)
                if max_transition < tmp_score:
                    max_transition = tmp_score
    return max_transition


def regex_tip(str):
    """
    extracts the number from a tip label in a phylogeny
    :param str: tip label
    :return: count of tip
    """
    tmp = re.search("(\d+)", str)
    if tmp:  # there is a number at the tip, and not X
        num = int(tmp.group(1))
    return num


def regex_internal(str):
    """
    extracts the number from the internal node's name from a phylogeny, in a NX-XX format. Used in get_max_transition(tree_file)
    :param str: internal node's label
    :return: count of node
    """
    tmp = re.search("N\d+\-(\d+)", str)
    if tmp:
        num = int(tmp.group(1))
    return num


def get_counts(counts_file):
    file = open(counts_file)
    content = file.read()
    file.close()
    pattern = re.compile(">[\S]+[\s]+([\d]+)", re.MULTILINE)
    counts_str = pattern.findall(content)
    counts = [int(count) for count in counts_str]
    return counts



def get_functions(model):
    dict_func_as_params = {}
    pattern = re.compile("([a-z]+)\+([A-Z]+)")
    dict_functions_names = {"g":"_gainFunc", "l":"_lossFunc", "du":"_duplFunc", "de":"_demiDuplFunc", "b":"_baseNumRFunc"}
    dict_fucntions_defs = {"L":"LINEAR", "C":"CONST", "E":"EXP", "I":"IGNORE"}
    rate_and_func = pattern.findall(model)
    for abb_rate, abb_func in rate_and_func:
        dict_func_as_params[dict_functions_names[abb_rate]] = dict_fucntions_defs[abb_func]
    return dict_func_as_params


def set_paths(counts_path, tree_path, results_path):
    parameters = ""
    parameters += "_dataFile = " + counts_path +"\n"
    parameters += "_treeFile = "+ tree_path +"\n"
    parameters += "_resultsPathDir = "+ results_path +"\n"
    return parameters


def set_fixed_parameters():
    parameters = ""
    parameters += "_maxChrNum = -10\n"
    parameters += "_minChrNum = -1\n"
    parameters += "_tolParamOptimization = 0.1\n"
    parameters += "_optimizePointsNum = 10,2,1\n"
    parameters += "_optimizeIterNum = 0,1,3\n"
    parameters += "_optimizeIterNumNextRounds = 0,1\n"
    parameters += "_optimizePointsNumNextRounds = 2,1\n"
    parameters += "_seed = 1\n"
    parameters += "_optimizationMethod = Brent\n"
    parameters += "_maxParsimonyBound = true\n"
    parameters += "_baseNumOptimizationMethod = Ranges\n"
    return parameters


def search_node_sizes(node):
    "Finds all node sizes"
    sizes = {}
    for n in node.traverse():
        size = len(n)
        if not (size in sizes):
            sizes[size] = [n]
        else:
            sizes[size].append(n)
    return sizes

#################################################################
# Simulations                                                   #
#################################################################


def get_params(res_path, freq_file, nodes_path_file, heterogeneous = "true"):
    """
    parses results file, writes a root frequency file and creates parameters dictionary
    :param res_path: chromEvol results file (chromEvol.res(=)
    :param freq_file: root frequencies file to be written to
    :param nodes_path_file: file which dictates the node splits in case of heterogeneous models
    :return: parameters list
    """
    parameters = {}
    res_file = open(res_path, "r")
    res_content = res_file.read()
    res_file.close()
    min_chr_num = get_min_chromosome_number(res_content)
    max_chr_num = get_max_chromosome_number(res_content)
    tree_scaling_factor = get_tree_scaling_factor(res_content)
    functions_params = get_functions_settings(res_content)
    parameters["_fixedFrequenciesFilePath"] = freq_file
    parameters["_minChrNum"] = str(min_chr_num)
    parameters["_maxChrInferred"] = str(max_chr_num)
    parameters["_maxChrNum"] = str(200)
    parameters["_branchMul"] = str(tree_scaling_factor)
    parameters["_heterogeneousModel"] = heterogeneous
    if heterogeneous == "true":
        parameters["_nodeIdsFilePath"] = str(nodes_path_file)
        num_of_models = 2
    else:
        parameters["_nodeIdsFilePath"] = "none"
        num_of_models = 1

    parameters["_numOfModels"] = str(num_of_models)
    parameters["_maxNumOfModels"] = str(num_of_models)
    parameters.update(functions_params)
    return parameters


def get_tree_scaling_factor(content):
    """
    extracts the tree scaling factor that should be used for the simulations
    :param content: chromEvol results file (chromEvol.res(=)
    :return: the tree scaling factor
    """
    pattern = re.compile("Tree scaling factor is:[\s]+([\S]+)")
    scaling_factor = float(pattern.findall(content)[0])
    return scaling_factor


def get_functions_settings(res_content):
    """
    extracts the functions set for each transition rate
    :param res_content: chromEvol results file content (chromEvol.res(=)
    :return: list of functions settings
    """
    max_num_of_functions = 5
    dict_of_functions = {}
    pattern_section = re.compile("Assigned functions for each rate parameter:[\s]+(.*)?#", re.DOTALL)
    pattern_functions = re.compile("([\S]+):[\s]+([\S]+)")
    section = pattern_section.findall(res_content)[0]
    functions_raw = pattern_functions.findall(section)
    dict_functions = {"baseNumR": "_baseNumRFunc", "dupl": "_duplFunc",
                      "demi": "_demiDuplFunc", "gain": "_gainFunc", "loss": "_lossFunc"}
    for i in range(min(max_num_of_functions, len(functions_raw))):
        transition_type, function = functions_raw[i]
        if transition_type in dict_functions:
            dict_of_functions[dict_functions[transition_type]] = function
    return dict_of_functions


def get_rate_parameters(res_file_path):
    """
    get the rate parameters from the results file
    :param res_file_path: path for the results chromEvol file
    :param freq_file: root frequencies file to be written to
    :return: a dictionary where the key is the model, and the value is a dictionary
            with a parameter as a key, and
    """
    res_file = open(res_file_path, 'r')
    file_content = res_file.read()
    res_file.close()
    pattern_section = re.compile("Best chosen model[\s]+[#]+[\s]+(.*?)(\*|$)", re.DOTALL)
    pattern_params = re.compile("Chromosome\.([\D]+)([\d]*)_([\d]+)[\s]+=[\s]+([\S]+)")
    section = pattern_section.findall(file_content)[0][0]
    params_index_models = pattern_params.findall(section)
    dict_res_param_name_input_names = {"dupl": "dupl", "gain": "gain", "loss": "loss", "demi": "demiPloidyR",
                                       "baseNum": "baseNum", "baseNumR": "baseNumR"}
    dict_model_parameters = {} # model: {param:[(index, value)}
    for rate_param, index_str, model_str, value_str in  params_index_models:
        param = dict_res_param_name_input_names[rate_param]
        model = int(model_str)
        if index_str == "":
            index = 0
        else:
            index = int(index_str)
        value = float(value_str)
        if not (model in dict_model_parameters):
            dict_model_parameters[model] = {param: [(index, value)]}
        else:
            if param in dict_model_parameters[model]:
                dict_model_parameters[model][param].append((index, value))
            else:
                dict_model_parameters[model][param] = [(index, value)]
    for model in dict_model_parameters:
        for param in dict_model_parameters[model]:
            dict_model_parameters[model][param].sort(key=lambda x: x[0])

    return dict_model_parameters


def create_freq_file(res_file_path, freq_file_path):
    """
    create root frequency file
    :param res_file_path: the file path of the res chromEvol file
    :param freq_file_path: root frequencies file to be written to
    :return:
    """
    res_file = open(res_file_path, 'r')
    content = res_file.read()
    res_file.close()
    pattern_section = re.compile("Best chosen model[\s]+[#]+[\s]+(.*?)(\*|$)", re.DOTALL)
    section_best_model = pattern_section.findall(content)[0][0]
    min_chr_num = get_min_chromosome_number(content)
    pattern_freqs = re.compile("F\[([\d]+)\][\s]+=[\s]([\S]+)")
    root_frequencies = pattern_freqs.findall(section_best_model)
    inidices_with_freqs = [(int(state)-min_chr_num, float(freq)) for state,freq in root_frequencies]
    freq_file = open(freq_file_path, 'w')
    for i in range(inidices_with_freqs[0][0]):
        freq_file.write("0\n")
    for i in range(len(inidices_with_freqs)):
        freq_file.write(str(inidices_with_freqs[i][1])+"\n")

    freq_file.close()


def get_min_chromosome_number(res_content, file_handler=None):
    """
    get min chromosome number from results file
    :param res_content: the content of the results file (str)
    :param file_handler: file handler of the results file
    :return: min chromosome number
    """
    if res_content is None:
        res_content = file_handler.read()
    pattern_min_chr_num = re.compile("Min[\s]+allowed[\s]+chromosome[\s]+number[\s]+=[\s]+([\d]+)")
    min_chr_num_sec = pattern_min_chr_num.findall(res_content)
    min_chr_num = int(min_chr_num_sec[0])
    return min_chr_num


def get_max_chromosome_number(res_content, file_handler=None):
    """
    get max chromosome number from results file
    :param res_content: the content of the results file (str)
    :param file_handler: file handler of the results file
    :return: max chromosome number
    """
    if res_content is None:
        res_content = file_handler.read()
    pattern_max_chr_num = re.compile("Max allowed chromosome number[\s]+=[\s]+([\d]+)")
    max_chr_num = int(pattern_max_chr_num.findall(res_content)[0])
    return max_chr_num


def create_rescaled_tree(manipulated_rates, multiplier, nodes_file_path, original_tree_path, new_tree_path, expectation_file_path):
    tree = Tree(original_tree_path, format=1)
    original_tree_length = get_tree_length(tree)
    if manipulated_rates == "all":
        factor = 1/multiplier
    else:
        factor = 1/get_scaling_factor_from_expectations(expectation_file_path, manipulated_rates, multiplier)
    mrca_node = get_mrca(nodes_file_path, tree)
    rescale_tree(mrca_node, factor)
    factor_to_original_length = original_tree_length/get_tree_length(tree)
    rescale_tree(tree, factor_to_original_length)
    print("tree length", str(get_tree_length(tree)))
    tree.write(format=1, outfile=new_tree_path)



def get_mrca(nodes_file_path, tree):
    nodes_file = open(nodes_file_path, 'r')
    content = nodes_file.read()
    nodes_file.close()
    pattern = re.compile("Model[\s]+#[\d]+[\s]+=[\s]+\(([\S]+),([\S]+)\);")
    partitions = pattern.findall(content)
    species1, species2 = partitions[-1]
    ancestor = tree.get_common_ancestor(species1, species2)
    return ancestor


def get_descendent_leaves(node):
    leaves = node.get_leaves()
    leaves_names = [n.name for n in leaves]
    return set(leaves_names)


def get_leaves_under_shifts_tree_node(tree_path):
    tree = Tree(tree_path, format=1)
    pattern = re.compile("([\S]+)-([\d]+)$")
    leaves_names = []
    for node in tree.traverse():
        name_parts = pattern.findall(node.name)[0]
        name, model = name_parts[0], int(name_parts[1])
        if (model == 2) and (node.is_leaf()):
            leaves_names.append(name)
    return set(leaves_names)



def get_scaling_factor_from_expectations(expectation_file_path, manipulated_rates, multiplier):
    if manipulated_rates == "poly":
        rates = ["DUPLICATION", "DEMI-DUPLICATION", "BASE-NUMBER"]
    elif manipulated_rates == "dys":
        rates = ["GAIN", "LOSS"]
    else:
        raise Exception("get_scaling_factor_from_expectations()" + " such rate "+manipulated_rates + " does not exists!")
    expectation_file = open(expectation_file_path, 'r')
    content = expectation_file.read()
    expectation_file.close()
    pattern_section = re.compile("#TOTAL EXPECTATIONS:(.*)", re.DOTALL)
    section = pattern_section.findall(content)[0]
    pattern_exp = re.compile("([\S]+):[\s]+([\S]+)")
    expectations = pattern_exp.findall(section)
    weight = 0
    sum_of_weights = 0
    for rate_type, expecation in expectations:
        if rate_type == "TOMAX":
            continue
        if rate_type in rates:
            weight += float(expecation)
        sum_of_weights += float(expecation)
    factor = (weight/sum_of_weights) * multiplier
    return factor










