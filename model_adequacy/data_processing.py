from defs import *
from utils import *
import random


def get_counts(filename):
    """
    reads the .counts_edit file and extracts the counts
    :param filename: counts file (original or simulated)
    :return: list of counts
    """
    with open(filename, "r") as counts_file:
        counts = []
        for line in counts_file:
            line = line.strip()
            if line.startswith('>'):  # taxa name
                continue
            else:
                if line=="x" or line=="X":  # discard counts with x
                    continue
                counts.append(int(line))
    return counts


def match_counts_to_tree(tree_file, out_dir):
    """
    Matches tree file to counts, in case of missing taxa or missing counts.
    :param tree_file: mlAncTree (containing counts in tips labels) in NEWICK format
    :param out_dir: where the processed tree and counts should be written to.
    :return:(1) tree_wo_counts without X taxa and without counts in their tip names
            (2) tree_with_counts without X taxa and with counts in the tips
    """
    t = Tree(tree_file, format=1)
    t = prune_x_from_tree(t, out_dir)
    produce_tree_without_counts(t, out_dir)
    remove_internal_labels(out_dir)


def prune_x_from_tree(t, out_dir):
    """
    Prune taxa with Xs as counts
    :param t: tree object
    :param out_dir: output directory
    :return: pruned tree (unpruned if not needed)
    """
    tips_to_prune = []
    all_tips = []
    for leaf in t:
        all_tips.append(leaf.name)
        name_with_x = re.search(".*\-X", leaf.name)
        if name_with_x:
            tips_to_prune.append(leaf.name)
    t.prune(list(set(all_tips) - set(tips_to_prune)))
    t.write(format=1, outfile=out_dir + tree_with_counts)
    return t


def produce_tree_without_counts(t, out_dir):
    """
    trim the counts digits from the tip labels
    :param t: tree object
    :param out_dir: output directory
    :return: tree without counts in tips labels
    """
    for leaf in t:
        name = re.search("(.*)\-[\d]", leaf.name)
        leaf.name = name.group(1)
    t.write(format=1, outfile=out_dir + tree_wo_counts)
    return t


def remove_internal_labels(out_dir):
    """
    remove internal nodes labels and add ":-1" at the end of the tree to make it rooted
    :param out_dir:
    :return: NA
    """
    with open(out_dir + tree_wo_counts, "r") as tree:
        t = tree.readline()
    # t = re.sub(";", ":-1;", t)  # remove ";" and add ":-1;"
    t = re.sub(r"N\d+-\d+", r"", t)  # replace all _N\d+-\d+_ with ""
    with open(out_dir + tree_wo_counts, "w+") as tree:
        tree.write(t)


def create_freq_file(content, freq_file_path):
    """
    create root frequency file
    :param content: the content of the res file
    :param freq_file_path: root frequencies file to be written to
    :return:
    """
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
    # The following segment of code is removed because chromEvol should be able to add zeros
    # whenever needed
    # for i in range(inidices_with_freqs[-1][0]+1, global_max_chromosome_value-min_chr_num+1):
    #     freq_file.write("0\n")
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


def get_rate_parameters(file_content, freq_file):
    """
    get the rate parameters from the results file
    :param file_content: file content of results file
    :param freq_file: root frequencies file to be written to
    :return: list of rate parameters with their values (in the format that they should appear)
    """
    pattern_section = re.compile("Best chosen model[\s]+[#]+[\s]+(.*?)(\*|$)", re.DOTALL)
    pattern_params = re.compile("Chromosome\.([\D]+)([\d]*)_([\d]+)[\s]+=[\s]+([\S]+)")
    section = pattern_section.findall(file_content)[0][0]
    params_index_models = pattern_params.findall(section)
    dict_res_param_name_input_names = {"dupl": "dupl", "gain": "gain", "loss": "loss", "demi": "demiPloidyR",
                                       "baseNum": "baseNum", "baseNumR": "baseNumR"}
    dict_model_parameters = {} # model: {param:[(index, value)}
    for rate_param, index_str, model_str, value_str in params_index_models:
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
    final_params = []
    category_index = 1
    for model in range(1, len(dict_model_parameters)+1):
        for param in dict_model_parameters[model]:
            dict_model_parameters[model][param].sort(key=lambda x: x[0])
            param_setting = "_"+param+"_"+str(model)+ " = "+ str(category_index)+";"
            for i in range(len(dict_model_parameters[model][param])):
                param_setting += str(dict_model_parameters[model][param][i][1])
                if i != len(dict_model_parameters[model][param])-1:
                    param_setting += ","
            final_params.append(param_setting+"\n")
            category_index += 1
    final_params.append("_numOfModels = "+ str(len(dict_model_parameters))+"\n")
    final_params.append("_maxNumOfModels = "+ str(len(dict_model_parameters))+"\n")
    create_freq_file(file_content, freq_file)
    return final_params


def get_tree_scaling_factor(content):
    """
    extracts the tree scaling factor that should be used for the simulations
    :param content: chromEvol results file (chromEvol.res(=)
    :return: the tree scaling factor
    """
    pattern = re.compile("Tree scaling factor is:[\s]+([\S]+)")
    scaling_factor = float(pattern.findall(content)[0])
    return scaling_factor


def is_heterogeneous(res_content):
    pattern = re.compile("Number of models in the best model[\s]+=[\s]+([\d]+)")
    if int(pattern.findall(res_content)[0]) > 1:
        return "true"
    return "false"


def create_nodes_splits_file(res_content, nodes_path_file, tree_file_path):
    """
    extracts the partitions of the shifts regimes of the heterogeneous model, and creates an input file,
    which specifies the nodes partitions to the corresponding models
    :param res_content: chromEvol results file content(chromEvol.res(=)
    :param nodes_path_file: file which dictates the node splits in case of heterogeneous models
    :param tree_file_path: tree file with shifts.
    :return: No return
    """
    pattern_section = re.compile("Best chosen model[\s]+[#]+[\s]+(.*?)(\*|$)", re.DOTALL)
    section_best_model = pattern_section.findall(res_content)[0][0]
    pattern_models = re.compile("#[\s]+Model[\s]+\$([\d]+)[\s]+=[\s]+([\S]+)")
    models_with_nodes = pattern_models.findall(section_best_model)
    tree = Tree(tree_file_path, 1)
    file = open(nodes_path_file, 'w')
    number_of_indicated_nodes = 0
    for node in tree.traverse():
        node_name = node.name
        for model, subtree_root in models_with_nodes:
            if subtree_root+"-" in node_name:
                number_of_indicated_nodes += 1
                line_to_write = "Model #"+model+ " = ("
                sons = node.get_children()
                for i in range(len(sons)):
                    son = sons[i]
                    if son.is_leaf():
                        line_to_write += re.sub(r'([\S]+)(-[\d]+)', r'\1', son.name)
                    else:
                        son_leaves = son.get_leaves()
                        line_to_write += re.sub(r'([\S]+)(-[\d]+)', r'\1', son_leaves[0].name)
                    if i < len(sons)-1:
                        line_to_write += ","
                    else:
                        line_to_write += ");"
                file.write(line_to_write+"\n")
                break

        if number_of_indicated_nodes == len(models_with_nodes):
            break
    file.close()


def get_functions_settings(res_content):
    """
    extracts the functions set for each transition rate
    :param res_content: chromEvol results file content (chromEvol.res(=)
    :return: list of functions settings
    """
    lst_of_functions = []
    pattern_section = re.compile("Assigned functions for each rate parameter:[\s]+(.*)?#", re.DOTALL)
    pattern_funtions = re.compile("([\S]+):[\s]+([\S]+)")
    section = pattern_section.findall(res_content)[0]
    functions_raw = pattern_funtions.findall(section)
    dict_functions = {"baseNumR": "_baseNumRFunc", "dupl": "_duplFunc",
                      "demi": "_demiDuplFunc", "gain": "_gainFunc", "loss": "_lossFunc"}
    for i in range(min(max_num_of_functions, len(functions_raw))):
        transition_type, function = functions_raw[i]
        if transition_type in dict_functions:
            lst_of_functions.append(dict_functions[transition_type] +" = "+ function+"\n")
    return lst_of_functions


def get_params(res_path, freq_file, nodes_path_file, tree_shifts_path):
    """
    parses results file, writes a root frequency file and creates parameters dictionary
    :param res_path: chromEvol results file (chromEvol.res(=)
    :param freq_file: root frequencies file to be written to
    :param nodes_path_file: file which dictates the node splits in case of heterogeneous models
    :param tree_shifts_path: tree file with shifts
    :return: parameters list
    """
    list_of_params = []
    res_file = open(res_path, "r")
    res_content = res_file.read()
    res_file.close()
    rate_parameters = get_rate_parameters(res_content, freq_file)
    min_chr_num = get_min_chromosome_number(res_content)
    max_chr_num = get_max_chromosome_number(res_content)
    tree_scaling_factor = get_tree_scaling_factor(res_content)
    heterogeneous_model = is_heterogeneous(res_content)
    functions_params_lst = get_functions_settings(res_content)
    list_of_params.append("_fixedFrequenciesFilePath = "+ freq_file.replace("//", "/")+"\n")
    list_of_params.append("_minChrNum = "+str(min_chr_num)+"\n")
    list_of_params.append("_maxChrInferred = "+str(max_chr_num)+"\n")
    list_of_params.append("_branchMul = 1.0\n")#+str(tree_scaling_factor)+"\n")
    list_of_params.append("_heterogeneousModel = "+heterogeneous_model+"\n")
    if (not (nodes_path_file is None)) and (heterogeneous_model == "true"):
        create_nodes_splits_file(res_content, nodes_path_file, tree_shifts_path)
        list_of_params.append("_nodeIdsFilePath = "+str(nodes_path_file)+"\n")
    else:
        list_of_params.append("_nodeIdsFilePath = none\n")
    list_of_params.extend(functions_params_lst)
    list_of_params.extend(rate_parameters)
    return list_of_params


#####################################################################################################
#####################################################################################################
#                                                                                                   #
#                                 SIMULATIONS FUNCTIONS                                             #
#                                                                                                   #
#####################################################################################################
#####################################################################################################

#
# def initialize_defaults(ma_output_dir, max_for_sim, tree_path, freq_file):
#     """
#     initialize several parameters in parameters dictionary, to be printed in the simulations parameters file
#     :param ma_output_dir: where the simulations will be written to
#     :param max_for_sim: current maximum allowed and initial maximum computed
#     :param tree_path: phylogeny
#     :param freq_file: frequency file path
#     :return: parameters dictionary initialized with initial values
#     """
#     params_list = []
#     params_list.append("_resultsPathDir = "+ma_output_dir+"\n")
#     params_list.append("_treeFile = " + tree_path + "\n")
#     d = dict()
#     d["_resultsPathDir"] = ma_output_dir
#     d["_treeFile"] = tree_path
#     d["_fixedFrequenciesFilePath"] = freq_file
#
#     return d


def create_control_file(control_file, ma_output_dir, max_for_sim, tree_path, freq_file, params_lst,
                        orig_counts):
    """
    create chromEvol parameters file for simulations, based on the results file
    :param control_file: name of file
    :param ma_output_dir: where the simulations will be written to
    :param max_for_sim: current maximum allowed and initial maximum computed
    :param tree_path: phylogeny
    :param freq_file: frequency file path
    :param params_dict: parameters dictionary
    :param orig_counts: original counts
    :return: NA
    """
    params_lst_for_sim = params_lst.copy()
    params_lst_for_sim.append("_simulateData = true\n")
    seed = str(100)#str(random.randint(1, 100000))
    # if max_for_sim == 360:
    #     seed = str(4014)
    params_lst_for_sim.append("_seed = " + seed + "\n")
    params_dict = create_dict_of_parameters(params_lst_for_sim)
    params_lst_for_sim.append("_maxChrNum = " + str(max_for_sim)+"\n")
    if params_dict["_baseNumRFunc"] != "IGNORE":
        # there will either be _maxBaseNumTransition after the second run. If not then the range is appropriate.
        if not ("_maxBaseNumTransition" in params_dict):
            params_dict["_maxBaseNumTransition"] = max(max(orig_counts) - min(orig_counts),
                                                    get_max_base_number(params_lst_for_sim))
            params_lst_for_sim.append("_maxBaseNumTransition = " + str(params_dict["_maxBaseNumTransition"]) +"\n")
    with open(control_file, "w+") as cf:
        paths_params(cf, ma_output_dir, ma_output_dir, tree_path)
        for param_entry in params_lst_for_sim:
            cf.write(param_entry)
        cf.write("_numOfSimulatedData = " + str(get_nsims())+"\n")
        cf.write("_numOfRequiredSimulatedData = " + str(get_nsims())+"\n")
        cf.write("_fracAllowedFailedSimulations = "+ str(frac_allowed_failed_simulations)+"\n")


def get_initial_max_allowed(orig_counts, results_file):
    """
    get the initial maximum chromosome number allowed from the results file and the current maximum allowed to be used in the first iteration
    :param orig_counts: original counts
    :param results_file: chromEvol results file
    :return: current maximum allowed and initial maximum computed
    """
    real_max = max(orig_counts)
    init_max_for_sim = max(real_max, min(real_max*10, 200))
    max_allowed = extract_line_from_file(results_file, "Max allowed chromosome number", True, True)
    return [max_allowed, init_max_for_sim]


def update_max_for_sim(m, init_max, max_allowed):
    """
    updates the current maximal number allowed by a factor
    :param m: multiplication factor
    :param init_max: initial maximum chromosome number allowed
    :param max_allowed: previous maximum chromosome number allowed
    :return: the current updated maximum for next iteration
    """
    max_for_sim = 200 * m + init_max
    if max_for_sim < max_allowed:
        max_for_sim = max_allowed
    return max_for_sim


def is_reached_upper_bound(ma_output_dir, working_dir, max_for_sim, m):
    """
    checks if simulations reached upper bound of maximal chromosome number allowed
    :param ma_output_dir: where the simulations will be written to
    :param working_dir: output directory supplied by the user
    :param max_for_sim: current allowed maximal number
    :param m: multiplication factor
    :return: did simulations reach upper bound? True or False
    """
    bound_reached = False
    bad_simulations = 0
    max_bad_simulations = frac_allowed_failed_simulations * get_nsims()
    for i in range(0, get_nsims()):
        if not os.path.exists(os.path.join(ma_output_dir, str(i))):
            bound_reached = True
            break

        sim_events_file = ma_output_dir + "/"+str(i)+"/simulatedEvolutionPaths.txt"
        if not os.path.exists(sim_events_file):
            bound_reached = True
            break
        bad = is_reached_max_bound(sim_events_file, max_for_sim)
        if bad:
            bad_simulations += 1
            if bad_simulations > max_bad_simulations:
                bound_reached = True
                break

    if bound_reached:
        with open(working_dir + "/increasing_max_chr.txt", "w") as fh:
            fh.write("Iteration number " + str(m+1) + ", max number is currently " + str(max_for_sim))
        return True
    return False


def is_reached_max_bound(sim_events_file, max_for_sim):
    events_file = open(sim_events_file, 'r')
    content = events_file.read()
    events_file.close()
    pattern_from = re.compile("from state:[\s]+([\d]+)")
    pattern_to = re.compile("to state[\s]+=[\s]+([\d]+)")
    lst_of_states = []
    from_states = pattern_from.findall(content)
    to_states = pattern_to.findall(content)
    lst_of_states.extend(from_states)
    lst_of_states.extend(to_states)
    set_of_states = set([int(i) for i in lst_of_states])
    if len(set_of_states) == 0:
        return False
    max_state = max(set_of_states)
    return max_state == max_for_sim


def update_max_base_transition(second_run_param_file_path, parameters, chr_range):
    if not os.path.exists(second_run_param_file_path):
        parameters.append("_maxBaseNumTransition  = " + str(chr_range) + "\n")
        return
    file = open(second_run_param_file_path, 'r')
    content = file.read()
    file.close()
    pattern = re.compile("(_maxBaseNumTransition)[\s]+=[\s]+([\S]+)")
    param, value = pattern.findall(content)[0]
    parameters.append(param + " = "+ value+"\n")


def run_simulations(results_file, ma_output_dir, orig_counts, tree_path, freq_file, params_dict, user_out_dir):
    """
    creates parameters file for simulations.
    checks upper bound of maximal chromosome number for simulations and updates accordingly.
    :param results_file: chromEvol results file
    :param ma_output_dir: where the simulations will be written to
    :param orig_counts: original counts
    :param tree_path: phylogeny
    :param freq_file: frequency file path
    :param params_dict: parameters dictionary
    :param user_out_dir: user output directory
    :return: NA
    """
    if not os.path.exists(ma_output_dir):
        os.makedirs(ma_output_dir)

    [max_allowed, init_max_for_sim] = get_initial_max_allowed(orig_counts, results_file)
    for mult in range(0, 5):  # mult is the factor increasing _maxChrNumForSimulations
        max_for_sim = update_max_for_sim(mult, init_max_for_sim, max_allowed)
        create_control_file(ma_output_dir + sim_control, ma_output_dir, max_for_sim, tree_path, freq_file,
                            params_dict, orig_counts)
        control_file_path = ma_output_dir+sim_control
        control_file_path = control_file_path.replace("//", "/")
        os.system(chromevol_path + ' param=' + control_file_path)
        # should be true if the bound is not reached!

        tmp = is_reached_upper_bound(ma_output_dir, user_out_dir, max_for_sim, mult)
        if not tmp:  # did not hit upper bound, no need to increase max_for_sim again
            break
        else:
            for i in range(get_nsims()):
                sim_dir = os.path.join(ma_output_dir, str(i))
                if not os.path.exists(sim_dir):
                    continue
                files = os.listdir(sim_dir)
                for file in files:
                    os.remove(os.path.join(sim_dir, file))


#####################################################################################################
#####################################################################################################
#                                                                                                   #
#                                 BASE NUMBER SECOND RUN FUNCTIONS                                  #
#                                                                                                   #
#####################################################################################################
#####################################################################################################


def test_max_on_tree(base_numbers, counts_file, tree_file):
    """
    calculates the maximal transition on the inferred tree
    :param base_numbers: base number as parsed from the results file
    :param counts_file: counts input file
    :param tree_file: ml inferred tree
    :return: if the maximal transitions if equal or larger than the original range -
                return 0
            otherwise
                return a list of the new base number and the max base transition
    """
    counts = get_counts(counts_file)
    counts_range = range_of_lst(counts)
    max_base_on_tree = get_max_transition(tree_file)
    if max_base_on_tree >= counts_range:
        return 0
    max_base = max(max_base_on_tree, 3)
    base_num = min(max_base, min(base_numbers))
    return [base_num, max_base]


def create_control_file_second_run(control_file, out_dir, counts_file, tree_file, params_list, mb):
    """
    creates control file for second chromEvol run (for base_num models only)
    :param control_file: parameters file for chromEvol run
    :param out_dir: output directory of the second chromEvol run
    :param counts_file: counts file for the second chromEvol run
    :param tree_file: tree file for the second chromEvol run
    :param params_list: list of model specific parameters
    :param mb: max base transition
    :return: NA
    """
    with open(control_file, "w+") as cf:
        paths_params(cf, out_dir, counts_file, tree_file)
        fixed_params(cf, params_list)
        write_model_specific_params(cf, params_list)
        cf.write("_maxBaseNumTransition = " + str(mb) + "\n")
        cf.write("_useMaxBaseTransitonNumForOpt = true\n")


def write_model_specific_params(cf, params_list):
    """
    Adds model specific parameters
    :param cf: file handler of control file
    :param params_list: list of model specific parameters
    :return: NA
    """
    pattern_function_line = re.compile("_[\S]+Func[\s]+=[\s]+[\S]+")
    pattern_function = re.compile("(_[\S]+)Func[\s]+=[\s]+([\S]+)")
    pattern_rates = re.compile("(_[\S]+)(_[\d]+[\s]+=[\s]+[\d]+);[\S]+")
    dict_functions = {}
    cf.write("_maxParsimonyBound = true\n")
    cf.write("_baseNumOptimizationMethod = Ranges\n")
    for param_setting in params_list:
        function_sec = pattern_function.findall(param_setting)
        if function_sec != []:
            dict_functions[function_sec[0][0]] = function_sec[0][1]
            function_line = pattern_function_line.findall(param_setting)[0]
            cf.write(function_line+"\n")
            continue
        rates_section = pattern_rates.findall(param_setting)

        if rates_section != []:
            cf.write(param_setting + "\n")
            #set_default_params(cf, param_setting, pattern_rates, dict_functions)
            continue
        if "_numOfModels" in param_setting:
            cf.write(param_setting+"\n")
        elif "_maxNumOfModels" in param_setting:
            cf.write(param_setting+"\n")
        elif "_minCladeSize" in param_setting:
            cf.write(param_setting+"\n")
        elif "_heterogeneousModel" in param_setting:
            cf.write(param_setting+"\n")
        elif "_nodeIdsFilePath" in param_setting:
            cf.write(param_setting+"\n")
        else:
            continue


def set_default_params(cf, param_line, pattern_rates, dict_functions):
    rate = pattern_rates.findall(param_line)[0][0]

    if rate == "_baseNum":
        cf.write(param_line)
    else:

        if rate == "_demiPloidyR":
            rate = "_demiDupl"
        function = dict_functions[rate]
        if function == "CONST":
            rate_to_repl = r'\1\2;1'
        elif function == "LINEAR":
            rate_to_repl = r'\1\2;1,0.1'
        elif function == "EXP":
            rate_to_repl = r'\1\2;1,0.01'
        elif function == "LINEAR_BD":
            rate_to_repl = r'\1\2;0.1'
        else:
            raise Exception("Unknown function!")
        replacement_with_default = re.sub(pattern_rates, rate_to_repl, param_line)
        cf.write(replacement_with_default+"\n")


def paths_params(control_file, out_dir, counts_file, tree_file):
    """
    Adds paths to the control file for the second chromEvol run
    :param control_file: file handler of control file
    :param out_dir: output directory for the second run
    :param counts_file: counts file for the second run
    :param tree_file: tree file for the second run
    :return: NA
    """
    control_file.write("_resultsPathDir = " + out_dir.replace("//","/") + "\n")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    control_file.write("_dataFile = " + counts_file.replace("//","/") + "\n")
    control_file.write("_treeFile = " + tree_file.replace("//","/") + "\n")


def fixed_params(control_file, params_list):
    """
    Adds fixed parameters to the control file for the second chromEvol run
    :param control_file: file handler of control file
    :param params_list: list of inferred parameters, to get the min chr number as used in the inference
    :return: NA
    """
    pattern_min_chr_num = re.compile("_minChrNum[\s]+=[\s]+([\d]+)")
    min_chr_num = pattern_min_chr_num.findall("\n".join(params_list))[0]
    control_file.write("_maxChrNum = -10\n")
    control_file.write("_minChrNum = " + min_chr_num+"\n")
    control_file.write("_branchMul = 999\n")
    control_file.write("_tolParamOptimization = 0.1\n")
    control_file.write("_optimizePointsNum = 1\n")
    control_file.write("_optimizeIterNum = 1\n")
    control_file.write("_optimizeIterNumNextRounds = 0,1\n")
    control_file.write("_optimizePointsNumNextRounds = 2,1\n")
    control_file.write("_seed = 1\n")
    control_file.write("_optimizationMethod = Brent\n")
    control_file.write("_backwardPhase = false\n")


def second_run(params_list, results_path, counts_file, tree_file):
    """
    checks if a second chromEvol run is needed based on the maximal transition on the phylogeny
    :param params_list: parameters list
    :param results_path: where to print results to
    :param counts_file: original coutns file
    :param tree_file: phylogeny
    :return: NA
    """
    base_numbers = get_base_numbers(params_list)
    res = test_max_on_tree(base_numbers, counts_file, results_path + tree_with_counts)
    if isinstance(res, list):  # need to re-run chromEvol with updated parameters
        create_control_file_second_run(results_path + "/second_run.params", results_path, counts_file, tree_file,
                                       params_list, max(max(base_numbers)+1,res[1]))
        param_file = results_path + "/second_run.params"
        param_file = param_file.replace("//", "/")
        os.system(chromevol_path + ' param=' + param_file)
    open(results_path + "/second_run_tested", 'a').close()
