from utils import *
from defs import *
import random
from create_job_files import create_job_format_for_all_jobs
from subprocess import Popen
import argparse


def simulate_models(tree_dir, sampling_fractions, model, multipliers, num_of_simulations, out_dir, manipul_rates, queue):
    manipulations_on = {"poly":["demiPloidyR", "dupl", "baseNumR"], "dys":["gain", "loss"],
                        "all":["demiPloidyR", "dupl", "baseNumR", "gain", "loss"]}
    model_results_dir = os.path.join(tree_dir, "Results_" + model)
    original_tree_path = os.path.join(tree_dir, "tree.newick")
    counts_path = os.path.join(tree_dir, "counts.fasta")
    freq_file_path = os.path.join(model_results_dir, "root_freq")
    chromevol_res_path = os.path.join(model_results_dir, "chromEvol.res")
    ml_tree = os.path.join(model_results_dir, "MLAncestralReconstruction.tree")
    expectation_file_path = os.path.join(model_results_dir, "expectations_second_round.txt")
    if not os.path.exists(freq_file_path):
        create_freq_file(chromevol_res_path, freq_file_path)

    for sampling_frac in sampling_fractions:
        sampling_frac_dir = os.path.join(out_dir, str(sampling_frac))
        if not os.path.exists(sampling_frac_dir):
            os.makedirs(sampling_frac_dir)
        nodes_file_path = os.path.join(sampling_frac_dir, "nodes.txt")
        if not os.path.exists(nodes_file_path):
            create_nodes_split_file(nodes_file_path, sampling_frac, original_tree_path)
        parameters = get_params(chromevol_res_path, freq_file_path, nodes_file_path)
        rate_parameters = get_rate_parameters(chromevol_res_path)
        for multiplier in multipliers:
            multiplier_dir_path = os.path.join(sampling_frac_dir, str(multiplier))
            if not os.path.exists(multiplier_dir_path):
                os.makedirs(multiplier_dir_path)
            for manipulation_model in manipul_rates:
                mult_manupulation_dir = os.path.join(multiplier_dir_path, manipulation_model)
                if not os.path.exists(mult_manupulation_dir):
                    os.makedirs(mult_manupulation_dir)

                param_file_path = os.path.join(mult_manupulation_dir, "sim_param.params")
                tree_path = os.path.join(mult_manupulation_dir, "tree.newick")
                create_rescaled_tree(manipulation_model, multiplier, nodes_file_path, original_tree_path, tree_path, expectation_file_path)
                create_simulation_param_file(parameters, multiplier, manipulations_on[manipulation_model], rate_parameters, mult_manupulation_dir,
                                            num_of_simulations, tree_path, counts_path, ml_tree, param_file_path)
                module_python = PYTHON_MODULE_COMMAND
                cmd = module_python
                chromevol_exe = CHROMEVOL_EXE
                cmd += chromevol_exe + " param=" + param_file_path + " > " \
                        + os.path.join(mult_manupulation_dir, "out.txt") + "\n"
                job_path = os.path.join(mult_manupulation_dir,
                                        "CE_sim_"+str(sampling_frac)+"_"+ str(multiplier)+"_"+str(manipulation_model)+".sh")
                job_file = open(job_path, 'w')
                job_content = create_job_format_for_all_jobs(mult_manupulation_dir, "CE_sim", "10gb", queue, 1, cmd)
                job_file.write(job_content)
                job_file.close()
                Popen(["qsub", job_path])


def simulate_homogeneous_model(tree_dir, out_dir, model, num_of_simulations, queue):
    model_results_dir = os.path.join(tree_dir, "Results_" + model)
    tree_path = os.path.join(tree_dir, "tree.newick")
    counts_path = os.path.join(tree_dir, "counts.fasta")
    freq_file_path = os.path.join(model_results_dir, "root_freq")
    chromevol_res_path = os.path.join(model_results_dir, "chromEvol.res")
    ml_tree = os.path.join(model_results_dir, "MLAncestralReconstruction.tree")
    create_freq_file(chromevol_res_path, freq_file_path)
    parameters = get_params(chromevol_res_path, freq_file_path, None, "false")
    rate_parameters = get_rate_parameters(chromevol_res_path)
    param_file_path = os.path.join(out_dir, "sim_param.params")
    create_simulation_param_file(parameters, None, {}, rate_parameters, out_dir,
                                 num_of_simulations, tree_path, counts_path, ml_tree, param_file_path)

    module_python = PYTHON_MODULE_COMMAND
    cmd = module_python
    chromevol_exe = CHROMEVOL_EXE
    cmd += chromevol_exe + " param=" + param_file_path + " > " \
           + os.path.join(out_dir, "out.txt") + "\n"
    job_path = os.path.join(out_dir,"CE_sim_homo.sh")
    job_file = open(job_path, 'w')
    job_content = create_job_format_for_all_jobs(out_dir, "CE_sim", "10gb", queue, 1, cmd)
    job_file.write(job_content)
    job_file.close()
    Popen(["qsub", job_path])


def create_simulation_param_file(parameters, multiplier, manipulated_rates, rate_parameters, multiplier_dir_path,
                                 num_of_simulations, tree_path, counts_path, ml_tree, param_file_path):
    dic_func_conversion = {"dupl":"_duplFunc", "demiPloidyR":"_demiDuplFunc", "gain":"_gainFunc", "loss":"_lossFunc",
                           "baseNumR":"_baseNumRFunc"}
    parameters["_seed"] = str(random.randint(1, 100000))
    parameters["_simulateData"] = "true"
    parameters["_numOfSimulatedData"] = str(num_of_simulations)
    parameters["_fracAllowedFailedSimulations"] = 0.1
    parameters["_numOfRequiredSimulatedData"] = str(num_of_simulations)
    if not (multiplier is None):
        rate_parameters[2] = {}
        for param in rate_parameters[1]:
            rate_parameters[2][param] = []
            for index, value in rate_parameters[1][param]:
                if param in manipulated_rates:
                    func = parameters[dic_func_conversion[param]]
                    if (func == "LINEAR") or (func == "CONST"):
                        rate_parameters[2][param].append((index, multiplier*value))
                    elif (func == "EXP") and (index == 0):
                        rate_parameters[2][param].append((index, multiplier*value))
                    else:
                        rate_parameters[2][param].append((index, value))
                else:
                    rate_parameters[2][param].append((index, value))
    if parameters["_baseNumRFunc"] != "IGNORE":
        parameters["_maxBaseNumTransition"] = max(test_max_on_tree(counts_path, ml_tree), int(max(rate_parameters[1]["baseNum"][0][1], rate_parameters[2]["baseNum"][0][1]))+1)
    param_file_content = set_paths(multiplier_dir_path, tree_path, multiplier_dir_path)
    param_file_content += write_general_params(parameters)
    param_file_content += write_rate_parameters(rate_parameters)
    param_file = open(param_file_path, 'w')
    param_file.write(param_file_content)
    param_file.close()


def create_nodes_split_file(nodes_file_path, sampling_frac, tree_path):
    tree = Tree(tree_path, format=1)
    dict_size_node = search_node_sizes(tree)
    subclade_size = find_the_closest_size(dict_size_node, sampling_frac, len(tree))
    nodes = dict_size_node[subclade_size]
    node = random.sample(nodes, 1)[0]
    mrcas_defs = []
    mrcas_defs.append(get_farthest_leaves(tree))
    mrcas_defs.append(get_farthest_leaves(node))
    num_of_models = 2
    node_splits = ""
    for i in range(1, num_of_models+1):
        node_splits += write_nodes_splits_line(mrcas_defs[i-1], i)
    nodes_file = open(nodes_file_path, 'w')
    nodes_file.write(node_splits)
    nodes_file.close()


def write_nodes_splits_line(farthest_leaves, model_num):
    line = "Model #"+str(model_num)+ " = ("+farthest_leaves[0]+","+farthest_leaves[1]+");\n"
    return line


def get_farthest_leaves(tree):
    sons = tree.get_children()
    farthest_leaves = []
    for i in range(len(sons)):
        son = sons[i]
        if son.is_leaf():
            farthest_leaves.append(son.name)
        else:
            son_leaves = son.get_leaves()
            farthest_leaves.append(son_leaves[0].name)

    return farthest_leaves


def find_the_closest_size(clade_sizes, sampling_frac, full_tree_size):
    min_tree_size = 400
    required_size = int(sampling_frac*full_tree_size)
    for clade_size in clade_sizes:
        if abs(required_size-clade_size) < abs(min_tree_size-required_size):
            min_tree_size = clade_size
    return min_tree_size

def simulate_with_empirical_params(tree_dir, sampling_fractions, num_of_simulations,
                                   out_dir, model_results_dir, queue):

    original_tree_path = os.path.join(tree_dir, "tree.newick")
    counts_path = os.path.join(tree_dir, "counts.fasta")
    freq_file_path = os.path.join(model_results_dir, "root_freq")
    chromevol_res_path = os.path.join(model_results_dir, "chromEvol.res")
    ml_tree = os.path.join(model_results_dir, "MLAncestralReconstruction.tree")
    if not os.path.exists(freq_file_path):
        create_freq_file(chromevol_res_path, freq_file_path)
    for sampling_frac in sampling_fractions:
        sampling_frac_dir = os.path.join(out_dir, str(sampling_frac))
        if not os.path.exists(sampling_frac_dir):
            os.makedirs(sampling_frac_dir)
        nodes_file_path = os.path.join(sampling_frac_dir, "nodes.txt")
        if not os.path.exists(nodes_file_path):
            create_nodes_split_file(nodes_file_path, sampling_frac, original_tree_path)
        parameters = get_params(chromevol_res_path, freq_file_path, nodes_file_path)
        rate_parameters = get_rate_parameters(chromevol_res_path)
        empirical_dir_path = os.path.join(sampling_frac_dir, "empirical_params")
        if not os.path.exists(empirical_dir_path):
            os.makedirs(empirical_dir_path)
        param_file_path = os.path.join(empirical_dir_path, "sim_param.params")
        create_simulation_param_file(parameters, None, None, rate_parameters, empirical_dir_path,
                                        num_of_simulations, original_tree_path, counts_path, ml_tree, param_file_path)
        module_python = PYTHON_MODULE_COMMAND
        cmd = module_python
        chromevol_exe = CHROMEVOL_EXE
        cmd += chromevol_exe + " param=" + param_file_path + " > " \
                + os.path.join(empirical_dir_path, "out.txt") + "\n"
        job_path = os.path.join(empirical_dir_path,
                                    "CE_sim_"+str(sampling_frac)+"_emp.sh")
        job_file = open(job_path, 'w')
        job_content = create_job_format_for_all_jobs(empirical_dir_path, "CE_sim", "10gb", queue, 1, cmd)
        job_file.write(job_content)
        job_file.close()
        Popen(["qsub", job_path])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='simulates homogeneous and heterogeneous models')
    parser.add_argument('--tree_dir', '-t', help='tree directory')
    parser.add_argument('--num_of_simulations', '-n', type=int, help='number of simulations')
    parser.add_argument('--sampling_fractions', '-s', type=float, nargs='+', help='clade size fractions')
    parser.add_argument('--model', '-m', help='adequate model')
    parser.add_argument('--multipliers', '-k', type=float, nargs = '+', help="rate multipliers")
    parser.add_argument('--out_dir', '-o', help='simulations directory')
    parser.add_argument('--manipulated_rate', '-r', nargs='+', help='manipulted rates')
    parser.add_argument('--empirical_hetero_params', '-e', default=None,
                        help="chrom res file to get empirical heterogeneous model data")
    parser.add_argument('--queue', '-q', help="queue")

    # parse arguments
    args = parser.parse_args()
    tree_dir = args.tree_dir
    sampling_fractions = args.sampling_fractions
    model = args.model
    multipliers = args.multipliers
    num_of_simulations = args.num_of_simulations
    out_dir = args.out_dir
    empirical_hetero_params_dir = args.empirical_hetero_params
    manipulated_rate = args.manipulated_rate
    queue = args.queue
    if 1 in multipliers:
        simulate_homogeneous_model(tree_dir, out_dir, model, num_of_simulations,queue)
    elif (0 in multipliers) and (not (empirical_hetero_params_dir is None)):
        simulate_with_empirical_params(tree_dir, sampling_fractions, num_of_simulations, out_dir, empirical_hetero_params_dir, queue)
    else:
        simulate_models(tree_dir, sampling_fractions, model, multipliers, num_of_simulations, out_dir, manipulated_rate, queue)


