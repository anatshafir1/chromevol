import os, re, argparse
from ete3 import Tree

#  fixed names
log = "/log.txt"
mlAncTree = "/MLAncestralReconstruction.tree"
posterior_tree = "/MarginalAncestralReconstruction.tree"
exp_tree = "/exp.tree"
expectation_file = "/expectations.txt"
anc_prob = "/ancestorsProbs.txt"
tree_with_counts = "/tree_with_counts.tree"
tree_wo_counts = "/tree_wo_counts.tree"
root_freq_filename = "/root_freq"
sim_control = "/param_sim"
adequacy_vec = "/adequacy_vec"
#module_loads = "module load R/4.1.0"
R_path = "Rscript"  # PLACE YOUR R EXECUTABLE FULL PATH HERE
chromevol_path = "/groups/itay_mayrose/anatshafir1/phd_project/scripts/chromevol_trait/chromevol/ChromEvol/chromEvol"  # PLACE YOUR CHROMEVOL EXECUTABLE FULL PATH HERE
max_num_of_functions = 5
nsims = None
frac_allowed_failed_simulations = 0.01


def get_nsims():
    return nsims


def set_nsims(n):
    global nsims
    nsims = n


def get_arguments():
    parser = argparse.ArgumentParser(description='Pipeline to test model adequacy of a chromEvol model')
    parser.add_argument('--counts', '-c', help='Counts file', required=True)
    parser.add_argument('--tree', '-t', help='Tree file', required=True)
    parser.add_argument('--results_file', '-r', help='chromEvol results file', required=True)
    parser.add_argument('--user_output_dir', '-out', help='Output directory', required=True)
    parser.add_argument('--sims_per_tree', '-n', help='Number of simulations for the adequacy test', required=False, default=1000)
    parser.add_argument('--nodes_file_path', '-s', help='file of model assignment if the model is heterogeneous', default=None)
    parser.add_argument('--shifts_tree_path', '-m', help = 'tree shifts file', default=None)

    # parse arguments
    args = parser.parse_args()
    counts_file = args.counts
    tree_file = args.tree
    results_file = args.results_file
    user_output_dir = args.user_output_dir
    nodes_file_path = args.nodes_file_path
    shifts_tree_path = args.shifts_tree_path
    set_nsims(int(args.sims_per_tree))

    return counts_file, tree_file, results_file, user_output_dir, nodes_file_path, shifts_tree_path
