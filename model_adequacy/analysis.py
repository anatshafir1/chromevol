from defs import *
from utils import *
from data_processing import *
from scipy import stats
from numpy import var
import subprocess


def calculate_statistics(counts, tree_file, out_file, run_flag, simulated_counts_file=False, emp_var=True, emp_ent=True, emp_pars=True, emp_pars_time=True):
    """
    :param counts: list of counts
    :param tree_file: for parsimony (fitch) and time-parsimony calculations
    :param out_file: output file to where the stats will be printed per simulation
    :param run_flag: if "orig" - calculation on empirical counts, single run of acctran().
    :param empirical_stats: the empirical statistics calculated, to identify which statistics should be calculated
    :param simulated_counts_file: if supplied - stats are calculated on simulations
    :return: list of statistics representing the counts
    """
    v = emp_var and var_calc(counts)  # if None in empirical - no calculation will be done on simulated
    e = emp_ent and entropy_calc(counts)
    p = emp_pars and fitch(tree_file, simulated_counts_file)
    if run_flag:
        a = emp_pars_time and acctran(tree_file)
    else:
        a = 0
    lst_of_stats = [v, e, p, a]
    round_stats = [round(x, 2) if x is not None else None for x in lst_of_stats]
    if run_flag:
        with open(out_file, "w+") as stats:
            stats.write(','.join([str(x) for x in round_stats]))

    return round_stats


def var_calc(counts_vec):
    try:
        v = round(var(counts_vec), 2)  # variance
    except Exception as ex:
        print_error("Unable to calculate variance", ex)
        v = None
    return v


def entropy_calc(counts_vec):
    """
    calculated Shannon's entropy on a vector of chromosome counts
    :param counts_vec: list of counts
    :return: entropy
    """
    try:
        d = {}
        for i in counts_vec:
            d[i] = counts_vec.count(i)
        prob_lst = [x / len(counts_vec) for x in list(d.values())]
        e = stats.entropy(prob_lst)
    except Exception as ex:
        print_error("Unable to calculate entropy", ex)
        e = None
    return e


def fitch(tree_file, c=False):
    """
    calculates the maximum parsimony score following Fitch algorithm, where the states are chromosome counts.
    :param tree_file: tree file
    :param c: simulated counts file
    :return: number of unions (parsimony score)
    """
    try:
        t = Tree(tree_file, format=1)
        score = 0
        if c:
            d = create_counts_hash(c)  # function found in utils.py

        for node in t.traverse("postorder"):
            if not node.is_leaf():  # internal node
                lst = []  # list version
                intersect, union = None, None
                for child in node.get_children():
                    if child.is_leaf():  # if the child is a tip - parse number from tip label
                        if c:  # dictionary exists if the tree is simulated --> take the number from it
                            name = re.search("(.*)\-\d+", child.name)
                            if name:
                                num = {int(d.get(name.group(1)))}
                        else:  # calculation on original counts
                            tmp = re.search("(\d+)", child.name)
                            if tmp:  # there is a number at the tip, and not X
                                num = {int(tmp.group(1))}
                    else:  # if the child is an internal node - take number
                        num = child.name
                    lst.append(num)
                intersect = lst[0] & lst[1]
                union = lst[0] | lst[1]
                if len(intersect) == 0:
                    result = union
                    score += 1
                else:
                    result = intersect
                node.name = result
    except Exception as ex:
        print_error("Unable to calculate parsimony score", ex)
        score = None
    return score


def acctran(tree_file, num_iter=1, sims_dir=0):
    """
    calculates parsimony score.
    :param tree_file: phylogeny to be analyzed
    :param num_iter: number of iterations. 1 is for the empirical data, any other positive number is for the simulations
    :param sims_dir: directory of simulations, if supplied. otherwise calculation is on the original counts
    :return: list of calculated statistic
    """
    try:
        if sims_dir != 0:  # simulated trees, need to add a semicolon at the end
            if tree_file is not None:  # simulations trees, add semicolon to at the end to enable downstream reading
                fix_simulated_tree_file(tree_file)
        else:
            sims_dir = ""
        #os.system(module_loads)
        script = os.path.dirname(os.path.realpath(__file__)) + "/calculate_acctran.R"
        arg = [tree_file, str(num_iter), sims_dir]
        cmd = [R_path, script] + arg
        res = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
        out, err = res.communicate()
        if out == "":
            return None
        else:
            res = [float(i) for i in out.split()]
            if len(res) == 1:
                res = res[0]
    except Exception as ex:
        print_error("Unable to calculate parsimony-time score", ex)
        res = None
    return res


def add_stat_to_lst(all_stats_lst, single_stat_lst):
    """
    replaces in the statistics list the a = 0 with the calculated values from the simulations
    :param all_stats_lst: [[x1,x2,x3,0],[y1,y2,y3,0]....]
    :param single_stat_lst: [a1,a2,....,an]
    :return: [[x1,x2,x3,a1],[y1,y2,y3,a2]....]
    """
    for i in range(len(single_stat_lst)):
        all_stats_lst[i][-1] = single_stat_lst[i]
    return all_stats_lst


#####################################################################################################
#####################################################################################################
#                                                                                                   #
#                                 MODEL ADEQUACY FUNCTIONS                                          #
#                                                                                                   #
#####################################################################################################
#####################################################################################################


def get_list_of_stats_per_sim(sim_dir, user_out_dir, empirical_stats):
    """
    returns a list of three statistics, each calculated on a single simulation
    :param sim_dir: directory of current simulation
    :param user_out_dir: directory with trees
    :param empirical_stats: list of empirical statistics that were calculated
    :return: list of statistics [a,b,c]
    """
    sim_counts_file = sim_dir + "/counts.fasta"
    sim_counts = get_counts(sim_counts_file)
    simulated_counts_statistics = calculate_statistics(sim_counts, user_out_dir + tree_with_counts, None, False, sim_counts_file, empirical_stats[0], empirical_stats[1], empirical_stats[2], empirical_stats[3])
    return simulated_counts_statistics


def create_simulated_stats_distribution(out_dir, user_out_dir, empirical_stats):
    """
    creates list of lists, each sub-list is one statistic over all simulations:
    [[a1,a2,...,an], [b1,b2,...,bn],...]
    to avoid long-time calculation of the pars-time statistic, it is called only once, for all simulations
    :param out_dir: where the simulations are located
    :param user_out_dir: directory with trees
    :param empirical_stats: list of empirical statistics that were calculated
    :return: list of lists of statistics
    """
    simulated_counts_stats_dist = []
    for i in range(get_nsims()):  # for each simulation
        sim_dir = os.path.join(out_dir, str(i))
        simulated_counts_statistics = get_list_of_stats_per_sim(sim_dir, user_out_dir, empirical_stats)
        simulated_counts_stats_dist.append(simulated_counts_statistics)  # creates list of lists by appending
    if empirical_stats[3] is not None:
        simulated_counts_stats_dist = add_stat_to_lst(simulated_counts_stats_dist, acctran("simulatedDataAncestors.tree", get_nsims(), out_dir))
    else:
        simulated_counts_stats_dist = add_stat_to_lst(simulated_counts_stats_dist, [None]*get_nsims())
    return simulated_counts_stats_dist


def remove_none_from_stats(orig_lst, sim_lst):
    statistics_names = ["Variance", "Entropy", "Parsimony", "Time_parsimony"]
    none_ind = [i for i, val in enumerate(orig_lst) if val == None]
    for i in sorted(none_ind, reverse=True):
        del orig_lst[i]
        del statistics_names[i]
    sim_lst = [[x for x in sub_lst if x is not None] for sub_lst in sim_lst]
    return orig_lst, sim_lst, statistics_names


def two_sided_midpoint_PV(dist, s, out_dir, PVs_file):
    """
    two-sided midpoint PV - Hohna et al. 2017
    :param dist: simulated distributions
    :param s: current original statistic to compare
    :param out_dir: where to write results to
    :param PVs_file: the calculated two-sided midpoint p-values
    :return: 0 or 1 - inadequate or adequate
    """
    sim_is_lower = len([x for x in dist if x < s])  # V
    sim_equals = len([x for x in dist if x == s])  # V
    nsims = len(dist)
    Pml = round((sim_is_lower / nsims) + (sim_equals / (2 * nsims)), 3)  # midpoint one-tailed lower p-value
    Pmu = 1 - Pml
    Pmt = 2 * min(Pml, Pmu)
    with open(out_dir + PVs_file, "a") as p_vals:
        p_vals.write(str(Pmt) + "\n")
    if Pmt <= 0.05:
        return 0
    else:
        return 1


def test_adequacy(sim_stats, orig_stats, out_dir, pv_file, stats_dist_file):
    """
    calculates two-sided midpoint p-values over simulated distributions per statistic.
    print the p-value per statistics to a separate file
    :param sim_stats: statistics distributions
    :param orig_stats: original statistics
    :param out_dir: simulations directory
    :param pv_file: output file to print the p-values to
    :param stats_dist_file: output file to print the statistics distribution to
    :return: NA
    """
    adequacy_lst = []
    orig_stats, sim_stats, statistics_names = remove_none_from_stats(orig_stats, sim_stats)
    for i in range(len(orig_stats)):
        sim_stat_dist = [x[i] for x in sim_stats]  # a single statistic distribution
        model = two_sided_midpoint_PV(sim_stat_dist, orig_stats[i], out_dir, pv_file)
        adequacy_lst.append(model)
        handle_distributions(sim_stat_dist, out_dir, stats_dist_file)
    write_final_result(out_dir, adequacy_lst, statistics_names)


def write_final_result(out_dir, adequacy_lst, statistics_names):
    """
    writes the final output of the model adequacy framework.
    :param out_dir: output directory provided by the user
    :param adequacy_lst: vector of 0/1 denoting each statistics adequacy
    :param statistics_names: statistics names
    :return: NA
    """
    suffix = {0:"Inadequate", 1:"Adequate"}
    with open(out_dir + "final_results", "w") as fh:
        for i in range(len(adequacy_lst)):
            fh.write(statistics_names[i] + " is " + suffix[adequacy_lst[i]] + "\n")


def handle_distributions(dist, out_dir, stats_dist_file):
    """
    writes the distributions of the statistics to a file
    :param dist: the distribution of a single statistic
    :param out_dir: the directory in which the final file will be printed to
    :param stats_dist_file: the file to print to
    :return: NA
    """
    with open(os.path.join(out_dir, stats_dist_file), "a") as distribution_file:
        sim_stat_dist = [round(x, 4) for x in dist]
        distribution_file.write(str(sim_stat_dist) + "\n")


def model_adequacy(ma_output_dir, orig_counts_stats, user_out_dir, empirical_stats):
    """
    calculates the distribution of each statistics and compares the original respective statistic to it, to determine adequacy of statistic of inadequacy. Prints output files of the results and targz all simulations to a single file.
    :param ma_output_dir: model adequacy results path
    :param orig_counts_stats: original statistic
    :param user_out_dir: the output directory specified by the user
    :param empirical_stats: list of empirical statistic to be calculated
    :return: NA
    """
    sim_dist = create_simulated_stats_distribution(ma_output_dir, user_out_dir, empirical_stats)
    test_adequacy(sim_dist, orig_counts_stats, user_out_dir, "PVs", "stats_dist_sims")

