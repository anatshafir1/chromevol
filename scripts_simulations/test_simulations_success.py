import re
import os
import pandas as pd
import argparse
import shutil
import random


def check_simulations(sim_dir, num_of_simulations, stats_file_path):
    absolute_max_state_evol = 0
    absolute_max_state_counts = 0
    dic_stats = {"simulation":[], "max_state_counts":[], "max_state_evol":[]}
    for i in range(num_of_simulations):
        sim_dir_i = os.path.join(sim_dir, str(i))
        if not os.path.exists(sim_dir_i):
            continue
        counts_path = os.path.join(sim_dir_i, "counts.fasta")
        evol_path = os.path.join(sim_dir_i, "simulatedEvolutionPaths.txt")
        max_state_evol = get_max_state_transition(evol_path)
        max_state_counts = get_max_state_in_counts(counts_path)
        if max_state_counts > absolute_max_state_counts:
            absolute_max_state_counts = max_state_counts
        if max_state_evol > absolute_max_state_evol:
            absolute_max_state_evol = max_state_evol
        if (max_state_evol == 200) or (max_state_counts == 200):
            shutil.rmtree(sim_dir_i)
            print("removing dir", str(i))
        dic_stats["simulation"].append(i)
        dic_stats["max_state_evol"].append(max_state_evol)
        dic_stats["max_state_counts"].append(max_state_counts)
    df = pd.DataFrame(dic_stats)
    df.to_csv(stats_file_path)
    print("max state in evolution path:", str(absolute_max_state_evol))
    print("max state in counts:", str(absolute_max_state_counts))


def get_max_state_in_counts(counts_path):
    file = open(counts_path, 'r')
    content = file.read()
    file.close()
    pattern = re.compile(">[\S]+[\s]+([\d]+)", re.MULTILINE)
    counts_str = pattern.findall(content)
    counts = [int(count_str) for count_str in counts_str]
    return max(counts)


def get_max_state_transition(evol_path):
    file = open(evol_path, 'r')
    content = file.read()
    file.close()
    pattern_state_to_state = re.compile("from state:[\s]+([\d]+)[\s]+t[\s]+=[\s]+[\S]+[\s]+to state[\s]+=[\s]+([\d]+)")
    states = pattern_state_to_state.findall(content)
    max_state = 0
    for state_from, state_to in states:
        if int(state_from) > max_state:
            max_state = int(state_from)
        elif int(state_to) > max_state:
            max_state = int(state_to)
        else:
            continue
    return max_state


def pick_successful_simualtions(sim_dir, num_of_simulations):
    dirs = os.listdir(sim_dir)
    lst_of_dirs = []
    for file in dirs:
        sim_i_dir = os.path.join(sim_dir, file)
        if not os.path.isdir(sim_i_dir):
            continue
        lst_of_dirs.append(int(file))
    sample_of_successful = random.sample(lst_of_dirs, num_of_simulations)
    set_of_dirs_to_remove = set(lst_of_dirs)-set(sample_of_successful)
    for dir_to_remove in set_of_dirs_to_remove:
        full_dir_to_remove = os.path.join(sim_dir, str(dir_to_remove))
        print("Removing", full_dir_to_remove)
        shutil.rmtree(full_dir_to_remove)


def rename_directories(sim_dir, num_of_simulations):
    dirs = os.listdir(sim_dir)
    lst_of_dirs = []
    for file in dirs:
        sim_i_dir = os.path.join(sim_dir, file)
        if not os.path.isdir(sim_i_dir):
            continue
        lst_of_dirs.append(int(file))
    lst_of_dirs.sort()
    for i in range(num_of_simulations):
        curr_dir = os.path.join(sim_dir, str(lst_of_dirs[i]))
        dst = os.path.join(sim_dir, str(i))
        if i != lst_of_dirs[i]:
            os.rename(curr_dir, dst)





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Retains only successful simulations')
    parser.add_argument('--num_of_simulations', '-n', type=int, help='number of simulations')
    parser.add_argument('--output_file_path', '-o', help='output file path of stats')
    parser.add_argument('--sim_dir', '-d', help='simulation directory')
    parser.add_argument('--initial_num_of_simulations', '-i', type=int, help='initial number of simulations')

    # parse arguments
    args = parser.parse_args()
    num_of_simulations = args.num_of_simulations
    stats_file_path = args.output_file_path
    sim_dir = args.sim_dir
    initial_num_simulations = args.initial_num_of_simulations
    check_simulations(sim_dir, initial_num_simulations, stats_file_path)
    pick_successful_simualtions(sim_dir, num_of_simulations)
    rename_directories(sim_dir, num_of_simulations)
    check_simulations(sim_dir, num_of_simulations, stats_file_path)


