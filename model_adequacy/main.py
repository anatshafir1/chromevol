from utils import *
from defs import *
from data_processing import *
from analysis import *

if __name__ == '__main__':
    # counts_file: original counts file on which was based the inference
    # tree_file: ml ancestors reconstruction tree
    # results_file: chromEvol.res file of the inferred model
    # user_output_dir: directory for outputs
    # nodes_to_models_path: nodes partition path (if the model is heterogeneous)
    # tree_shifts_path: tree shifts path (if the model is heterogeneous)
    counts_file, tree_file, results_file, user_output_dir, nodes_to_models_path, tree_shifts_path = get_arguments()
    model_adequacy_out_dir = user_output_dir + "/simulations/"
    if not os.path.exists(model_adequacy_out_dir):
        os.makedirs(model_adequacy_out_dir)
    parameters = get_params(results_file, user_output_dir + root_freq_filename, nodes_to_models_path, tree_shifts_path)
    match_counts_to_tree(tree_file, user_output_dir)
    original_counts = get_counts(counts_file)
    counts_range = range_of_lst(original_counts)
    original_counts_statistics = calculate_statistics(original_counts, user_output_dir + tree_with_counts, os.path.join(user_output_dir, "orig_stats"), True)
    has_baseNumber = contains_base_number(parameters)
    if has_baseNumber:  # base_num flag - execute second chromEvol run, if needed
        second_run_param_file_path = user_output_dir + "/second_run.params"
        if not os.path.exists(user_output_dir + "/second_run_tested"):
            second_run(parameters, user_output_dir, counts_file, user_output_dir + tree_wo_counts)
        update_max_base_transition(second_run_param_file_path, parameters, counts_range)


    nsim = get_nsims()
    run_simulations(results_file, model_adequacy_out_dir, original_counts, user_output_dir + tree_wo_counts, user_output_dir + root_freq_filename, parameters, user_output_dir)
    model_adequacy(model_adequacy_out_dir, original_counts_statistics, user_output_dir, original_counts_statistics)