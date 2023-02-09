import random
import argparse
from utils import *
from defs import *
from subprocess import Popen
from create_job_files import create_job_format_for_all_jobs


def prepare_data(dir_of_trees, out_dir, frac_sizes, size_to_sample, models, queue):
    original_trees_dir = os.listdir(dir_of_trees)
    for original_tree_dir in original_trees_dir:
        src_tree_dir = os.path.join(dir_of_trees, original_tree_dir)
        if not os.path.isdir(src_tree_dir):
            continue
        dst_tree_dir = os.path.join(out_dir, original_tree_dir)
        if not os.path.exists(dst_tree_dir):
            os.makedirs(dst_tree_dir)
        src_tree_path = os.path.join(src_tree_dir, "tree.newick")
        dst_tree_path = os.path.join(dst_tree_dir, "tree.newick")
        src_counts_path = os.path.join(src_tree_dir, "counts.fasta")
        dst_counts_path = os.path.join(dst_tree_dir, "counts.fasta")
        sample_tree_of_size(src_tree_path, size_to_sample, dst_tree_path, frac_sizes)
        sample_counts_file(dst_tree_path, src_counts_path, dst_counts_path)
        for model in models:
            control_file_path = os.path.join(dst_tree_dir, "paramFile_"+model)
            results_path = os.path.join(dst_tree_dir, "Results_"+model)
            write_homogeneous_model_control_file(dst_tree_dir, control_file_path, dst_counts_path,
                                             dst_tree_path, results_path, model)
            module_python = PYTHON_MODULE_COMMAND
            cmd = module_python
            chromevol_exe = CHROMEVOL_EXE
            cmd += chromevol_exe + " param=" + control_file_path + " > " \
                   + os.path.join(dst_tree_dir, "out_"+model+".txt") + "\n"
            job_path = os.path.join(dst_tree_dir, "inference_"+model+".sh")
            job_file = open(job_path, 'w')
            job_content = create_job_format_for_all_jobs(dst_tree_dir, "CE_"+model, "4gb", queue, 1, cmd)
            job_file.write(job_content)
            job_file.close()
            Popen(["qsub", job_path])


def sample_counts_file(tree_path, src_counts_path, dst_counts_path):
    tree = Tree(tree_path, format=1)
    leaves = tree.get_leaves()
    leaves_names = set([leaf.name for leaf in leaves])
    species_with_counts = read_counts_data(src_counts_path)
    sampled_counts = ""
    for species, count in species_with_counts:
        if species in leaves_names:
            sampled_counts += ">"+species+"\n"
            sampled_counts += count +"\n"
    counts_file = open(dst_counts_path, 'w')
    counts_file.write(sampled_counts)
    counts_file.close()


def sample_tree_of_size(tree_path, size, sampled_tree_path, frac_sizes_orig):
    tree = Tree(tree_path, format=1)
    leaves = tree.get_leaves()
    leaves_names = [leaf.name for leaf in leaves]
    max_num_of_iterations = 1000
    success = False
    for i in range(max_num_of_iterations):
        frac_sizes = frac_sizes_orig.copy()
        print("number of iteration is:", str(i))
        leaves_to_sample = random.sample(leaves_names, size)
        tree_to_sample = tree.copy("newick")
        tree_to_sample.prune(leaves_to_sample, preserve_branch_length=True)
        all_clade_sizes = search_node_sizes(tree_to_sample)
        for clade_size in all_clade_sizes:
            for frac in frac_sizes:
                frac_current = clade_size/size
                if (frac_current > frac + 0.02) or (frac_current < frac-0.02):
                    continue
                print(str(frac_current))
                frac_sizes.remove(frac)
                break
            if frac_sizes == []:
                success = True
                break
        if success:
            tree_to_sample.write(format=1, outfile=sampled_tree_path)
            return



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='samples trees and chromosome counts from empirical trees')
    parser.add_argument('--dir_of_trees', '-t', help='list of tree_paths')
    parser.add_argument('--out_dir', '-o', help='output_dir')
    parser.add_argument('--frac_sizes', '-f', nargs='+', type=float, help="fractions of tree sizes to sample")
    parser.add_argument('--size', '-s', type=int, help='size to sample')
    parser.add_argument("--models", "-m", help="models", default=["g+L_l+L_du+E_de+C_b+C"])
    parser.add_argument('--queue', '-q', help="queue")

    # parse arguments
    args = parser.parse_args()
    dir_of_trees = args.dir_of_trees
    out_dir = args.out_dir
    frac_sizes = args.frac_sizes
    size_to_sample = args.size
    models = args.models
    queue = args.queue
    prepare_data(dir_of_trees, out_dir, frac_sizes, size_to_sample, models, queue)











