from defs import *


def average(lst):
    """
    Calculates the average of a given list
    :param lst: list of numbers
    :return: average value
    """
    return sum(lst) / len(lst)


def range_of_lst(lst):
    """
    calculates the range of a given list
    :param lst: list of numbers
    :return: max-min
    """
    return max(lst) - min(lst)


def fix_simulated_tree_file(tree_file):
    """
    a patch that adds a semicolon at the end of .phr tree (simulated tree) so the can tree can be read
    :param tree_file:
    :return: NA
    """
    with open(tree_file, "a") as add:
        add.write(";")


def create_counts_hash(counts_file):
    """
    puts all counts from a counts file in a dictionary: taxa_name: count. Taxa with an X count are skipped.
    If there are counts that are X the function returns two hashes and a set of the taxa to prune.
    :param counts_file in FASTA format
    :return:(2) dictionary of counts
    """
    d = {}
    with open(counts_file, "r") as counts_handler:
        for line in counts_handler:
            line = line.strip()
            if line.startswith('>'):  # taxon name
                name = line[1:]
            else:
                if line != "x":
                    num = int(line)
                    d[name] = num
    return d


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


def get_base_numbers(params_list):
    parameters_joined = "".join(params_list)
    pattern = re.compile("_baseNum_[\d]+[\s]+=[\s]+[\d]+;([\d]+)")
    base_numbers_str = pattern.findall(parameters_joined)
    base_numbers = [int(i) for i in base_numbers_str]
    return base_numbers


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


def contains_base_number(parameters_lst):
    """
    checks whether the model contains base number
    :param: list of parameters (as lines)
    :return: boolean. true if base number is included in the model. False otherwise
    """
    pattern = re.compile("_baseNumRFunc[\s]+=[\s]+CONST")
    for param_line in parameters_lst:
        if pattern.findall(param_line) != []:
            return True
    return False


def create_dict_of_parameters(params_lst):
    """
    creates a dictionary from a list of parameters
    :param: list of parameters (as lines)
    :return: a dictionary with parameters and their values
    """
    pattern = re.compile("([\S]+)[\s]+=[\s]+([\S]+)")
    params_dict = {}
    for param_entry in params_lst:
        param, value = pattern.findall(param_entry)[0]
        params_dict[param] = value
    return params_dict


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


def extract_line_from_file(filename, str_search, num=False, integer=False):
    """
    uses regular expression to search a string in a file
    :param filename: the file to be searched in
    :param str_search: the string to search for
    :param num: should a number be extracted?
    :param integer: should an integer be returned? if False and num=True then a float will be returned
    :return: True/False for a sentence
                integer if num=True and integer=True
                float if num=True and integer=False
    """
    with open(filename, "r") as fh:
        for line in fh:
            line = line.strip()
            if num:
                tmp = re.search(str_search + ".*?(\d+)", line)
                if tmp:
                    return int(tmp.group(1)) if integer else float(tmp.group(1))
            else:
                tmp = re.search(str_search, line)
                if tmp:
                    return True
        return False


def print_error(msg, exception="", sep="*****"):
    print(sep)
    print(msg, "\n")
    print(sep)
    if exception is not "":
        print(exception)
        print(sep)


def get_max_base_number(params_lst):
    pattern = re.compile("_baseNum_[\d]+[\s]+=[\s]+[\d]+;([\d]+)")
    joined_params_lst = "\n".join(params_lst)
    base_numbers_str = pattern.findall(joined_params_lst)
    base_numbers = [int(i) for i in base_numbers_str]
    return max(base_numbers)


