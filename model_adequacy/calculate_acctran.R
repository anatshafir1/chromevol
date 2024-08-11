rm(list = ls())
library(phytools)
library(phangorn)

# arguments: (1) tree_file (2) simulations directory (3) number of iterations - 1 for original counts, n for simulations)
# output: list of calculated statistic N times

args = commandArgs(trailingOnly=TRUE)
tree_file = args[1]
iter = args[2]
sims_dir = args[3]


fix_tree = function(tree_path){
  fileConn <- file(tree_path ,open="a")
  cat(";",file = fileConn)
  close(fileConn)
}

calculate_stat = function(tree){
  counts = c()
  
  for (i in 1:(length(tree$tip.label))){
    counts = c(counts, as.numeric(gsub("[^0-9]", "",  tree$tip.label[i])))
    tree$tip.label[i] = gsub('-[0-9]+', '', tree$tip.label[i])
  }
  mat_data = matrix(counts,dimnames = list(tree$tip.label,NULL), nrow=length(counts), byrow=TRUE)
  data = phyDat(mat_data, type = "USER", levels = unique(counts))
  tmp = acctran(tree,data)
  
  tmp_mat = as.data.frame(cbind(tmp$edge,tmp$edge.length))
  names(tmp_mat) = c("begin","end","parsimony")
  tree_mat = as.data.frame(cbind(tree$edge,tree$edge.length))
  names(tree_mat) = c("begin","end","time")
  mat = merge(tree_mat, tmp_mat, by=c("begin","end"))
  mat$time = 0
  distances = dist.nodes(tree)
  root_ind = which(node.depth.edgelength(tree)==0)
  for (i in 1:nrow(mat)){
    end = mat$end[i]
    mat$time[i] = distances[root_ind,end]
  }
  
  x = mat$time
  y = mat$parsimony
  
  lm = lm(y~x)
  result = round(lm$coefficients[2],2)
  names(result) = NULL
  return(result)
}

results_lst = c()

if (iter==1){
  tree = read.newick(tree_file)
  results_lst = calculate_stat(tree)
} else {
  tree_file = "simulatedDataAncestors.tree"
  for (n in (1:iter)-1){
    #fix_tree(paste(sims_dir,n,tree_file,sep="/")) # fix .phy file
    tree = read.newick(paste(sims_dir,n,tree_file,sep="/"))
    stat = try(calculate_stat(tree))
    if (inherits(stat,'try-error')){
      next
    }
    results_lst = append(results_lst,stat)
  }
}
cat(results_lst)


