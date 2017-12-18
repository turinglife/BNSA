# Copyright (c) 2014 by
#
# Data Analytics Group
# Advanced Computing Research Centre (ACRC)
# School of Information Technology & Mathematical Sciences
# University of South Australia (UniSA)
#
# Implemented by
# Xin Zhu (A visiting research student at UniSA from USTC)
# zhuxy021@mymail.unisa.edu.au (UniSA)
# sa612114@mail.ustc.edu.cn (USTC)
#
# This is the main function for BNSAs.
# Citation: 
# Bing Liu, Jiuyong Li, Anna Tsykin, Lin Liu, Arti B. Gaur and Gregory J. Goodall, 
# Exploring Complex miRNA-mRNA Interactions with Bayesian Networks by Splitting-Averaging Strategy, 
# BMC Bioinformatics 2009, 10:408
#
# BNSA is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either this version of the License,
# or (at your option) any later version.
#
# BNSA is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.


###################################################################
# Description: BNSA user inferface

# @p_exp: matrix of expression values, for parent node set (data.frame)
# @c_exp: matrix of expression values, for children node set (data.frame)
# @class: class lable vector for the sample matched expression profiles in p_exp and c_exp. 
#         For example, class <- c(rep('Cancer', 6), rep('Normal', 6)) (vector)
# @t_db: target database used for initial network structure
# @p_nodes: parent nodes used for regulatory networks
# @c_nodes: children nodes used for regulatory networks
# @numiter: number of iterations
# @seed: the default random seed is 12345. It's used to test if the results are reproducible.
# @path: the path to store the results produced by BNSA, default path is current one.

# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
bnsa <- function (p_exp, c_exp, class, t_db, p_nodes, c_nodes, numiter, seed = 12345, path = "./")
{
	# Check parameters
	if(nargs() < 7)
    {
		usage()        

	 	return()       
    }

	if(!is.double(numiter))
	{
		print("the type of last parameter is not correct!")
		usage()

		return()
	}

	if(length(p_exp) < length(class) || length(c_exp) < length(class))
	{
		print("set number of samples is larger than actual number of samples")
		usage()

		return()	
	}

    # Install basic packages
	install_packages(t_db)

	# set random seed
	set.seed(seed)
	
	# Retrieve the dimension names of p_exp and c_exp
	p_dimensionnames <- dimnames(p_exp)
	p_rownames = p_dimensionnames[[1]]
	p_colnames = p_dimensionnames[[2]]

	c_dimensionnames <- dimnames(c_exp)
	c_rownames = c_dimensionnames[[1]]
	c_colnames = c_dimensionnames[[2]]

	p_exp_m <- matrix(unlist(p_exp), ncol = length(p_colnames))
	c_exp_m <- matrix(unlist(c_exp), ncol = length(c_colnames))
	
	# Retrieve the category information
	len_class1 <- length(grep((class[1]), class, ignore.case = FALSE, fixed = TRUE))
	len_class2 <- length(class) - len_class1

	# Compute the index of specific nodes(p_nodes, c_nodes) in node sets(p_rownames, c_rownames)
	index_of_pnodes_in_prownames <- match(p_nodes, p_rownames)
	index_of_cnodes_in_crownames <- match(c_nodes, c_rownames)
	    
	# Generate matrix of parent nodes and children nodes
	p_nodes_c1 <- p_exp_m[index_of_pnodes_in_prownames, 1:len_class1]
	p_nodes_c2 <- p_exp_m[index_of_pnodes_in_prownames, (len_class1+1):length(p_colnames)]

	c_nodes_c1 <- c_exp_m[index_of_cnodes_in_crownames, 1:len_class1]
	c_nodes_c2 <- c_exp_m[index_of_cnodes_in_crownames, (len_class1+1):length(c_colnames)]
	
	dimnames <- list(p_rownames[index_of_pnodes_in_prownames], c_rownames[index_of_cnodes_in_crownames])
	
	pnodes_c1 <- discretize(p_nodes_c1)
	pnodes_c2 <- discretize(p_nodes_c2)
	cnodes_c1 <- discretize(c_nodes_c1)
	cnodes_c2 <- discretize(c_nodes_c2)

	expression_matrices <- list(pnodes_c1, pnodes_c2, cnodes_c1, cnodes_c2)
	
	db <- hgu133aENTREZID
	mRNA.probe <- retrieveEntrezIDwithProbe(mRNA.DE.list, db)

	# Given miRNA, mRNA DE list, retrieve miRNA targets using gene probe names
	predtargets <- retrievemiRtargetsfromHsDB(miR.DE.list, mRNA.probe)

	# format path
	len <- nchar(path)
	last_character <- substr(path, len, len)
	print(last_character)
	
	if(last_character == "/")
	{
		savepath <- sprintf("%s", path)
	}
	else
	{
		savepath <- sprintf("%s/", path)	
	}

	outdatapath = savepath
	midoutputfile = 'results.RData'
	finaloutputfile = 'results_sig.RData'

	outputfilename <- construct_path(outdatapath, midoutputfile)

	# Invoke BNSA interface
	retval <- main_bnsa(dimnames, expression_matrices, predtargets, len_class1, len_class2, numiter, outputfilename)
  	if(retval)
	{
    	print("main_bnsa() retval is not correct!\n")
       	return()
	}
	
	# caculate pvalue
	main_cal_pvalue(outdatapath, midoutputfile, numiter, finaloutputfile)	

}


###################################################################
# Description: Usage information
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
usage <- function()
{
 	print("At least former seven input paramaters must be specified!")
    print("------------")
    print("Usage:")
    print("bnsa(p_exp, c_exp, class, t_db, p_nodes, c_nodes, numiter, seed, path)")
    print("------------")
    print("Arguments:")
    print("p_exp -- matrix of expression values, for parent node set (data.frame).")
    print("c_exp -- matrix of expression values, for children node set (data.frame).")
    print("class -- class lable vector for the sample matched expression profiles in p_exp and c_exp.")
	print("         For example, class <- c(rep('Cancer', 29), rep('Normal', 30)) (vector)")
    print("t_db -- target database used for initial network structure")
    print("p_nodes -- parent nodes used for regulatory networks")
    print("c_nodes -- children nodes used for regulatory networks")  
	print("numiter -- number of iterations")  
	print("seed (optional) -- the default random seed is 12345. It's used to test if the results are reproducible.")  
    print("path (optional) -- the path for saving the results. the default path is current path.")
    print("------------")
    print("Examples:")
	print("rd <- load('rBNSAObjects.RData')")
	print("rd")
	print("class <- c(rep('Cancer', 29), rep('Normal', 30))")
	print("bnsa(miRNA, mRNA, class, 'targetscan.Hs.eg.db', miR.DE.list, mRNA.DE.list, 10)")
	print("bnsa(miRNA, mRNA, class, 'targetscan.Hs.eg.db', miR.DE.list, mRNA.DE.list, 10, 80)")
	print("bnsa(miRNA, mRNA, class, 'targetscan.Hs.eg.db', miR.DE.list, mRNA.DE.list, 10, , 'd:/') for Windows platform")
    print("bnsa(miRNA, mRNA, class, 'targetscan.Hs.eg.db', miR.DE.list, mRNA.DE.list, 10, , '/') for Mac and Linux platform")
}

