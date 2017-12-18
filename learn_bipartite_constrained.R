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
# SA-BNs is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either this version of the License,
# or (at your option) any later version.
#
# SA-BNs is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

###################################################################
# Description: Learn on Bayesian Network
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
learn_bipartite_constrained <- function(es_miR_E, es_mR_E, target)
{
	r_miR_E <- nrow(es_miR_E)
	c_miR_E <- ncol(es_miR_E)
	
	r_mR_E <- nrow(es_mR_E)
	c_mR_E <- ncol(es_mR_E)
	
	if(c_miR_E != c_mR_E)
	{
		print("Datasets do not have same column size!")

		stop()
	}
	
	# create a class that encapsulates return information 
	#setClass(Class="bipartite_learned_with_constraint", representation(dag = "matrix", dag_score = "numeric"))
	
	# merge es_miR_E and es_mR_E into one matrix
	data = rbind(es_miR_E, es_mR_E)
	
	dag = matrix(data = 0, nrow = r_miR_E, ncol = (r_miR_E + r_mR_E))
	dagscore = matrix(nrow = r_mR_E, ncol = 1)

	t_miR_probe = target[[1]]
	t_mR_probe = target[[2]]
	
	for(index in 1:r_mR_E)
	{
		print(sprintf("Searching parents for node %2d (out of %d) ...", index, r_mR_E))
		pattern <- sprintf("\\b%d\\b", index)
		
		location = grep(pattern, t_mR_probe)
		
		num_miR = length(location)
		index_t_miR_probe = t_miR_probe[location]

		subdata = rbind(es_miR_E[index_t_miR_probe, ], es_mR_E[index, ])
		subdags <- enumerate_bipartite_graphs(num_miR, 1)
		
		subscore = score_bipartite_dags(subdata, subdags)	
		dagscore[index] = max(subscore)	
		index_maxscore = which.max(subscore)

		print(sprintf("done, dag No.%d, score is %6.4f", index_maxscore, dagscore[index]))
    	#fprintf('\n dag %d is the best one for node %d, score is %6.4f\n', best_p, i, dagscore(i))

		# transform the subgraph to graph
		dag_with_maxscore = subdags[[index_maxscore]]
		ind = which(dag_with_maxscore[1:num_miR, num_miR + 1] == 1)
		trans_index = index_t_miR_probe[ind]

		dag[trans_index, (r_miR_E + index)] = 1
	}
	
	all_dag_score = sum(dagscore)
	print(sprintf("The final score of the graph is %6.8f", all_dag_score))
	
	#return(new("bipartite_learned_with_constraint", dag = dag, dag_score = all_dag_score))
	
	ret <- list(dag, all_dag_score)

	return(ret)

}
