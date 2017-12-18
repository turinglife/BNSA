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
# Description: Merge two local dags into one global dag
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
average_dags_bootstrp <- function (dags, num_miR_node, num_mR_node)
{
	numdags <- length(dags)
	
	# cut out bipartite from dags
	l_b1 = dags[1][[1]]
	l_b2 = dags[2][[1]]
	
	# transform from list to matrix type
	m_dag1 = l_b1[1][[1]]
	m_dag2 = l_b2[1][[1]]
	
	# obtain dag socre
	dag1_score = l_b1[2][[1]]
	dag2_score = l_b2[2][[1]]

	# choose the corresponding strategy according to numdags
	if(numdags == 1)
	{
		return (dags)
	}
	else if(numdags == 2)
	{
		# average on the dags
	}
	else
	{
		print("The number of data set is not eligibile for this algorithm!")
		
		return
	}
	
	# create two matrix for storing average values
	dag = matrix(data = 0, nrow = num_miR_node, ncol = (num_miR_node + num_mR_node))
	score = matrix(data = 0, nrow = num_miR_node, ncol = (num_miR_node + num_mR_node))

	for(i in 1:num_miR_node)
	{
		for(j in (num_miR_node + 1): (num_miR_node + num_mR_node))
		{
			if(m_dag1[i, j] == 0 && m_dag2[i ,j] == 0)
			{
				# do nothing
			}
			else if(m_dag1[i, j] != 0 && m_dag2[i ,j] == 0)
			{
				dag[i, j] = 1
				score[i, j] = dag1_score
			}
			else if(m_dag1[i, j] == 0 && m_dag2[i ,j] != 0)
			{
				dag[i, j] = 1
				score[i, j] = dag2_score
			}
			else
			{
				dag[i, j] = 2
				score[i, j] = max(dag1_score, dag2_score)
			}
		}
	}

	return(list(dag, score))
}



