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
# Description: Generate a string representing the corresponding graph
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
gen_graph <- function(dag)
{
	matrix_dag = (dag)[[1]]
	
	# number of miR node
	num_miR = nrow(matrix_dag) - 1
	
	# the index of the last column in dag matrix
	index_last_col = nrow(matrix_dag)

	# concatenate miR nodes string
	miR_node_sequence = "[miR1]"
	index = 2
	
	while(index <= num_miR)
	{
		miR_node_sequence <- paste(miR_node_sequence, sprintf("[miR%d]", index), sep = "")
	
		index = index + 1
	}

	# generate mR node string
	num_miR_regulate_mR = 0
	for(index in 1:num_miR)
	{
		if(matrix_dag[index, index_last_col] == 1)
		{
			num_miR_regulate_mR = num_miR_regulate_mR + 1	
		}
	}
	
	need_vertical_bar = TRUE

	mR_node_sequence = "[mR"
	for(index in 1:num_miR)
	{
		if(matrix_dag[index, index_last_col] == 1)
		{
			# add a vertical bar after the first miR regulating mR
			# for instance, it looks like "[mR|miR1]"
			if(need_vertical_bar)
			{
				mR_node_sequence <- paste(mR_node_sequence, sprintf("|", index), sep = "")
				need_vertical_bar = FALSE 
			}
			
			mR_node_sequence <- paste(mR_node_sequence, sprintf("miR%d", index), sep = "")
			
			# add a colon between each miR regulating mR
			# for instance, it looks like "[mR|miR1:miR2:miR3]"
			num_miR_regulate_mR = num_miR_regulate_mR - 1
			if(num_miR_regulate_mR > 0)
			{
				mR_node_sequence <- paste(mR_node_sequence, sprintf(":"), sep = "")
			}
		}
	}
	mR_node_sequence <- paste(mR_node_sequence, sprintf("]"), sep = "")

	# merge miR and mR into node graph
	# for instance, the final graph may look like "[miR1][miR2][miR3][mR|miR1:miR2]"
	graph <- paste(miR_node_sequence, mR_node_sequence, sep = "")

	return (graph)
}
