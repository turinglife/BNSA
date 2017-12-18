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
# Description: Enumerate all possible graph
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
enumerate_bipartite_graphs <- function(num_miR, num_mR)
{
	sub_num = num_miR + num_mR
	num_dags = 2^num_miR

	dags <- list()

	# canditate directed acyclic graph
	for(i in 1:num_dags)
	{
		dags[i] <- list(matrix(data = 0, nrow = sub_num, ncol = sub_num))
			
		#matrix_dags = dags[[i]] 
		for(j in 1:num_miR)
		{
			dags[[i]][j, sub_num] = as.integer(intToBits(i-1))[j] 
		}
	}

	return(dags)
}
