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
# Description: Score for every possible candidate graph
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
score_bipartite_dags <- function(subdata, subdags)
{	
	LV3 = c(1,2)

	num_node = nrow(subdata)
        
	# generate a data frame
	fdata = data.frame(miR = factor(subdata[1, ], levels = LV3))

	for(index in 2:num_node)
	{
		if(index != num_node)
		{
			fdata = data.frame(fdata, miR = factor(subdata[index, ], levels = LV3))
		}
		else
	    {
			mR = factor(subdata[num_node, ], levels = LV3)		
		}
	}

	fdata <- data.frame(fdata, mR)
	
	# assign new column name(node name)
	for(index in 1:num_node)
	{
		if(index != num_node)
		{
			# the column will be named from miR1 to miRx
			colnames(fdata)[index] <- c(sprintf("miR%d", index))
		}
		else
		{
			# the last column will be named mR
			colnames(fdata)[index] <- c("mR")
		}
			
	}

	num_candidate_dags = 2^(num_node - 1)
	
	score_vector <- vector()
	for(index in 1:num_candidate_dags)
	{
		res = empty.graph(names(fdata))
		modelstring(res) = gen_graph(subdags[index])
		#plot(res)
		
		score_vector[index] <- score(res, fdata, type = "bde", debug = F, iss = 1)
	}
	
	return (score_vector)
}
