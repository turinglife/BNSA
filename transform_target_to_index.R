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
# Description: Find the target index in the expression values
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
transform_target_to_index <- function(nodenames, target)
{
    node_miR = nodenames[[1]]
    node_mR = nodenames[[2]]
                
    list_miR_probe = target[1]
    list_mR_probe = target[3]
 
    len_miR_probe <- length(unlist(list_miR_probe))
    len_mR_probe <- length(unlist(list_mR_probe))

    vec_miR_probe <- vector()
	vec_mR_probe <- vector()
    for(index in 1:len_miR_probe)
    {
        vec_miR_probe[index] <- c(tolower(list_miR_probe[[1]][index]))
		vec_mR_probe[index] <- c(tolower(list_mR_probe[[1]][index]))
    }

	# target_index_miR_p and target_index_mR_p are one-dimension column vector
    target_index_miR_p <- array(0, dim = c(len_miR_probe, 1))
    target_index_mR_p <- array(0, dim = c(len_mR_probe, 1))

    for(index in 1:length(node_miR))
    {
		# return the location of node_miR[i] in vec_miR_probe 
    	location <- grep(node_miR[index], vec_miR_probe, fixed = TRUE)

		target_index_miR_p[location] = index		
    }

	for(index in 1:length(node_mR))
	{
		location <- grep(node_mR[index], vec_mR_probe, fixed = TRUE)

		target_index_mR_p[location] = index
	}

	target_index <- list(target_index_miR_p, target_index_mR_p)		
	
	return (target_index)
}
