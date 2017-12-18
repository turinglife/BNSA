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
# Description: Calculate the pvalue
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
main_cal_pvalue <- function(outdatapath, outputfile, numiter, finaloutput)
{
	## Step 1. 	

	## Step 2. Set the paramaters used for significace estimation
	alpha = 0.05;
	p = 0.5

	## Step 3. Read the results produced by main_BNSA.m
	filename <- construct_path(outdatapath, outputfile)
	load(filename)
	
	o_miR_probe = results[1] 
	o_miR = results[2] 
	o_mR_probe = results[3]
	o_gene_sym = results[4]
	o_edge_score = results[5] 
	o_p_score = results[6] 

    vec_o_miR_probe <- convert_l_to_v_with_lowercase(o_miR_probe)
	vec_o_miR <- convert_l_to_v_with_lowercase(o_miR)
	vec_o_mR_probe <- convert_l_to_v_with_lowercase(o_mR_probe)
	vec_o_gene_sym <- convert_l_to_v_with_lowercase(o_gene_sym)
	vec_o_edge_score <- convert_l_to_v_with_lowercase(o_edge_score)
	vec_o_p_score <- convert_l_to_v_with_lowercase(o_p_score)


	## Step 4. Estimate the confidence of the results
	num_pair <- length(vec_o_miR_probe)

	score <- vector()
	pvalue <- vector()
	comment <- vector()

	for(index in 1:length(vec_o_edge_score))
	{
		score[index] = (as.double(vec_o_edge_score))[index] * (2 * numiter)		
	}	
	
	for(index in 1:num_pair)
	{
		pvalue[index] = 1 - pbinom(score[index], (2 * numiter), p)

		if(pvalue[index] < alpha)
		{
			comment[index] = "High Confidence"
		}
		else
		{
			comment[index] = "Low Confidence"
		}

		#print("pvalue -- comment")
		#print(pvalue[index])
		#print(comment[index])
	}
	
	s <- sort(score, decreasing = TRUE, index.return = TRUE)
	
	vec_index <- vector()

	for(index in 1:length(s$i))
	{
		vec_index[index] = s$i[index]
	}
		
	## Step 5. Save the results
	outputfilename <- construct_path(outdatapath, finaloutput)
	
	results_sig <- vector()
	for(index in 1:num_pair)
	{
		ind = vec_index[index]

		miR_probe <- vec_o_miR_probe[ind]
		miR <- vec_o_miR[ind]
		mr_probe <- vec_o_mR_probe[ind]
		mr <- vec_o_gene_sym[ind]
		edge_score <- as.double(vec_o_edge_score[ind])
		p_score <- as.double(vec_o_p_score[ind])
		p_value <- pvalue[ind]
		comments <- comment[ind]

		c <- data.frame(miR_probe, miR, mr_probe, mr, edge_score, p_score, p_value, comments)
		results_sig <- rbind(results_sig, c)
	
	}
	
	save(results_sig, file = outputfilename)

	print("Complete Successfully...!")
}


