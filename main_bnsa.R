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
# Description: This is a main function of BNSA, which contains the 
#              most important flowchart.
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
main_bnsa <- function (dimnames, expression_datasets, target, len_class1, len_class2, numiter, outputfilename)
{  
    ## Step 1. Unify the node names for the data set. 
	
	# Step 1.1 data from target databases
	# convert list to vector and uppercase to lowercase
	vec_t_col1 <- convert_l_to_v_with_lowercase(target[1])
	vec_t_col3 <- convert_l_to_v_with_lowercase(target[3])
	vec_t_col2 <- convert_l_to_v_with_lowercase(target[2])
	vec_t_col4 <- convert_l_to_v_with_lowercase(target[4])
	
	# remove redundant data
	t_pnodes <- unique(vec_t_col1, fromLast = TRUE)               
    t_cnodes <- unique(vec_t_col3, fromLast = TRUE)                 
	
    # Step 1.2 expression data from file or datafame
    e_pnodes <- vector()
    e_cnodes <- vector()

	# covert list to vector and uppercase to lowercase
	for(index in 1:length((unlist(dimnames[1]))))
        e_pnodes[index] <- c(tolower((unlist(dimnames[1]))[index]))
	
	for(index in 1:length((unlist(dimnames[2]))))
        e_cnodes[index] <- c(tolower((unlist(dimnames[2]))[index]))
    
	# Step 1.3 unify nodes between nodes from target and nodes from dataframe.
    nodename_miR <- intersect(t_pnodes, e_pnodes)
	nodename_miR <- sort(nodename_miR, decreasing = FALSE)
    index_miRnode_in_es_miR_p <- match(nodename_miR, e_pnodes)
	
	nodename_mR <- intersect(t_cnodes, e_cnodes)
    nodename_mR <- sort(nodename_mR, decreasing = FALSE)
    index_mRnode_in_es_mR_p = match(nodename_mR, e_cnodes)

	# Step 1.4 integrate nodename_miR with nodename_mR into nodenames
    nodenames <- list(nodename_miR, nodename_mR)
	
	num_miR_node <- length(nodename_miR)
    num_mR_node <- length(nodename_mR)


    ## Step 2. Transform the target to match the node index
    trans_target <- transform_target_to_index(nodenames, target);


    ## Step 3. Resampling the data set using bootstrapping 
    bootInd_E <- matrix(nrow = len_class1, ncol = numiter)
    bootInd_M <- matrix(nrow = len_class2, ncol = numiter)
    
    for(index in 1:numiter)
    {
        bootInd_E[, index] <- sample(len_class1, replace = TRUE)
        bootInd_M[, index] <- sample(len_class2, replace = TRUE)    
	}

    
    ## Step 4. Learn on the data sets with contraints
    l_dags <- list()     # local dags for storing the result of every iterations.
	g_dags <- list()     # global dags for storing the final result of all iterations.

    print(sprintf("Iterating start ... num is %d\n", numiter))
    for(index in 1:numiter)
    {
        print(sprintf("Round %2d out of %d \n", index, numiter))   
        
        # Learn on one category(class) sample data 
		value <- matrix(unlist(expression_datasets[1]), ncol = len_class1)
		miR_E = value[, bootInd_E[, index]]
		
        value <- matrix(unlist(expression_datasets[3]), ncol = len_class1)
        mR_E = value[, bootInd_E[, index]]

        bipartite = learn_bipartite_constrained(miR_E[index_miRnode_in_es_miR_p, ], mR_E[index_mRnode_in_es_mR_p, ], trans_target)
		b1 <- list(bipartite[[1]], bipartite[[2]])
		
        # Learn on the other category(class) sample data	
		value <- matrix(unlist(expression_datasets[2]), ncol = len_class2)
        miR_M = value[, bootInd_M[, index]]
		
        value <- matrix(unlist(expression_datasets[4]), ncol = len_class2)
        mR_M = value[, bootInd_M[, index]]

		bipartite = learn_bipartite_constrained(miR_M[index_miRnode_in_es_miR_p, ], mR_M[index_mRnode_in_es_mR_p, ], trans_target)
		b2 <- list(bipartite[[1]], bipartite[[2]])
		
		l_dags <- list(b1, b2)
		
        # Integrate results learning on two categories respectively. 
        ave_dag <- average_dags_bootstrp(l_dags, num_miR_node, num_mR_node)
		g_dags[index] <- list(ave_dag)
    }
    print("Iterating finished ...!\n")
	

    ## Step 5. Get the final DAG and score
    bootstrap_dag.dag <- matrix(data = 0, nrow = num_miR_node, ncol = num_miR_node + num_mR_node)
    bootstrap_dag.score <- matrix(data = 0, nrow = num_miR_node, ncol = num_miR_node + num_mR_node)
    
    for(i in 1:num_miR_node)
    {
        for(j in (num_miR_node + 1):(num_miR_node + num_mR_node))
        {
            for(k in 1:numiter)
            {	
				dag_with_score = g_dags[k][[1]]
				dag = dag_with_score[1][[1]]
				score = dag_with_score[2][[1]]

                bootstrap_dag.dag[i, j] = bootstrap_dag.dag[i, j] + dag[i, j]
				bootstrap_dag.score[i, j] = bootstrap_dag.score[i, j] + score[i, j]
            }
        }
    }
	
	# average the bootstrap dags
	bootstrap_dag.dag = bootstrap_dag.dag / (2 * numiter)
	bootstrap_dag.score = bootstrap_dag.score / numiter


    ## Step 6. Output results
	results <- vector()
	for(i in 1:num_miR_node)
	{
		for(j in 1:num_mR_node)
		{
			if(bootstrap_dag.dag[i, j + num_miR_node] != 0)
			{
				miR = nodenames[1][[1]]
				mR = nodenames[2][[1]]

				ind_miR <- match_index(miR[i], vec_t_col1)
				ind_mR <- match_index(mR[j], vec_t_col3)
				
				# generate output objects
				miR <- miR[i]
				mir_name <- vec_t_col2[ind_miR[1]]
				mr_probe <- mR[j]
    			genesym <- vec_t_col4[ind_mR[1]]
				edgescore <- bootstrap_dag.dag[i, j + num_miR_node]
				pscore <- bootstrap_dag.score[i, j + num_miR_node]

				c <- data.frame(miR, mir_name, mr_probe, genesym, edgescore, pscore)

				results <- rbind(results, c)				
			}
		}
	}

	save(results, file = outputfilename)

    print("Complete Successfully...!\n")
	
    return (0)
}
