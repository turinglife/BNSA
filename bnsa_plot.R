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
# Description: Plot the running results
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
bnsa_plot <- function(data, level = "gene", path = "./")
{
	if(nargs() < 1)
    {
		plot_usage()        

	 	return()       
    }

	# format 
	v_col1_of_d <- convert_l_to_v_with_lowercase(data[1])
	v_col2_of_d <- convert_l_to_v_with_lowercase(data[2])
	v_col3_of_d <- convert_l_to_v_with_lowercase(data[3])
	
	if(level == "probe")
	{
		# plot in the probe level
		plot_in_probe_level(v_col1_of_d, v_col3_of_d, path)
	}
	else if(level == "gene")
	{
		# plot in the gene level
		plot_in_gene_level(v_col2_of_d, v_col3_of_d, path)
	}
	else
	{
		print("The default plotting levels consist of gene(default) and probe!")
	}
}

###################################################################
# Description: Plot in the probe level
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
plot_in_probe_level <- function(v_col1_of_data, v_col3_of_data, path)
{
	# Total number of rows	
	i_total_rows <- length(v_col1_of_data)

	# Extract unique elements and return a vector 
	v_pnodes <- unique(v_col1_of_data, fromLast = TRUE)
	v_cnodes <- unique(v_col3_of_data, fromLast = TRUE)
	
	# Length of vector
	i_num_pnodes <- length(v_pnodes)
	i_num_cnodes <- length(v_cnodes)

	# Build a matrix for representing the relatioships between parent nodes and children nodes
   	m_relation_matrix <- matrix(data = 0, nrow = (i_num_pnodes + i_num_cnodes), ncol = (i_num_pnodes + i_num_cnodes))
	
 	# Vector for storing the name of nodes 
	row_names <- c(v_pnodes, v_cnodes)
	col_names <- c(v_pnodes, v_cnodes)
	
	# Assign name of rows and columns for matrix
	rownames(m_relation_matrix) <- c(row_names)
	colnames(m_relation_matrix) <- c(col_names)
	
	for(index in 1:i_total_rows)
	{
		# Retrieve the position of relations(parents[miR] -> children[mR]) in matrix
    	pnode_index <- grep(v_col1_of_data[index], row_names, fixed = TRUE)
        cnode_index <- grep(v_col3_of_data[index], col_names, fixed = TRUE)

		# Store the value representing relation into matrix
        m_relation_matrix[pnode_index, cnode_index] = 1
	}   
	
	# Generate graph according to m_relation_matrix 
	m_graph = graph.adjacency(m_relation_matrix, mode = "directed", weighted = NULL, diag = FALSE)

    # Parent nodes coloring as red
	for(index in 1:i_num_pnodes)
	{
		V(m_graph)[index]$color <- c("red")
    }

	# Children nodes coloring as green
	for(index in (i_num_pnodes + 1):(i_num_pnodes + i_num_cnodes))
	{
	    V(m_graph)[index]$color <- c("green")
	}
	
	# Figure will be saved as pdf file.
	len <- nchar(path)
	last_character <- substr(path, len, len)
	
	if(last_character == "/")
	{
		savepath <- sprintf("%srelation.pdf", path)
	}
	else
	{
		savepath <- sprintf("%s/relation.pdf", path)	
	}
	
	pdf(savepath, width = 30, height = 20)
  	
	# Plotting with some parameters 
	plot(m_graph, layout = layout.fruchterman.reingold, vertex.size = 4, vertex.label.dist = 0, edge.arrow.size = 0.3)
	
    # After pdf was built, closing needs to be invoked.
	dev.off()

	print("Plot in probe level successfully!")
}

###################################################################
# Description: Plot in the gene level
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
plot_in_gene_level <- function(v_col2_of_data, v_col3_of_data, path)
{
	# total rows
	i_total_rows <- length(v_col2_of_data)
	
	# hgu133aSYMBOL is an R object that provides mappings between manufacturer identifiers and gene abbreviations
	db <- hgu133aSYMBOL
	symdb <- as.list(db[as.character(v_col3_of_data)])
	
	# produce the gene symbol vector corresponding to mRNA probes
	v_symbol <- vector()
	for(index in 1:length(v_col3_of_data))
	{
		v_symbol[index] <- symdb[as.character(v_col3_of_data[index])]	
	}
	
	# miRNA family vector
	v_mir_family <- v_col2_of_data

	# Extract unique elements and return a vector 
	v_pnodes <- unique(v_mir_family, fromLast = TRUE)
	v_cnodes <- unique(v_symbol, fromLast = TRUE)
	
	# Length of vector
	i_num_pnodes <- length(v_pnodes)
	i_num_cnodes <- length(v_cnodes)

	# Build a matrix for representing the relatioships between parent nodes and children nodes
   	m_relation_matrix <- matrix(data = 0, nrow = (i_num_pnodes + i_num_cnodes), ncol = (i_num_pnodes + i_num_cnodes))
	
 	# Vector for storing the name of nodes 
	row_names <- c(v_pnodes, v_cnodes)
	col_names <- c(v_pnodes, v_cnodes)
	
	# Assign name of rows and columns for matrix
	rownames(m_relation_matrix) <- c(row_names)
	colnames(m_relation_matrix) <- c(col_names)
	
	for(index in 1:i_total_rows)
	{
		# Retrieve the position of relations(parents[miR] -> children[mR]) in matrix
    	pnode_index <- grep(v_mir_family[index], row_names, fixed = TRUE)
        cnode_index <- grep(v_symbol[index], col_names, fixed = TRUE)

		# Store the value representing relation into matrix
        m_relation_matrix[pnode_index, cnode_index] = 1
	}   
	
	# Generate graph according to m_relation_matrix 
	m_graph = graph.adjacency(m_relation_matrix, mode = "directed", weighted = NULL, diag = FALSE)

    # Parent nodes coloring as red
	for(index in 1:i_num_pnodes)
	{
		V(m_graph)[index]$color <- c("red")
    }

	# Children nodes coloring as green
	for(index in (i_num_pnodes + 1):(i_num_pnodes + i_num_cnodes))
	{
	    V(m_graph)[index]$color <- c("green")
	}
	
	# Figure will be saved as pdf file.	
	len <- nchar(path)
	last_character <- substr(path, len, len)
	
	if(last_character == "/")
	{
		savepath <- sprintf("%srelation.pdf", path)
	}
	else
	{
		savepath <- sprintf("%s/relation.pdf", path)	
	}
	
	pdf(savepath, width = 30, height = 20)
  	
	# Plotting with some parameters 
	plot(m_graph, layout = layout.fruchterman.reingold, vertex.size = 4, vertex.label.dist = 0, edge.arrow.size = 0.3)
	
    # After pdf was built, closing needs to be invoked.
	dev.off()
	
	print("Plot in gene level successfully!")
}

###################################################################
# Description: Usage
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
plot_usage <- function()
{
	print("Plot the result generated from BNSA")
	print("------------")
    print("Usage:")
	print("bnsa_plot(data, level, path)")
	print("data -- to be analyzed")
	print("level -- probe or gene(default)")
	print("path -- address in which plotting will be stored(default path is current address)")
	print("Examples:")
	print("bnsa_plot(results)")
	print("bnsa_plot(results, probe)")
	print("bnsa_plot(results, ,'d:/') for Windows platform")
    print("bnsa_plot(results, ,'/') for Mac and Linux platform")

}

