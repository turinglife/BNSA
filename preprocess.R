# Copyright (c) 2014 by
#
# Data Analytics Group
# Advanced Computing Research Centre (ACRC)
# School of Information Technology & Mathematical Sciences
# University of South Australia (UniSA)
#
# Implemented by
# Bing Liu
# BLiu@ccia.unsw.edu.au
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
# Description: Calculate the cut-off of all samples according to median
# Author: Bing Liu
# Date: 2014-01-21
# Comments: 
###################################################################
BaselinetoAllSamplesMedian <- function(x)
{
	# The median is sample wise ( The margin 1 indicates rows)
	md <- apply(x, 1, median)
	x.centralized <- x - md
	x.centralized
}


###################################################################
# Description: Calculate the cut-off of specified samples according to median
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
BaselinetoControlSamplesMedian <- function(x, ind_ctrl)
{
	ctrlset <- x[,ind_ctrl]

	if(is.vector(ctrlset))
	{
  		x.centralized <- x - ctrlset 
	}
	else  
	{
		md <- apply(ctrlset, 1, median); 
		
		x.centralized <- x - md;
	}

	x.centralized  
}


###################################################################
# Description: Calculate the cut-off of all samples according to mean
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
BaselinetoAllSamplesMean <- function(x)
{
	# The median is sample wise ( The margin 1 indicates rows)
	md <- apply(x, 1, mean)
	x.centralized <- x - md
	x.centralized;
}


###################################################################
# Description: Calculate the cut-off of specified samples according to mean
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
BaselinetoControlSamplesMean <- function(x, ind_ctrl)
{
	ctrlset <- x[,ind_ctrl]

	if(is.vector(ctrlset)) 
		x.centralized <- x - ctrlset else  {md <- apply(ctrlset, 1, mean); x.centralized <- x - md;}
	x.centralized;  
}


#For Bayesian network
###################################################################
# Description: Discretization
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
dichotomy <- function(x, dval)
{
	x.dis <- x;
	up_index <- (x >= 0)
	down_index <- !up_index;

	x.dis[up_index] <- dval[2]
	x.dis[down_index] <- dval[1]

	x.dis
}


###################################################################
# Description: Discretization function based on column (sample) median or mean

# @expression_matrix: microarray expression data, column is samples, row is probes/gene/miRNA
# @method: using either median or mean
# @dval: what values are assigned to discretized data, default is 0 for negative values, and 1 for positive value based on cutoff of median or mean.

# Author: Bing Liu
# Date: 2014-01-21
# Modifier: Xin Zhu
# Comments: 
###################################################################
discretize <- function(pExprs, ind_ctrl, method = c('median', 'mean'), dval = c(1, 2))
{
	method <- match.arg(method)
	pExprs.scale <- switch(method,
						   median = BaselinetoControlSamplesMedian(pExprs, ind_ctrl),
						   mean = BaselinetoControlSamplesMean(pExprs, ind_ctrl),
						   default = BaselinetoAllSamplesMedian(pExprs))
    
	pExprs.dis <- dichotomy(pExprs.scale, dval)
	
	pExprs.dis
}


###################################################################
# Description: Convert list to table
# Author: Bing Liu
# Date: 2014-01-21
# Comments: 
###################################################################
list2table <- function(plist, cnames)
{
	library(plyr)
  	tb <- ldply(plist, cbind)
  	colnames(tb) <- cnames
  	tb <- unique(tb)
  	ind <- complete.cases(tb)
  	tb <- tb[ind, ]
  	tb
}


###################################################################
# Description: retrieve EntrezID with Probe
# Author: Bing Liu
# Date: 2014-01-21
# Comments: 
###################################################################
retrieveEntrezIDwithProbe <- function(probelist, db)
{
	entrezID <- as.list(db[as.character(probelist)])

	probe_entrezID <- data.frame(mR = probelist, EntrezID = unlist(entrezID))
	ind <- complete.cases(probe_entrezID)
  
  	probe_entrezID <- probe_entrezID[ind,]
	probe_entrezID <- unique(probe_entrezID)
	
    probe_entrezID
}


###################################################################
# Description: retrieve EntrezID with GeneSymbol 
# Author: Bing Liu
# Date: 2014-01-21
# Comments: 
###################################################################
retrieveEntrezIDwithGeneSymbol <- function(genelist, db)
{
	entrezID <- as.list(mget(genelist, db, ifnotfound = NA))
 	gene_entrezID <- data.frame(mR = genelist, EntrezID = unlist(entrezID))
 	ind <- complete.cases(gene_entrezID)
 	gene_entrezID <- gene_entrezID[ind,]
 	gene_entrezID <- unique(gene_entrezID)
 	gene_entrezID
}


###################################################################
# Description: Generate a target info according to rules 
# Author: Bing Liu
# Date: 2014-01-21
# Comments: 
###################################################################
retrievemiRtargetsfromHsDB <- function(miR.list, mRNA.list)
{
	library(targetscan.Hs.eg.db)
  	#retrieve miRNA-mRNA target from target database

  	targets <- mget(as.character(mRNA.list$EntrezID), envir = targetscan.Hs.egTARGETS, ifnotfound = NA)
  	miR.family <- mget(miR.DE.list, envir = targetscan.Hs.egMIRBASE2FAMILY, ifnotfound = NA)
  
  	target.table <- list2table(targets, c('EntrezID', 'miR.family'))
  	miR.family <- list2table(miR.family, c('miR', 'miR.family'))
  	target.pred <- merge(miR.family, target.table)
  	target.pred <- merge(target.pred, mRNA.list)
  	target.pred <- target.pred[, c('miR', 'miR.family', 'mR', 'EntrezID')]  
}


