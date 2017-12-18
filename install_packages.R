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
# Description: Install and load packages BNSA needs
# Author: Xin Zhu
# Date: 2014-01-21
# Comments: 
###################################################################
install_packages <- function(package)
{
	# Retrieve the installed packages information
	installed_packages <- installed.packages()
	
	source("http://bioconductor.org/biocLite.R")

	# If hgu133a.db has not been installed, install it.
	#packages_info <- grep("hgu133a.db", installed_packages)
	#if(length(packages_info) == 0)
	#{
		# Install package from Bioconductor
	#	biocLite("hgu133a.db")
	#}
	
	# If targetscan.Hs.eg.db has not been installed, install it.
	packages_info <- grep(package, installed_packages)
	if(length(packages_info) == 0)
	{
		# Install package from Bioconductor
		biocLite(package)
	}
	
	#packages_info <- grep("Rgraphviz", installed_packages)
	#if(length(packages_info) == 0)
	#1{
	#	# Install package frobiobiom Bioconductor
	#	biocLite("Rgraphviz")
	#}

	#If plyr has not been installed, install it.
	packages_info <- grep("plyr", installed_packages)
	if(length(packages_info) == 0)
	{
		# Install package
		install.packages("plyr")
	}
	
	# If bnlearn has not been installed, install it.
	packages_info <- grep("bnlearn", installed_packages)
	if(length(packages_info) == 0)
	{
		install.packages("bnlearn")
	}

	#packages_info <- grep("base", installed_packages)
	#if(length(packages_info) == 0)
	#{
	#	install.packages("base")
	#}

	#packages_info <- grep("Biobase", installed_packages)
	#if(length(packages_info) == 0)
	#{
	#	install.packages("Biobase")
	#}

	#packages_info <- grep("igraph", installed_packages)
	#if(length(packages_info) == 0)
	#{
	#	install.packages("igraph")
	#}

	# Load packages
	library(hgu133a.db)
	library(package, character.only = TRUE)
	library(plyr)
	#library(Rgraphviz)
	library(bnlearn)
	#library(base)
	#library(Biobase)
	library(graph)
	library(igraph)

}
