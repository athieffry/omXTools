#!/bin/bash
# Axel Thieffry - August 2020
# v0.1

# get number of cpus from system
CPUS=$(nproc)
echo "${CPUS}"
	# see if user provided cpu number
	# read -@ parameter
	# if some -@: override CPU number

# check if samtools is found in PATH
	# if not throw error "Samtools not found"

# read input SAM file
	# check if equivalent file in BAM already exists
	# if yes: throw error "BAM version of $input already exists"

# check if verbose mode is activated (by default nothing is said except errors)
	# verbose=0
	# if -v flag is found:
		# verbose=1

# execute SAM to BAM conversion
# -----------------------------
# sam to bam
samtools view
# sort and rename
samtools sort
# index

### EOF ###
