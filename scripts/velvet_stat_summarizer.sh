#!/usr/bin/env bash
# Axel Thieffry - August 2020
set -e
set -u
VERSION="1.0"
DEBUG=0
# USAGE
usage() {
		echo "NAME"
		echo "    velvet_stat_summarizer.sh - gathers main Velvet de novo assembly stats into one convenient table."
		echo "SYNOPSIS"
		echo "    velvet_stat_summarizer.sh [OPTIONS] directory_1 [directory_2 directory_3 ...]"
		echo " DESCRIPTION"
		echo "    Takes one or several velvet denovo output directory and create a single table"
		echo "    of main assembly statistics gathered from the Log files. Assumes that directory"
		echo "    names contain the sample and the k-mer length separated by underscore."
		echo "    For example: CEN.PK2_67."
		echo "OPTIONS"
		echo "    -h  Show this help."
		echo ""
		echo "Author: Axel Thieffry"
		echo "version: $VERSION"
		}

# fail function
function fail {
	echo "$@" >&2
	exit 1
}

# parse options
while getopts ":h" opt; do
	case $opt in
		h)
          usage
          exit 0
          ;;
		\?)
		  echo -e "Invalid option: -$OPTARG\n" >&2
		  usage >&2
		  ;;
    esac
done

# check that all directories exist
for f in "$@"; do
	[[ -e $f ]] || fail "Could not find directory: '${f}'"
done
[[ DEBUG -eq 0 ]] || echo "INPUT FILES: $@"

# check that all directories have a "Log" file
for f in "$@"; do
	[[ -e $f/Log ]] || fail "Could not find Log file: '${f}/Log'"
done

# make column names
echo -ne "sample\tkmer\tmed_cov_depth\tnodes\tn50\tmax\ttotal_len\treads_used\treads_tot"

# compiles Log outputs
# this is very bad coding, probably error prone (~shell & ~OS)
# improve this parsing
for f in "$@" ; do
	echo -ne "\n${f/_/\t}\t"
	tail -n 2 "$f/Log" | head -1 | grep -Eo '[0-9]+([.][0-9]+)?' | tr "\n" "\t "
	tail -n 1 "$f/Log" | sed "s/n50//" | grep -Eo '[0-9]+([.][0-9]+)?' | tr "\n" "\t"
done