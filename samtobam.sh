#!/usr/bin/env bash
# Axel Thieffry - August 2020
set -e
set -u
VERSION="1.1"
VERBOSE=0
KEEP=0
DEBUG=0

# detect if Mac or Linux OS
if [[ "$(uname -s)" == "Darwin" ]]
then
	CPU_SYSTEM=$(sysctl -n hw.logicalcpu)
else
	CPU_SYSTEM=$(nproc)
fi


# USAGE
usage() {
	echo "NAME"
	echo "    samtobam.sh - converts SAM into sorted and indexed BAM."
	echo "SYNOPSIS"
	echo "    samtobam.sh [-hvk] [-p <int>] -i <input.sam>"
	echo "DESCRIPTION"
	echo "    Takes a SAM alignment file and converts it into BAM format."
	echo "    The resulting BAM is then sorted and indexed."
	echo "OPTIONS"
	echo "    -i  SAM input file."
	echo "    -p  Specify the number of CPUs to be used (integer)."
	echo "        Default: detected from system."
	echo "    -v  Enable verbose output to stdout."
	echo "        Default: silent."
	echo "    -k  Keep original SAM file."
	echo "        Default: removes SAM file."
	echo "    -h  Show this help."
	echo ""
	echo "Note: this script will create (temporary) intermediate files."
	echo "Be sure to have enough disk space."
	echo ""
	echo "Author: Axel Thieffry"
	echo "version: $VERSION"
 	}

# check if any argument has been provided
if [[ $# -eq 0 ]]
then
	usage
	echo "samtobam.sh requires an input file" >&2
	exit 1
fi

# read and parse arguments
while getopts ":i:p:kvh" opts
	do
	  case "$opts" in
		i)
		  INPUT=$OPTARG

		  # check that SAM input file exists
		  if ! [[ -f "$INPUT" ]]
		  then
		  	echo "SAM input file not found ($INPUT)." >&2
			exit 1
		  fi

		  # check that input_file.bam does not already exist
		  if [[ -f ${INPUT%.*}.bam ]]
		  then
			echo "BAM file seems to already exist for this input ($INPUT)." >&2
			exit 1
		  fi
		  ;;

		p)
		  CPU_USER=${OPTARG%.*}

		  # check that number of CPU is greater than 1
		  if [[ $CPU_USER -lt "1" ]]
		  then
			echo "Provided CPU number ($CPU_USER) must be greater than 0." >&2
			echo ""
			usage
			exit 1
		  fi

		  # check that number of CPU is not greater than what the system has
		  if [[ $CPU_USER -gt $CPU_SYSTEM ]]
		  then
			echo "Provided CPU number ($CPU_USER) cannot exceed available CPUs ($CPU_SYSTEM)." >&2
			usage
			exit 1
		  fi
		  ;;
		v)
		  VERBOSE=1 ;;
		k)
		  KEEP=1 ;;
		h)
		  usage
		  exit 0 ;;
		\?)
		  echo "Invalid input/option/argument." >&2
		  echo ""
		  usage
		  exit 1 ;;
		:)
		  echo "Invalid option: $OPTARG" >&2
		  ;;
	  esac
	done

# if user provided senseful number of CPUs, use it, otherwise fall back on system
if [[ -z ${CPU_USER+x} ]]
then
	CPU=$CPU_SYSTEM
else
	CPU=$CPU_USER
fi


# low-tech DEBUG
if [[ $DEBUG -eq 1 ]]
then
	echo ""
	echo "//DEBUG//"
	echo "INPUT: ${INPUT:-}"
	echo "CPU_SYSTEM: ${CPU_SYSTEM:-}"
	echo "CPU_USER: ${CPU_USER:-}"
	echo "CPU: ${CPU:-}"
	echo "VERBOSE: ${VERBOSE:-}"
	echo "KEEP: ${KEEP:-}"
fi

# check if samtools can be found in PATH
if ! [[ -x "$(command -v samtools)" ]]
then
	echo "SAMTOOLS command could not be found." >&2
	echo "Make sure samtools is installed and accessible within PATH." >&2
	exit 1
fi


# execute SAM to BAM conversion
# -----------------------------
# 0. make filenames
SAM=$INPUT
BAM=${INPUT%.*}.bam
SORTEDBAM=${INPUT%.*}.sorted.bam

# 1. sam to bam
if [[ $VERBOSE -eq 1 ]]
	then echo "Converting $SAM into $BAM using $CPU cores..."
fi

samtools view -b $SAM -o $BAM -@ $CPU 2> /dev/null

# 2. remove input.sam
if [[ $KEEP -eq 0 ]]
then
	
	if [[ $VERBOSE -eq 1 ]]
	then
		echo "Removing input SAM..."
	fi
	
	echo "Deleting original $SAM file ..."
	rm $INPUT
fi

# 3. sort and rename
if [[ $VERBOSE -eq 1 ]]
then
	echo "Sorting and renaming $BAM ..."
fi

samtools sort $BAM -l 9 -o $SORTEDBAM -@ $CPU  2> /dev/null
mv $SORTEDBAM $BAM  2> /dev/null

# 4. index
if [[ $VERBOSE -eq 1 ]]
then
	echo "Indexing $BAM ..."
fi

samtools index $BAM -@ $CPU  2> /dev/null

if [[ $VERBOSE -eq 1 ]]
then
	echo "~DONE~"
fi

exit 0
### EOF ###
