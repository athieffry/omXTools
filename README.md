# omXTools
A collection of scripts for recurring omics-related analyses.

### List of tools:
   - samtobam.sh
   - velvet_stat_summarizer.sh

### Usage and details
   - [samtobam.sh](https://github.com/athieffry/omXTools/blob/master/samtobam.sh)<br>
Converts a SAM alignment file into its binary sorted and indexed BAM version.<br>
**Requires:** samtools
```
USAGE: $ samtobam.sh [-hvk] [-p <int>] -i <input.sam>
  -i     SAM input file.
  -p     Specify the number of CPUs to be used (integer).
         Default: detected by calling nproc.
  -v     Enable verbose output to stdout.
         Default: silent.
  -k     Keep original SAM file.
         Default: removes SAM file.
  -h     Show this help.
```

- [velvet_stat_summarizer.sh](https://github.com/athieffry/omXTools/blob/master/velvet_stat_summarizer.sh)<br>
Gathers main [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/) de novo assembly stats into one convenient table.<br>
**Requires:** nproc (for Mac users, install nproc by running '_brew install coreutils_')
```
SYNOPSIS
    velvet_stat_summarizer.sh [OPTIONS] directory_1 [directory_2 directory_3 ...]
DESCRIPTION
    Takes one or several velvet denovo output directory and create a single table
    of main assembly statistics gathered from the Log files. Assumes that directory
    names contain the sample and the k-mer length separated by underscore.
    For example: CEN.PK2_67.
OPTIONS
    -h   Show this help.
```
