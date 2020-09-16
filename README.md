# omXTools
A collection of scripts for recurring omics-related analyses.

### List of tools:
   - samtobam.sh
   - velvet_stat_summarizer.sh

### Usage and details
   - [contig_stats.pl](https://github.com/athieffry/omXTools/blob/master/contig_stats.pl)<br>
Summarize a collection of descriptive statistics for a FASTA file, typically for initial assessment of novo assemblies (contigs).
Output to stdout is tab-delimited with columns for: input file, total number of sequences, total number of bases, a collection of sequence length stats (min/max/average/median/N50), GC% content, and N% content.<br>
**Requires:** perl<br>
(Initial script was found somewhere on the internet (source lost) and modified to handle several input files)<br>
```
Usage incoming soon
```
   - [samtobam.sh](https://github.com/athieffry/omXTools/blob/master/samtobam.sh)<br>
Converts a SAM alignment file into its binary sorted and indexed BAM version.<br>
**Requires:** samtools, nproc (for Mac users, install nproc by running '_brew install coreutils_')
```
NAME
    samtobam.sh - converts SAM into sorted and indexed BAM.
SYNOPSIS"
    samtobam.sh [-hvk] [-p <int>] -i <input.sam>
DESCRIPTION
    Takes a SAM alignment file and converts it into BAM format.
    The resulting BAM is then sorted and indexed.
OPTIONS
    -i  SAM input file.
    -p  Specify the number of CPUs to be used (integer).
        Default: detected from system.
    -v  Enable verbose output to stdout.
        Default: silent.
    -k  Keep original SAM file.
        Default: removes SAM file.
    -h  Show this help.
Note: this script will create (temporary) intermediate files.
Be sure to have enough disk space.
```

- [velvet_stat_summarizer.sh](https://github.com/athieffry/omXTools/blob/master/velvet_stat_summarizer.sh)<br>
Gathers main [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/) de novo assembly stats into one convenient table.<br>
```
NAME
    velvet_stat_summarizer.sh - summarizes main velvet de novo assembly stats
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
