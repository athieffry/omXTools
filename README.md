# omXTools
Collection of scripts for recurring omics-related analyses.

- [samtobam.sh](omXTools/samtobam.sh)<br>
Converts a SAM alignment file into its binary sorted and indexed BAM version.
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
