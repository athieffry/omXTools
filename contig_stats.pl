#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum min max);
use Getopt::Long;
use File::Basename;

my $As = 0;
my $Ts = 0;
my $Gs = 0;
my $Cs = 0;
my $Ns = 0;

# Parameter variables
my $file;
my $helpAsked;

GetOptions("i=s" => \$file,
	   "h|help" => \$helpAsked );

if(defined($helpAsked)) { prtUsage(); exit; }
if(!defined($file)) { prtError("No input provided"); }

my ($fileName, $filePath) = fileparse($file);
printf "%s\t%s\n", $fileName,$filePath ;


open(I, "<$file") or die "Can not open file: $file\n";

my @len = ();
my $prevFastaSeqId = "";
my $fastaSeqId = "";
my $fastaSeq = "";

while(my $line = <I>) {
	chomp $line;
	if($line =~ /^>/) {
		$prevFastaSeqId = $fastaSeqId;
		$fastaSeqId = $line;
		if($fastaSeq ne "") {
			push(@len, length $fastaSeq);
			baseCount($fastaSeq);
		}
		$fastaSeq = "";
	}
	else {
		$fastaSeq .= $line;
	}
}
if($fastaSeq ne "") {
	$prevFastaSeqId = $fastaSeqId;
	push(@len, length $fastaSeq);
	baseCount($fastaSeq);
}

my $totalReads = scalar @len;
my $bases = sum(@len);
my $minReadLen = min(@len);
my $maxReadLen = max(@len);
my $avgReadLen = sprintf "%0.2f", $bases/$totalReads;
my $gc_content = sprintf "%0.2f", (($Gs+$Cs)/$bases*100);
my $n_content = sprintf "%0.2f", ($Ns/$bases*100);
my $medianLen = calcMedian(@len);
my $n25 = calcN50(\@len, 25);
my $n50 = calcN50(\@len, 50);
my $n75 = calcN50(\@len, 75);
my $n90 = calcN50(\@len, 90);
my $n95 = calcN50(\@len, 95);

printf "input_file\tcontigs\tbases\tminlen\tmaxlen\tmeanlen\tmedianlen\tn50\tgc\tNs\n";
printf "%s\t%d\t%d\t%d\t%d\t%0.2f\t%d\t%d\t%0.2f\t%0.2f\n",
$file,$totalReads,$bases,$minReadLen,$maxReadLen,$avgReadLen,$medianLen,$n50,$gc_content,$n_content ;

exit;

sub calcN50 {
	my @x = @{$_[0]};
	my $n = $_[1];
	@x=sort{$b<=>$a} @x;
	my $total = sum(@x);
	my ($count, $n50)=(0,0);
	for (my $j=0; $j<@x; $j++){
        $count+=$x[$j];
        if(($count>=$total*$n/100)){
            $n50=$x[$j];
            last;
        }
	}
	return $n50;
}

sub calcMedian {
	my @arr = @_;
	my @sArr = sort{$a<=>$b} @arr;
	my $arrLen = @arr;
	my $median;
	if($arrLen % 2 == 0) {
		$median = ($sArr[$arrLen/2-1] + $sArr[$arrLen/2])/2;
	}
	else {
		$median = $sArr[$arrLen/2];
	}
	return $median;
}

sub baseCount {
	my $seq = $_[0];
	my $tAs += $seq =~ s/A/A/gi;
	my $tTs += $seq =~ s/T/T/gi;
	my $tGs += $seq =~ s/G/G/gi;
	my $tCs += $seq =~ s/C/C/gi;
	$Ns += (length $seq) - $tAs - $tTs - $tGs - $tCs;
	$As += $tAs;
	$Ts += $tTs;
	$Gs += $tGs;
	$Cs += $tCs;
}

sub prtHelp {
	print "\n$0 options:\n\n";
	print "### Input reads/sequences (FASTA) (Required)\n";
	print "  -i <Read/Sequence file>\n";
	print "    Read/Sequence in fasta format\n";
	print "\n";
	print "### Other options [Optional]\n";
	print "  -h | -help\n";
	print "    Prints this help\n";
	print "\n";
}

sub prtError {
	my $msg = $_[0];
	print STDERR "+======================================================================+\n";
	printf STDERR "|%-70s|\n", "  Error:";
	printf STDERR "|%-70s|\n", "       $msg";
	print STDERR "+======================================================================+\n";
	prtUsage();
	exit;
}

sub prtUsage {
	print "\nUsage: perl $0 <options>\n";
	prtHelp();
}
