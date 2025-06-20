#!/usr/bin/env perl
use strict;
use warnings;

#files must be in same order (rows)
if (!@ARGV or scalar @ARGV != 2) {
   print "usage: makeTable.pl fileWithList.txt output.tab\n";
   print "List should contain tab separated directory peaks.bed file.bw \n";
   exit;
}

#files in list should be directory peak bam
my $flist = shift @ARGV; #file with list of files to use
my $outtab = shift @ARGV;

open(FH, $flist) or die "Couldn't open $flist, $!\n";
my $bams = '';
my @pks;
while (<FH>) {
   chomp;
   if (/^\s*$/) { next; } #skip blanks
   my @f = split(/\s+/);
   push(@pks, "$f[0]/$f[1]");
   $bams = "$bams $f[0]/$f[2]";
}
close FH or die "Couldn't close $flist, $!\n";

print "BAMS $bams\n\n";
my $comm = "cat " . join(" ", @pks) . " | sort -k1,1 -k2,2n > combinedPeaks.bed";
#system($comm) == 0 or die "Couldn't cat peak files, $!\n";

#$comm = "bedtools merge -i combinedPeaks.bed |";
#open(FH, $comm) or die "Couldn't merge peak files, $!\n";
#open(OUT, ">", "combinedPeaks.merged.bed") or die "Couldn't open outfile, $!\n";
#my $i = 1;
#while (<FH>) {
#chomp; 
#print OUT "$_\t$i\n";
#$i++;
#} 
#close OUT or die "Couldn't close outfile, $!\n";
#close FH or die "Couldn't finish merge peak files, $!\n";

#use bams, bigwigs are fold-change or p-value not counts
#bedtools multicov will generate table from list of bams
$comm = "bedtools multicov -bams $bams -bed combinedPeaks.merged.bed > $outtab";
system($comm) == 0 or die "Couldn't get signals, $!\n";

#outtab has bed fields plus scores in order of inputs
exit;
