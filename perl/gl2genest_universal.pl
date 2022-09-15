#!/usr/bin/perl

## This script returns a genotype matirx (locus by ind) with genotype
## means or modes (point estimates) from a genotype likelihood file
## now returns NA if all three genotype are equally likely

## USAGE: perl gl2genest.pl file.gl [mean/mode]
use warnings;

unless(scalar @ARGV ==2){
    die "USAGE: perl gl2genest.pl file.gl [mean/mode]";
}

$in = shift (@ARGV);
$summarytype = shift (@ARGV);

open (IN, $in) or die;
$out = $in;
$out =~ s/mpgl$/txt/;
open (OUT, ">", "pntest_"."$summarytype"."_$out") or die "Failed to open";

while (<IN>){
    chomp;
    if (s/^([\w.\d]+:\d+)\s+//){ ## this line has genotype data, get rid of locus id
	$locusid = $1;
	@line = split /\s+/;

	if(scalar @line % 3){  ## don't have triplets for some reason: would happen if more than two alleles at a locus
	    m/^(.{40})/;
	    print scalar @line, " >>> $locusid\n";
	    print "$1 ... ", $line[$#line - 2], ",", $line[$#line - 1], ",", $line[$#line], "\n";
	}
	else{
	    foreach $i (0..$#line){
		$line[$i] = 10 ** ($line[$i]/-10);
	    }
	    
	    for($i=0; $i<scalar @line; $i +=3){
		if($line[$i] == $line[$i+1] && $line[$i] == $line[$i+2]){
		    $gest[$i/3] = 'NA';	
		}
		elsif($summarytype eq 'mean'){
		    $sum = $line[$i] + $line[$i+1] + $line[$i+2];
		    $gest[$i/3] = sprintf("%.5f", $line[$i+1]/$sum + 2*($line[$i+2]/$sum));  ## genotype on scale of 0 to 2
		}
		else{
		    $gest[$i/3] = &whichmax($line[$i], $line[$i+1], 2*$line[$i+2]);
		}
	    }
	    print OUT "@gest\n";
	}
    }
}
close (IN);
close (OUT);


sub whichmax{
    @n = @_;
    my $max = 0; 
    my $index = -1;
    foreach $i (0..$#n){
	if($n[$i] > $max){
	    $max = $n[$i];
	    $index = $i;
	}
    }
    return($index);
}
