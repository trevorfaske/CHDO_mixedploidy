#!/usr/bin/perl
#
# Edited by Trevor to change regex to accept contig annotation of the GM genome


# Time-stamp: <Friday, 25 March 2016, 13:21 MDT -- monia>

# This is a copy of Alex's script for use with vcf-files having extra FORMAT fields and "./." uncalled
# genotypes. 

# This scripts converts a vcf file to a simpler format for downstream
# analysis. I am calling this format multiple population genetoype
# likelihood (mpgl). The first line lists: the number of individuals,
# loci.  There is a following line, which we do not use.  Instead,
# look in the indiv_ids.txt for the order and identifiers for inds.
# with one entry per individual. This is followed by one line per SNP
# that gives the SNP id (scaffold, position) and the phred-scaled
# genotype likelihoods, three per individual.

#
# USAGE: vcf2mpgl.pl in.vcf
# Follow with cat header_mpglvariants_0.9.txt variants_0.9.mpgl > entropy_in_variants_0.9_all.mpgl

use warnings;

my @line = ();
my $word;
my %keepind;
my $nind = 0;
my $nloc = 0;
my $nline = 0;


# my $idin = shift @ARGV;
# open (ID, "$idin") or die "Could not read the id file\n";
# ## individual names we want to keep
# chomp(<ID>);
# @line = split / /;
# foreach $word (@line){
#     $keepind{$word} = 1;
# }
# push (@ids,$_);
# while (<ID>){
#     chomp;
#     push (@ids,$_);
#     $nline++;
# }
# close (ID);


my $in = shift (@ARGV);
open (IN, $in) or die "Could not read the vcf file\n";
my $out = $in;
if ($out =~ s/vcf/mpgl/){
    open (OUT, "> $out") or die "Could not write $out\n";
}
else {
    die "File extension was not vcf\n";
}

while (<IN>){
    chomp;

    ## read genetic data lines, write gl
    ## match both scaffold111 and C113 and exclude loci with multiple alternative alleles
    if (m/^(\w+\d+.\d+)\s+(\d+)\s+[.\w]+\s+[AGCT]\s+[AGCT]\s+/){ 
	print OUT "$1:$2 ";   
	@line = split /\s+/;
	$i = 0;
	foreach $word (@line){
	    if($word =~ m#^([.\d])/([.\d])\:(\d+),(\d+),(\d+)#){
		$a1 = $1;
		$a2 = $2;
		$g1 = $3;
		$g2 = $4;
		$g3 = $5;
		## (16 Feb 2016) updated both patterns above to use #
		## for better readibility of the regexp, added . to
		## [.\d] set to capture genotypes of ./. in vcf, which
		## is a new output in samtools/bcftools 1.x relative
		## to 0.1.19
		## Monia changed the end of the matching pattern to match more FORMAT fields in the
		##vcf.
		if( $a1 eq "." && $a2 eq "."){
		    $g1 = 0;
		    $g2 = 0;
		    $g3 = 0;
		}
		
		print OUT " $g1 $g2 $g3";
		$i++; #This was originally in the wrong place
	    }

	}
	if($i > 0){
	    print OUT "\n";
	    $nloc++;
	}
    }
    ## get individual ids
    elsif (m/^#CHROM/){
	   @line = split /\s+/;
	   foreach $word (@line){
	       if ($word =~ m/[A-Z]+_\d+/){  ### delicate, depends on pattern of ind names
		   push (@inds, $word);
	       }
	   }
	   $nind = scalar @inds;
       }
}
close (OUT);

open (INDS, "> indiv_ids.txt") or die "Could not write the id file\n";
print INDS "ind,pop,taxon\n";
foreach $word (@inds){
    $word =~ m/([A-Za-z\d-]+).([A-Za-z\d-]+)/;
    print INDS "$word,$1,$2\n";
}
close (INDS);

$out =~ s/mpgl$/txt/;
open (INDS, "> header_mpgl$out") or die "Could not write the id file\n";
print INDS "$nind $nloc\n";
print INDS "discarded ind line\n";
close (INDS);

print "Number of loci: $nloc; number of individuals $nind\n";
