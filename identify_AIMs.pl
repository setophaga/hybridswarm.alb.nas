use strict ; 
use warnings ; 
use Getopt::ArgParse ;

### parser for importing options
my $ap = Getopt::ArgParse->new_parser(
       prog        => 'Identify AIMs',
       description => 'This program takes angsd output and an mpileup to produce AHMM input files'
) ;

### import ANGSD file 
$ap->add_arg('--ANGSD', '-A', required => 1) ;

### import mpileup file
$ap->add_arg('--mpileup',required => 1) ; 

### import mpileup file
$ap->add_arg('--output',required => 1) ;

### minimum distance between snps
$ap->add_arg('-m', default => 100, dest => 'min_distance' ) ; 

### mean recombination rate per basepair in morgans/bp
$ap->add_arg('-r', default => 3e-8 ) ; 

### min samples of pop1 
$ap->add_arg('--min_p1', default => 1 ) ; 

### min samples of pop2 
$ap->add_arg('--min_p2', default => 1 ) ;

### min allele frequency difference
$ap->add_arg('--freq_diff', default => 1 ) ;

## get args
my $ns = $ap->parse_args();

### data object to store ancestry informative sites
my %aim ; 

### import AIMs from ANGSD File
print STDERR "Reading Reference Panel Input\t", $ns->ANGSD, "\n" ; 
my $total_aims = 0 ; 

open IN, "<", $ns->ANGSD ;
<IN> ; 
while (<IN>) { 
	chomp ; 
	$_ =~ s/\"//g ;
	my @split = split ( /\,/, $_ ) ; 

	### require a minimum of 6 chromosomes and a fixation 
	if ( $split[2] > $ns->min_p1 && $split[7] > $ns->min_p2 && abs( $split[4] - $split[9] ) >= $ns->freq_diff ) { 
		if ( $split[4] == 1 ) { 
			$aim{$split[0]}{$split[1]} = $split[3]."\t".$split[5]."\t".$split[2]*$split[4]."\t".$split[2]*(1-$split[4])."\t".$split[7]*$split[9]."\t".$split[7]*(1-$split[9]) ;
		}
		else {
			$aim{$split[0]}{$split[1]} = $split[5]."\t".$split[3]."\t".$split[2]*(1-$split[4])."\t".$split[2]*$split[4]."\t".$split[7]*(1-$split[9])."\t".$split[7]*$split[9] ;
		}
		$total_aims ++ ; 
	}
}
close IN ; 
print STDERR "Found ", $total_aims, " AIMs\n" ; 


### read from mpileup file, do not give samtools a reference when you make this
print STDERR "Populating SNP matrix from\t", $ns->mpileup, "\nOutputting to\t", $ns->output, "\n" ;
my $last = -1 * $ns->min_distance ; 
my $chrom = "NA" ; 
open IN, "<", $ns->mpileup ; 
open OUT, ">", $ns->output ;
while (<IN>) { 

	chomp $_ ; 
	my @split = split ( /\t/, $_ ) ; 

	### update to the next chromosome
	if ( $split[0] ne $chrom ) { 
		$chrom = $split[0] ; 
		$last = -1 * $ns->min_distance ; 
	}

	### check if the site exists and is sufficiently far from the last site
	if ( !exists( $aim{$split[0]}{$split[1]} ) || $split[1] - $last < $ns->min_distance ) { 
		next ;
	}

	### get the site information 
	my @a = split ( /\t/, $aim{$split[0]}{$split[1]} ) ;

	### this uses a uniform recombination rate, 1e-8 per site, or 1cm/mb 
	print OUT $split[0], "\t", $split[1], "\t", $a[2], "\t", $a[3], "\t", $a[4], "\t", $a[5], "\t", ( $split[1] - $last ) * $ns->r ; 
	$last = $split[1] ;

	### iterate through the mpileup and add sites to our data input file 
	for ( my $i = 4 ; $i < $#split + 1 ; $i += 3 ) { 
		my $c1 = 0 ; my $c2 = 0 ; 
		for ( my $p = 0 ; $p < length( $split[$i] ) ; $p ++ ) { 
			if ( substr( $split[$i], $p, 1 ) =~ m/$a[0]/i ) { 
				$c1 ++ ; 
			}
			elsif ( substr( $split[$i], $p, 1 ) =~ m/$a[1]/i ) { 
        	                $c2 ++ ; 
			}
		}
		print OUT "\t", $c1, "\t", $c2 ;		
	}
	print OUT "\n" ;
}
close OUT ; 
